#!/usr/bin/env python
# Writer/reader registry of the /input/potential nodes (upside3.md 2.2).
#
# Static pass: every node type registered in src/ (RegisterNodeType, RegisterRotamerSidechain) with the datasets and
# attributes its constructor reads, following the members and helpers it hands the group to; and every function in
# py/ that creates or extends a group under /input/potential, with what it writes there.
# Dynamic pass: the engine library is constructed on each fixture under an HDF5 trace (h5trace.dylib, interposed on
# H5Dopen2, H5Aopen_by_name, H5Aexists_by_name, ...), so the names the engine really asked for are compared with the
# file's own contents.  Files are opened read-only.  The trace library is macOS-only (DYLD interposing), built with
#   clang -dynamiclib -I$HDF5_ROOT/include -L$HDF5_ROOT/lib -lhdf5 -o h5trace.dylib h5trace.c
# in scratchpad/upside3_audit/registry/; without it the dynamic pass is skipped.
#
# Writes scratchpad/upside3_audit/registry_<date>.md and, beside it, registry_<date>.json.
import os, re, ast, sys, copy, glob, json, difflib, datetime, subprocess
from collections import defaultdict

import numpy as np
import h5py

UPSIDE_HOME = os.environ['UPSIDE_HOME']
SRC_DIR = os.path.join(UPSIDE_HOME, 'src')
PY_DIR = os.path.join(UPSIDE_HOME, 'py')
AUDIT_DIR = os.path.join(UPSIDE_HOME, 'scratchpad', 'upside3_audit')
WORK_DIR = os.path.join(AUDIT_DIR, 'registry')
TRACE_LIB = os.path.join(WORK_DIR, 'h5trace.dylib')
ENGINE_LIB = os.path.join(UPSIDE_HOME, 'obj', 'libupside.dylib')
STAMP = datetime.date.today().strftime('%Y%m%d')
OUT_MD = os.path.join(AUDIT_DIR, 'registry_%s.md' % STAMP)
OUT_JSON = os.path.join(AUDIT_DIR, 'registry_%s.json' % STAMP)

# The fixtures named by the audit plan, then older single-purpose configs that carry the advanced_config and
# hybrid nodes the four fixtures lack.  Those were written by the 2026-10-19/20 writers, so they are reader
# evidence first; a writer-side difference against them may be the writer's later change.
FIXTURES = (sorted(glob.glob(os.path.join(AUDIT_DIR, 'fixtures_w3', '*.up'))) +
            [os.path.join(UPSIDE_HOME, 'scratchpad/bilayer_ref/upside_from_gmx_eq/1rkl_gmxeq/run/checkpoints/'
                          '1rkl.stage_7.0.up')])
SUPPLEMENTARY = [os.path.join(AUDIT_DIR, p) for p in
                 ('writers_fix/tc.up', 'writers_fix/rg.up', 'lateral/1qhj_memb_lat.up', 'remd_check/nail.up')]

ENGINE_ATTRS = {'@arguments': True, '@integrator_level': False}   # read by initialize_engine_from_hdf5 for every node
PYTABLES_ATTRS = {'CLASS', 'TITLE', 'VERSION', 'EXTDIM', 'FLAVOR', 'FILTERS', 'PYTABLES_FORMAT_VERSION'}


# ----------------------------------------------------------------------------------------------------------------
# C++ readers
# ----------------------------------------------------------------------------------------------------------------

def strip_comments(text):
    out, i, n = [], 0, len(text)
    while i < n:
        if text.startswith('//', i):
            j = text.find('\n', i)
            j = n if j < 0 else j
            out.append(' '*(j-i)); i = j
        elif text.startswith('/*', i):
            j = text.find('*/', i+2)
            j = n if j < 0 else j+2
            out.append(re.sub(r'[^\n]', ' ', text[i:j])); i = j
        elif text[i] in '"\'':
            q, j = text[i], i+1
            while j < n and text[j] != q:
                j += 2 if text[j] == '\\' else 1
            out.append(text[i:j+1]); i = j+1
        else:
            out.append(text[i]); i += 1
    return ''.join(out)

def match_close(text, i, angle=False):
    # index of the bracket closing text[i]; angle brackets nest only when asked for
    opens = '([{<' if angle else '([{'
    closes = ')]}>' if angle else ')]}'
    depth, j, n = 0, i, len(text)
    while j < n:
        c = text[j]
        if c in '"\'':
            k = j+1
            while k < n and text[k] != c:
                k += 2 if text[k] == '\\' else 1
            j = k+1
            continue
        if c in opens: depth += 1
        elif c in closes:
            depth -= 1
            if depth == 0: return j
        j += 1
    raise ValueError('unbalanced bracket at %d' % i)

def split_args(s, angle=False):
    out, depth, cur, j = [], 0, [], 0
    opens = '([{<' if angle else '([{'
    closes = ')]}>' if angle else ')]}'
    while j < len(s):
        c = s[j]
        if c in '"\'':
            k = j+1
            while k < len(s) and s[k] != c:
                k += 2 if s[k] == '\\' else 1
            cur.append(s[j:k+1]); j = k+1
            continue
        if c in opens: depth += 1
        elif c in closes: depth -= 1
        if c == ',' and depth == 0:
            out.append(''.join(cur).strip()); cur = []
        else:
            cur.append(c)
        j += 1
    if ''.join(cur).strip(): out.append(''.join(cur).strip())
    return out

def line_of(text, pos):
    return text.count('\n', 0, pos) + 1

def base_name(type_text):
    t = re.sub(r'\b(const|typename|struct|class|static|mutable)\b', '', type_text).strip()
    m = re.match(r'(?:\w+::)*(\w+)', t)
    return m.group(1) if m else t


class CppSource:
    def __init__(self):
        self.files = {}
        for path in sorted(glob.glob(os.path.join(SRC_DIR, '*.cpp')) + glob.glob(os.path.join(SRC_DIR, '*.h'))):
            self.files[os.path.relpath(path, UPSIDE_HOME)] = strip_comments(open(path).read())
        self.classes = {}                # name -> [dict(file, line, start, end, tparams, members)]
        self.ctors = defaultdict(list)   # class name -> [dict(file, line, params, gparams, code, pos, outside)]
        self.funcs = defaultdict(list)   # free function or method name -> [...]
        self.index_classes()
        self.index_callables()

    def index_classes(self):
        # name -> [definitions]; one name may be defined in several translation units
        for fn, t in self.files.items():
            for m in re.finditer(r'\b(struct|class)\s+(\w+)\s*(final\s*)?(:[^;{]*)?\{', t):
                name = m.group(2)
                tparams = []
                hm = re.search(r'template\s*<([^;{}]*)>\s*$', t[max(0, m.start()-300):m.start()])
                if hm:
                    for p in split_args(hm.group(1), angle=True):
                        ids = re.findall(r'\w+', p.split('=')[0])
                        if ids: tparams.append(ids[-1])
                start = m.end()-1
                end = match_close(t, start)
                self.classes.setdefault(name, []).append(dict(
                    name=name, file=fn, line=line_of(t, m.start()), start=start, end=end, tparams=tparams,
                    members=self.members(t[start+1:end])))

    def class_def(self, name, near_file=None):
        defs = self.classes.get(name, [])
        same = [d for d in defs if d['file'] == near_file]
        return (same or defs or [None])[0]

    def ctors_of(self, name, near_file=None):
        d = self.class_def(name, near_file)
        if d is None: return []
        return [r for r in self.ctors.get(name, []) if r['file'] == d['file'] and
                (d['start'] < r['pos'] < d['end'] or r['outside'])]

    @staticmethod
    def members(body):
        out = {}
        for m in re.finditer(r'(?m)^\s*(?:mutable\s+|const\s+|static\s+)*'
                             r'((?:\w+::)*\w+\s*(?:<[^;{}()]*>)?[\s\*&]*)\s+(\w+)\s*(?:\[[^\]]*\])?\s*(?:=[^;]*)?;', body):
            typ, name = m.group(1).strip(), m.group(2)
            if base_name(typ) in ('return', 'delete', 'throw', 'goto', 'using', 'typedef'): continue
            out.setdefault(name, typ)
        return out

    def body_after_params(self, t, close):
        # from the parameter list's closing paren, past an optional initialiser list, to the body's closing brace
        j = close+1
        n = len(t)
        while j < n and t[j].isspace(): j += 1
        if t.startswith('const', j): j += 5
        while j < n and t[j].isspace(): j += 1
        if t.startswith('try', j): j += 3
        while j < n and t[j].isspace(): j += 1
        if j >= n or t[j] not in ':{': return None
        if t[j] == ':':
            j += 1
            while j < n:
                c = t[j]
                if c in '([':
                    j = match_close(t, j)+1; continue
                if c == ';': return None
                if c == '{':
                    k = j-1
                    while k > 0 and t[k].isspace(): k -= 1
                    if re.match(r'[\w>]', t[k]):
                        j = match_close(t, j)+1; continue
                    break
                j += 1
        if j >= n or t[j] != '{': return None
        return match_close(t, j)

    def index_callables(self):
        for fn, t in self.files.items():
            for m in re.finditer(r'(?<![\w.>~])(\w+)\s*(?:::\s*(\w+)\s*)?\(', t):
                name = m.group(2) or m.group(1)
                if name in ('if', 'for', 'while', 'switch', 'return', 'sizeof', 'catch', 'decltype'): continue
                op = m.end()-1
                try: close = match_close(t, op)
                except ValueError: continue
                params = t[op+1:close]
                if not re.search(r'\bhid_t\b', params): continue
                if re.search(r'"', params): continue     # a call, not a definition
                end = self.body_after_params(t, close)
                if end is None: continue
                plist = split_args(params)
                gparams = [re.search(r'hid_t\s+(\w+)', p).group(1) if re.search(r'hid_t\s+(\w+)', p) else None
                           if 'hid_t' in p else False for p in plist]
                outside = m.group(2) is not None
                rec = dict(file=fn, line=line_of(t, m.start()), params=plist, gparams=gparams,
                           code=t[m.start():end+1], pos=m.start(), outside=outside)
                owner = m.group(1) if outside else self.enclosing_class(fn, m.start())
                if name in self.classes and owner == name:
                    self.ctors[name].append(rec)
                else:
                    self.funcs[name].append(rec)

    def enclosing_class(self, fn, pos):
        best = None
        for defs in self.classes.values():
            for c in defs:
                if c['file'] == fn and c['start'] < pos < c['end']:
                    if best is None or c['end']-c['start'] < best['end']-best['start']: best = c
        return best['name'] if best else None

    def registrations(self):
        regs = []
        for fn, t in self.files.items():
            for m in re.finditer(r'\bRegister(NodeType|RotamerSidechain)\s*<', t):
                lt = m.end()-1
                gt = match_close(t, lt, angle=True)
                mm = re.match(r'\s*\w+\s*\(\s*"([^"]+)"\s*\)', t[gt+1:])
                if not mm: continue
                targs = split_args(t[lt+1:gt], angle=True)
                if m.group(1) == 'NodeType':
                    cls_expr, nargs = ','.join(targs[:-1]).strip(), int(targs[-1])
                else:
                    cls_expr, nargs = 'RotamerSidechain<%s>' % targs[0], -1
                regs.append(dict(prefix=mm.group(1), cls_expr=cls_expr, cls=base_name(cls_expr), nargs=nargs,
                                 reg_file=fn, reg_line=line_of(t, m.start())))
        return regs


H5_HELPERS = ('traverse_dset', 'traverse_string_dset', 'check_size', 'get_dset_size', 'read_attribute',
              'h5_exists', 'attribute_exists', 'attribute_exists_hybrid', 'read_string_attribute_or_default',
              'open_group', 'read_int_dataset', 'H5Aexists', 'H5Dopen2')
NOT_READERS = set(H5_HELPERS) | {'if', 'for', 'while', 'switch', 'return', 'sizeof', 'catch', 'group_name',
                                 'add_logger', 'string', 'to_string', 'h5_obj', 'logging', 'throw', 'get',
                                 'c_str', 'H5Iget_name', 'H5Gclose', 'H5Dclose', 'decltype'}

def c_literal(expr):
    # the names a C++ argument can stand for: a literal, a ternary of literals, or suffix1("x")
    e = expr.strip()
    m = re.fullmatch(r'"([^"]*)"', e)
    if m: return [m.group(1)]
    m = re.fullmatch(r'suffix1\(\s*"([^"]*)"\s*\)\.c_str\(\)', e)
    if m: return [m.group(1), m.group(1)+'1']
    lits = re.findall(r'"([^"]*)"', e)
    if lits and re.fullmatch(r'[\w\s:]*\?\s*"[^"]*"\s*:\s*"[^"]*"', e): return lits
    return None


class Reads:
    def __init__(self):
        self.dsets = {}    # name -> dict(ranks=set, dims=[], required=bool, lines=set)
        self.attrs = {}    # key -> dict(required, default, lines)
        self.probes = set()
        self.unresolved = []
        self.via = set()

    def dset(self, name, rank, dims, where, required=True):
        d = self.dsets.setdefault(name, dict(ranks=set(), dims=[], required=False, lines=set()))
        if rank is not None: d['ranks'].add(rank)
        if dims: d['dims'].append(dims)
        d['required'] |= required
        d['lines'].add(where)

    def attr(self, key, default, where):
        a = self.attrs.setdefault(key, dict(required=False, default=None, lines=set()))
        if default is None: a['required'] = True
        else: a['default'] = default
        a['lines'].add(where)

    def finish(self):
        # a dataset or attribute guarded by an existence probe is optional
        for name, d in self.dsets.items():
            if any(n in self.probes for n in name.split('|')): d['required'] = False
        for key, a in self.attrs.items():
            if key in self.probes: a['required'] = False


def resolve_group(expr, gmap):
    e = expr.strip()
    e = re.sub(r'\.get\(\)$', '', e).strip()
    while e.startswith('(') and e.endswith(')') and match_close(e, 0) == len(e)-1: e = e[1:-1].strip()
    if e in gmap: return gmap[e]
    m = re.fullmatch(r'(?:h5::)?open_group\((.*)\)', e, re.S)
    if m:
        a = split_args(m.group(1))
        if len(a) == 2:
            parent = resolve_group(a[0], gmap)
            names = c_literal(a[1])
            if parent is not None: return parent + (names[0] if names else '<%s>' % a[1]) + '/'
    return None


def collect_reads(cpp, code, gmap, fn, line0, reads, cls, tsubst, depth, seen):
    # reads made in `code` on the group variables of gmap (name -> path prefix inside the node's group)
    gmap = dict(gmap)
    for m in re.finditer(r'(?:auto|H5Obj|h5::H5Obj)\s+(\w+)\s*=\s*(?:h5::)?open_group\(\s*(\w+)\s*,\s*([^)]*)\)', code):
        parent = gmap.get(m.group(2))
        names = c_literal(m.group(3))
        if parent is not None: gmap[m.group(1)] = parent + (names[0] if names else '<%s>' % m.group(3)) + '/'
    for m in re.finditer(r'hid_t\s+(\w+)\s*=\s*(\w+)\.get\(\)', code):
        if m.group(2) in gmap: gmap[m.group(1)] = gmap[m.group(2)]

    for m in re.finditer(r'(?<![\w.])(?:h5::)?(%s)\s*(<[^;(){}]*?>)?\s*\(' % '|'.join(H5_HELPERS), code):
        helper, targs = m.group(1), (m.group(2) or '')
        op = m.end()-1
        try: args = split_args(code[op+1:match_close(code, op)])
        except ValueError: continue
        where = '%s:%d' % (fn, line0 + code.count('\n', 0, m.start()))
        tnums = re.findall(r'\d+', targs)
        if helper == 'get_dset_size':
            if len(args) < 3: continue
            g, names, rank = args[1], c_literal(args[2]), (int(args[0]) if args[0].isdigit() else None)
        elif helper in ('traverse_dset', 'traverse_string_dset'):
            g, names, rank = args[0], c_literal(args[1]), (int(tnums[0]) if tnums else None)
        elif helper == 'check_size':
            g, names = args[0], c_literal(args[1])
            dims = args[2:]
            if len(dims) == 1 and dims[0].lstrip().startswith(('{', 'std::vector', 'vector')):
                inner = re.search(r'\{(.*)\}', dims[0], re.S)
                dims = split_args(inner.group(1).strip('{}')) if inner else [None]
            rank = len(dims)
        elif helper == 'read_int_dataset':
            g, names, rank = args[0], c_literal(args[1]), 1
        elif helper == 'H5Dopen2':
            g, names, rank = args[0], c_literal(args[1]), None
        else:
            g = args[0]
        prefix = resolve_group(g, gmap)
        if prefix is None: continue
        if helper == 'read_attribute':
            if len(args) < 3: continue
            path, attr = c_literal(args[1]), c_literal(args[2])
            if not path or not attr:
                reads.unresolved.append('%s(%s) at %s' % (helper, ', '.join(args[1:3]), where)); continue
            owner = '' if path[0] == '.' else path[0]
            key = (prefix + owner).rstrip('/') + '@' + attr[0]
            reads.attr(key, args[3] if len(args) > 3 else None, where)
            continue
        if helper in ('attribute_exists', 'attribute_exists_hybrid'):
            path, attr = c_literal(args[1]), c_literal(args[2])
            if path and attr: reads.probes.add((prefix + ('' if path[0] == '.' else path[0])).rstrip('/') + '@' + attr[0])
            continue
        if helper == 'H5Aexists':
            attr = c_literal(args[1])
            if attr: reads.probes.add(prefix.rstrip('/') + '@' + attr[0])
            continue
        if helper == 'read_string_attribute_or_default':
            attr = c_literal(args[1])
            if attr: reads.attr(prefix.rstrip('/') + '@' + attr[0], args[2] if len(args) > 2 else '""', where)
            continue
        if helper == 'h5_exists':
            names = c_literal(args[1])
            if names: reads.probes.update(prefix + n for n in names)
            continue
        if helper == 'open_group':
            names = c_literal(args[1])
            if names: reads.dset(prefix + names[0] + '/', None, None, where)
            continue
        if not names:
            reads.unresolved.append('%s(%s) at %s' % (helper, ', '.join(args[:3]), where)); continue
        full = '|'.join((n if n.startswith('/') else prefix + n) for n in names)
        reads.dset(full, rank, args[2:] if helper == 'check_size' else None, where)

    if depth > 4: return
    cdef = cpp.class_def(cls, fn) if cls else None
    members = cdef['members'] if cdef else {}
    for m in re.finditer(r'(?<![\w.>])(\w+)\s*(<[^;(){}]*?>)?\s*\(', code):
        name = m.group(1)
        if name in NOT_READERS: continue
        op = m.end()-1
        try: args = split_args(code[op+1:match_close(code, op)])
        except ValueError: continue
        gargs = [(k, resolve_group(a, gmap)) for k, a in enumerate(args)]
        gargs = [(k, p) for k, p in gargs if p is not None]
        if not gargs: continue
        k, prefix = gargs[0]
        callee_type = None
        if name in members:
            callee_type = base_name(members[name])
        elif name in cpp.classes:
            callee_type = name
        if callee_type is not None:
            callee_type = tsubst.get(callee_type, callee_type)
            targets = [(callee_type, r) for r in cpp.ctors_of(callee_type, fn)]
        else:
            targets = [(cpp.enclosing_class(r['file'], r['pos']), r) for r in cpp.funcs.get(name, [])]
        for tcls, r in targets:
            if len(r['params']) <= k or not r['gparams'][k]: continue
            key = (r['file'], r['line'], prefix)
            if key in seen: continue
            seen.add(key)
            reads.via.add('%s (%s:%d)' % (callee_type or name, r['file'], r['line']))
            collect_reads(cpp, r['code'], {r['gparams'][k]: prefix}, r['file'], r['line'], reads, tcls, tsubst,
                          depth+1, seen)


def reader_registry(cpp):
    out = {}
    for reg in cpp.registrations():
        cls = reg['cls']
        info = cpp.class_def(cls, reg['reg_file'])
        tsubst = {}
        if info and info['tparams']:
            lt = reg['cls_expr'].find('<')
            if lt >= 0:
                actual = split_args(reg['cls_expr'][lt+1:reg['cls_expr'].rfind('>')], angle=True)
                for p, a in zip(info['tparams'], actual): tsubst[p] = base_name(a)
        reads = Reads()
        ctor_lines, argnames = [], []
        for r in cpp.ctors_of(cls, reg['reg_file']):
            gp = [g for g in r['gparams'] if g]
            ctor_lines.append('%s:%d' % (r['file'], r['line']))
            if not argnames:
                argnames = [re.findall(r'\w+', p.split('=')[0])[-1].rstrip('_') for p in r['params']
                            if 'hid_t' not in p and re.findall(r'\w+', p)]
            if not gp: continue
            collect_reads(cpp, r['code'], {gp[0]: ''}, r['file'], r['line'], reads, cls, tsubst, 0, set())
        reads.finish()
        reg.update(ctor=ctor_lines, reads=reads, tsubst=tsubst, argnames=argnames)
        out[reg['prefix']] = reg
    return out


# ----------------------------------------------------------------------------------------------------------------
# Python writers
# ----------------------------------------------------------------------------------------------------------------

POT = '/input/potential'
FILE_OPENERS = ('File', 'open_file')
ARRAY_WRITERS = ('create_array', 'create_earray', 'create_carray')

def flatten_attr(node):
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr); node = node.value
    if isinstance(node, ast.Name): parts.append(node.id)
    else: return None
    return parts[::-1]


class Group:
    def __init__(self, path, created):
        self.path, self.created = path, created


class FuncScan(ast.NodeVisitor):
    """Writes made by one function: to groups it creates or opens, and to its own parameters."""
    def __init__(self, module, func, env, param_writes):
        self.module, self.func = module, func
        self.env = env                      # local name -> list of value nodes
        self.param_writes = param_writes    # (module-agnostic) func name -> {param index: [writes]}
        self.vars = {}                      # local name -> Group or ('param', i) or ('dset', Group, name)
        self.groups = {}                    # path -> dict(line, created, writes=[...])
        self.pwrites = defaultdict(list)    # param index -> [(kind, name, rank, line, cond)]
        self.cond = []                      # (test source, branch) of each enclosing if
        self.files = {'t', 'tf', 'h5', 'up', 'f'}   # HDF5 file handles by the names the writers use
        params = [a.arg for a in func.args.args] if func is not None else []
        for i, p in enumerate(params): self.vars[p] = ('param', i)

    # string values ---------------------------------------------------------------------------------------------
    def strvals(self, node, depth=0):
        if depth > 6: return {'*'}
        if isinstance(node, ast.Constant):
            if isinstance(node.value, str): return {node.value}
            if isinstance(node.value, bytes): return {node.value.decode()}
            return {'*'}
        if isinstance(node, ast.JoinedStr):
            res = {''}
            for v in node.values:
                part = self.strvals(v.value, depth+1) if isinstance(v, ast.FormattedValue) else self.strvals(v, depth+1)
                res = {a+b for a in res for b in part}
            return set(list(res)[:16])
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            l, r = self.strvals(node.left, depth+1), self.strvals(node.right, depth+1)
            return set(list({a+b for a in l for b in r})[:16])
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mod):
            fmts = self.strvals(node.left, depth+1)
            items = node.right.elts if isinstance(node.right, ast.Tuple) else [node.right]
            vals = [self.strvals(x, depth+1) for x in items]
            return self.subst(fmts, vals, r'%[-+ #0]*\d*(?:\.\d+)?[sdifr]')
        if isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Attribute) and f.attr == 'format':
                fmts = self.strvals(f.value, depth+1)
                return self.subst(fmts, [self.strvals(a, depth+1) for a in node.args], r'\{[^{}]*\}')
            if isinstance(f, ast.Attribute) and f.attr in ('encode', 'decode') : return self.strvals(f.value, depth+1)
            fname = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else '')
            if fname in ('bstring', 'str', 'bytes', 'bytes_', 'h5_as_text') and node.args:
                return self.strvals(node.args[0], depth+1)
            return {'*'}
        if isinstance(node, ast.IfExp):
            return self.strvals(node.body, depth+1) | self.strvals(node.orelse, depth+1)
        if isinstance(node, ast.Name) and node.id in self.env:
            res = set()
            for v in self.env[node.id]: res |= self.strvals(v, depth+1)
            return res or {'*'}
        return {'*'}

    @staticmethod
    def subst(fmts, vals, pat):
        out = set()
        for fmt in fmts:
            pieces = re.split(pat, fmt)
            res = {pieces[0]}
            for i, p in enumerate(pieces[1:]):
                v = vals[i] if i < len(vals) else {'*'}
                res = {a+b+p for a in res for b in v}
            out |= res
        return set(list(out)[:16])

    # ranks and argument counts --------------------------------------------------------------------------------
    def rank_of(self, node, depth=0):
        if node is None or depth > 6: return None
        if isinstance(node, ast.Call):
            f = node.func
            fname = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else '')
            if fname == 'arange': return 1
            if fname in ('zeros', 'ones', 'empty', 'full') and node.args:
                s = node.args[0]
                return len(s.elts) if isinstance(s, (ast.Tuple, ast.List)) else 1
            if fname in ('array', 'asarray', 'ascontiguousarray', 'require') and node.args:
                return self.rank_of(node.args[0], depth+1)
            if fname in ('astype', 'copy'):
                return self.rank_of(f.value, depth+1) if isinstance(f, ast.Attribute) else None
            if fname in ('concatenate',) and node.args and isinstance(node.args[0], (ast.List, ast.Tuple)) and node.args[0].elts:
                return self.rank_of(node.args[0].elts[0], depth+1)
            if fname in ('ravel', 'flatten'): return 1
            if fname == 'repeat' and node.args:
                return self.rank_of(node.args[0], depth+1)
            return None
        if isinstance(node, (ast.List, ast.Tuple)):
            if not node.elts: return 1
            inner = [self.rank_of(e, depth+1) if isinstance(e, (ast.List, ast.Tuple, ast.ListComp)) else 0
                     for e in node.elts]
            return None if any(r is None for r in inner) else 1 + max(inner)
        if isinstance(node, ast.ListComp):
            e = node.elt
            if isinstance(e, (ast.List, ast.Tuple)):
                r = self.rank_of(e, depth+1)
                return None if r is None else r+1
            if isinstance(e, (ast.Constant, ast.Compare, ast.IfExp, ast.BinOp, ast.Subscript)): return 1   # scalars
            return None
        if isinstance(node, ast.BinOp):
            rs = [self.rank_of(x, depth+1) for x in (node.left, node.right)]
            rs = [r for r in rs if r is not None and r > 0]
            return max(rs) if rs else None
        if isinstance(node, ast.Constant): return 0
        if isinstance(node, ast.Name) and node.id in self.env and len(self.env[node.id]) == 1:
            return self.rank_of(self.env[node.id][0], depth+1)
        return None

    def arglists(self, node, depth=0):
        # the argument lists an `arguments` value can hold, as tuples of node names
        if node is None or depth > 6: return None
        if isinstance(node, ast.Call) and node.args:
            f = node.func
            fname = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else '')
            return self.arglists(node.args[0], depth+1) if fname in ('array', 'asarray') else None
        if isinstance(node, (ast.List, ast.Tuple)):
            if any(isinstance(e, ast.Starred) for e in node.elts): return None
            return {tuple('|'.join(sorted(self.strvals(e))) for e in node.elts)}
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            l, r = self.arglists(node.left, depth+1), self.arglists(node.right, depth+1)
            return None if l is None or r is None else {a+b for a in l for b in r}
        if isinstance(node, ast.IfExp):
            l, r = self.arglists(node.body, depth+1), self.arglists(node.orelse, depth+1)
            return None if l is None or r is None else l | r
        if isinstance(node, ast.Name) and node.id in self.env:
            res = set()
            for v in self.env[node.id]:
                a = self.arglists(v, depth+1)
                if a is None: return None
                res |= a
            return res
        return None

    # group paths -------------------------------------------------------------------------------------------------
    def target(self, node):
        # Group / ('param', i) / ('dset', owner, name) denoted by an expression, or None
        if isinstance(node, ast.Name):
            v = self.vars.get(node.id)
            if v is not None: return v
            if node.id in ('potential', 'pot', 'potential_grp', 'pot_grp'): return Group(POT, False)
            if node.id in ('input', 'inp', 'input_grp'): return Group('/input', False)
            return None
        if isinstance(node, ast.Attribute):
            parts = flatten_attr(node)
            if parts and 'root' in parts:
                return Group('/' + '/'.join(parts[parts.index('root')+1:]), False)
            if parts and len(parts) >= 2:
                owner = self.target(node.value)
                if isinstance(owner, Group) and node.attr not in ('_v_attrs', 'attrs'):
                    return ('dset', owner, node.attr)
            return None
        if isinstance(node, ast.Subscript):
            owner = self.target(node.value)
            keys = self.strvals(node.slice) if not isinstance(node.slice, ast.Slice) else {'*'}
            if isinstance(owner, Group) and len(keys) == 1 and '*' not in next(iter(keys)):
                k = next(iter(keys))
                return Group(k if k.startswith('/') else owner.path.rstrip('/') + '/' + k, False)
            return None
        if isinstance(node, ast.Call):
            return self.created_group(node)
        return None

    def created_group(self, call):
        f = call.func
        fname = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else None)
        if fname in ('create_group', 'require_group') and isinstance(f, ast.Attribute):
            if len(call.args) >= 2:
                parent, name = self.target(call.args[0]), call.args[1]
            elif len(call.args) == 1:
                parent, name = self.target(f.value), call.args[0]
            else:
                return None
        elif fname in ('recreate_group', 'require_group') and len(call.args) >= 2:
            parent, name = self.target(call.args[0]), call.args[1]
        else:
            return None
        names = self.strvals(name)
        if not isinstance(parent, Group):
            if any(n.startswith('/') for n in names): parent = Group('', False)
            elif isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name) and len(call.args) == 1:
                parent = Group('', False) if f.value.id in self.files else None
            if parent is None: return None
        paths = sorted(n if n.startswith('/') else parent.path.rstrip('/') + '/' + n for n in names)
        g = Group(paths[0] if len(paths) == 1 else paths, True)
        g.line = call.lineno
        return g

    # recording ---------------------------------------------------------------------------------------------------
    def branch_paths(self, paths):
        # A group named from `a if T else b` (anywhere in the name's expression) and written under `if T:` (or
        # its else) is only the a (or b) group.
        for test, branch in self.cond:
            for vals in self.env.values():
                for v in (n for val in vals for n in ast.walk(val)):
                    if not (isinstance(v, ast.IfExp) and ast.unparse(v.test) == test): continue
                    keep = {k for k in self.strvals(v.body if branch else v.orelse) if k and k != '*'}
                    drop = {d for d in self.strvals(v.orelse if branch else v.body) if d and d != '*'} - keep
                    sel = [p for p in paths if (not keep or any(k in p for k in keep)) and
                           not any(d in p for d in drop)]
                    if sel: paths = sel
        return paths

    def record(self, owner, kind, name, rank, line, via=None):
        if isinstance(owner, tuple) and owner[0] == 'param':
            self.pwrites[owner[1]].append((kind, name, rank, line, bool(self.cond), via))
            return
        if isinstance(owner, tuple) and owner[0] == 'dset':
            owner, name = owner[1], owner[2] + '@' + name.lstrip('@')
            kind = 'attr'
        if not isinstance(owner, Group): return
        paths = self.branch_paths(owner.path) if isinstance(owner.path, list) else [owner.path]
        for path in paths:
            g = self.groups.setdefault(path, dict(line=getattr(owner, 'line', line), created=owner.created,
                                                  writes=[], arglists=None))
            g['writes'].append((kind, name, rank, line, bool(self.cond), via))

    def touch(self, g):
        for path in (g.path if isinstance(g.path, list) else [g.path]):
            self.groups.setdefault(path, dict(line=g.line, created=True, writes=[], arglists=None))

    def visit_With(self, node):
        for item in node.items:
            c = item.context_expr
            if isinstance(c, ast.Call) and isinstance(item.optional_vars, ast.Name):
                f = c.func
                fname = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else '')
                if fname in FILE_OPENERS:
                    self.vars[item.optional_vars.id] = Group('', False)
                    self.files.add(item.optional_vars.id)
        for s in node.body: self.visit(s)

    def visit_If(self, node):
        self.visit(node.test)
        test = ast.unparse(node.test)
        for branch, body in ((True, node.body), (False, node.orelse)):
            self.cond.append((test, branch))
            for s in body: self.visit(s)
            self.cond.pop()

    def visit_Assign(self, node):
        self.visit(node.value)
        value = node.value
        for tgt in node.targets:
            elts = tgt.elts if isinstance(tgt, ast.Tuple) else [tgt]
            for e in elts: self.assign_target(e, value if len(elts) == 1 else None, node)

    def assign_target(self, tgt, value, node):
        if isinstance(tgt, ast.Name) and value is not None:
            if isinstance(value, ast.Call):
                f = value.func
                fname = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else '')
                if fname in FILE_OPENERS:
                    self.vars[tgt.id] = Group('', False); self.files.add(tgt.id); return
                if fname == 'create_dataset' and isinstance(f, ast.Attribute):
                    owner = self.target(f.value)
                    names = self.strvals(value.args[0]) if value.args else {'*'}
                    if isinstance(owner, Group) and len(names) == 1:
                        self.vars[tgt.id] = ('dset', owner, next(iter(names))); return
                if fname in ARRAY_WRITERS:
                    owner, names = self.array_call(value)
                    if owner is not None and len(names) == 1:
                        self.vars[tgt.id] = ('dset', owner, next(iter(names))); return
            g = self.target(value)
            if isinstance(g, (Group, tuple)):
                self.vars[tgt.id] = g
                if isinstance(g, Group) and g.created: self.touch(g)
            elif tgt.id in self.vars:
                del self.vars[tgt.id]
            return
        line = node.lineno
        # G._v_attrs.x = v ; G.attrs['x'] = v ; G._v_attrs['x'] = v ; G['name'] = v
        if isinstance(tgt, ast.Attribute) and isinstance(tgt.value, ast.Attribute) and tgt.value.attr == '_v_attrs':
            owner = self.target(tgt.value.value)
            self.record(owner, 'attr', '@' + tgt.attr, None, line)
            if tgt.attr == 'arguments' and isinstance(owner, Group) and value is not None:
                self.set_arglists(owner, value, line)
            return
        if isinstance(tgt, ast.Subscript):
            if isinstance(tgt.value, ast.Attribute) and tgt.value.attr in ('attrs', '_v_attrs'):
                owner = self.target(tgt.value.value)
                for n in self.strvals(tgt.slice):
                    self.record(owner, 'attr', '@' + n, None, line)
                    if n == 'arguments' and isinstance(owner, Group) and value is not None:
                        self.set_arglists(owner, value, line)
                return
            owner = self.target(tgt.value)
            if isinstance(owner, (Group, tuple)) and not (isinstance(owner, tuple) and owner[0] == 'dset'):
                for n in self.strvals(tgt.slice):
                    if n.startswith('/'):
                        parent, _, leaf = n.rpartition('/')
                        self.record(Group(parent, False), 'dset', leaf, self.rank_of(value), line)
                    else:
                        self.record(owner, 'dset', n, self.rank_of(value), line)

    def set_arglists(self, owner, value, line):
        # `arguments` as a set of argument-name tuples per group.  When the group's name and its arguments both
        # come from `a if T else b` on the same test, each branch's name gets that branch's arguments.
        paths = owner.path if isinstance(owner.path, list) else [owner.path]
        per_path = {p: self.arglists(value) for p in paths}
        tests = {ast.unparse(n.test) for n in ast.walk(value) if isinstance(n, ast.IfExp)}
        if len(paths) > 1:
            for test in tests:
                for branch in (True, False):
                    class Pick(ast.NodeTransformer):
                        def visit_IfExp(s, n):
                            s.generic_visit(n)
                            return (n.body if branch else n.orelse) if ast.unparse(n.test) == test else n
                    self.cond.append((test, branch))
                    sel = self.branch_paths(paths)
                    self.cond.pop()
                    if len(sel) < len(paths):
                        lists = self.arglists(Pick().visit(copy.deepcopy(value)))
                        for p in sel: per_path[p] = lists
        for path, lists in per_path.items():
            g = self.groups.setdefault(path, dict(line=line, created=owner.created, writes=[], arglists=None))
            g['arglists'] = lists if lists is not None else 'unknown'

    def array_call(self, call):
        # (owner, names) of create_array(grp, nm, obj) helpers and of t.create_array(where, name, obj)
        groups, params = [], []
        for i, a in enumerate(call.args[:3]):
            if isinstance(a, ast.Name) and a.id in self.files: continue
            t = self.target(a)
            if isinstance(t, Group): groups.append(i)
            elif isinstance(t, tuple) and t[0] == 'param': params.append(i)
        cands = groups or params
        if not cands: return None, set()
        i = cands[0]
        if not groups and len(params) > 1 and len(call.args) > params[1] + 1:
            i = params[1]    # create_array(h5, grp, nm, obj): the file comes first
        if len(call.args) <= i+1: return None, set()
        return self.target(call.args[i]), self.strvals(call.args[i+1])

    def array_obj(self, call, i):
        for kw in call.keywords:
            if kw.arg in ('obj', 'data'): return kw.value
        return call.args[i+2] if len(call.args) > i+2 else None

    def visit_Call(self, node):
        self.generic_visit(node)
        f = node.func
        fname = f.attr if isinstance(f, ast.Attribute) else (f.id if isinstance(f, ast.Name) else '')
        line = node.lineno
        g = self.created_group(node)
        if g is not None:
            self.touch(g)
            return
        if fname in ARRAY_WRITERS:
            owner, names = self.array_call(node)
            if owner is not None:
                i = [k for k, a in enumerate(node.args[:3]) if self.target(a) is owner or
                     (isinstance(a, ast.Name) and self.vars.get(a.id) is owner)]
                obj = self.array_obj(node, i[0] if i else 0)
                for n in names: self.record(owner, 'dset', n, self.rank_of(obj), line)
            return
        if fname == 'create_dataset' and isinstance(f, ast.Attribute):
            owner = self.target(f.value)
            names = self.strvals(node.args[0]) if node.args else self.strvals(
                next((k.value for k in node.keywords if k.arg == 'name'), ast.Constant('*')))
            obj = next((k.value for k in node.keywords if k.arg == 'data'), node.args[1] if len(node.args) > 1 else None)
            for n in names:
                if n.startswith('/'):
                    parent, _, leaf = n.rpartition('/')
                    self.record(Group(parent, False), 'dset', leaf, self.rank_of(obj), line)
                else:
                    self.record(owner, 'dset', n, self.rank_of(obj), line)
            return
        if fname == 'create' and isinstance(f, ast.Attribute) and isinstance(f.value, ast.Attribute) and \
                f.value.attr == 'attrs' and node.args:
            owner = self.target(f.value.value)
            for n in self.strvals(node.args[0]): self.record(owner, 'attr', '@' + n, None, line)
            return
        if fname == 'setattr' and len(node.args) >= 2 and isinstance(node.args[0], ast.Attribute) and \
                node.args[0].attr == '_v_attrs':
            owner = self.target(node.args[0].value)
            for n in self.strvals(node.args[1]): self.record(owner, 'attr', '@' + n, None, line)
            return
        # a group handed to another function: that function's writes to the parameter count here
        callee = self.param_writes.get(fname)
        if callee:
            for i, a in enumerate(node.args):
                t = self.target(a)
                if not isinstance(t, (Group, tuple)) or (isinstance(t, tuple) and t[0] == 'dset'): continue
                for (kind, name, rank, l, cond, via) in callee.get(i, []):
                    self.record(t, kind, name, rank, line, via=via or '%s:%d' % (fname, l))

    def visit_FunctionDef(self, node):
        for s in node.body: self.visit(s)   # a nested helper shares the enclosing function's scope


def local_env(func):
    env = defaultdict(list)
    for n in ast.walk(func):
        if isinstance(n, ast.Assign):
            for t in n.targets:
                if isinstance(t, ast.Name): env[t.id].append(n.value)
    return env


def scan_writers():
    modules = {}
    for path in sorted(glob.glob(os.path.join(PY_DIR, '*.py'))):
        try: modules[os.path.relpath(path, UPSIDE_HOME)] = ast.parse(open(path).read())
        except SyntaxError: continue
    funcs = []
    for mod, tree in modules.items():
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)): funcs.append((mod, node))
    param_writes = {}
    # each round carries a helper's writes to its parameters one call further up
    for rnd in range(3):
        results = []
        for mod, fn in funcs:
            s = FuncScan(mod, fn, local_env(fn), param_writes)
            s.visit(fn)
            results.append((mod, fn, s))
        param_writes = {}
        for mod, fn, s in results:
            if s.pwrites: param_writes.setdefault(fn.name, {}).update(s.pwrites)
    writers = []
    for mod, fn, s in results:
        for path, g in s.groups.items():
            writers.append(dict(file=mod, func=fn.name, func_line=fn.lineno, path=path, line=g['line'],
                                created=g['created'], writes=g['writes'], arglists=g['arglists']))
    return writers


def node_of(path):
    # (node name, sub-path inside the node) of a path under /input/potential
    if not path.startswith(POT + '/'): return None, None
    rest = path[len(POT)+1:]
    node, _, sub = rest.partition('/')
    return node, sub


def prefix_for(name, prefixes):
    # registered prefixes the engine could pick for a node name (registrations are prefix-free, so a concrete
    # name has at most one); a '*' left in a writer's pattern is an unresolved part of the name
    head = name.split('*')[0]
    exact = [p for p in prefixes if head.startswith(p)]
    if exact or '*' not in name or not head: return exact
    return [p for p in prefixes if p.startswith(head)]


# ----------------------------------------------------------------------------------------------------------------
# Dynamic pass
# ----------------------------------------------------------------------------------------------------------------

def file_contents(path):
    nodes = {}
    with h5py.File(path, 'r', locking=False) as h:
        pot = h['/input/potential']
        for name, g in pot.items():
            items = {'@' + k: 'attr' for k in g.attrs if k not in PYTABLES_ATTRS}
            def visit(sub, obj):
                items[sub + ('' if isinstance(obj, h5py.Dataset) else '/')] = (
                    '%s %s' % (tuple(obj.shape), obj.dtype) if isinstance(obj, h5py.Dataset) else 'group')
                for k in obj.attrs:
                    if k not in PYTABLES_ATTRS: items[sub + '@' + k] = 'attr'
            g.visititems(visit)
            args = [a.decode() if isinstance(a, bytes) else str(a) for a in g.attrs.get('arguments', [])]
            nodes[name] = dict(items=items, nargs=len(args))
        n_atom = h['/input/pos'].shape[0]
    return nodes, n_atom


def trace_engine(path, n_atom, out):
    code = ('import ctypes as ct\n'
            'lib = ct.cdll.LoadLibrary(%r)\n'
            'lib.construct_deriv_engine.restype = ct.c_void_p\n'
            'lib.construct_deriv_engine.argtypes = [ct.c_int, ct.c_char_p, ct.c_bool]\n'
            'e = lib.construct_deriv_engine(%d, %r, True)\n'
            'print("ENGINE_OK" if e else "ENGINE_FAIL")\n') % (ENGINE_LIB, n_atom, path.encode())
    env = dict(os.environ, DYLD_INSERT_LIBRARIES=TRACE_LIB, H5TRACE_OUT=out, HDF5_USE_FILE_LOCKING='FALSE')
    r = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True)
    ok = 'ENGINE_OK' in r.stdout
    reads = defaultdict(lambda: defaultdict(set))   # node -> name -> set(op)
    external = defaultdict(set)
    for ln in open(out):
        op, loc, name, attr = (ln.rstrip('\n').split('\t') + ['', '', '', ''])[:4]
        if loc == POT:
            if attr: reads[name]['@' + attr].add(op)
            continue
        if not loc.startswith(POT + '/'): continue
        node, sub = node_of(loc)
        sub = sub + '/' if sub else ''
        if attr:
            owner = '' if name in ('.', '') else name
            if owner.startswith('/'): external[node].add(owner + '@' + attr); continue
            reads[node][(sub + owner).rstrip('/') + '@' + attr].add(op)
        else:
            if name.startswith('/'): external[node].add(name); continue
            reads[node][sub + name].add(op)
    return ok, reads, external, r.stderr.strip()[-400:]


def dynamic_pass(paths):
    os.makedirs(WORK_DIR, exist_ok=True)
    out = {}
    for path in paths:
        if not os.path.exists(path): continue
        nodes, n_atom = file_contents(path)
        tr = os.path.join(WORK_DIR, os.path.basename(os.path.dirname(path)) + '__' + os.path.basename(path) + '.trace')
        ok, reads, external, err = trace_engine(path, n_atom, tr)
        for node, r in reads.items():     # a group the engine opened or probed is named as the file names it
            items = nodes.get(node, {}).get('items', {})
            for k in [k for k in r if k + '/' in items]: r[k + '/'] |= r.pop(k)
        out[path] = dict(nodes=nodes, ok=ok, reads=reads, external=external, err=err, trace=tr)
    return out


# ----------------------------------------------------------------------------------------------------------------
# Comparison and report
# ----------------------------------------------------------------------------------------------------------------

def names_of(key):
    return key.split('|')

def written(name, wset):
    return any(n in wset for n in names_of(name))

def fmt_reads(reads):
    items = []
    for name, d in sorted(reads.dsets.items()):
        r = '/'.join(str(x) for x in sorted(d['ranks'])) if d['ranks'] else ''
        items.append('%s%s%s' % (name, '[%s]' % r if r else '', '' if d['required'] else '?'))
    for key, a in sorted(reads.attrs.items()):
        items.append('%s%s' % (key, '' if a['required'] else '?'))
    return items

def fmt_writes(ws):
    seen = {}
    for kind, name, rank, line, cond, via in ws:
        key = name
        s = seen.setdefault(key, [None, False])
        if rank is not None: s[0] = rank
        s[1] |= cond
    out = []
    for name, (rank, cond) in sorted(seen.items()):
        out.append('%s%s%s' % (name, '[%d]' % rank if rank is not None else '', '*' if cond else ''))
    return out


STATUS_ORDER = ['argument-count disagreement', 'rank/shape disagreement', 'read-but-never-written',
                'written-but-never-read']


def analyse(readers, writers, dyn):
    prefixes = sorted(readers)
    rows = {p: dict(reader=readers[p], writers=[], dyn=[], flags=[]) for p in prefixes}
    orphans = []
    for w in writers:
        if not isinstance(w['path'], str): continue
        node, sub = node_of(w['path'])
        if not node: continue
        cands = prefix_for(node, prefixes)
        if not cands: orphans.append(w)
        for p in cands: rows[p]['writers'].append(dict(w, node=node, sub=sub))

    for path, d in dyn.items():
        for node, info in d['nodes'].items():
            cands = prefix_for(node, prefixes)
            if len(cands) == 1:
                rows[cands[0]]['dyn'].append(dict(file=path, node=node, items=info['items'], nargs=info['nargs'],
                                                  reads=d['reads'].get(node, {}), ok=d['ok'],
                                                  external=d['external'].get(node, set())))

    for p, row in rows.items():
        rd = row['reader']['reads']
        nargs = row['reader']['nargs']
        args_ok = lambda a: a >= 1 if nargs == -1 else a == nargs
        req = {n for n, d in rd.dsets.items() if d['required'] and not n.startswith('/')}
        req |= {k for k, a in rd.attrs.items() if a['required'] and not k.startswith('/')}
        read_names = {n for k in set(rd.dsets) | set(rd.attrs) | set(ENGINE_ATTRS) for n in names_of(k)}

        # one entry per writing function; a subgroup's writes are folded into its node with the sub-path
        merged = {}
        for w in row['writers']:
            sub = (w['sub'] + '/') if w['sub'] else ''
            m = merged.setdefault((w['file'], w['func']), dict(file=w['file'], func=w['func'], nodes=set(),
                                                               writes=[], args={}, lines=set(), creates=False))
            m['nodes'].add(w['node'])
            m['lines'].add(w['line'])
            for (kind, name, rank, line, cond, via) in w['writes']:
                full = (sub.rstrip('/') + name) if name.startswith('@') else sub + name
                m['writes'].append((kind, full, rank, line, cond, via))
            if sub: m['writes'].append(('group', sub, None, w['line'], False, None))
            if w['arglists'] is not None: m['args'][w['node']] = w['arglists']
            if w['created'] and not sub: m['creates'] = True
        row['merged'] = sorted(merged.values(), key=lambda m: (m['file'], min(m['lines'])))

        flags, notes = [], []
        for dd in row['dyn']:
            traced = set(dd['reads'])
            dd['unread'] = sorted(k for k, v in dd['items'].items() if v != 'group' and k not in traced)
            dd['absent'] = sorted(k for k, ops in dd['reads'].items()
                                  if ops & {'attr_no', 'exists_no'} and k != '@integrator_level')
            dd['static_missed'] = sorted(k for k in traced if k not in read_names and not k.endswith('/'))
            if dd['nargs'] and not args_ok(dd['nargs']):
                flags.append(('argument-count disagreement', '%s in %s has %d arguments, reader takes %s'
                              % (dd['node'], os.path.basename(dd['file']), dd['nargs'], '1+' if nargs == -1 else nargs)))
        # a static read the engine never asked for while building a traced file is conditional
        conditional = {n for n in req for dd in row['dyn']
                       if dd['ok'] and not any(x in dd['reads'] or x + '/' in dd['reads'] for x in names_of(n))}

        def two_of_one(n, wset):
            # InteractionGraph reads `X2` only for an asymmetric graph, one written as `X1` and `X2`
            alts = [k for k in rd.dsets if n[:-1] in names_of(k) and n[:-1] + '1' in names_of(k)]
            return n.endswith('2') and alts and n[:-1] in wset and n[:-1] + '1' not in wset

        all_written = {name for m in row['merged'] for (k, name, r, l, c, v) in m['writes']}
        extenders = [m for m in row['merged'] if not m['creates']]
        for m in row['merged']:
            if not m['creates']: continue
            # a function that only adds to a group made elsewhere counts with its creator
            pool = [m] + [e for e in extenders if e['nodes'] & m['nodes']]
            wset = {name for x in pool for (k, name, r, l, c, v) in x['writes']}
            where = '%s:%d %s' % (m['file'], min(m['lines']), m['func'])
            missing = sorted(n for n in req if not written(n, wset) and not two_of_one(n, wset))
            hard = [n for n in missing if n not in conditional]
            if hard: flags.append(('read-but-never-written', '%s lacks %s' % (where, ', '.join(hard))))
            optional = [k for k in set(rd.dsets) | set(rd.attrs) if k not in req]
            unread = []
            for n in sorted({name for (k, name, r, l, c, v) in m['writes'] if k != 'group'} - read_names):
                close = difflib.get_close_matches(n, [x for k in optional + sorted(req) for x in names_of(k)
                                                      if x not in wset], n=1, cutoff=0.75)
                unread.append(n + (' (close to the read name %s)' % close[0] if close else ''))
            if unread: flags.append(('written-but-never-read', '%s writes %s' % (where, ', '.join(unread))))
            for (k, name, r, l, c, v) in m['writes']:
                if k != 'dset' or r is None: continue
                for rn, d in rd.dsets.items():
                    if name in names_of(rn) and d['ranks'] and r not in d['ranks']:
                        flags.append(('rank/shape disagreement', '%s:%d writes %s with rank %d, reader expects %s'
                                      % (m['file'], l, name, r, '/'.join(map(str, sorted(d['ranks']))))))
            for node, lists in m['args'].items():
                if isinstance(lists, set) and not any(args_ok(len(a)) for a in lists):
                    flags.append(('argument-count disagreement', '%s gives `%s` %s arguments, reader takes %s'
                                  % (where, node, '/'.join(sorted({str(len(a)) for a in lists})),
                                     '1+' if nargs == -1 else nargs)))
        if row['merged']:
            unserved = sorted(n for n in conditional if not written(n, all_written) and not two_of_one(n, all_written))
            if unserved:
                flags.append(('read-but-never-written', 'conditional reads no writer provides (not reached in the '
                              'traced files): %s' % ', '.join(unserved)))
            defaulted = sorted(k for k in set(rd.dsets) | set(rd.attrs)
                               if k not in req and not k.startswith('/') and not written(k, all_written))
            if defaulted: notes.append('optional, never written: ' + ', '.join(defaulted))
        for name, d in rd.dsets.items():
            if len(d['ranks']) > 1:
                flags.append(('rank/shape disagreement', 'the reader itself reads %s with ranks %s (%s)'
                              % (name, '/'.join(map(str, sorted(d['ranks']))), ', '.join(sorted(d['lines'])))))
        row['notes'] = notes
        row['flags'] = flags
        kinds = {f[0] for f in flags}
        row['status'] = ('no writer' if not row['merged'] else
                         next((k for k in STATUS_ORDER if k in kinds), 'matched'))
    return rows, orphans


def md_escape(s):
    return s.replace('|', '\\|')


def report(rows, orphans, dyn):
    counts = defaultdict(int)
    for r in rows.values(): counts[r['status']] += 1
    L = ['# Upside 3 writer/reader registry, %s' % STAMP, '',
         'Generated by `tests/registry.py`: a static parse of the node constructors in `src/` and of the group '
         'writers in `py/`, and an HDF5 read trace of the engine library constructing each fixture. '
         '%d registered node types: %s; %d written groups match no registered prefix.'
         % (len(rows), ', '.join('%s %d' % (k, v) for k, v in sorted(counts.items())), len(orphans)), '',
         'Reads: `name[rank]` is a dataset, its rank from `traverse_dset`, `check_size` or `get_dset_size`; `@attr` an '
         'attribute of the node; `sub/` a subgroup; `dset@attr` an attribute of a dataset; `a|b` either name '
         '(template-dependent); `?` optional (default value or existence probe). Every node also has `@arguments` '
         '(required) and `@integrator_level` (default 2), read by `initialize_engine_from_hdf5`. Writes: `[rank]` '
         'where the writer makes it evident; `*` written only under a condition. Dynamic: for each traced file '
         'holding the node, `unread` is what the file holds and the engine never opened, `absent` what the engine '
         'probed for and did not find, `traced-not-in-static` a read the static pass missed.', '',
         'Traced: ' + ', '.join('`%s`%s' % (os.path.relpath(p, UPSIDE_HOME), '' if d['ok'] else ' (engine refused it)')
                               for p, d in dyn.items()) +
         '. The first five are the audit fixtures; the others are earlier single-purpose configs that carry nodes '
         'the fixtures lack.', '',
         '## Registry', '',
         '| node prefix | reader (ctor `class`; arg count (arg names)) | reads | writers | writes | dynamic | status |',
         '|---|---|---|---|---|---|---|']
    for p in sorted(rows, key=str.lower):
        row = rows[p]
        rg, rd = row['reader'], row['reader']['reads']
        reader = '%s `%s`; %s (%s)' % (', '.join(rg['ctor']) or '%s:%d' % (rg['reg_file'], rg['reg_line']),
                                       rg['cls'], '1+' if rg['nargs'] == -1 else rg['nargs'], ', '.join(rg['argnames']))
        reads = ', '.join(fmt_reads(rd)) or '(none)'
        if rd.unresolved: reads += '; unresolved: ' + '; '.join(rd.unresolved)
        wcells, wrcells = [], []
        for m in row['merged']:
            fmt = lambda al: (' or '.join('(%s)' % ', '.join(a) for a in sorted(al)) if isinstance(al, set) else '?')
            byargs = defaultdict(list)
            for node in sorted(m['nodes']): byargs[fmt(m['args'][node]) if node in m['args'] else ''].append(node)
            nodes = '; '.join('`%s`%s' % (', '.join(ns), ' args ' + a if a else '') for a, ns in byargs.items())
            wcells.append('%s:%d `%s` -> %s%s' % (m['file'], min(m['lines']), m['func'], nodes,
                                                   '' if m['creates'] else ' (extends)'))
            wrcells.append('%s: %s' % (m['func'], ', '.join(fmt_writes([w for w in m['writes'] if w[0] != 'group']))
                                       or '(none)'))
        dcells = defaultdict(list)    # identical results from several files share one entry
        for dd in row['dyn']:
            parts = []
            if dd['unread']: parts.append('unread ' + ', '.join(dd['unread']))
            if dd['absent']: parts.append('absent ' + ', '.join(dd['absent']))
            if dd['static_missed']: parts.append('traced-not-in-static ' + ', '.join(dd['static_missed']))
            if dd['external']: parts.append('outside the node ' + ', '.join(sorted(dd['external'])))
            dcells['`%s`: %s' % (dd['node'], '; '.join(parts) or 'all read')].append(os.path.basename(dd['file'])[:-3])
        dcells = ['%s %s' % (','.join(sorted(f)), k) for k, f in dcells.items()]
        status = '<br>'.join([row['status']] + [v if k == row['status'] else '%s: %s' % (k, v) for k, v in row['flags']]
                             + row['notes'])
        L.append('| `%s` | %s | %s | %s | %s | %s | %s |' % (
            p, md_escape(reader), md_escape(reads), md_escape('<br>'.join(wcells) or '-'),
            md_escape('<br>'.join(wrcells) or '-'), md_escape('<br>'.join(sorted(dcells)) or '-'), md_escape(status)))
    for w in orphans:
        L.append('| `%s` | (none) | - | %s:%d `%s` | %s | - | no reader |' % (
            node_of(w['path'])[0], w['file'], w['line'], w['func'], md_escape(', '.join(fmt_writes(w['writes'])))))
    return L


def main():
    cpp = CppSource()
    readers = reader_registry(cpp)
    writers = scan_writers()
    dyn = dynamic_pass(FIXTURES + SUPPLEMENTARY) if os.path.exists(TRACE_LIB) else {}
    rows, orphans = analyse(readers, writers, dyn)
    lines = report(rows, orphans, dyn)
    # the hand-written Findings section of an earlier run survives a regeneration of the table
    old = open(OUT_MD).read() if os.path.exists(OUT_MD) else ''
    if '## Findings' in old and '## Registry' in old:
        kept = old[old.index('## Findings'):old.index('## Registry')].rstrip('\n').split('\n') + ['']
        at = lines.index('## Registry')
        lines = lines[:at] + kept + lines[at:]
    with open(OUT_MD, 'w') as f: f.write('\n'.join(lines) + '\n')
    dump = {}
    for p, row in rows.items():
        rd = row['reader']['reads']
        dump[p] = dict(status=row['status'], flags=row['flags'], ctor=row['reader']['ctor'],
                       nargs=row['reader']['nargs'], reads=fmt_reads(rd), via=sorted(rd.via),
                       unresolved=rd.unresolved,
                       writers=[dict(file=m['file'], func=m['func'], nodes=sorted(m['nodes']), line=min(m['lines']),
                                     writes=fmt_writes(m['writes'])) for m in row['merged']],
                       dynamic=[dict(file=os.path.relpath(d['file'], UPSIDE_HOME), node=d['node'],
                                     unread=d['unread'], absent=d['absent'], static_missed=d['static_missed'],
                                     items=d['items']) for d in row['dyn']])
    with open(OUT_JSON, 'w') as f: json.dump(dict(rows=dump, orphans=[dict(path=w['path'], file=w['file'],
                                                  func=w['func'], line=w['line']) for w in orphans]), f, indent=1,
                                            default=str)
    print(OUT_MD)
    for p in sorted(rows): print('%-36s %s' % (p, rows[p]['status']))


if __name__ == '__main__':
    main()
