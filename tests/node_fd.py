#!/usr/bin/env python
# Finite-difference checks of a config's potential terms, one term at a time: the config's other potential
# leaves are removed, then the force and the parameter derivative (get_param_deriv) of every node feeding the
# term are compared with central differences of the term's energy.  BP is converged to tol 1e-5 in the copies,
# so its early stop does not enter the differences.
# usage: node_fd.py WORK_DIR CONFIG [CONFIG ...]
import sys, os, shutil
import numpy as np
import tables as tb

sys.path.insert(0, os.environ['UPSIDE_HOME'] + '/py')
import upside_engine as ue

rng = np.random.default_rng(1)

def potential_leaves(t):
    pot = t.root.input.potential
    used = set(a.decode() for g in pot for a in g._v_attrs.arguments)
    return [g._v_name for g in pot if g._v_name not in used]

def single_term_config(config, term, path):
    # returns the nodes the term depends on, with each node's attributes
    shutil.copyfile(config, path)
    with tb.open_file(path, 'a') as t:
        pot = t.root.input.potential
        for node in potential_leaves(t):
            if node != term: t.remove_node(pot, node, recursive=True)
        if 'rotamer' in pot:
            pot.rotamer._v_attrs.tol = 1e-5
            pot.rotamer._v_attrs.max_iter = 20000
        attrs = dict((g._v_name, dict((k, g._v_attrs[k]) for k in g._v_attrs._v_attrnamesuser)) for g in pot)
    upstream, todo = {}, [term]
    while todo:
        node = todo.pop()
        if node in upstream or node == 'pos': continue
        upstream[node] = attrs[node]
        todo.extend(a.decode() for a in attrs[node]['arguments'])
    return upstream

def tie(node, attrs, v):
    # A tied parameter set reports in every slot the derivative with respect to the shared parameter, so a
    # direction is made tied and the analytic derivative is divided by the tie multiplicity.
    if node == 'rotamer':
        # entry [b,a] is [a,b] with the beads' angular splines exchanged, and every slot reads 2x (findings 3.12)
        n_type = int(round(np.sqrt(v.size/54)))
        p = v.reshape(n_type, n_type, 54)
        q = p.transpose(1,0,2)[:,:,np.r_[15:30, 0:15, 30:54]]
        return (0.5*(p+q)).ravel(), 0.5
    if attrs.get('number_independent_weights') == 20:
        # environment coupling: the last 20x20 parameters are the weights (after 3 per type for the sigmoid, the
        # spline coefficients for the nonlinear coupling), one per environment residue type shared by every coupling type
        w = v[-400:].reshape(20, 20)
        scale = np.ones_like(v); scale[-400:] = 1./20
        return np.concatenate([v[:-400], np.tile(w.mean(axis=0), 20)]), scale
    return v, 1.

def stencil(f, v, h):
    # fourth-order central difference along v.  f returns the energy and the step it really took, so the
    # analytic side is applied to the same steps and float32 rounding of the perturbed values cancels.
    e, d = zip(*[f(s*h*v) for s in (1., -1., 2., -2.)])
    c = np.array([8., -8., -1., 1.]) / (12.*h)
    return np.dot(c, e), np.tensordot(c, np.array(d), 1)

def energy_noise(f, v, h):
    # the energy's own noise along v: for a smooth energy a fourth difference at spacing h/8 is noise alone,
    # with variance 70 sigma^2; five of them from nine points
    e = np.array([f(k*h/8.*v)[0] for k in np.arange(-4., 5.)])
    d4 = e[:-4] - 4.*e[1:-3] + 6.*e[2:-2] - 4.*e[3:-1] + e[4:]
    return np.sqrt(np.mean(d4**2) / 70.)

def compare(label, analytic, f, directions, e0, step):
    # A mismatch at step h/2 counts when it is above 1% and above the difference formula's own error: three
    # times the stencil's response to the measured energy noise (0.95 sigma/step, sigma at least one ulp) and
    # three times its change to step h and to step h/4.  A wrong derivative leaves the same mismatch at every
    # step; a stencil that straddles a kink (a cutoff edge, a clamped cosine) changes it from one step to the
    # next.  A kink much closer than h/4 is not resolved: a flagged direction is settled by the line integral of
    # the analytic derivative against the energy (line_integral.py).
    ulp = float(np.spacing(np.float32(abs(e0) + 1.)))
    rows = []
    for name, v in directions:
        h = step(float(np.sum(analytic*v)), ulp)
        f1, d1 = stencil(f, v, h)
        f2, d2 = stencil(f, v, h/2)
        f4, d4 = stencil(f, v, h/4)
        a1, a2, a4 = np.sum(analytic*d1), np.sum(analytic*d2), np.sum(analytic*d4)
        res = 0.95*max(energy_noise(f, v, h), ulp)/(h/2)   # never below the energy's float32 rounding
        noise = max(3.*res, 3.*abs((f1-a1) - (f2-a2)), 3.*abs((f4-a4) - (f2-a2)))
        rows.append((abs(a2-f2)/(abs(f2)+res), abs(a2-f2) > noise, name, a2, f2, noise))
    rows.sort(key=lambda r: -r[0])
    bad = [r for r in rows if r[0] > 1e-2 and r[1]]
    print('  %-56s n %4d  median err %.1e  bad %3d' % (label, len(rows), np.median([r[0] for r in rows]), len(bad)))
    for err, _, name, a, f2, noise in bad[:5]:
        print('      %-18s analytic % .6e  fd % .6e  err %.1e  fd noise %.1e' % (name, a, f2, err, noise))

def check_term(config, term, path, frames):
    nodes = single_term_config(config, term, path)
    engine = ue.Upside(path, quiet=True)
    for i, x in enumerate(frames):
        e0 = float(engine.energy(x))
        g = engine.deriv(x).astype('f8')
        def energy_at_pos(d):
            y = (x+d).astype('f4')
            return float(engine.energy(y)), y.astype('f8')-x
        n = x.size
        coords = range(n) if n <= 300 else rng.choice(n, 300, replace=False)
        dirs  = [('atom %d dim %d' % (k//3, k%3), np.eye(1, n, k).reshape(x.shape)) for k in coords]
        dirs += [('random %d' % j, v/np.sqrt(np.sum(v*v))) for j, v in enumerate(rng.normal(size=(20,)+x.shape))]
        compare('%s frame %d force' % (term, i), g, energy_at_pos, dirs, e0, lambda a, ulp: 1e-2)

        for node, attrs in sorted(nodes.items()):
            n_param = engine.get_param_size(node)
            if not n_param: continue
            p0 = engine.get_param((n_param,), node).astype('f8')
            engine.energy(x)
            try:
                an = engine.get_param_deriv((n_param,), node).astype('f8')
            except RuntimeError:
                print('  %-56s not implemented' % ('%s frame %d param %s [%d]' % (term, i, node, n_param)))
                continue
            scale = tie(node, attrs, an)[1]
            w = np.minimum(np.abs(p0), 1.) + 1e-2   # step scale; capped so coordinates move at most ~0.2 A
            def energy_at_param(dp):
                p = (p0+dp).astype('f4')
                engine.set_param(p, node)
                e = float(engine.energy(x))
                engine.set_param(p0.astype('f4'), node)
                return e, p.astype('f8')-p0
            slots = range(n_param) if n_param <= 200 else rng.choice(n_param, 200, replace=False)
            dirs  = [('param %d' % k, tie(node, attrs, w*np.eye(1, n_param, k)[0])[0]) for k in slots]
            dirs += [('random %d' % j, tie(node, attrs, w*v/np.sqrt(np.sum(v*v)))[0])
                     for j, v in enumerate(rng.normal(size=(20, n_param)))]
            compare('%s frame %d param %s [%d]' % (term, i, node, n_param), scale*an, energy_at_param, dirs, e0,
                    lambda a, ulp: float(np.clip(1e3*ulp/(abs(a)+1e-12), 2e-3, 0.2)))
    del engine

def main():
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    for config in sys.argv[2:]:
        print(config)
        with tb.open_file(config) as t:
            x0 = t.root.input.pos[:,:,0].astype('f4').astype('f8')
            terms = potential_leaves(t)
        frames = [x0, (x0 + rng.normal(0., 0.2, x0.shape)).astype('f4').astype('f8')]
        for term in terms:
            check_term(config, term, os.path.join(work, 'term_%s.up' % term), frames)

if __name__ == '__main__':
    main()
