#!/usr/bin/env python
# List every /input dataset and attribute that differs between two configs, so a writer change can be shown to
# touch only what it means to.
# usage: config_diff.py CONFIG_A CONFIG_B
import sys
import numpy as np
import tables as tb

def contents(path):
    out = {}
    with tb.open_file(path) as t:
        for node in t.walk_nodes('/input'):
            if isinstance(node, tb.Leaf):
                out[node._v_pathname] = node[:]
            for k in node._v_attrs._v_attrnamesuser:
                out[node._v_pathname + '@' + k] = np.asarray(node._v_attrs[k])
    return out

def main():
    a, b = contents(sys.argv[1]), contents(sys.argv[2])
    n_diff = 0
    for k in sorted(set(a) | set(b)):
        if k not in a or k not in b:
            print('%-70s only in %s' % (k, 'A' if k in a else 'B')); n_diff += 1
        elif a[k].shape != b[k].shape:
            print('%-70s shape %s vs %s' % (k, a[k].shape, b[k].shape)); n_diff += 1
        elif a[k].tobytes() != b[k].tobytes():
            if a[k].dtype.kind in 'fiub':
                d = np.abs(a[k].astype('f8') - b[k].astype('f8'))
                print('%-70s max |diff| %.3e  (%d of %d differ)' % (k, d.max(), np.sum(d != 0), d.size))
            else:
                print('%-70s differs' % k)
            n_diff += 1
    print('%d of %d entries differ' % (n_diff, len(set(a) | set(b))))

if __name__ == '__main__':
    main()
