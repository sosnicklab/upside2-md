
def get_cb_points(param1, param2):

    assert(param1.shape == param2.shape)

    n_aa, n_bl, n_sp = param1.shape
    cb_param = np.zeros((n_aa, n_bl, 3))

    xx = np.arange(1, n_sp-1)
    for i in range(n_aa):
        for j in range(n_bl):
            p1 = param1[i,j]
            p2 = param2[i,j]

            y1 = ue.clamped_spline_value(p1, xx)
            y2 = ue.clamped_spline_value(p2, xx)

            cb_param[i,j,0] = y1[0]
            cb_param[i,j,1] = y1[-1]
            cb_param[i,j,2] = y2[-1]
            assert(y1[-1] == y2[0])

    return cb_param

def get_hb_points(param1, param2):

    assert(param1.shape == param2.shape)

    n_aa, n_sp = param1.shape
    hb_param = np.zeros((n_aa, 3))

    xx = np.arange(1, n_sp-1)
    for i in range(n_aa):
        p1 = param1[i]
        p2 = param2[i]
        y1 = ue.clamped_spline_value(p1, xx)
        y2 = ue.clamped_spline_value(p2, xx)

        hb_param[i,0] = y1[0]
        hb_param[i,1] = y1[-1]
        hb_param[i,2] = y2[-1]
        assert(y1[-1] == y2[0])

    return hb_param

def set_cb_points(param1, param2, cb_param):
    assert(param1.shape == param2.shape)
    n_aa, n_bl, n_sp = param1.shape

    old_cb_param = get_cb_points(param1, param2)
    dy = cb_param - old_cb_param

    n_param1 = param1 * 1.0
    n_param2 = param2 * 1.0
    for i in range(n_aa):
        for j in range(n_bl):
            p1 = param1[i,j]*1.
            p2 = param2[i,j]*1.
            np1 = p1*1.
            np2 = p2*1.

            if dy[i,j, 0] != 0.0:
                np1 = shift_left(p1, dy[i,j,0])
            if dy[i,j, 2] != 0.0:
                np2 = shift_right(p2, dy[i,j,2])
            if dy[i,j, 1] != 0.0:
                np1,np2 = shift_middle(np1, np2, dy[i,j,1])

            n_param1[i,j] = np1
            n_param2[i,j] = np2
    return n_param1, n_param2

def set_hb_points(param1, param2, hb_param):
    assert(param1.shape == param2.shape)
    n_aa, n_sp = param1.shape

    old_hb_param = get_hb_points(param1, param2)
    dy = hb_param - old_hb_param

    n_param1 = param1 * 1.0
    n_param2 = param2 * 1.0
    for i in range(n_aa):
        p1 = param1[i]*1.
        p2 = param2[i]*1.
        np1 = p1*1.
        np2 = p2*1.
        if dy[i, 0] != 0.0:
            np1 = shift_left(p1, dy[i,0])
        if dy[i, 2] != 0.0:
            np2 = shift_right(p2, dy[i,2])
        if dy[i, 1] != 0.0:
            np1,np2 = shift_middle(np1, np2, dy[i,1])
        n_param1[i] = np1
        n_param2[i] = np2
    return n_param1, n_param2

def shift_left(p1, dy):
    n_node = p1.size
    xx     = np.arange(1, n_node-1)
    y      = ue.clamped_spline_value(p1, xx)
    s      = y[0] - y[-1]
    ns     = s + dy
    ny     = y-y[-1]
    ny     = ns/s*ny+y[-1]
    np     = ue.clamped_spline_solve(ny)
    return np

def shift_right(p2, dy):
    n_node = p2.size
    xx     = np.arange(1, n_node-1)
    y      = ue.clamped_spline_value(p2, xx)
    s      = y[-1] - y[0]
    ns     = s + dy
    ny     = y-y[0]
    ny     = ns/s*ny+y[0]
    np     = ue.clamped_spline_solve(ny)
    return np

def shift_middle(p1, p2, dy):
    n_node = p1.size
    xx     = np.arange(1, n_node-1)
    assert(n_node == p2.size)

    y1 = ue.clamped_spline_value(p1, xx)
    y2 = ue.clamped_spline_value(p2, xx)
    assert(y1[-1] == y2[0])

    np1 = shift_right(p1, dy)
    np2 = shift_left(p2, dy)
    return np1, np2

def compute_divergence(config_base, pos):

    restype        = t.root.input.potential.cb_membrane_potential.res_type[:]
    uniq_restype   = np.unique(restype)
    n_aa,n_bl,n_sp = t.root.input.potential.cb_membrane_potential.coeff_left.shape
    cbz_param_size = (n_aa, n_bl, 3)
    cbz_size       = n_aa*n_bl*n_sp

    n_hb,n_sp2     = t.root.input.potential.hb_membrane_potential.coeff_left.shape
    hbz_param_size = (n_hb, 3)
    cbz_size       = n_hb*n_sp2

    engine    = ue.Upside(config_base)
    contrast  = Update([],[])

    cb_param = engine.get_param_deriv((cbz_size*2,), 'cb_membrane_potential')
    hb_param = engine.get_param_deriv((hbz_size*2,), 'hb_membrane_potential')

    dp = 0.1
    dp_scale = dp*2

    for i in range(pos.shape[0]):

        cb_dparam = np.zeros(cbz_param_size)
        hb_dparam = np.zeros(hbz_param_size)

        for r in uniq_restype:
            for l in range(n_bl):

                st1 = (r*n_bl + l)*n_sp
                st2 = st1 + cbz_size
                p1 = cb_param[st1:st1+n_sp]*1.
                p2 = cb_param[st2:st2+n_sp]*1.

                new_cb_param = cb_param*1.0
                np1 = shift_left(p1, dy)
                new_cb_param[st1:st1+n_sp] = np1*1.
                engine.set_param(new_cb_param, 'cb_membrane_potential')
                cb_param[r,l,0] = engine.get_output('cb_membrane_potential')[0,0]*dp_scale
                new_cb_param = cb_param*1.0
                np1 = shift_left(p1, -dy)
                new_cb_param[st1:st1+n_sp] = np1*1.
                engine.set_param(new_cb_param, 'cb_membrane_potential')
                cb_param[r,l,0] -= engine.get_output('cb_membrane_potential')[0,0]*dp_scale

                new_cb_param = cb_param*1.0
                np2 = shift_right(p2, dy)
                new_cb_param[st2:st2+n_sp] = np2*1.
                engine.set_param(new_cb_param, 'cb_membrane_potential')
                cb_param[r,l,2] = engine.get_output('cb_membrane_potential')[0,0]*dp_scale
                new_cb_param = cb_param*1.0
                np2 = shift_right(p2, -dy)
                new_cb_param[st2:st2+n_sp] = np2*1.
                engine.set_param(new_cb_param, 'cb_membrane_potential')
                cb_param[r,l,2] -= engine.get_output('cb_membrane_potential')[0,0]*dp_scale

                new_cb_param = cb_param*1.0
                np1,np2 = shift_middle(p1, p2, dy)
                new_cb_param[st1:st1+n_sp] = np1*1.
                new_cb_param[st2:st2+n_sp] = np2*1.
                engine.set_param(new_cb_param, 'cb_membrane_potential')
                cb_param[r,l,1] = engine.get_output('cb_membrane_potential')[0,0]*dp_scale
                new_cb_param = cb_param*1.0
                np1,np2 = shift_middle(p1, p2, -dy)
                new_cb_param[st1:st1+n_sp] = np1*1.
                new_cb_param[st2:st2+n_sp] = np2*1.
                engine.set_param(new_cb_param, 'cb_membrane_potential')
                cb_param[r,l,1] -= engine.get_output('cb_membrane_potential')[0,0]*dp_scale

            st1 = r*n_sp2
            st2 = st1 + hbz_size
            p1 = hb_param[st1:st1+n_sp2]*1.
            p2 = hb_param[st2:st2+n_sp2]*1.

            np1 = shift_left(p1, dy)
            new_hb_param = hb_param*1.0
            new_hb_param[st1:st1+n_sp2] = np1*1.
            engine.set_param(new_hb_param, 'hb_membrane_potential')
            hb_param[r,0] = engine.get_output('hb_membrane_potential')[0,0]*dp_scale
            np1 = shift_left(p1, -dy)
            new_hb_param = hb_param*1.0
            new_hb_param[st1:st1+n_sp2] = np1*1.
            engine.set_param(new_hb_param, 'hb_membrane_potential')
            hb_param[r,0] -= engine.get_output('hb_membrane_potential')[0,0]*dp_scale

            np2 = shift_right(p2, dy)
            new_hb_param = hb_param*1.0
            new_hb_param[st2:st2+n_sp2] = np2*1.
            engine.set_param(new_hb_param, 'hb_membrane_potential')
            hb_param[r,2] = engine.get_output('hb_membrane_potential')[0,0]*dp_scale
            np2 = shift_right(p2, -dy)
            new_hb_param = hb_param*1.0
            new_hb_param[st2:st2+n_sp2] = np2*1.
            engine.set_param(new_hb_param, 'hb_membrane_potential')
            hb_param[r,2] -= engine.get_output('hb_membrane_potential')[0,0]*dp_scale

            np1,np2 = shift_middle(p1, p2, dy)
            new_hb_param = hb_param*1.0
            new_hb_param[st1:st1+n_sp2] = np1*1.
            new_hb_param[st2:st2+n_sp2] = np2*1.
            engine.set_param(new_hb_param, 'hb_membrane_potential')
            hb_param[r,1] = engine.get_output('hb_membrane_potential')[0,0]*dp_scale
            np1,np2 = shift_middle(p1, p2, -dy)
            new_hb_param = hb_param*1.0
            new_hb_param[st1:st1+n_sp2] = np1*1.
            new_hb_param[st2:st2+n_sp2] = np2*1.
            engine.set_param(new_hb_param, 'hb_membrane_potential')
            hb_param[r,1] -= engine.get_output('hb_membrane_potential')[0,0]*dp_scale

        contrast.cbz.append(cb_param)
        contrast.hbz.append(hb_param)

    # convert to numpy arrays
    contrast = Update(*[np.array(x) for x in contrast])
    return contrast
