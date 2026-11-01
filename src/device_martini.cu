// Device implementations of martini_hybrid.cpp's martini_hybrid_position and martini_potential.cpp's
// martini_potential, martini_sc_table_1body, dist_spring and angle_spring.  Their periodic box is fixed (the device
// engine refuses the barostat).  As on the CPU, a node's cache of pairs within its cutoff plus cache_buffer is
// rebuilt, per system, when one of its particles has moved more than half the buffer since the last rebuild; each
// cache is a slab per particle of all its candidate partners, so it has no capacity to exceed.
#include "device_nodes.h"
#include "martini_interaction.h"

using namespace std;

namespace {

// The lines of martini_potential.cpp's MartiniScTableOneBody, DistSpring and AngleSpring that compute their terms,
// copied: the CPU nodes keep their own (martini_interaction.h)

// martini_sc_table_1body's tables: B-spline coefficients, (n_grid+2) x (n_theta+2) row-major per layer, on a radial
// grid of step grid_step_ang from 0 and a theta grid of step theta_step from 0 to pi
struct ScTable {
    const float* spline_coeff;
    int n_grid, n_theta;
    float grid_step_ang, theta_step, cutoff_ang;

    // Value of a layer and its derivatives per Angstrom and per radian, at r < cutoff_ang and theta in [0, pi].
    UPSIDE_HD void evaluate(float& value, float& dVdr, float& dVdtheta, int layer, float dist, float theta) const {
        float x = dist / grid_step_ang;
        float y = theta / theta_step;
        int ix = n_grid - 2 < int(x) ? n_grid - 2 : int(x);
        int iy = n_theta - 2 < int(y) ? n_theta - 2 : int(y);
        float wx[4], dwx[4], wy[4], dwy[4];
        cubic_bspline_weights(wx, dwx, x - float(ix));
        cubic_bspline_weights(wy, dwy, y - float(iy));
        const float* c = spline_coeff + (size_t(layer)*(n_grid + 2) + ix)*(n_theta + 2) + iy;
        float v = 0.f, dx = 0.f, dy = 0.f;
        for(int a = 0; a < 4; ++a) {
            const float* row = c + a*(n_theta + 2);
            float rv = wy[0]*row[0] + wy[1]*row[1] + wy[2]*row[2] + wy[3]*row[3];
            float rd = dwy[0]*row[0] + dwy[1]*row[1] + dwy[2]*row[2] + dwy[3]*row[3];
            v  +=  wx[a]*rv;
            dx += dwx[a]*rv;
            dy +=  wx[a]*rd;
        }
        value = v;
        dVdr = dx / grid_step_ang;
        dVdtheta = dy / theta_step;
    }

    // std::max(lo, std::min(hi, x)), NaN going to hi as there
    UPSIDE_HD static float clamp_cos(float x) {
        float m = x < 1.f ? x : 1.f;
        return -1.f < m ? m : -1.f;
    }

    // The energy of layer `layer` for an environment atom at minimum-image displacement dr from a CB at whose side
    // chain vector is cbv; false, with no energy, at or past the cutoff
    UPSIDE_HD bool contact_value(float& value, int layer, const Vec<3>& dr, const Vec<3>& cbv) const {
        const float kMinDistance = 1.0e-8f;
        float d2 = mag2(dr);
        float dist = sqrtf(d2 < kMinDistance ? kMinDistance : d2);
        if(dist >= cutoff_ang) return false;

        float cos_theta = clamp_cos(dot(dr, cbv) / (dist * mag(cbv)));
        float dVdr, dVdtheta;
        evaluate(value, dVdr, dVdtheta, layer, dist, acosf(cos_theta));
        return true;
    }

    // The derivative of row_scale times that energy by the environment atom (grad_env; by the CB point its negative)
    // and by the side chain vector; false, with neither, at or past the cutoff and at the CB
    UPSIDE_HD bool contact_deriv(Vec<3>& grad_env, Vec<3>& vector_grad, float row_scale, int layer, const Vec<3>& dr,
            const Vec<3>& cbv) const {
        const float kMinDistance = 1.0e-8f;
        float d2 = mag2(dr);
        float dist = sqrtf(d2 < kMinDistance ? kMinDistance : d2);
        if(dist >= cutoff_ang || dist <= 1.0e-6f) return false;

        Vec<3> u = (1.f / dist) * dr;
        float inv_cbv = 1.f / mag(cbv);
        Vec<3> v = inv_cbv * cbv;
        float cos_theta = clamp_cos(dot(u, v));
        float value, dVdr, dVdtheta;
        evaluate(value, dVdr, dVdtheta, layer, dist, acosf(cos_theta));

        // theta is the angle between the directions u and v = cbv/|cbv|: dtheta/d(target) = -w/dist and
        // dtheta/d(cbv) = -w_v/|cbv|, with w and w_v the unit vectors perpendicular to u (toward v) and to
        // v (toward u). On the CB axis they are undefined and dV/dtheta is exactly zero (clamped spline),
        // so the angular term vanishes there.
        Vec<3> perp_u = v - cos_theta * u;
        Vec<3> perp_v = u - cos_theta * v;
        float norm_u = mag(perp_u);
        float norm_v = mag(perp_v);
        grad_env = row_scale * dVdr * u;
        vector_grad = make_zero<3>();
        if(norm_u > 0.f) grad_env += (-row_scale * dVdtheta / (dist * norm_u)) * perp_u;
        if(norm_v > 0.f) vector_grad = (-row_scale * dVdtheta * inv_cbv / norm_v) * perp_v;
        return true;
    }
};

// dist_spring's bond between atoms at x1 and x2: spring_constant times the canonical spline of delta^2/2 (a
// LayeredClampedSpline1D<1>'s arrays on n_knot points from delta_min, coord_scale knots per unit delta) at
// delta = r - equil_dist, r their minimum-image distance, and in deriv its derivative by x1 (by x2 its negative)
UPSIDE_HD inline float dist_spring_term(Vec<3>& deriv, const Vec<3>& x1, const Vec<3>& x2, float equil_dist,
        float spring_constant, const float* coefficients, const float* left_clamped_value,
        const float* right_clamped_value, int n_knot, float delta_min, float coord_scale,
        float box_x, float box_y, float box_z) {
    auto disp = simulation_box::minimum_image(x1 - x2, box_x, box_y, box_z);
    float dist = mag(disp);

    // Canonical spline in delta-r: value delta^2/2 and its derivative per unit delta
    float result[2];
    clamped_spline_value_and_deriv<1>(result, coefficients, left_clamped_value, right_clamped_value, n_knot, 0,
            (dist - equil_dist - delta_min) * coord_scale);
    deriv = (spring_constant * result[0] * coord_scale / dist) * disp;
    return spring_constant * result[1];
}

// angle_spring's angle at x2 between x1 and x3: spring_constant times the canonical spline of delta^2/2 (arrays as
// dist_spring_term's) at delta = cos(theta) - cos_theta0, and its derivatives by the three atoms
UPSIDE_HD inline float angle_spring_term(Vec<3>& grad_a, Vec<3>& grad_b, Vec<3>& grad_c, const Vec<3>& x1,
        const Vec<3>& x2, const Vec<3>& x3, float cos_theta0, float spring_constant, const float* coefficients,
        const float* left_clamped_value, const float* right_clamped_value, int n_knot, float delta_min,
        float coord_scale, float box_x, float box_y, float box_z) {
    auto disp1 = simulation_box::minimum_image(x1 - x2, box_x, box_y, box_z);
    auto disp2 = simulation_box::minimum_image(x3 - x2, box_x, box_y, box_z);

    float norm1_sq = mag2(disp1);
    float norm2_sq = mag2(disp2);
    float norm1 = sqrtf(norm1_sq);
    float norm2 = sqrtf(norm2_sq);
    float dot_p = dot(disp1, disp2);
    float dp = dot_p / (norm1 * norm2); // cos(theta); rounding can take it past +-1, inside the spline's range
    float delta_cos = dp - cos_theta0;

    // Canonical spline in delta_cos: value delta_cos^2/2 and its derivative per unit delta_cos
    float result[2];
    clamped_spline_value_and_deriv<1>(result, coefficients, left_clamped_value, right_clamped_value, n_knot, 0,
            (delta_cos - delta_min) * coord_scale);

    // Forces via chain rule: F = -dV/dx = -(dV/d(cos(theta))) * (d(cos(theta))/dx)
    float dV_dcos = spring_constant * result[0] * coord_scale;
    float inv_norm1 = 1.0f / norm1;
    float inv_norm2 = 1.0f / norm2;
    Vec<3> dcos_dr1 = (disp2 * inv_norm1 - disp1 * (dot_p * inv_norm1 / norm1_sq)) * inv_norm2;
    Vec<3> dcos_dr2 = (disp1 * inv_norm2 - disp2 * (dot_p * inv_norm2 / norm2_sq)) * inv_norm1;

    grad_a = dV_dcos * dcos_dr1;
    grad_c = dV_dcos * dcos_dr2;
    grad_b = -grad_a - grad_c;
    return spring_constant * result[1];
}


struct DeviceHybridPosition : public DeviceNode {
    DeviceCoord &out, &source, &hbond;
    DeviceArray<unsigned char> derived;
    DeviceArray<int> bb, o, o_elem, src;
    DeviceArray<float> weight, back_weight, rescale;

    DeviceHybridPosition(DeviceCoord& out_, DeviceCoord& source_, DeviceCoord& hbond_,
            const vector<unsigned char>& derived_, const vector<int>& bb_, const vector<int>& o_,
            const vector<int>& o_elem_, const vector<int>& src_, const vector<float>& weight_,
            const vector<float>& back_weight_, const vector<float>& rescale_):
        out(out_), source(source_), hbond(hbond_), derived(derived_), bb(bb_), o(o_), o_elem(o_elem_), src(src_),
        weight(weight_), back_weight(back_weight_), rescale(rescale_) {}

    void compute_value(ComputeMode) override {
        VecArray x = source.value(), oh = hbond.value();
        float* y = out.output.ptr;
        int w = out.row_width;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int na) {
            for(int d=0; d<3; ++d) y[na*w+d] = x(d,na);
        });
        const int *b = bb.ptr, *oa = o.ptr, *oe = o_elem.ptr, *s = src.ptr;
        const float *wt = weight.ptr, *rs = rescale.ptr;
        parallel_for(int(bb.n), DEVICE_LAMBDA (int k) {
            const bool have_o = oe[k] >= 0;
            auto ov = have_o ? load_vec<3>(oh, oe[k]) : make_vec3(0.f, 0.f, 0.f);
            auto v = make_vec3(0.f, 0.f, 0.f);
            for(int d = 0; d < 4; ++d) {
                float wd = wt[4*k+d];
                if(wd == 0.f) continue;
                v += wd * ((d == 3) ? ov : load_vec<3>(x, s[3*k+d]));
            }
            v *= rs[k];
            if(have_o) for(int d=0; d<3; ++d) y[oa[k]*w+d] = ov[d];
            for(int d=0; d<3; ++d) y[b[k]*w+d] = v[d];
        });
    }

    void propagate_deriv() override {
        auto s = out.sens.view(), ss = source.sens.view(), hs = hbond.sens.view();
        int w = out.row_width, ws = source.row_width, wh = hbond.row_width;
        const unsigned char* der = derived.ptr;
        parallel_for(out.n_elem, DEVICE_LAMBDA (int na) {
            if(!der[na]) for(int d=0; d<3; ++d) ss.add(na*ws+d, s.get(na*w+d));
        });
        const int *b = bb.ptr, *oa = o.ptr, *oe = o_elem.ptr, *sr = src.ptr;
        const float* bw = back_weight.ptr;
        parallel_for(int(bb.n), DEVICE_LAMBDA (int k) {
            const bool have_o = oe[k] >= 0;
            auto bb_s = make_vec3(s.get(b[k]*w+0), s.get(b[k]*w+1), s.get(b[k]*w+2));
            auto o_s  = make_vec3(s.get(oa[k]*w+0), s.get(oa[k]*w+1), s.get(oa[k]*w+2));
            // O carries its own sensitivity plus BB's O-weighted share
            auto o_total = have_o ? o_s : make_vec3(0.f, 0.f, 0.f);
            for(int d = 0; d < 4; ++d) {
                float wd = bw[4*k+d];
                if(wd == 0.f) continue;
                if(d == 3) o_total += wd*bb_s;
                else for(int c=0; c<3; ++c) ss.add(sr[3*k+d]*ws+c, wd*bb_s[c]);
            }
            if(have_o) for(int c=0; c<3; ++c) hs.add(oe[k]*wh+c, o_total[c]);
            else       for(int c=0; c<3; ++c) ss.add(oa[k]*ws+c, o_s[c]);
        });
    }
};


// Starts the per-system rebuild flags of a skin cache: every system flagged at the first evaluation and while the
// buffer is not positive (true), otherwise none, for flag_moved to add to (false)
bool flag_all(DeviceArray<int>& rebuild, bool& first, float cache_buffer) {
    bool all = first || !(cache_buffer > 0.f);
    first = false;
    if(all) rebuild.upload(vector<int>(rebuild.n, 1));
    else    rebuild.zero();
    return all;
}

// Flags each system with a particle (n of them, at x with stride w, particle i at element at(i)) more than half the
// buffer from its position at the last rebuild, in the minimum image
template <typename At>
void flag_moved(DeviceArray<int>& rebuild, float cache_buffer, int n, const float* x, int w, const At& at,
        const int* system, const float* cached, float box_x, float box_y, float box_z) {
    int* rb = rebuild.ptr;
    float max_cache_dist2 = sqr(0.5f * cache_buffer);
    parallel_for(n, DEVICE_LAMBDA (int i) {
        int e = at(i);
        float dx = x[e*w+0] - cached[3*i+0];
        float dy = x[e*w+1] - cached[3*i+1];
        float dz = x[e*w+2] - cached[3*i+2];
        simulation_box::minimum_image_scalar(dx, dy, dz, box_x, box_y, box_z);
        if(dx*dx + dy*dy + dz*dz > max_cache_dist2) rb[system[i]] = 1;
    });
}

// The positions of the particles of the systems that rebuilt
template <typename At>
void store_cached(const DeviceArray<int>& rebuild, int n, const float* x, int w, const At& at, const int* system,
        float* cached) {
    const int* rb = rebuild.ptr;
    parallel_for(n, DEVICE_LAMBDA (int i) {
        if(rb[system[i]]) for(int d=0; d<3; ++d) cached[3*i+d] = x[at(i)*w+d];
    });
}


// The pair table, shared by every system, holds each atom's partners of higher index (local to the system) and
// their splines; the cache of atom i of system s is the slab of the table's entries at s*n_entry + row_start[i]
struct DeviceMartiniPotential : public DeviceNode {
    AccumView potential;
    DeviceCoord& pos;
    int n_atom, n_entry, n_grid, spline_stride;
    float r_min, r_max, cutoff, cache_buffer, box_x, box_y, box_z;
    DeviceArray<int> row_start, neighbor, spline;
    DeviceArray<float> coefficients, left, right;
    DeviceArray<int> atom_system, rebuild, n_cached, cached;
    DeviceArray<float> cached_pos;
    bool first;

    DeviceMartiniPotential(AccumView potential_, DeviceCoord& pos_, const vector<int>& row_start_,
            const vector<int>& neighbor_, const vector<int>& spline_, const vector<LayeredClampedSpline1D<1>>& splines,
            float r_min_, float r_max_, float cutoff_, float cache_buffer_, float box_x_, float box_y_, float box_z_):
        potential(potential_), pos(pos_), n_atom(int(row_start_.size())-1), n_entry(neighbor_.size()),
        n_grid(splines.size() ? splines[0].nx : 0), spline_stride(splines.size() ? splines[0].coefficients.size() : 0),
        r_min(r_min_), r_max(r_max_), cutoff(cutoff_), cache_buffer(cache_buffer_),
        box_x(box_x_), box_y(box_y_), box_z(box_z_),
        row_start(row_start_), neighbor(neighbor_), spline(spline_),
        atom_system(pos_.system_of_element()), rebuild(pos_.offset.size()-1), n_cached(pos_.n_elem),
        cached(size_t(pos_.offset.size()-1)*neighbor_.size()), cached_pos(3*size_t(pos_.n_elem)), first(true)
    {
        vector<float> c, l, r;
        for(auto& s: splines) {
            c.insert(c.end(), s.coefficients.begin(), s.coefficients.end());
            l.push_back(s.left_clamped_value[0]);
            r.push_back(s.right_clamped_value[0]);
        }
        coefficients = DeviceArray<float>(c);
        left = DeviceArray<float>(l);
        right = DeviceArray<float>(r);
        n_cached.zero();
    }

    void compute_value(ComputeMode mode) override {
        if(!n_entry) return;
        const float* x = pos.output.ptr;
        int w = pos.row_width, na_sys = n_atom, ne = n_entry;
        const int *sys = atom_system.ptr, *rs = row_start.ptr, *nb = neighbor.ptr, *sp = spline.ptr;
        float bx = box_x, by = box_y, bz = box_z;
        if(!flag_all(rebuild, first, cache_buffer))
            flag_moved(rebuild, cache_buffer, pos.n_elem, x, w, DEVICE_LAMBDA (int i) {return i;}, sys,
                    cached_pos.ptr, bx, by, bz);

        const int* rb = rebuild.ptr;
        int *nc = n_cached.ptr, *ca = cached.ptr;
        float active_cutoff = cutoff + (cache_buffer > 0.f ? cache_buffer : 0.f);
        float active_cutoff2 = sqr(active_cutoff);
        lanes_for(pos.n_elem, DEVICE_LAMBDA (int g, Lanes l) {
            int s = sys[g];
            if(!rb[s]) return;
            int a = g - s*na_sys, count = 0;
            size_t base = size_t(s)*ne + rs[a];
            l.chunks(rs[a], rs[a+1], [&](int k, bool in) {
                bool keep = false;
                if(in) {
                    int h = s*na_sys + nb[k];
                    float dx = x[g*w+0]-x[h*w+0], dy = x[g*w+1]-x[h*w+1], dz = x[g*w+2]-x[h*w+2];
                    simulation_box::minimum_image_scalar(dx, dy, dz, bx, by, bz);
                    keep = dx*dx + dy*dy + dz*dz < active_cutoff2;
                }
                int slot = l.place(keep, count);
                if(keep) ca[base + slot] = k;
            });
            nc[g] = count;
        });
        store_cached(rebuild, pos.n_elem, x, w, DEVICE_LAMBDA (int i) {return i;}, sys, cached_pos.ptr);

        auto ps = pos.sens.view();
        auto pot = potential;
        bool want_pot = mode == PotentialAndDerivMode;
        const float *co = coefficients.ptr, *lv = left.ptr, *rv = right.ptr;
        float cut = cutoff, rmin = r_min, rmax = r_max;
        int ng = n_grid, stride = spline_stride;
        lanes_for(pos.n_elem, DEVICE_LAMBDA (int g, Lanes l) {
            int s = sys[g], a = g - s*na_sys;
            size_t base = size_t(s)*ne + rs[a];
            auto pa = make_vec3(x[g*w+0], x[g*w+1], x[g*w+2]);
            AccumLocal pot_g, sens_g[3];
            l.chunks(0, nc[g], [&](int m, bool in) {
                if(!in) return;
                int k = ca[base + m];
                int h = s*na_sys + nb[k], ip = sp[k];
                float e;
                Vec<3> f;
                if(!martini_pair_force(e, f, pa, make_vec3(x[h*w+0], x[h*w+1], x[h*w+2]), co + size_t(ip)*stride,
                            lv + ip, rv + ip, cut, rmin, rmax, ng, bx, by, bz)) return;
                if(want_pot) pot_g.add(e);
                for(int d=0; d<3; ++d) {
                    sens_g[d].add(-f[d]);
                    ps.add(h*w+d, f[d]);
                }
            });
            if(want_pot) l.flush(pot_g, pot, s);
            for(int d=0; d<3; ++d) l.flush(sens_g[d], ps, g*w+d);
        });
    }
};


// The cache of row r is the slab at row_base[r] of its system's environment entries (local indices)
struct DeviceScTable : public DeviceNode {
    DeviceCoord &out, &pos, &cb;
    int n_grid, n_theta;
    float grid_step_ang, theta_step, cutoff_ang, cache_buffer, box_x, box_y, box_z;
    DeviceArray<float> spline_coeff;
    DeviceArray<int> row_cb, row_layer, row_group, row_system, row_base;
    DeviceArray<int> env_atom, env_target, env_start, env_system, cb_system;
    DeviceArray<int> rebuild, n_cached, cached;
    DeviceArray<float> cached_cb, cached_env;
    bool first;

    DeviceScTable(DeviceCoord& out_, DeviceCoord& pos_, DeviceCoord& cb_, const vector<float>& spline_coeff_,
            int n_grid_, int n_theta_, float grid_step_ang_, float theta_step_, float cutoff_ang_,
            float cache_buffer_, float box_x_, float box_y_, float box_z_, const vector<int>& row_cb_,
            const vector<int>& row_layer_, const vector<int>& row_group_, const vector<int>& row_system_,
            const vector<int>& env_atom_, const vector<int>& env_target_, const vector<int>& env_start_):
        out(out_), pos(pos_), cb(cb_), n_grid(n_grid_), n_theta(n_theta_), grid_step_ang(grid_step_ang_),
        theta_step(theta_step_), cutoff_ang(cutoff_ang_), cache_buffer(cache_buffer_),
        box_x(box_x_), box_y(box_y_), box_z(box_z_), spline_coeff(spline_coeff_),
        row_cb(row_cb_), row_layer(row_layer_), row_group(row_group_), row_system(row_system_),
        env_atom(env_atom_), env_target(env_target_), env_start(env_start_), cb_system(cb_.system_of_element()),
        rebuild(env_start_.size()-1), n_cached(row_cb_.size()), cached_cb(3*size_t(cb_.n_elem)),
        cached_env(3*env_atom_.size()), first(true)
    {
        vector<int> base, es;
        size_t total = 0;
        for(int s: row_system_) {base.push_back(int(total)); total += env_start_[s+1]-env_start_[s];}
        for(int s=0; s+1<int(env_start_.size()); ++s) es.insert(es.end(), env_start_[s+1]-env_start_[s], s);
        row_base = DeviceArray<int>(base);
        env_system = DeviceArray<int>(es);
        cached = DeviceArray<int>(total);
        n_cached.zero();
    }

    ScTable table() const {
        return ScTable{spline_coeff.ptr, n_grid, n_theta, grid_step_ang, theta_step, cutoff_ang};
    }

    void find_contacts() {
        const float *xp = pos.output.ptr, *xc = cb.output.ptr;
        int wp = pos.row_width, wc = cb.row_width;
        const int *ea = env_atom.ptr, *es = env_start.ptr, *rsy = row_system.ptr, *rcb = row_cb.ptr;
        const int *rl = row_layer.ptr, *rbs = row_base.ptr;
        float bx = box_x, by = box_y, bz = box_z;
        if(!flag_all(rebuild, first, cache_buffer)) {
            flag_moved(rebuild, cache_buffer, cb.n_elem, xc, wc, DEVICE_LAMBDA (int i) {return i;}, cb_system.ptr,
                    cached_cb.ptr, bx, by, bz);
            flag_moved(rebuild, cache_buffer, int(env_atom.n), xp, wp, DEVICE_LAMBDA (int i) {return ea[i];},
                    env_system.ptr, cached_env.ptr, bx, by, bz);
        }

        const int* rb = rebuild.ptr;
        int *nc = n_cached.ptr, *ca = cached.ptr;
        float active_cutoff = cutoff_ang + (cache_buffer > 0.f ? cache_buffer : 0.f);
        float active_cutoff2 = sqr(active_cutoff);
        lanes_for(int(row_cb.n), DEVICE_LAMBDA (int r, Lanes l) {
            int s = rsy[r];
            if(!rb[s]) return;
            int count = 0;
            if(rl[r] >= 0) {
                auto cbp = make_vec3(xc[rcb[r]*wc+0], xc[rcb[r]*wc+1], xc[rcb[r]*wc+2]);
                l.chunks(es[s], es[s+1], [&](int i, bool in) {
                    bool keep = false;
                    if(in) {
                        Vec<3> dr = cbp - make_vec3(xp[ea[i]*wp+0], xp[ea[i]*wp+1], xp[ea[i]*wp+2]);
                        dr = simulation_box::minimum_image(dr, bx, by, bz);
                        keep = mag2(dr) < active_cutoff2;
                    }
                    int slot = l.place(keep, count);
                    if(keep) ca[rbs[r] + slot] = i - es[s];
                });
            }
            nc[r] = count;
        });
        store_cached(rebuild, cb.n_elem, xc, wc, DEVICE_LAMBDA (int i) {return i;}, cb_system.ptr, cached_cb.ptr);
        store_cached(rebuild, int(env_atom.n), xp, wp, DEVICE_LAMBDA (int i) {return ea[i];}, env_system.ptr,
                cached_env.ptr);
    }

    void compute_value(ComputeMode) override {
        find_contacts();
        auto tab = table();
        const float *xp = pos.output.ptr, *xc = cb.output.ptr;
        float* y = out.output.ptr;
        int wp = pos.row_width, wc = cb.row_width, wy = out.row_width;
        const int *ea = env_atom.ptr, *et = env_target.ptr, *es = env_start.ptr, *rsy = row_system.ptr;
        const int *rcb = row_cb.ptr, *rl = row_layer.ptr, *rg = row_group.ptr, *rbs = row_base.ptr;
        const int *nc = n_cached.ptr, *ca = cached.ptr;
        float bx = box_x, by = box_y, bz = box_z;
        lanes_for(int(row_cb.n), DEVICE_LAMBDA (int r, Lanes l) {
            float sum = 0.f;
            int s = rsy[r];
            auto cbp = make_vec3(xc[rcb[r]*wc+0], xc[rcb[r]*wc+1], xc[rcb[r]*wc+2]);
            auto cbv = make_vec3(xc[rcb[r]*wc+3], xc[rcb[r]*wc+4], xc[rcb[r]*wc+5]);
            l.chunks(0, nc[r], [&](int m, bool in) {
                // a contact past the table's reach adds +0, which leaves the sum's bits as skipping it does
                float value = 0.f;
                if(in) {
                    int i = es[s] + ca[rbs[r]+m];
                    Vec<3> dr = make_vec3(xp[ea[i]*wp+0], xp[ea[i]*wp+1], xp[ea[i]*wp+2]) - cbp;
                    dr = simulation_box::minimum_image(dr, bx, by, bz);
                    if(!tab.contact_value(value, rl[r] + et[i], dr, cbv)) value = 0.f;
                }
                sum = l.add_in_order(sum, value);
            });
            y[r*wy] = sum / float(rg[r]);
        });
    }

    void propagate_deriv() override {
        auto tab = table();
        const float *xp = pos.output.ptr, *xc = cb.output.ptr;
        auto s_out = out.sens.view(), ps = pos.sens.view(), cs = cb.sens.view();
        int wp = pos.row_width, wc = cb.row_width, wy = out.row_width;
        const int *ea = env_atom.ptr, *et = env_target.ptr, *es = env_start.ptr, *rsy = row_system.ptr;
        const int *rcb = row_cb.ptr, *rl = row_layer.ptr, *rg = row_group.ptr, *rbs = row_base.ptr;
        const int *nc = n_cached.ptr, *ca = cached.ptr;
        float bx = box_x, by = box_y, bz = box_z;
        lanes_for(int(row_cb.n), DEVICE_LAMBDA (int r, Lanes l) {
            float row_scale = s_out.get(r*wy);
            if(row_scale == 0.f) return;
            row_scale = row_scale / float(rg[r]);
            int s = rsy[r], c = rcb[r];
            auto cbp = make_vec3(xc[c*wc+0], xc[c*wc+1], xc[c*wc+2]);
            auto cbv = make_vec3(xc[c*wc+3], xc[c*wc+4], xc[c*wc+5]);
            AccumLocal sens_cb[6];
            l.chunks(0, nc[r], [&](int m, bool in) {
                if(!in) return;
                int i = es[s] + ca[rbs[r]+m];
                Vec<3> dr = make_vec3(xp[ea[i]*wp+0], xp[ea[i]*wp+1], xp[ea[i]*wp+2]) - cbp;
                dr = simulation_box::minimum_image(dr, bx, by, bz);
                Vec<3> grad_env, vector_grad;
                if(!tab.contact_deriv(grad_env, vector_grad, row_scale, rl[r] + et[i], dr, cbv)) return;
                for(int d=0; d<3; ++d) {
                    sens_cb[d].add(-grad_env[d]);
                    sens_cb[3+d].add(vector_grad[d]);
                    ps.add(ea[i]*wp+d, grad_env[d]);
                }
            });
            for(int d=0; d<6; ++d) l.flush(sens_cb[d], cs, c*wc+d);
        });
    }
};


// One of the two canonical springs of martini_potential.cpp: n_dep atoms per term
template <int n_dep>
struct DeviceMartiniSpring : public DeviceNode {
    AccumView potential;
    DeviceCoord& pos;
    DeviceArray<int> atom, system;
    DeviceArray<float> equil, spring_const, coefficients;
    float left, right, delta_min, coord_scale, box_x, box_y, box_z;
    int n_knot;

    DeviceMartiniSpring(AccumView potential_, DeviceCoord& pos_, const vector<int>& atom_,
            const vector<float>& equil_, const vector<float>& spring_const_, const vector<int>& system_,
            const LayeredClampedSpline1D<1>& spline, float delta_min_, float coord_scale_,
            float box_x_, float box_y_, float box_z_):
        potential(potential_), pos(pos_), atom(atom_), system(system_), equil(equil_), spring_const(spring_const_),
        coefficients(spline.coefficients), left(spline.left_clamped_value[0]), right(spline.right_clamped_value[0]),
        delta_min(delta_min_), coord_scale(coord_scale_), box_x(box_x_), box_y(box_y_), box_z(box_z_),
        n_knot(spline.nx) {}

    void compute_value(ComputeMode mode) override {
        VecArray x = pos.value();
        auto ps = pos.sens.view();
        auto pot = potential;
        bool want_pot = mode == PotentialAndDerivMode;
        int w = pos.row_width, nk = n_knot;
        const int *at = atom.ptr, *sys = system.ptr;
        const float *eq = equil.ptr, *k = spring_const.ptr, *co = coefficients.ptr;
        float lv = left, rv = right, dmin = delta_min, scale = coord_scale, bx = box_x, by = box_y, bz = box_z;
        parallel_for(int(system.n), DEVICE_LAMBDA (int nt) {
            Vec<3> g[n_dep];
            float energy;
            if(n_dep == 2)
                energy = dist_spring_term(g[0], load_vec<3>(x, at[n_dep*nt+0]), load_vec<3>(x, at[n_dep*nt+1]),
                        eq[nt], k[nt], co, &lv, &rv, nk, dmin, scale, bx, by, bz);
            else
                energy = angle_spring_term(g[0], g[1], g[n_dep-1], load_vec<3>(x, at[n_dep*nt+0]),
                        load_vec<3>(x, at[n_dep*nt+1]), load_vec<3>(x, at[n_dep*nt+n_dep-1]), eq[nt], k[nt],
                        co, &lv, &rv, nk, dmin, scale, bx, by, bz);
            if(n_dep == 2) g[1] = -g[0];
            if(want_pot) pot.add(sys[nt], energy);
            for(int j=0; j<n_dep; ++j)
                for(int d=0; d<3; ++d) ps.add(at[n_dep*nt+j]*w+d, g[j][d]);
        });
    }
};

}

DeviceNode* device_hybrid_position(DeviceCoord& out, DeviceCoord& source, DeviceCoord& hbond,
        const vector<unsigned char>& derived, const vector<int>& bb, const vector<int>& o, const vector<int>& o_elem,
        const vector<int>& src, const vector<float>& weight, const vector<float>& back_weight,
        const vector<float>& rescale) {
    return new DeviceHybridPosition(out, source, hbond, derived, bb, o, o_elem, src, weight, back_weight, rescale);
}
DeviceNode* device_martini_potential(AccumView potential, DeviceCoord& pos, const vector<int>& row_start,
        const vector<int>& neighbor, const vector<int>& spline, const vector<LayeredClampedSpline1D<1>>& splines,
        float r_min, float r_max, float cutoff, float cache_buffer, float box_x, float box_y, float box_z) {
    return new DeviceMartiniPotential(potential, pos, row_start, neighbor, spline, splines, r_min, r_max, cutoff,
            cache_buffer, box_x, box_y, box_z);
}
DeviceNode* device_sc_table_1body(DeviceCoord& out, DeviceCoord& pos, DeviceCoord& cb,
        const vector<float>& spline_coeff, int n_grid, int n_theta, float grid_step_ang, float theta_step,
        float cutoff_ang, float cache_buffer, float box_x, float box_y, float box_z, const vector<int>& row_cb,
        const vector<int>& row_layer, const vector<int>& row_group, const vector<int>& row_system,
        const vector<int>& env_atom, const vector<int>& env_target, const vector<int>& env_start) {
    return new DeviceScTable(out, pos, cb, spline_coeff, n_grid, n_theta, grid_step_ang, theta_step, cutoff_ang,
            cache_buffer, box_x, box_y, box_z, row_cb, row_layer, row_group, row_system, env_atom, env_target,
            env_start);
}
DeviceNode* device_dist_spring(AccumView potential, DeviceCoord& pos, const vector<int>& atom,
        const vector<float>& equil_dist, const vector<float>& spring_const, const vector<int>& system,
        const LayeredClampedSpline1D<1>& spline, float delta_min, float coord_scale,
        float box_x, float box_y, float box_z) {
    return new DeviceMartiniSpring<2>(potential, pos, atom, equil_dist, spring_const, system, spline, delta_min,
            coord_scale, box_x, box_y, box_z);
}
DeviceNode* device_angle_spring(AccumView potential, DeviceCoord& pos, const vector<int>& atom,
        const vector<float>& cos_theta0, const vector<float>& spring_const, const vector<int>& system,
        const LayeredClampedSpline1D<1>& spline, float delta_min, float coord_scale,
        float box_x, float box_y, float box_z) {
    return new DeviceMartiniSpring<3>(potential, pos, atom, cos_theta0, spring_const, system, spline, delta_min,
            coord_scale, box_x, box_y, box_z);
}
