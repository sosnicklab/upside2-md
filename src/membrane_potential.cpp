#include "deriv_engine.h"
#include <string>
#include "timing.h"
#include "affine.h"
#include <cmath>
#include "h5_support.h"
#include <vector>
#include "spline.h"
#include <iostream>

#define N_AMINO_ACID_TYPE 20

using namespace h5;
using namespace std;

struct MembranePotential : public PotentialNode
{
    struct ResidueParams {
        // Logically the two residues are the same, but they may have different
        // indices in the CB and Env outputs respectively
        index_t cb_index;
        index_t env_index;
        int restype;
    };

    struct PotentialCBParams {
        float cov_midpoint;
        float cov_sharpness;
    };

    int n_elem;        // number of residues to process
    int n_restype;
    int n_donor, n_acceptor;

    CoordNode& res_pos;  // CB atom
    CoordNode& environment_coverage;
    CoordNode& protein_hbond;

    vector<ResidueParams>     res_params;
    vector<PotentialCBParams> pot_params;
    vector<float> weights;

    LayeredClampedSpline1D<1> membrane_energy_cb_spline;
    LayeredClampedSpline1D<1> membrane_energy_uhb_spline;

    // shift and scale to convert z coordinates to spline coordinates
    float cb_z_shift, cb_z_scale;
    float uhb_z_shift, uhb_z_scale;

    MembranePotential(hid_t grp, CoordNode& res_pos_,
                                 CoordNode& environment_coverage_,
                                 CoordNode& protein_hbond_):
        PotentialNode(),

        n_elem    (get_dset_size(1, grp,             "cb_index")[0]),
        n_restype (get_dset_size(2, grp,            "cb_energy")[0]),
        n_donor   (get_dset_size(1, grp,    "donor_residue_ids")[0]),
        n_acceptor(get_dset_size(1, grp, "acceptor_residue_ids")[0]),

        res_pos(res_pos_),
        environment_coverage(environment_coverage_),
        protein_hbond(protein_hbond_),

        res_params(n_elem),
        pot_params(n_restype),
        weights(get_dset_size(1, grp, "weights")[0]),

        membrane_energy_cb_spline(
                get_dset_size(2, grp, "cb_energy")[0],
                get_dset_size(2, grp, "cb_energy")[1]),

        membrane_energy_uhb_spline(
                get_dset_size(2, grp, "uhb_energy")[0],
                get_dset_size(2, grp, "uhb_energy")[1]),

        cb_z_shift(-read_attribute<float>(grp, "cb_energy", "z_min")),
        cb_z_scale((membrane_energy_cb_spline.nx-1)/(read_attribute<float>(grp, "cb_energy", "z_max")+cb_z_shift)),

        uhb_z_shift(-read_attribute<float>(grp, "uhb_energy", "z_min")),
        uhb_z_scale((membrane_energy_uhb_spline.nx-1)/(read_attribute<float>(grp, "uhb_energy", "z_max")+uhb_z_shift))
    {
        check_elem_width_lower_bound(res_pos, 3);
        check_elem_width_lower_bound(environment_coverage, 1);

        check_size(grp,      "cb_index",    n_elem);
        check_size(grp,     "env_index",    n_elem);
        check_size(grp,  "residue_type",    n_elem);
        check_size(grp,  "cov_midpoint", n_restype);
        check_size(grp, "cov_sharpness", n_restype);
        check_size(grp,     "cb_energy", n_restype, membrane_energy_cb_spline.nx);
        check_size(grp,    "uhb_energy",         2, membrane_energy_uhb_spline.nx); // type 0 for unpaired donor, type 1 for unpaired acceptor

        traverse_dset<1,  int>(grp,      "cb_index", [&](size_t nr,   int  x) {res_params[nr].cb_index  = x;});
        traverse_dset<1,  int>(grp,     "env_index", [&](size_t nr,   int  x) {res_params[nr].env_index = x;});
        traverse_dset<1,  int>(grp,  "residue_type", [&](size_t nr,   int rt) {res_params[nr].restype   = rt;});
        traverse_dset<1,float>(grp,  "cov_midpoint", [&](size_t rt, float bc) {pot_params[rt].cov_midpoint  = bc;});
        traverse_dset<1,float>(grp, "cov_sharpness", [&](size_t rt, float bw) {pot_params[rt].cov_sharpness = bw;});
        traverse_dset<1,float>(grp,       "weights", [&](size_t n,  float x)  {weights[n]=x;});

        vector<double> cb_energy_data;
        traverse_dset<2,double>(grp, "cb_energy", [&](size_t rt, size_t z_index, double value) {
                cb_energy_data.push_back(value);});
        membrane_energy_cb_spline.fit_spline(cb_energy_data.data());

        vector<double> uhb_energy_data;
        traverse_dset<2,double>(grp, "uhb_energy", [&](size_t rt, size_t z_index, double value) {
                uhb_energy_data.push_back(value);});
        membrane_energy_uhb_spline.fit_spline(uhb_energy_data.data());
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("membrane_potential"));

        VecArray cb_pos       = res_pos.output;
        VecArray cb_pos_sens  = res_pos.sens;
        VecArray env_cov      = environment_coverage.output;
        VecArray env_cov_sens = environment_coverage.sens;
        VecArray hb_pos       = protein_hbond.output;
        VecArray hb_sens      = protein_hbond.sens;

        potential = 0.f;

        for(int nr=0; nr<n_elem; ++nr) {
            auto &p = res_params[nr];
            float cb_z = cb_pos(2, p.cb_index);

            float result[2];    // deriv then value
            membrane_energy_cb_spline.evaluate_value_and_deriv(result, p.restype,
                    (cb_z + cb_z_shift) * cb_z_scale);
            float spline_value = result[1];
            float spline_deriv = result[0]*cb_z_scale; // scale to get derivative in *unnormalized* coordinates

            auto& pp = pot_params[p.restype];

            float bl = 0.0;
            for(int aa: range(20)) 
                bl += weights[aa] * env_cov(0, p.env_index*20+aa);

            auto cover_sig = compact_sigmoid(bl-pp.cov_midpoint, pp.cov_sharpness);

            potential                    += spline_value*cover_sig.x();
            cb_pos_sens (2, p.cb_index)  += spline_deriv*cover_sig.x();

            auto d_pot  = spline_value*cover_sig.y();
            for(int aa: range(20)) 
                env_cov_sens(0, p.env_index*20+aa) +=  weights[aa] * d_pot;
        }

        int n_virtual = n_donor+n_acceptor;
        for(int nv=0; nv<n_virtual; ++nv) {
            float hb_z    = hb_pos(2, nv);
            float hb_prob = hb_pos(6, nv);  // probability that his virtual participates in any HBond

            float result[2];
            membrane_energy_uhb_spline.evaluate_value_and_deriv(result, int(nv>=n_donor),
                    (hb_z + uhb_z_shift) * uhb_z_scale);
            float spline_value = result[1];
            float spline_deriv = result[0]*uhb_z_scale;

            float uhb_prob  = 1.f-hb_prob;

            potential      += spline_value*sqr(uhb_prob);
            hb_sens(2, nv) += spline_deriv*sqr(uhb_prob);
            hb_sens(6, nv) += -2.f*spline_value*uhb_prob; 
        }
    }

    virtual std::vector<float> get_param() const override {
        vector<float> params(40, 0.f);
        for(int n=0;n<20;n++) {
            auto& pp = pot_params[n];
            params[n] = pp.cov_midpoint;
            params[20+n] = pp.cov_sharpness;
        }
        return params;
    }

    virtual void set_param(const std::vector<float>& new_params) override {
        for(int n=0;n<20;n++) {
            auto& pp = pot_params[n];
            pp.cov_midpoint = new_params[n];
            pp.cov_sharpness = new_params[20+n];
        }
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        vector<float> deriv(40, 0.f);

        VecArray cb_pos  = res_pos.output;
        VecArray env_cov = environment_coverage.output;

        for(int nr=0; nr<n_elem; ++nr) {
            auto &p = res_params[nr];
            float cb_z = cb_pos(2, p.cb_index);

            float result[2]; // deriv then value
            membrane_energy_cb_spline.evaluate_value_and_deriv(result, p.restype, (cb_z + cb_z_shift) * cb_z_scale);
            float spline_value = result[1];

            float bl = 0.0;
            for(int aa: range(20)) 
                bl += weights[aa] * env_cov(0, p.env_index*20+aa);

            auto& pp        = pot_params[p.restype];
            auto  env_coord = bl-pp.cov_midpoint;
            auto cover_sig  = compact_sigmoid(env_coord, pp.cov_sharpness);

            deriv[p.restype]    += -spline_value * cover_sig.y();
            deriv[p.restype+20] +=  spline_value * env_coord/pp.cov_sharpness * cover_sig.y();
        }

        return deriv;
    }
#endif

};

static RegisterNodeType<MembranePotential, 3> membrane_potential_node("membrane_potential");


struct MembraneCBPotential : public PotentialNode
{
    int n_elem;        // number of residues to process
    int n_restype;

    vector<int> cb_index;
    vector<int> res_type;

    CoordNode& res_pos;  // CB atom
    CoordNode& environment_coverage;
    CoordNode& center_of_curvature;
    vector<float> weights;

    int n_bl;
    vector<float> bl_nodes;
    float bl_range_sharpness;
    float left_right_node;
    float left_right_sharpness;

    bool  use_curvature;
    float curvature_radius;
    float curvature_sign;

    int n_node;
    float zstart;
    float zscale;

    vector<float> coeff;      // length n_restype*n_coeff

    MembraneCBPotential(hid_t grp, CoordNode& res_pos_,
                                   CoordNode& environment_coverage_,
                                   CoordNode& centre_of_curvature_ ):
        PotentialNode(),
        n_elem              (get_dset_size(1, grp, "cb_index")[0]),
        n_restype           (get_dset_size(3, grp, "coeff")[0]),
        cb_index            (n_elem),
        res_type            (n_elem),
        res_pos             (res_pos_),
        environment_coverage(environment_coverage_),
        center_of_curvature (centre_of_curvature_),
        weights             (get_dset_size(1, grp, "weights")[0]),
        n_bl                (get_dset_size(3, grp, "coeff")[1]),
        bl_nodes            (n_bl-1),
        bl_range_sharpness  (read_attribute<float>(grp, ".", "bl_range_sharpness" )),
        left_right_node     (read_attribute<float>(grp, ".", "left_right_node" )),
        left_right_sharpness(read_attribute<float>(grp, ".", "left_right_sharpness" )),
        use_curvature       (read_attribute<int>  (grp, ".", "use_curvature" )),
        curvature_radius    (read_attribute<float>(grp, ".", "curvature_radius" )),
        curvature_sign      (read_attribute<float>(grp, ".", "curvature_sign" )),
        n_node              (get_dset_size(3, grp, "coeff")[2]),
        zstart              (read_attribute<float>(grp, ".", "z_start" )),
        zscale              (read_attribute<float>(grp, ".", "z_scale" )),
        coeff               (n_restype*n_bl*n_node)
    {
        check_elem_width_lower_bound(res_pos, 3);
        check_elem_width_lower_bound(environment_coverage, 1);

        check_size(grp, "res_type", n_elem);
        check_size(grp, "bl_nodes", n_bl-1);
        check_size(grp, "coeff",    n_restype, n_bl, n_node);

        traverse_dset<1, float>(grp, "cb_index",  [&](size_t n, float x) {cb_index[n] = x;});
        traverse_dset<1, int>  (grp, "res_type",  [&](size_t n, int x)   {res_type[n] = x;});
        traverse_dset<1, float>(grp, "weights",   [&](size_t n, float x) {weights[n]  = x;});
        traverse_dset<1, float>(grp, "bl_nodes",  [&](size_t n, float x) {bl_nodes[n] = x;});

        traverse_dset<3, float>(grp, "coeff",[&](size_t n, size_t l, size_t c, float x) 
                        { coeff[n*n_node*n_bl + l*n_node + c] = x;} );
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("cb_membrane_potential"));

        VecArray cb_pos       = res_pos.output;
        VecArray cb_pos_sens  = res_pos.sens;
        VecArray env_cov      = environment_coverage.output;
        VecArray env_cov_sens = environment_coverage.sens;

        potential = 0.f;

        vector<float> bl_score1(n_bl-1);
        vector<float> bl_score2(n_bl-1);
        vector<float> d_bl_score1(n_bl-1);
        vector<float> d_bl_score2(n_bl-1);
        vector<float> bl_score(n_bl);
        vector<float> d_bl_score(n_bl);

        vector<float> fs(n_bl);
        vector<float> dfs(n_bl);

        auto ccenter = load_vec<3>(center_of_curvature.output, 0);
        curvature_radius = ccenter.z()*-1.f*curvature_sign;

        for(int nr=0; nr<n_elem; ++nr) {
            int rt = res_type[nr];
            int ri = cb_index[nr];
        
            // burial level
            float bl = 0.0;
            for(int aa: range(n_restype)) 
                bl += weights[aa] * env_cov(0, ri*n_restype+aa);

            for (int i: range(n_bl-1)) {
                auto s1 = compact_sigmoid(bl-bl_nodes[i], bl_range_sharpness);
                auto s2 = compact_sigmoid(bl_nodes[i]-bl, bl_range_sharpness);
                bl_score1[i]   = s1.x();
                d_bl_score1[i] = s1.y();
                bl_score2[i]   = s2.x();
                d_bl_score2[i] = -s2.y();
            }

            bl_score[0]   = bl_score1[0];
            d_bl_score[0] = d_bl_score1[0];
            bl_score[n_bl-1]   = bl_score2[n_bl-2];
            d_bl_score[n_bl-1] = d_bl_score2[n_bl-2];
            for (int i: range(n_bl-2)) {
                bl_score[i+1]   = bl_score2[i]*bl_score1[i+1];
                d_bl_score[i+1] = d_bl_score2[i]*bl_score1[i+1] + bl_score2[i]*d_bl_score1[i+1];
            }

            // z
            float cb_z = cb_pos(2, ri);
            float dist = 0.f;
            float r_dist = 1.f;
            auto xyz = load_vec<3>(cb_pos, ri);

            if (use_curvature) {
                dist = mag( xyz - ccenter );
                r_dist = 1.f/dist;
                cb_z = curvature_sign* ( dist - curvature_radius );
            }

            auto sig_left = compact_sigmoid(cb_z-left_right_node, left_right_sharpness);
            float ll  = sig_left.x();
            float lr  = 1.f-sig_left.x();
            float dll = sig_left.y();
            float dlr = -sig_left.y();

            auto coord_left  = ( cb_z-zstart)*zscale;
            auto coord_right = (-cb_z-zstart)*zscale;

            float df_dx  = 0.0;
            float df_dy  = 0.0;
            float df_dz  = 0.0;
            float df_dbl = 0.0;
            for (int i: range(n_bl)) {
                int shift = (rt*n_bl + i)*n_node;
                auto v = clamped_deBoor_value_and_deriv(coeff.data() + shift, coord_left, n_node);
                float f_left  = v.x();
                float df_left = v.y()*zscale;

                v = clamped_deBoor_value_and_deriv(coeff.data() + shift, coord_right, n_node);
                float f_right  = v.x();
                float df_right = -v.y()*zscale;

                float fs  = ll*f_left + lr*f_right;
                float dfs = ll*df_left + dll*f_left + lr*df_right + dlr*f_right;

                potential +=   bl_score[i] *  fs;
                df_dbl    += d_bl_score[i] *  fs;
                if (use_curvature) {
                    float sen_d = curvature_sign * bl_score[i] * dfs * r_dist; 
                    df_dx += sen_d * xyz.x();
                    df_dy += sen_d * xyz.y();
                    df_dz += sen_d * (xyz.z()-curvature_sign*curvature_radius*-1.f);
                }
                else
                    df_dz += bl_score[i] * dfs;
            }

            if (use_curvature) {
                cb_pos_sens (0, ri)  += df_dx;
                cb_pos_sens (1, ri)  += df_dy;
            }
            cb_pos_sens (2, ri)  += df_dz;
            for(int aa: range(n_restype)) 
                env_cov_sens(0, ri*n_restype+aa) += weights[aa] * df_dbl;
        }
    }

    virtual std::vector<float> get_param() const override {
        int np = n_restype*n_bl*n_node;
        vector<float> params(np, 0.f);
        for(int n=0;n<np;n++) {
            params[n] = coeff[n];
        }
        return params;
    }

    virtual void set_param(const std::vector<float>& new_params) override {
        int np = n_restype*n_bl*n_node;
        for(int n=0;n<np;n++) {
            coeff[n] = new_params[n];
        }
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        int np = n_restype*n_bl*n_node;
        vector<float> deriv(np, 0.f);

        VecArray cb_pos  = res_pos.output;
        VecArray env_cov = environment_coverage.output;

        vector<float> bl_score1(n_bl-1);
        vector<float> bl_score2(n_bl-1);
        vector<float> d_bl_score1(n_bl-1);
        vector<float> d_bl_score2(n_bl-1);
        vector<float> bl_score(n_bl);
        vector<float> d_bl_score(n_bl);

        auto ccenter = load_vec<3>(center_of_curvature.output, 0);
        curvature_radius = ccenter.z()*-1.f*curvature_sign;

        for(int nr=0; nr<n_elem; ++nr) {
            int rt = res_type[nr];
            int ri = cb_index[nr];

            float bl = 0.0;
            for(int aa: range(n_restype)) 
                bl += weights[aa] * env_cov(0, ri*n_restype+aa);
            for (int i: range(n_bl-1)) {
                auto s1 = compact_sigmoid(bl-bl_nodes[i], bl_range_sharpness);
                auto s2 = compact_sigmoid(bl_nodes[i]-bl, bl_range_sharpness);
                bl_score1[i]   = s1.x();
                bl_score2[i]   = s2.x();
            }
            bl_score[0]       = bl_score1[0];
            for (int i: range(n_bl-2)) 
                bl_score[i+1] = bl_score2[i]*bl_score1[i+1];
            bl_score[n_bl-1]  = bl_score2[n_bl-2];

            float cb_z        = cb_pos(2, ri);
            if (use_curvature) {
                auto xyz = load_vec<3>(cb_pos, ri);
                float dist = mag( xyz - ccenter );
                cb_z = curvature_sign* ( dist - curvature_radius );
            }
            auto  sig_left    = compact_sigmoid(cb_z-left_right_node, left_right_sharpness);
            float ll          = sig_left.x();
            float lr          = 1.f-sig_left.x();

            auto  coord_left  = ( cb_z-zstart)*zscale;
            auto  coord_right = (-cb_z-zstart)*zscale;

            int starting_bin;
            float result[4];
            clamped_deBoor_coeff_deriv(&starting_bin, result, coord_left, n_node);
            for (int bi: range(n_bl)) {
                int shift = (rt*n_bl + bi)*n_node;
                for(int i: range(4)) deriv[shift+starting_bin+i] += result[i]*bl_score[bi]*ll;
            }

            clamped_deBoor_coeff_deriv(&starting_bin, result, coord_right, n_node);
            for (int bi: range(n_bl)) {
                int shift = (rt*n_bl + bi)*n_node + np;
                for(int i: range(4)) deriv[shift+starting_bin+i] += result[i]*bl_score[bi]*lr;
            }
        }
        return deriv;
    }
#endif
};

static RegisterNodeType<MembraneCBPotential, 3> membrane_cb_potential_node("cb_membrane_potential");

struct MembraneHBPotential : public PotentialNode
{
    int n_elem;        // number of residues to process

    vector<int> hb_index;
    vector<int> hb_type;

    CoordNode& protein_hbond;
    CoordNode& center_of_curvature;

    int   n_node;
    float zstart;
    float zscale;
    float left_right_node;
    float left_right_sharpness;

    bool  use_curvature;
    float curvature_radius;
    float curvature_sign;

    vector<float> coeff;      // length n_restype*n_coeff

    MembraneHBPotential(hid_t grp, CoordNode& protein_hbond_, CoordNode& center_of_curvature_):
        PotentialNode(),
        n_elem              (get_dset_size(1, grp, "hb_index")[0]),
        hb_index            (n_elem),
        hb_type             (n_elem),
        protein_hbond       (protein_hbond_),
        center_of_curvature (center_of_curvature_),
        n_node              (get_dset_size(3, grp, "coeff")[2]),
        zstart              (read_attribute<float>(grp, ".", "z_start" )),
        zscale              (read_attribute<float>(grp, ".", "z_scale" )),
        left_right_node     (read_attribute<float>(grp, ".", "left_right_node" )),
        left_right_sharpness(read_attribute<float>(grp, ".", "left_right_sharpness" )),
        use_curvature       (read_attribute<int>  (grp, ".", "use_curvature" )),
        curvature_radius    (read_attribute<float>(grp, ".", "curvature_radius" )),
        curvature_sign      (read_attribute<float>(grp, ".", "curvature_sign" )),
        coeff               (2*2*n_node)
    {
        check_size(grp, "hb_type",  n_elem);
        check_size(grp, "coeff",  2, 2, n_node);
        traverse_dset<1, float>( grp, "hb_index", [&](size_t n, float x) {hb_index[n] = x;} );
        traverse_dset<1, int>  ( grp, "hb_type",  [&](size_t n, int x)   {hb_type[n]  = x;} );

        traverse_dset<3, float>(grp, "coeff",[&](size_t n, size_t l, size_t c, float x) 
                        { coeff[(n*2+l)*n_node + c] = x;} );
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("hb_membrane_potential"));

        VecArray hb_pos  = protein_hbond.output;
        VecArray hb_sens = protein_hbond.sens;

        potential = 0.f;

        auto ccenter = load_vec<3>(center_of_curvature.output, 0);
        curvature_radius = ccenter.z()*-1.f*curvature_sign;

        for(int nv=0; nv<n_elem; ++nv) {

            int   rt          = hb_type[nv];
            int   ri          = hb_index[nv];
            float hb_z        = hb_pos(2, ri);
            float hb_prob     = hb_pos(6, ri);  // probability that his virtual participates in any HBond

            float dist = 0.f;
            float r_dist = 1.f;
            auto xyz = load_vec<3>(hb_pos, ri);
            if (use_curvature) {
                dist = mag( xyz - ccenter );
                r_dist = 1.f/dist;
                hb_z = curvature_sign* ( dist - curvature_radius );
            }

            auto  coord_left  = ( hb_z-zstart)*zscale;
            auto  coord_right = (-hb_z-zstart)*zscale;
        
            auto sig_left = compact_sigmoid(hb_z-left_right_node, left_right_sharpness);
            float ll      = sig_left.x();
            float lr      = 1.f-sig_left.x();
            float dll     = sig_left.y();
            float dlr     = -sig_left.y();

            auto f1l   = clamped_deBoor_value_and_deriv(coeff.data() + (rt*2)*n_node, coord_left,  n_node);
            auto f1r   = clamped_deBoor_value_and_deriv(coeff.data() + (rt*2)*n_node, coord_right, n_node);
            float fs1  = ll*f1l.x() + lr*f1r.x();
            float dfs1 = ll*f1l.y()*zscale + dll*f1l.x() - lr*f1r.y()*zscale + dlr*f1r.x();

            float uhb_prob      = 1.f-hb_prob;
            float sqr_uhb_prob  = sqr(uhb_prob);
            potential          += sqr_uhb_prob  *  fs1;
            hb_sens(6, ri)     += -2.f*uhb_prob *  fs1;


            auto f2l   = clamped_deBoor_value_and_deriv(coeff.data() + (rt*2+1)*n_node, coord_left,  n_node);
            auto f2r   = clamped_deBoor_value_and_deriv(coeff.data() + (rt*2+1)*n_node, coord_right, n_node);
            float fs2  = ll*f2l.x() + lr*f2r.x();
            float dfs2 = ll*f2l.y()*zscale + dll*f2l.x() - lr*f2r.y()*zscale + dlr*f2r.x();

            float sqr_hb_prob   = 1.f - sqr(1.f-hb_prob);
            potential          += sqr_hb_prob *  fs2;
            hb_sens(6, ri)     += 2.f*uhb_prob *  fs2;

            if (use_curvature) {
                float sen_d = curvature_sign * sqr_uhb_prob  * dfs1 * r_dist; 
                hb_sens(0, ri) += sen_d * xyz.x();
                hb_sens(1, ri) += sen_d * xyz.y();
                hb_sens(2, ri) += sen_d * (xyz.z()-curvature_sign*curvature_radius*-1.f);
                sen_d = curvature_sign * sqr_hb_prob  * dfs2 * r_dist; 
                hb_sens(0, ri) += sen_d * xyz.x();
                hb_sens(1, ri) += sen_d * xyz.y();
                hb_sens(2, ri) += sen_d * (xyz.z()-curvature_sign*curvature_radius*-1.f);
            }
            else {
                hb_sens(2, ri) += sqr_uhb_prob  * dfs1;
                hb_sens(2, ri) += sqr_hb_prob * dfs2;
            }
        }
    }

    virtual std::vector<float> get_param() const override {
        int np = 2*2*n_node;
        vector<float> params(np, 0.f);
        for(int n=0;n<np;n++) {
            params[n] = coeff[n];
        }
        return params;
    }

    virtual void set_param(const std::vector<float>& new_params) override {
        int np = 2*2*n_node;
        for(int n=0;n<np;n++) {
            coeff[n]  = new_params[n];
        }
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        int np = 2*2*n_node;
        vector<float> deriv(np, 0.f);
        VecArray hb_pos = protein_hbond.output;

        auto ccenter = load_vec<3>(center_of_curvature.output, 0);
        curvature_radius = ccenter.z()*-1.f*curvature_sign;

        for(int nv=0; nv<n_elem; ++nv) {

            int   rt           = hb_type[nv];
            int   ri           = hb_index[nv];

            float hb_z         = hb_pos(2, ri);
            if (use_curvature) {
                auto xyz = load_vec<3>(hb_pos, ri);
                float dist = mag( xyz - ccenter );
                hb_z = curvature_sign* ( dist - curvature_radius );
            }

            auto  coord_left   = ( hb_z-zstart)*zscale;
            auto  coord_right  = (-hb_z-zstart)*zscale;
            auto  sig_left     = compact_sigmoid(hb_z-left_right_node, left_right_sharpness);
            float ll           = sig_left.x();
            float lr           = 1.f-sig_left.x();

            float hb_prob      = hb_pos(6, ri);  // probability that his virtual participates in any HBond
            float uhb_prob     = 1.f-hb_prob;
            float sqr_uhb_prob = sqr(uhb_prob);
            float sqr_hb_prob  = 1.f - sqr_uhb_prob;

            int starting_bin;
            float result[4];
            clamped_deBoor_coeff_deriv(&starting_bin, result, coord_left, n_node);
            for(int i: range(4)) deriv[(rt*2+0)*n_node+starting_bin+i] += result[i]*ll*sqr_uhb_prob;
            for(int i: range(4)) deriv[(rt*2+1)*n_node+starting_bin+i] += result[i]*ll*sqr_hb_prob;

            clamped_deBoor_coeff_deriv(&starting_bin, result, coord_right, n_node);
            for(int i: range(4)) deriv[(rt*2+0)*n_node+starting_bin+i] += result[i]*lr*sqr_uhb_prob;
            for(int i: range(4)) deriv[(rt*2+1)*n_node+starting_bin+i] += result[i]*lr*sqr_hb_prob;
        }
        return deriv;
    }
#endif
};

static RegisterNodeType<MembraneHBPotential, 2> membrane_hb_potential_node("hb_membrane_potential");

struct MembraneSurfCBPotential : public PotentialNode
{
    int n_elem;        // number of residues to process
    int n_restype;

    vector<int> cb_index;
    vector<int> res_type;

    CoordNode& res_pos;  // CB atom
    CoordNode& environment_coverage;
    CoordNode& surface;
    CoordNode& center_of_curvature;

    vector<float> weights;

    int n_bl;
    vector<float> bl_nodes;
    float bl_range_sharpness;
    float left_right_node;
    float left_right_sharpness;

    bool  use_curvature;
    float curvature_radius;
    float curvature_sign;

    int n_node;
    float zstart;
    float zscale;
    float half_thickness;

    vector<float> coeff;
    vector<float> coeff_inner;

    MembraneSurfCBPotential(hid_t grp, CoordNode& res_pos_,
                                       CoordNode& environment_coverage_,
                                       CoordNode& surface_,
                                       CoordNode& center_of_curvature_):
        PotentialNode(),
        n_elem               (get_dset_size(1, grp, "cb_index")[0]),
        n_restype            (get_dset_size(3, grp, "coeff")[0]),
        cb_index             (n_elem),
        res_type             (n_elem),
        res_pos              (res_pos_),
        environment_coverage (environment_coverage_),
        surface              (surface_),
        center_of_curvature  (center_of_curvature_),
        weights              (get_dset_size(1, grp, "weights")[0]),
        n_bl                 (get_dset_size(3, grp, "coeff")[1]),
        bl_nodes             (n_bl-1),
        bl_range_sharpness   (read_attribute<float>(grp, ".", "bl_range_sharpness" )),
        left_right_node      (read_attribute<float>(grp, ".", "left_right_node" )),
        left_right_sharpness (read_attribute<float>(grp, ".", "left_right_sharpness" )),
        use_curvature        (read_attribute<int>  (grp, ".", "use_curvature" )),
        curvature_radius     (read_attribute<float>(grp, ".", "curvature_radius" )),
        curvature_sign       (read_attribute<float>(grp, ".", "curvature_sign" )),
        n_node               (get_dset_size(3, grp, "coeff")[2]),
        zstart               (read_attribute<float>(grp, ".", "z_start" )),
        zscale               (read_attribute<float>(grp, ".", "z_scale" )),
        half_thickness       (read_attribute<float>(grp, ".", "half_thickness" )),
        coeff                (n_restype*n_bl*n_node),
        coeff_inner          (n_restype*n_node)
    {
        check_elem_width_lower_bound(res_pos, 3);
        check_elem_width_lower_bound(environment_coverage, 1);

        check_size(grp, "res_type",    n_elem);
        check_size(grp, "bl_nodes",    n_bl-1);
        check_size(grp, "coeff",       n_restype, n_bl, n_node);
        check_size(grp, "coeff_inner", n_restype, n_node);

        traverse_dset<1, float>(grp, "cb_index", [&](size_t n, float x) {cb_index[n] = x;});
        traverse_dset<1, int>  (grp, "res_type", [&](size_t n, int x)   {res_type[n] = x;});
        traverse_dset<1, float>(grp, "weights",  [&](size_t n, float x) {weights[n]  = x;});
        traverse_dset<1, float>(grp, "bl_nodes", [&](size_t n, float x) {bl_nodes[n] = x;});

        traverse_dset<3, float>(grp, "coeff",    [&](size_t n, size_t l, size_t c, float x) 
                        { coeff[n*n_node*n_bl + l*n_node + c] = x;} );

        traverse_dset<2, float>(grp, "coeff_inner", [&](size_t n, size_t c, float x) 
                        { coeff_inner[n*n_node + c] = x;} );
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("cb_surf_membrane_potential"));

        VecArray cb_pos       = res_pos.output;
        VecArray cb_pos_sens  = res_pos.sens;
        VecArray env_cov      = environment_coverage.output;
        VecArray env_cov_sens = environment_coverage.sens;

        VecArray surf         = surface.output;

        potential = 0.f;

        vector<float> bl_score1(n_bl-1);
        vector<float> bl_score2(n_bl-1);
        vector<float> d_bl_score1(n_bl-1);
        vector<float> d_bl_score2(n_bl-1);
        vector<float> bl_score(n_bl);
        vector<float> d_bl_score(n_bl);

        vector<float> fs(n_bl);
        vector<float> dfs(n_bl);

        auto ccenter = load_vec<3>(center_of_curvature.output, 0);
        curvature_radius = ccenter.z()*-1.f*curvature_sign;

        for(int nr=0; nr<n_elem; ++nr) {
            int rt = res_type[nr];
            int ri = cb_index[nr];
        
            // burial level
            float bl = 0.0;
            for(int aa: range(n_restype)) 
                bl += weights[aa] * env_cov(0, ri*n_restype+aa);

            for (int i: range(n_bl-1)) {
                auto s1 = compact_sigmoid(bl-bl_nodes[i], bl_range_sharpness);
                auto s2 = compact_sigmoid(bl_nodes[i]-bl, bl_range_sharpness);
                bl_score1[i]   = s1.x();
                d_bl_score1[i] = s1.y();
                bl_score2[i]   = s2.x();
                d_bl_score2[i] = -s2.y();
            }

            bl_score[0]        = bl_score1[0];
            d_bl_score[0]      = d_bl_score1[0];
            bl_score[n_bl-1]   = bl_score2[n_bl-2];
            d_bl_score[n_bl-1] = d_bl_score2[n_bl-2];
            for (int i: range(n_bl-2)) {
                bl_score[i+1]   = bl_score2[i]*bl_score1[i+1];
                d_bl_score[i+1] = d_bl_score2[i]*bl_score1[i+1] + bl_score2[i]*d_bl_score1[i+1];
            }

            // z
            float cb_z = cb_pos(2, ri);
            float dist = 0.f;
            float r_dist = 1.f;
            auto xyz = load_vec<3>(cb_pos, ri);

            if (use_curvature) {
                dist = mag( xyz - ccenter );
                r_dist = 1.f/dist;
                cb_z = curvature_sign* ( dist - curvature_radius );
            }

            auto sig_left = compact_sigmoid(cb_z-left_right_node, left_right_sharpness);
            float ll  = sig_left.x();
            float lr  = 1.f-sig_left.x();
            float dll = sig_left.y();
            float dlr = -sig_left.y();

            auto coord_left  = ( cb_z-zstart)*zscale;
            auto coord_right = (-cb_z-zstart)*zscale;

            float df_dx  = 0.0;
            float df_dy  = 0.0;
            float df_dz  = 0.0;
            float df_dbl = 0.0;
            float surfv  = 0.0;

            for (int i: range(n_bl)) {

                if (i>0 and (cb_z > half_thickness or cb_z < -half_thickness)) 
                    surfv = 1.0;
                else
                    surfv = surf(0, ri);

                int shift = (rt*n_bl + i)*n_node;
                auto v = clamped_deBoor_value_and_deriv(coeff.data() + shift, coord_left, n_node);
                float f_left  = v.x() * surfv;
                float df_left = v.y() * zscale * surfv;

                v = clamped_deBoor_value_and_deriv(coeff.data() + shift, coord_right, n_node);
                float f_right  =  v.x() * surfv;
                float df_right = -v.y() * surfv * zscale;

                if ( surfv < 1.f) {
                    int shift = rt*n_node;
                    auto v    = clamped_deBoor_value_and_deriv(coeff_inner.data() + shift, coord_left, n_node);
                    f_left   += v.x() * (1.f-surfv);
                    df_left  += v.y() * (1.f-surfv) * zscale;
        
                    v         = clamped_deBoor_value_and_deriv(coeff_inner.data() + shift, coord_right, n_node);
                    f_right  += v.x() * (1.f-surfv);
                    df_right -= v.y() * (1.f-surfv) * zscale;
                }

                float fs  = ll*f_left + lr*f_right;
                float dfs = ll*df_left + dll*f_left + lr*df_right + dlr*f_right;

                potential +=   bl_score[i] *  fs;
                df_dbl    += d_bl_score[i] *  fs;
                if (use_curvature) {
                    float sen_d = curvature_sign * bl_score[i] * dfs * r_dist; 
                    df_dx += sen_d * xyz.x();
                    df_dy += sen_d * xyz.y();
                    df_dz += sen_d * (xyz.z()-curvature_sign*curvature_radius*-1.f);
                }
                else
                    df_dz += bl_score[i] * dfs;
            }

            if (use_curvature) {
                cb_pos_sens (0, ri)  += df_dx;
                cb_pos_sens (1, ri)  += df_dy;
            }
            cb_pos_sens (2, ri)  += df_dz;
            for(int aa: range(n_restype)) 
                env_cov_sens(0, ri*n_restype+aa) += weights[aa] * df_dbl;
        }
    }

    virtual std::vector<float> get_param() const override {
        int np = n_restype*n_bl*n_node;
        vector<float> params(np, 0.f);
        for(int n=0;n<np;n++) {
            params[n] = coeff[n];
        }
        return params;
    }

    virtual void set_param(const std::vector<float>& new_params) override {
        int np = n_restype*n_bl*n_node;
        for(int n=0;n<np;n++) {
            coeff[n] = new_params[n];
        }
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        int np = n_restype*n_bl*n_node;
        vector<float> deriv(np, 0.f);

        VecArray cb_pos  = res_pos.output;
        VecArray env_cov = environment_coverage.output;
        VecArray surf    = surface.output;

        vector<float> bl_score1(n_bl-1);
        vector<float> bl_score2(n_bl-1);
        vector<float> d_bl_score1(n_bl-1);
        vector<float> d_bl_score2(n_bl-1);
        vector<float> bl_score(n_bl);
        vector<float> d_bl_score(n_bl);

        auto ccenter = load_vec<3>(center_of_curvature.output, 0);
        curvature_radius = ccenter.z()*-1.f*curvature_sign;

        for(int nr=0; nr<n_elem; ++nr) {
            int rt = res_type[nr];
            int ri = cb_index[nr];

            float bl = 0.0;
            for(int aa: range(n_restype)) 
                bl += weights[aa] * env_cov(0, ri*n_restype+aa);
            for (int i: range(n_bl-1)) {
                auto s1 = compact_sigmoid(bl-bl_nodes[i], bl_range_sharpness);
                auto s2 = compact_sigmoid(bl_nodes[i]-bl, bl_range_sharpness);
                bl_score1[i]   = s1.x();
                bl_score2[i]   = s2.x();
            }
            bl_score[0]       = bl_score1[0];
            for (int i: range(n_bl-2)) 
                bl_score[i+1] = bl_score2[i]*bl_score1[i+1];
            bl_score[n_bl-1]  = bl_score2[n_bl-2];

            float cb_z     = cb_pos(2, ri);
            if (use_curvature) {
                auto xyz = load_vec<3>(cb_pos, ri);
                float dist = mag( xyz - ccenter );
                cb_z = curvature_sign* ( dist - curvature_radius );
            }

            auto  sig_left = compact_sigmoid(cb_z-left_right_node, left_right_sharpness);
            float ll       = sig_left.x();
            float lr       = 1.f-sig_left.x();

            float surfv  = 0.0;

            auto  coord_left  = ( cb_z-zstart)*zscale;
            auto  coord_right = (-cb_z-zstart)*zscale;

            int starting_bin;
            float result[4];
            clamped_deBoor_coeff_deriv(&starting_bin, result, coord_left, n_node);
            for (int bi: range(n_bl)) {

                if (bi>0 and (cb_z > half_thickness or cb_z < -half_thickness)) 
                    surfv = 1.0;
                else
                    surfv = surf(0, ri);

                int shift = (rt*n_bl + bi)*n_node;
                for(int i: range(4)) deriv[shift+starting_bin+i] += result[i]*bl_score[bi]*ll*surfv;
            }

            clamped_deBoor_coeff_deriv(&starting_bin, result, coord_right, n_node);
            for (int bi: range(n_bl)) {

                if (bi>0 and (cb_z > half_thickness or cb_z < -half_thickness)) 
                    surfv = 1.0;
            else
                    surfv = surf(0, ri);

                int shift = (rt*n_bl + bi)*n_node + np;
                for(int i: range(4)) deriv[shift+starting_bin+i] += result[i]*bl_score[bi]*lr*surfv;
            }
        }
        return deriv;
    }
#endif
};

static RegisterNodeType<MembraneSurfCBPotential, 4> membrane_surf_cb_potential_node("cb_surf_membrane_potential");

struct MembraneSurfHBPotential : public PotentialNode
{
    int n_elem;        // number of residues to process
    int n_restype;

    vector<int> hb_index;
    vector<int> hb_type;
    vector<int> hb_res;

    CoordNode& protein_hbond;
    CoordNode& environment_coverage;
    CoordNode& surface;
    CoordNode& center_of_curvature;
    vector<float> weights;

    int   n_node;
    float zstart;
    float zscale;
    float half_thickness;
    float left_right_node;
    float left_right_sharpness;

    bool  use_curvature;
    float curvature_radius;
    float curvature_sign;

    vector<float> coeff;
    vector<float> coeff_inner;

    MembraneSurfHBPotential(hid_t grp, CoordNode& protein_hbond_, 
                               CoordNode& environment_coverage_,
                               CoordNode& surface_,
                               CoordNode& center_of_curvature_):
        PotentialNode(),
        n_elem               (get_dset_size(1, grp, "hb_index")[0]),
        n_restype            (get_dset_size(1, grp, "weights")[0]),
        hb_index             (n_elem),
        hb_type              (n_elem),
        hb_res               (n_elem),
        protein_hbond        (protein_hbond_),
        environment_coverage (environment_coverage_),
        surface              (surface_),
        center_of_curvature  (center_of_curvature_),
        weights              (n_restype),
        n_node               (get_dset_size(3, grp, "coeff")[2]),
        zstart               (read_attribute<float>(grp, ".", "z_start" )),
        zscale               (read_attribute<float>(grp, ".", "z_scale" )),
        half_thickness       (read_attribute<float>(grp, ".", "half_thickness" )),
        left_right_node      (read_attribute<float>(grp, ".", "left_right_node" )),
        left_right_sharpness (read_attribute<float>(grp, ".", "left_right_sharpness" )),
        use_curvature        (read_attribute<int>  (grp, ".", "use_curvature" )),
        curvature_radius     (read_attribute<float>(grp, ".", "curvature_radius" )),
        curvature_sign       (read_attribute<float>(grp, ".", "curvature_sign" )),
        coeff                (2*2*n_node),
        coeff_inner          (2*2*n_node)
    {
        check_size(grp, "hb_type", n_elem);
        check_size(grp, "hb_res",  n_elem);
        check_size(grp, "coeff",       2, 2, n_node);
        check_size(grp, "coeff_inner", 2, 2, n_node);

        traverse_dset<1, int> ( grp, "hb_index", [&](size_t n, int x) {hb_index[n] = x;} );
        traverse_dset<1, int> ( grp, "hb_type",  [&](size_t n, int x) {hb_type[n]  = x;} );
        traverse_dset<1, int> ( grp, "hb_res",   [&](size_t n, int x) {hb_res[n]   = x;} );

        traverse_dset<1, float> ( grp, "weights", [&](size_t n, float x) {weights[n] = x;} );

        traverse_dset<3, float>(grp, "coeff", [&](size_t n, size_t l, size_t c, float x) 
                                                   { coeff[(n*2+l)*n_node + c] = x;} );

        traverse_dset<3, float>(grp, "coeff_inner", [&](size_t n, size_t l, size_t c, float x) 
                                                   { coeff_inner[(n*2+l)*n_node + c] = x;} );

    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("hb_membrane_potential"));

        VecArray hb_pos  = protein_hbond.output;
        VecArray hb_sens = protein_hbond.sens;
        VecArray surf    = surface.output;

        VecArray env_cov      = environment_coverage.output;
        VecArray env_cov_sens = environment_coverage.sens;

        float pot     = 0.f;
        float comb_f  = 0.f;
        float comb_df = 0.f;

        auto ccenter = load_vec<3>(center_of_curvature.output, 0);
        curvature_radius = ccenter.z()*-1.f*curvature_sign;

        for(int nv=0; nv<n_elem; ++nv) {

            int   rt          = hb_type[nv];
            int   ri          = hb_index[nv];
            int   rr          = hb_res[nv];

            float hb_z        = hb_pos(2, ri);
            float hb_prob     = hb_pos(6, ri);  // probability that his virtual participates in any HBond

            float dist = 0.f;
            float r_dist = 1.f;
            auto xyz = load_vec<3>(hb_pos, ri);
            if (use_curvature) {
                dist = mag( xyz - ccenter );
                r_dist = 1.f/dist;
                hb_z = curvature_sign* ( dist - curvature_radius );
            }

            // burial level
            float bl = 0.0;
            for(int aa: range(n_restype)) 
                bl += weights[aa] * env_cov(0, rr*n_restype+aa);
            auto sig_bl = compact_sigmoid(bl-2.5, 2.0);

            float surfv = surf(0, rr);
            if (hb_z > half_thickness or hb_z < -half_thickness)
                surfv = 1.0;

            float outer   = surfv*sig_bl.x();
            float d_outer = surfv*sig_bl.y();
            float inner   = 1.f - outer;
            float d_inner = -d_outer;

            float uhb_prob      = 1.f-hb_prob;
            float sqr_uhb_prob  = sqr(uhb_prob);
            float sqr_hb_prob   = 1.f - sqr_uhb_prob;

            auto sig_left = compact_sigmoid(hb_z-left_right_node, left_right_sharpness);
            float ll  =     sig_left.x();
            float lr  = 1.f-sig_left.x();
            float dll =     sig_left.y();
            float dlr =    -sig_left.y();

            auto  coord_left  = ( hb_z-zstart)*zscale;
            auto  coord_right = (-hb_z-zstart)*zscale;

            float dZ      = 0.f;
            float dHB     = 0.f;
            float dBL     = 0.f;

            if (outer > 0.f) {
                auto uhb_l = clamped_deBoor_value_and_deriv(coeff.data() + (rt*2)*n_node, coord_left,  n_node);
                auto uhb_r = clamped_deBoor_value_and_deriv(coeff.data() + (rt*2)*n_node, coord_right, n_node);
                comb_f     = ll*uhb_l.x() + lr*uhb_r.x();
                comb_df    = ll*uhb_l.y()*zscale + dll*uhb_l.x() - lr*uhb_r.y()*zscale + dlr*uhb_r.x();
                pot       += sqr_uhb_prob *   outer * comb_f;
                dHB       -= 2.f*uhb_prob *   outer * comb_f;
                dBL       += sqr_uhb_prob * d_outer * comb_f;
                dZ        += sqr_uhb_prob *   outer * comb_df;

                auto hb_l = clamped_deBoor_value_and_deriv(coeff.data() + (rt*2+1)*n_node, coord_left,  n_node);
                auto hb_r = clamped_deBoor_value_and_deriv(coeff.data() + (rt*2+1)*n_node, coord_right, n_node);
                comb_f    = ll*hb_l.x() + lr*hb_r.x();
                comb_df   = ll*hb_l.y()*zscale + dll*hb_l.x() - lr*hb_r.y()*zscale + dlr*hb_r.x();
                pot      +=  sqr_hb_prob *   outer * comb_f;
                dHB      += 2.f*uhb_prob *   outer * comb_f;
                dBL      +=  sqr_hb_prob * d_outer * comb_f;
                dZ       +=  sqr_hb_prob *   outer * comb_df;
            }

            if (inner > 0.f) {
                auto uhb_l_i = clamped_deBoor_value_and_deriv(coeff_inner.data() + (rt*2)*n_node, coord_left,  n_node);
                auto uhb_r_i = clamped_deBoor_value_and_deriv(coeff_inner.data() + (rt*2)*n_node, coord_right, n_node);
                comb_f       = ll*uhb_l_i.x() + lr*uhb_r_i.x();
                comb_df      = (ll*uhb_l_i.y() - lr*uhb_r_i.y())*zscale + dll*uhb_l_i.x() + dlr*uhb_r_i.x();
                pot         += sqr_uhb_prob *   inner * comb_f;
                dHB         -= 2.f*uhb_prob *   inner * comb_f;
                dBL         += sqr_uhb_prob * d_inner * comb_f;
                dZ          += sqr_uhb_prob *   inner * comb_df;

                auto hb_l_i = clamped_deBoor_value_and_deriv(coeff_inner.data() + (rt*2+1)*n_node, coord_left,  n_node);
                auto hb_r_i = clamped_deBoor_value_and_deriv(coeff_inner.data() + (rt*2+1)*n_node, coord_right, n_node);
                comb_f      = ll*hb_l_i.x() + lr*hb_r_i.x();
                comb_df     = ll*hb_l_i.y()*zscale + dll*hb_l_i.x() - lr*hb_r_i.y()*zscale + dlr*hb_r_i.x();
                pot        +=  sqr_hb_prob *   inner * comb_f;
                dHB        += 2.f*uhb_prob *   inner * comb_f;
                dBL        +=  sqr_hb_prob * d_inner * comb_f;
                dZ         +=  sqr_hb_prob *   inner * comb_df;
            }

            if (use_curvature) {
                float sen_d = curvature_sign * dZ * r_dist; 
                hb_sens(0, ri) += sen_d * xyz.x();
                hb_sens(1, ri) += sen_d * xyz.y();
                hb_sens(2, ri) += sen_d * (xyz.z()-curvature_sign*curvature_radius*-1.f);
            }
            else
                hb_sens(2, ri) += dZ;

            hb_sens(6, ri) += dHB;

            for(int aa: range(n_restype)) 
                env_cov_sens(0, rr*n_restype+aa) += weights[aa] * dBL;

        }

        potential = pot;
    }

    virtual std::vector<float> get_param() const override {
        int np = 2*2*n_node;
        vector<float> params(np, 0.f);
        for(int n=0;n<np;n++) {
            params[n] = coeff[n];
        }
        return params;
    }

    virtual void set_param(const std::vector<float>& new_params) override {
        int np = 2*2*n_node;
        for(int n=0;n<np;n++) {
            coeff[n]  = new_params[n];
        }
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        int np = 2*2*n_node;
        vector<float> deriv(np, 0.f);
        VecArray hb_pos = protein_hbond.output;
        VecArray surf   = surface.output;

        auto ccenter = load_vec<3>(center_of_curvature.output, 0);
        curvature_radius = ccenter.z()*-1.f*curvature_sign;

        for(int nv=0; nv<n_elem; ++nv) {

            int   rt          = hb_type[nv];
            int   ri          = hb_index[nv];
            int   rr          = hb_res[nv];

            float hb_z        = hb_pos(2, ri);
            float dist = 0.f;
            if (use_curvature) {
                auto xyz = load_vec<3>(hb_pos, ri);
                dist = mag( xyz - ccenter );
                hb_z = curvature_sign* ( dist - curvature_radius );
            }

            float surfv       = surf(0, rr);
            if (hb_z > half_thickness or hb_z < -half_thickness)
                surfv = 1.0;

            auto  coord_left   = ( hb_z-zstart)*zscale;
            auto  coord_right  = (-hb_z-zstart)*zscale;
            auto  sig_left     = compact_sigmoid(hb_z-left_right_node, left_right_sharpness);
            float ll           = sig_left.x();
            float lr           = 1.f-sig_left.x();

            float hb_prob      = hb_pos(6, ri);  // probability that his virtual participates in any HBond
            float uhb_prob     = 1.f-hb_prob;
            float sqr_uhb_prob = sqr(uhb_prob);
            float sqr_hb_prob  = 1.f - sqr_uhb_prob;

            int starting_bin;
            float result[4];
            clamped_deBoor_coeff_deriv(&starting_bin, result, coord_left, n_node);
            for(int i: range(4)) deriv[(rt*2+0)*n_node+starting_bin+i] += result[i]*ll*sqr_uhb_prob*surfv;
            for(int i: range(4)) deriv[(rt*2+1)*n_node+starting_bin+i] += result[i]*ll*sqr_hb_prob*surfv;

            clamped_deBoor_coeff_deriv(&starting_bin, result, coord_right, n_node);
            for(int i: range(4)) deriv[(rt*2+0)*n_node+starting_bin+i] += result[i]*lr*sqr_uhb_prob*surfv;
            for(int i: range(4)) deriv[(rt*2+1)*n_node+starting_bin+i] += result[i]*lr*sqr_hb_prob*surfv;
        }
        return deriv;
    }
#endif
};

static RegisterNodeType<MembraneSurfHBPotential, 4> membrane_surf_hb_potential_node("hb_surf_membrane_potential");

struct MembraneLateralPotential : public PotentialNode
{
    int n_residue;        // number of residues to process
    int n_restype;        // number of residues to process
    CoordNode& pos;     // CB atom
    CoordNode& BL;      // burial level
    CoordNode& surface; // surface score
    vector<int> cb_index;
    vector<int> cb_restype;
    vector<float> weights;
    vector<float> cov_midpoint;
    vector<float> cov_sharpness;
    float min_z, max_z, inv_dz;
    int n_knot;
    float half_thickness;
    vector<float> params;

    MembraneLateralPotential(hid_t grp, CoordNode& pos_,
                                        CoordNode& environment_coverage_,
                                        CoordNode& surface_):
        PotentialNode(),
        n_residue (get_dset_size(1, grp,  "cb_index")[0]),
        n_restype (N_AMINO_ACID_TYPE),
        pos(pos_),
        BL(environment_coverage_),
        surface(surface_),
        cb_index(n_residue),
        cb_restype(n_residue),
        weights(N_AMINO_ACID_TYPE),
        cov_midpoint(n_restype),
        cov_sharpness(n_restype),
        min_z(read_attribute<float>(grp, ".","min_z")),
        max_z(read_attribute<float>(grp, ".","max_z")),
        inv_dz(read_attribute<float>(grp, ".","inv_dz")),
        n_knot(read_attribute<int>(grp, ".", "n_knot")),
        half_thickness(read_attribute<float>(grp, ".", "half_thickness" )),
        params(n_knot)
    {
        check_elem_width_lower_bound(pos, 3);
        check_elem_width_lower_bound(BL, 1);
        check_elem_width_lower_bound(surface, 3);

        check_size(grp,      "restypes", n_residue);
        check_size(grp,  "cov_midpoint", n_restype);
        check_size(grp, "cov_sharpness", n_restype);
        check_size(grp,       "weights", N_AMINO_ACID_TYPE);
        check_size(grp,        "params", n_knot);
        traverse_dset<1,float>(grp,      "cb_index", [&](size_t n,   float x) {cb_index[n] = x;});
        traverse_dset<1,  int>(grp,      "restypes", [&](size_t nr,   int rt) {cb_restype[nr]    = rt;});
        traverse_dset<1,float>(grp,       "weights", [&](size_t rt, float bw) {weights[rt]       = bw;});
        traverse_dset<1,float>(grp,  "cov_midpoint", [&](size_t rt, float bc) {cov_midpoint[rt]  = bc;});
        traverse_dset<1,float>(grp, "cov_sharpness", [&](size_t rt, float bs) {cov_sharpness[rt] = bs;});
        traverse_dset<1,float>(grp,        "params", [&](size_t rt, float pv) {params[rt]        = pv;});
    }

    virtual void compute_value(ComputeMode mode) {
        Timer timer(string("membrane_lateral_potential"));

        VecArray surf    = surface.output;
        VecArray bl      = BL.output;
        //VecArray bl_sens = BL.sens; // FIXME calcu the sens of bl
        VecArray cb_pos  = pos.output;
        VecArray cb_sens = pos.sens;

        potential = 0.f;

        vector<float> weighted_bl(n_residue, 0.);
        for(int nr=0; nr<n_residue; ++nr)
            for(int aa: range(N_AMINO_ACID_TYPE))
                weighted_bl[nr] += weights[aa] * bl(0, nr*N_AMINO_ACID_TYPE+aa);

        float* pp = &params[0];
        for(int i=0; i<n_residue; ++i) {
            int nr = cb_index[i];

            float surfv = surf(0, nr);
            if (cb_pos(2, nr) > half_thickness or cb_pos(2, nr) < -half_thickness)
                surfv = 1.0;

            int rtype = cb_restype[nr];
            auto cover_sig = compact_sigmoid(weighted_bl[nr]-cov_midpoint[rtype], cov_sharpness[rtype]);

            if(surfv<=0 or cover_sig.x() ==0 ) continue;

            float           coord = (cb_pos(2, nr)-min_z)*inv_dz;
            auto   lateral_spline = clamped_deBoor_value_and_deriv(pp, coord, n_knot);
            float force_magnitude = surfv*lateral_spline.x();
            cb_sens (0, nr)      -= force_magnitude*cover_sig.x()*surf(1, nr);
            cb_sens (1, nr)      -= force_magnitude*cover_sig.x()*surf(2, nr);
        }
    }
};

static RegisterNodeType<MembraneLateralPotential, 3> membrane_lateral_potential_node("membranelateral_potential");

