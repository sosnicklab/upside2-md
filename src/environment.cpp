#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "spline.h"
#include "interaction_graph.h"
#include <algorithm>

using namespace std;
using namespace h5;

namespace {
    template <bool is_symmetric>
    struct SphereEnvironmentCoverageInteraction {
        // parameters are r0,r_sharpness, dot0,dot_sharpness

        constexpr static bool  symmetric = is_symmetric;
        constexpr static int   n_param=2, n_dim1=3, n_dim2=3, simd_width=1;

        static float cutoff(const float* p) {
            return p[0] + compact_sigmoid_cutoff(p[1]);
        }

        static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
            auto sequence_exclude = Int4(1);  // exclude i,i, i,i+1, and i,i+2
            return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
        }

        static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4],
                const Vec<n_dim1,Float4> &cb_pos, const Vec<n_dim2,Float4> &sc_pos) {

            Float4 one(1.f);
            auto displace = extract<0,3>(sc_pos)-extract<0,3>(cb_pos);

            auto dist2            = mag2(displace);
            auto inv_dist         = rsqrt(dist2);
            auto dist             = dist2*inv_dist;
            auto displace_unitvec = inv_dist*displace;

            // read parameters then transpose
            Float4 r0(p[0]);
            Float4 s0(p[1]);
            Float4 null0(0.f);
            Float4 null1(0.f);
            transpose4(r0, s0, null0, null1);

            auto radial_sig = compact_sigmoid(dist-r0, s0);
            auto score = radial_sig.x();

            // now we compute derivatives (minus sign is from the derivative of angular_sig)
            auto d_displace = radial_sig.y() * displace_unitvec;
            store<0,3>(d1, -d_displace);
            store<0,3>(d2,  d_displace);
            return score;
        }
        static void param_deriv(Vec<n_param> &d_param, const float* p,
                const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
            d_param = make_zero<n_param>();   // not implemented currently
        }
        static bool is_compatible(const float* p1, const float* p2) {return true;};
    };
    
    template <bool is_symmetric>
    struct WeightedSphereEnvironmentCoverageInteraction {
        // parameters are r0,r_sharpness, dot0,dot_sharpness

        constexpr static bool  symmetric = is_symmetric;
        constexpr static int   n_param=2, n_dim1=3, n_dim2=4, simd_width=1;

        static float cutoff(const float* p) {
            return p[0] + compact_sigmoid_cutoff(p[1]);
        }

        static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
            auto sequence_exclude = Int4(2);  // exclude i,i, i,i+1, and i,i+2
            return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
        }

        static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4],
                const Vec<n_dim1,Float4> &cb_pos, const Vec<n_dim2,Float4> &sc_pos) {

            Float4 one(1.f);
            auto displace = extract<0,3>(sc_pos)-extract<0,3>(cb_pos);
            auto prob     = sc_pos[3];

            auto dist2            = mag2(displace);
            auto inv_dist         = rsqrt(dist2);
            auto dist             = dist2*inv_dist;
            auto displace_unitvec = inv_dist*displace;

            // read parameters then transpose
            Float4 r0(p[0]);
            Float4 s0(p[1]);
            Float4 null0(0.f);
            Float4 null1(0.f);
            transpose4(r0, s0, null0, null1);

            auto radial_sig = compact_sigmoid(dist-r0, s0);
            auto score = radial_sig.x();

            // now we compute derivatives (minus sign is from the derivative of angular_sig)
            auto d_displace = prob * radial_sig.y() * displace_unitvec;
            store<0,3>(d1, -d_displace);
            store<0,3>(d2,  d_displace);
            d2[3] = score;
            return score*prob;
        }
        static void param_deriv(Vec<n_param> &d_param, const float* p,
                const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
            d_param = make_zero<n_param>();   // not implemented currently
        }
        static bool is_compatible(const float* p1, const float* p2) {return true;};
    };

    struct EnvironmentCoverageInteraction {
        // parameters are r0,r_sharpness, dot0,dot_sharpness

        constexpr static bool  symmetric = false;
        constexpr static int   n_param=4, n_dim1=6, n_dim2=4, simd_width=1;

        static float cutoff(const float* p) {
            return p[0] + compact_sigmoid_cutoff(p[1]);
        }

        static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
            auto sequence_exclude = Int4(2);  // exclude i,i, i,i+1, and i,i+2
            return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
        }

        static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4],
                const Vec<n_dim1,Float4> &cb_pos, const Vec<n_dim2,Float4> &sc_pos) {
            Float4 one(1.f);
            auto displace = extract<0,3>(sc_pos)-extract<0,3>(cb_pos);
            auto rvec1    = extract<3,6>(cb_pos);
            auto prob     = sc_pos[3];

            auto dist2 = mag2(displace);
            auto inv_dist = rsqrt(dist2);
            auto dist = dist2*inv_dist;
            auto displace_unitvec = inv_dist*displace;

            // read parameters then transpose
            Float4 r0(p[0]);
            Float4 r_sharpness(p[1]);
            Float4 dot0(p[2]);
            Float4 dot_sharpness(p[3]);
            transpose4(r0,r_sharpness,dot0,dot_sharpness);

            auto dp = dot(displace_unitvec,rvec1);
            auto radial_sig  = compact_sigmoid(dist-r0, r_sharpness);
            auto angular_sig = compact_sigmoid(dot0-dp, dot_sharpness);

            // now we compute derivatives (minus sign is from the derivative of angular_sig)
            auto d_displace = prob*(radial_sig.y()*angular_sig.x() * displace_unitvec -
                                    radial_sig.x()*angular_sig.y()* inv_dist*(rvec1 - dp*displace_unitvec));

            store<3,6>(d1, -prob*radial_sig.x()*angular_sig.y()*displace_unitvec);
            store<0,3>(d1, -d_displace);
            store<0,3>(d2,  d_displace);
            auto score = radial_sig.x() * angular_sig.x();
            d2[3] = score;
            return prob * score;
        }

        static void param_deriv(Vec<n_param> &d_param, const float* p,
                const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
            d_param = make_zero<n_param>();   // not implemented currently
        }

        static bool is_compatible(const float* p1, const float* p2) {return true;};
    };

    struct HbondEnvironmentCoverageInteraction {
        // parameters are r0,r_sharpness, dot0,dot_sharpness
    
        constexpr static bool  symmetric = false;
        constexpr static int   n_param=4, n_dim1=6, n_dim2=6, simd_width=1;
    
        static float cutoff(const float* p) {
    	return p[0] + compact_sigmoid_cutoff(p[1]);
        }
    
        static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
            auto sequence_exclude = Int4(2);  // exclude i,i, i,i+1, and i,i+2
            return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
            //return Int4() == Int4();  // No exclusions (all true)
        }
    
        static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4],
    	    const Vec<n_dim1,Float4> &x1, const Vec<n_dim2,Float4> &x2) {

    	Float4 one(1.f);
    	auto displace = extract<0,3>(x2)-extract<0,3>(x1);
    	auto rvec1 = extract<3,6>(x1);
    
    	auto dist2 = mag2(displace);
    	auto inv_dist = rsqrt(dist2);
    	auto dist = dist2*inv_dist;
    	auto displace_unitvec = inv_dist*displace;
    
    	// read parameters then transpose
    	Float4 r0(p[0]);
    	Float4 r_sharpness(p[1]);
    	Float4 dot0(p[2]);
    	Float4 dot_sharpness(p[3]);
    	transpose4(r0,r_sharpness,dot0,dot_sharpness);
    
    	auto dp = dot(displace_unitvec,rvec1);
    	auto radial_sig  = compact_sigmoid(dist-r0, r_sharpness);
    	auto angular_sig = compact_sigmoid(dot0-dp, dot_sharpness);
    
    	// now we compute derivatives (minus sign is from the derivative of angular_sig)
    	auto d_displace = radial_sig.y()*angular_sig.x() * displace_unitvec -
    	       		radial_sig.x()*angular_sig.y()* inv_dist*(rvec1 - dp*displace_unitvec);
    
    	store<3,6>(d1, -radial_sig.x()*angular_sig.y()*displace_unitvec);
    	store<0,3>(d1, -d_displace);
    	store<0,3>(d2,  d_displace);
    	auto score = radial_sig.x() * angular_sig.x();
    	return score;
        }
    
        static void param_deriv(Vec<n_param> &d_param, const float* p,
    	    const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
    	d_param = make_zero<n_param>();   // not implemented currently
        }
    
        static bool is_compatible(const float* p1, const float* p2) {return true;};
    };

    struct HbondEnvironmentCoverageInteraction2 {
        // parameters are r0,r_sharpness, dot0,dot_sharpness
    
        constexpr static bool  symmetric = false;
        constexpr static int   n_param=4, n_dim1=6, n_dim2=3, simd_width=1;
    
        static float cutoff(const float* p) {
    	return p[0] + compact_sigmoid_cutoff(p[1]);
        }
    
        static Int4 acceptable_id_pair(const Int4& id1, const Int4& id2) {
            auto sequence_exclude = Int4(2);  // exclude i,i, i,i+1, and i,i+2
            return (sequence_exclude < id1-id2) | (sequence_exclude < id2-id1);
            //return Int4() == Int4();  // No exclusions (all true)
        }
    
        static Float4 compute_edge(Vec<n_dim1,Float4> &d1, Vec<n_dim2,Float4> &d2, const float* p[4],
    	    const Vec<n_dim1,Float4> &x1, const Vec<n_dim2,Float4> &x2) {

    	Float4 one(1.f);
    	auto displace = extract<0,3>(x2)-extract<0,3>(x1);
    	auto rvec1 = extract<3,6>(x1);
    
    	auto dist2 = mag2(displace);
    	auto inv_dist = rsqrt(dist2);
    	auto dist = dist2*inv_dist;
    	auto displace_unitvec = inv_dist*displace;
    
    	// read parameters then transpose
    	Float4 r0(p[0]);
    	Float4 r_sharpness(p[1]);
    	Float4 dot0(p[2]);
    	Float4 dot_sharpness(p[3]);
    	transpose4(r0,r_sharpness,dot0,dot_sharpness);
    
    	auto dp = dot(displace_unitvec,rvec1);
    	auto radial_sig  = compact_sigmoid(dist-r0, r_sharpness);
    	auto angular_sig = compact_sigmoid(dot0-dp, dot_sharpness);
    
    	// now we compute derivatives (minus sign is from the derivative of angular_sig)
    	auto d_displace = radial_sig.y()*angular_sig.x() * displace_unitvec -
    	       		radial_sig.x()*angular_sig.y()* inv_dist*(rvec1 - dp*displace_unitvec);
    
    	store<3,6>(d1, -radial_sig.x()*angular_sig.y()*displace_unitvec);
    	store<0,3>(d1, -d_displace);
    	store<0,3>(d2,  d_displace);
    	auto score = radial_sig.x() * angular_sig.x();
    	return score;
        }
    
        static void param_deriv(Vec<n_param> &d_param, const float* p,
    	    const Vec<n_dim1> &hb_pos, const Vec<n_dim2> &sc_pos) {
    	d_param = make_zero<n_param>();   // not implemented currently
        }
    
        static bool is_compatible(const float* p1, const float* p2) {return true;};
    };

struct EnvironmentCoverage : public CoordNode {

    int n_aa;
    int n_res;
    std::vector<int> aa_types;
    InteractionGraph<EnvironmentCoverageInteraction> igraph;

    EnvironmentCoverage(hid_t grp, CoordNode& cb_pos_, CoordNode& weighted_sidechains_):
        CoordNode(get_dset_size(1,grp,"index1")[0]*read_attribute<int>(grp,".","num_aa_types") , 1),
        n_aa(read_attribute<int>(grp,".","num_aa_types")), 
        n_res(get_dset_size(1,grp,"index1")[0]),
        aa_types(n_res), // FIXME: aa_types should be initialized from (non bead rotamer redundant) n_res of SC beads. index2 is multiplied by rotamers, so can't use that directly
        igraph(grp, &cb_pos_, &weighted_sidechains_)
    {
        check_size(grp, "aa_types", n_res);
        traverse_dset<1,int>  (grp,"aa_types",[&](size_t ne, int x){aa_types[ne]=x;});

        if(logging(LOG_EXTENSIVE)) {
            default_logger->add_logger<float>("environment_coverage", {n_elem}, [&](float* buffer) {
                    for(int ne: range(n_elem))
                            buffer[ne] = output(0,ne);});
        }

        for(int i=0;i<n_aa*n_res; ++i)
            output(0, i) = 0.0;
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("environment_coverage"));

        igraph.compute_edges();
        fill(output, 0.f);

        // accumulate for each cb
        for(int ne=0; ne<igraph.n_edge; ++ne) {
            int indices1 = igraph.edge_indices1[ne];
            int id2      = igraph.edge_id2[ne];
            int type2    = aa_types[id2];
            output(0, indices1*n_aa+type2) += igraph.edge_value[ne];
        }
    }

    virtual void propagate_deriv() override {
        Timer timer(string("d_environment_coverage"));

        for(int ne: range(igraph.n_edge)) {
            int indices1 = igraph.edge_indices1[ne];
            int id2      = igraph.edge_id2[ne];
            int type2    = aa_types[id2];
            igraph.edge_sensitivity[ne] = sens(0, indices1*n_aa+type2);
        }
        igraph.propagate_derivatives();
    }

    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}
};
static RegisterNodeType<EnvironmentCoverage,2> environment_coverage_node("environment_coverage");

struct HbondEnvironmentCoverage : public CoordNode {
    InteractionGraph<HbondEnvironmentCoverageInteraction> igraph;

    HbondEnvironmentCoverage(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(get_dset_size(1,grp,"index1")[0], 1),
        igraph(grp, &pos1_, &pos2_) {
        if(logging(LOG_EXTENSIVE)) {
            default_logger->add_logger<float>("hbond_environment_coverage", {n_elem}, [&](float* buffer) {
                for(int ne: range(n_elem))
                    buffer[ne] = output(0,ne);});
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("hbond_environment_coverage"));
        igraph.compute_edges();
        fill(output, 0.f);
        for(int ne=0; ne<igraph.n_edge; ++ne)  // accumulate for each cb
            output(0, igraph.edge_indices1[ne]) += igraph.edge_value[ne];
    }

    virtual void propagate_deriv() override {
        Timer timer(string("d_hbond_environment_coverage"));

        for(int ne: range(igraph.n_edge))
            igraph.edge_sensitivity[ne] = sens(0,igraph.edge_indices1[ne]);
        igraph.propagate_derivatives();
    }

    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}
};
static RegisterNodeType<HbondEnvironmentCoverage,2> hbond_environment_coverage_node("hb_environment_coverage");

struct HbondBackBoneCoverage : public CoordNode {
    InteractionGraph<HbondEnvironmentCoverageInteraction2> igraph;

    HbondBackBoneCoverage(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(get_dset_size(1,grp,"index1")[0], 1),
        igraph(grp, &pos1_, &pos2_) {
        if(logging(LOG_EXTENSIVE)) {
            default_logger->add_logger<float>("hbond_backbone_coverage", {n_elem}, [&](float* buffer) {
                for(int ne: range(n_elem))
                    buffer[ne] = output(0,ne);});
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("hbond_backbone_coverage"));
        igraph.compute_edges();
        fill(output, 0.f);
        for(int ne=0; ne<igraph.n_edge; ++ne)  // accumulate for each cb
            output(0, igraph.edge_indices1[ne]) += igraph.edge_value[ne];
    }

    virtual void propagate_deriv() override {
        Timer timer(string("d_hbond_backbone_coverage"));

        for(int ne: range(igraph.n_edge))
            igraph.edge_sensitivity[ne] = sens(0,igraph.edge_indices1[ne]);
        igraph.propagate_derivatives();
    }

    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}
};
static RegisterNodeType<HbondBackBoneCoverage,2> hbond_backbone_coverage_node("hbbb_coverage");

struct SphereEnvironmentCoverage : public CoordNode {
    InteractionGraph<SphereEnvironmentCoverageInteraction<true>> igraph;

    SphereEnvironmentCoverage(hid_t grp, CoordNode& pos1_ ):
        CoordNode(get_dset_size(1,grp,"index")[0], 1),
        igraph(grp, &pos1_) {
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("sphere_environment_coverage"));
        igraph.compute_edges();
        fill(output, 0.f);
        for(int ne=0; ne<igraph.n_edge; ++ne)  // accumulate for each cb
            output(0, igraph.edge_indices1[ne]) += igraph.edge_value[ne];
    }

    virtual void propagate_deriv() override {
        Timer timer(string("d_hbond_environment_coverage"));

        for(int ne: range(igraph.n_edge))
            igraph.edge_sensitivity[ne] = sens(0,igraph.edge_indices1[ne]);
        igraph.propagate_derivatives();
    }

    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}
};
static RegisterNodeType<SphereEnvironmentCoverage,1> sphere_environment_coverage_node("SphereCoverage");

struct WeightedSphereEnvironmentCoverage : public CoordNode {
    InteractionGraph<WeightedSphereEnvironmentCoverageInteraction<false>> igraph;

    WeightedSphereEnvironmentCoverage(hid_t grp, CoordNode& pos1_):
        CoordNode(get_dset_size(1,grp,"index")[0], 1),
        igraph(grp, &pos1_, &pos1_) {
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("sphere_environment_coverage"));
        igraph.compute_edges();
        fill(output, 0.f);
        for(int ne=0; ne<igraph.n_edge; ++ne)  // accumulate for each cb
            output(0, igraph.edge_indices1[ne]) += igraph.edge_value[ne];
    }

    virtual void propagate_deriv() override {
        Timer timer(string("d_hbond_environment_coverage"));

        for(int ne: range(igraph.n_edge))
            igraph.edge_sensitivity[ne] = sens(0,igraph.edge_indices1[ne]);
        igraph.propagate_derivatives();
    }

    virtual std::vector<float> get_param() const override {return igraph.get_param();}
#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {return igraph.get_param_deriv();}
#endif
    virtual void set_param(const std::vector<float>& new_param) override {igraph.set_param(new_param);}
};
static RegisterNodeType<WeightedSphereEnvironmentCoverage,1> weighted_sphere_environment_coverage_node("WeightedSphereCoverage");


struct WeightedPos : public CoordNode {
    CoordNode &pos;
    CoordNode &energy;

    struct Param {
        index_t index_pos;
        index_t index_weight;
    };
    vector<Param> params;

    WeightedPos(hid_t grp, CoordNode& pos_, CoordNode& energy_):
        CoordNode(get_dset_size(1,grp,"index_pos")[0], 4),
        pos(pos_),
        energy(energy_),
        params(n_elem)
    {
        traverse_dset<1,int>(grp,"index_pos"   ,[&](size_t ne, int x){params[ne].index_pos   =x;});
        traverse_dset<1,int>(grp,"index_weight",[&](size_t ne, int x){params[ne].index_weight=x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("weighted_pos");

        for(int ne=0; ne<n_elem; ++ne) {
            auto p = params[ne];
            output(0,ne) = pos.output(0,p.index_pos);
            output(1,ne) = pos.output(1,p.index_pos);
            output(2,ne) = pos.output(2,p.index_pos);
            output(3,ne) = expf(-energy.output(0,p.index_weight));
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_weighted_pos");

        for(int ne=0; ne<n_elem; ++ne) {
            auto p = params[ne];
            pos.sens(0,p.index_pos) += sens(0,ne);
            pos.sens(1,p.index_pos) += sens(1,ne);
            pos.sens(2,p.index_pos) += sens(2,ne);
            energy.sens(0,p.index_weight) -= output(3,ne)*sens(3,ne); // exponential derivative
        }
    }
};
static RegisterNodeType<WeightedPos,2> weighted_pos_node("weighted_pos");

struct UniformTransform : public CoordNode {
    CoordNode& input;
    int n_coeff;
    float spline_offset;
    float spline_inv_dx;
    unique_ptr<float[]> bspline_coeff;
    unique_ptr<float[]> jac;


    UniformTransform(hid_t grp, CoordNode& input_):
        CoordNode(input_.n_elem, 1),
        input(input_),
        n_coeff(get_dset_size(1,grp,"bspline_coeff")[0]),
        spline_offset(read_attribute<float>(grp,"bspline_coeff","spline_offset")),
        spline_inv_dx(read_attribute<float>(grp,"bspline_coeff","spline_inv_dx")),
        bspline_coeff(new_aligned<float>(n_coeff,4)),
        jac          (new_aligned<float>(input.n_elem,4))
    {
        check_elem_width(input,1); // this restriction could be lifted
        fill_n(bspline_coeff.get(), round_up(n_coeff,4), 0.f);
        traverse_dset<1,float>(grp,"bspline_coeff",[&](size_t ne, float x){bspline_coeff[ne]=x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("uniform_transform");
        for(int ne=0; ne<n_elem; ++ne) {
            auto coord = (input.output(0,ne)-spline_offset)*spline_inv_dx;
            auto v = clamped_deBoor_value_and_deriv(bspline_coeff.get(), coord, n_coeff);
            output(0,ne) = v[0];
            jac[ne]      = v[1]*spline_inv_dx;
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_uniform_transform");
        for(int ne=0; ne<n_elem; ++ne)
            input.sens(0,ne) += jac[ne]*sens(0,ne);
    }

    virtual std::vector<float> get_param() const override{
        vector<float> ret(2+n_coeff);
        ret[0] = spline_offset;
        ret[1] = spline_inv_dx;
        copy_n(bspline_coeff.get(), n_coeff, &ret[2]);
        return ret;
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        vector<float> ret(2+n_coeff, 0.f);
        int starting_bin;
        float d[4];
        for(int ne=0; ne<n_elem; ++ne) {
            auto coord = (input.output(0,ne)-spline_offset)*spline_inv_dx;
            auto v = clamped_deBoor_value_and_deriv(bspline_coeff.get(), coord, n_coeff);
            clamped_deBoor_coeff_deriv(&starting_bin, d, coord, n_coeff);

            ret[0] += v[1];                                    // derivative of offset
            ret[1] += v[1]*(input.output(0,ne)-spline_offset); // derivative of inv_dx
            for(int i: range(4)) ret[2+starting_bin+i] += d[i];
        }
        return ret;
    }
#endif

    virtual void set_param(const std::vector<float>& new_param) override {
        if(new_param.size() < size_t(2+4)) throw string("too small of size for spline");
        if(int(new_param.size())-2 != n_coeff) {
            n_coeff = int(new_param.size())-2;
            bspline_coeff = new_aligned<float>(n_coeff,4);
            fill_n(bspline_coeff.get(), round_up(n_coeff,4), 0.f);
        }
        spline_offset = new_param[0];
        spline_inv_dx = new_param[1];
        copy_n(begin(new_param)+2, n_coeff, bspline_coeff.get());
    }
};
static RegisterNodeType<UniformTransform,1> uniform_transform_node("uniform_transform");

struct LinearCoupling : public PotentialNode {
    CoordNode& input;
    vector<float> couplings;      // length n_restype
    vector<int>   coupling_types; // length input n_elem
    CoordNode* inactivation;  // 0 to 1 variable to inactivate energy
    int inactivation_dim;

    LinearCoupling(hid_t grp, CoordNode& input_, CoordNode& inactivation_):
        LinearCoupling(grp, &input_, &inactivation_) {}

    LinearCoupling(hid_t grp, CoordNode& input_):
        LinearCoupling(grp, &input_, nullptr) {}

    LinearCoupling(hid_t grp, CoordNode* input_, CoordNode* inactivation_):
        PotentialNode(),
        input(*input_),
        couplings(get_dset_size(1,grp,"couplings")[0]),
        coupling_types(input.n_elem),
        inactivation(inactivation_),
        inactivation_dim(inactivation ? read_attribute<int>(grp, ".", "inactivation_dim") : 0)
    {
        check_elem_width(input, 1);  // could be generalized

        if(inactivation) {
            if(input.n_elem != inactivation->n_elem)
                throw string("Inactivation size must match input size");
            check_elem_width_lower_bound(*inactivation, inactivation_dim+1);
        }

        check_size(grp, "coupling_types", input.n_elem);
        traverse_dset<1,float>(grp,"couplings",[&](size_t nt, float x){couplings[nt]=x;});
        traverse_dset<1,int>(grp,"coupling_types",[&](size_t ne, int x){coupling_types[ne]=x;});
        for(int i: coupling_types) if(i<0 || i>=int(couplings.size())) throw string("invalid coupling type");

        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>(
                    (inactivation ? "linear_coupling_with_inactivation" : "linear_coupling_uniform"),
                    {input.n_elem}, [&](float* buffer) {
                    for(int ne: range(input.n_elem)) {
                        float c = couplings[coupling_types[ne]];
                        buffer[ne] = c*input.output(0,ne);
                    }});
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("linear_coupling");
        int n_elem = input.n_elem;
        float pot = 0.f;
        for(int ne=0; ne<n_elem; ++ne) {
            float c = couplings[coupling_types[ne]];
            float act = inactivation ? sqr(1.f-inactivation->output(inactivation_dim,ne)) : 1.f;
            float val = input.output(0,ne);
            pot += c * val * act;
            input.sens(0,ne) += c*act;
            if(inactivation) inactivation->sens(inactivation_dim,ne) -= c*val;
        }
        potential = pot;
    }

    virtual std::vector<float> get_param() const override {
        return couplings;
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        vector<float> deriv(couplings.size(), 0.f);

        int n_elem = input.n_elem;
        for(int ne=0; ne<n_elem; ++ne) {
            float act = inactivation ? 1.f - inactivation->output(inactivation_dim,ne) : 1.f;
            deriv[coupling_types[ne]] += input.output(0,ne) * act;
        }
        return deriv;
    }
#endif

    virtual void set_param(const std::vector<float>& new_param) override {
        if(new_param.size() != couplings.size())
            throw string("attempting to change size of couplings vector on set_param");
        copy(begin(new_param),end(new_param), begin(couplings));
    }
};
static RegisterNodeType<LinearCoupling,1> linear_coupling_node1("linear_coupling_uniform");
static RegisterNodeType<LinearCoupling,2> linear_coupling_node2("linear_coupling_with_inactivation");


struct NonlinearCoupling : public PotentialNode {
    CoordNode& input;
    int n_restype, n_coeff, n_res;
    float spline_offset, spline_inv_dx;
    vector<float> coeff;      // length n_restype*n_coeff
    vector<int>   coupling_types; // length input n_elem
    vector<float> weights;
    vector<float> wnumber;
    int num_independent_weight;

    NonlinearCoupling(hid_t grp, CoordNode& input_):
        PotentialNode(),
        input(input_),
        n_restype(get_dset_size(2,grp,"coeff")[0]),
        n_coeff  (get_dset_size(2,grp,"coeff")[1]),
        n_res(input.n_elem/n_restype),
        spline_offset(read_attribute<float>(grp,"coeff","spline_offset")),
        spline_inv_dx(read_attribute<float>(grp,"coeff","spline_inv_dx")),
        coeff(n_restype*n_coeff),
        coupling_types(n_res),
        weights(n_restype*n_restype),
        wnumber(n_res),
        num_independent_weight(read_attribute<int>(grp, ".", "number_independent_weights"))
    {
        check_elem_width(input, 1);  // could be generalized

        check_size(grp, "coupling_types", n_res);
        check_size(grp, "weights", n_restype*n_restype);

        traverse_dset<2,float>(grp,"coeff",[&](size_t nt, size_t nc, float x){coeff[nt*n_coeff+nc]=x;});
        traverse_dset<1,int>  (grp,"coupling_types",[&](size_t ne, int x){coupling_types[ne]=x;});
        traverse_dset<1,float>(grp,"weights",[&](size_t nw, float x){weights[nw]=x;});

        for(int i: coupling_types) if(i<0 || i>=n_restype) throw string("invalid coupling type");

        if (num_independent_weight != 1 and num_independent_weight != 20 and num_independent_weight != 400) 
            throw string("the number of independent weights should be 1, 20 or 400");

        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>("nonlinear_coupling", {n_res}, [&](float* buffer) {
                for(int nr: range(n_res)) {
                    wnumber[nr] = 0.0;
                    int ctype = coupling_types[nr];
                    for(int aa: range(n_restype)) 
                        wnumber[nr] += weights[ctype*n_restype+aa] * input.output(0, nr*n_restype+aa);
                    auto coord = (wnumber[nr]-spline_offset)*spline_inv_dx;
                    buffer[nr] = clamped_deBoor_value_and_deriv(
                            coeff.data() + ctype*n_coeff, coord, n_coeff)[0];
                }});
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("nonlinear_coupling");

        float pot = 0.f;
        for(int nr: range(n_res)) {
            wnumber[nr] = 0.0;
            int ctype = coupling_types[nr];
            for(int aa: range(n_restype)) 
                wnumber[nr] += weights[ctype*n_restype+aa] * input.output(0, nr*n_restype+aa);

            auto coord = (wnumber[nr]-spline_offset)*spline_inv_dx;
            auto v = clamped_deBoor_value_and_deriv(coeff.data() + ctype*n_coeff, coord, n_coeff);

            pot += v[0];
            for(int aa: range(n_restype)) 
                input.sens(0, nr*n_restype+aa) = weights[ctype*n_restype+aa] * spline_inv_dx * v[1];
        }

        potential = pot;
    }

    virtual std::vector<float> get_param() const override {
        int csize = coeff.size();
        int wsize = n_restype*n_restype;
        vector<float> params(csize+wsize, 0.f);

        for(int ct: range(csize)) 
            params[ct] = coeff[ct];

        for(int ct: range(wsize)) 
            params[ct+csize] = weights[ct];

        return params;
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {

        int csize = coeff.size();
        int wsize = n_restype*n_restype;

        vector<float> deriv(csize+wsize, 0.f);
        vector<float> sens(n_res, 0.f);

        for(int nr: range(n_res)) {
            int starting_bin;
            float result[4];
            wnumber[nr] = 0.0;
            int ctype = coupling_types[nr];

            for(int aa: range(n_restype)) 
                wnumber[nr] += weights[ctype*n_restype+aa] * input.output(0, nr*n_restype+aa);

            auto coord = (wnumber[nr]-spline_offset)*spline_inv_dx;

            clamped_deBoor_coeff_deriv(&starting_bin, result, coord, n_coeff);
            for(int i: range(4)) deriv[ctype*n_coeff+starting_bin+i] += result[i];

            if ( num_independent_weight > 1) {
                auto v = clamped_deBoor_value_and_deriv(coeff.data() + ctype*n_coeff, coord, n_coeff);
                sens[nr] = spline_inv_dx * v[1];
            }
        }

        if ( num_independent_weight == 1 )
            return deriv;

        for(int nr: range(n_res)) {
            int ctype = coupling_types[nr];
            for(int aa: range(n_restype)) {
                auto deriv_aa = sens[nr] * input.output(0, nr*n_restype+aa);
                if ( num_independent_weight == 400)
                    deriv[csize+ctype*n_restype+aa] += deriv_aa;
                else if ( num_independent_weight == 20 ) 
                    for (int aa2: range(n_restype)) 
                        deriv[csize+aa2*n_restype+aa] += deriv_aa;
            }
        }
        
        return deriv;
    }
#endif

    virtual void set_param(const std::vector<float>& new_param) override {
        int csize = coeff.size();
        int wsize = n_restype*n_restype;
        if(new_param.size() != size_t(csize+wsize))
            throw string("the size of parameters should be the size of coeff + (the number of amino acid type)^2");

        std::vector<float> params(csize+wsize);
        copy(begin(new_param),end(new_param), begin(params));

        for( int i : range(csize)) coeff[i] = params[i];
        for(int aa: range(wsize)) weights[aa] = params[csize+aa];
    }
};
static RegisterNodeType<NonlinearCoupling,1> nonlinear_coupling_node("nonlinear_coupling");

struct SigmoidCoupling : public PotentialNode {
    CoordNode& input;
    int n_restype, n_residue;
    vector<float> scale;      // length n_residue
    vector<float> center;     // length n_residue
    vector<float> sharpness;  // length n_residue
    vector<float> weights;
    vector<int>   coupling_types; // length input n_residue
    vector<float> wnumber;
    int num_independent_weight;

    SigmoidCoupling(hid_t grp, CoordNode& input_):
        PotentialNode(),
        input(input_),
        n_restype(get_dset_size(1,grp, "scale")[0]),
        n_residue(input.n_elem/n_restype),
        scale(n_restype),
        center(n_restype),
        sharpness(n_restype),
        weights(400),
        coupling_types(n_residue),
        wnumber(n_residue),
        num_independent_weight(read_attribute<int>(grp, ".", "number_independent_weights"))
    {
        check_elem_width(input, 1);  // could be generalized

        check_size(grp, "center", n_restype);
        check_size(grp, "sharpness", n_restype);
        check_size(grp, "coupling_types", n_residue);
        check_size(grp, "weights", n_restype*n_restype);

        traverse_dset<1,float>(grp,"scale",[&](size_t n, float x){scale[n]=x;});
        traverse_dset<1,float>(grp,"center",[&](size_t n, float x){center[n]=x;});
        traverse_dset<1,float>(grp,"sharpness",[&](size_t n, float x){sharpness[n]=x;});
        traverse_dset<1,float>(grp,"weights",[&](size_t n, float x){weights[n]=x;});
        traverse_dset<1,int>(grp,"coupling_types",[&](size_t n, int x){coupling_types[n]=x;});

        for(int i: coupling_types) if(i<0 || i>=n_restype) throw string("invalid coupling type");

        if (num_independent_weight != 1 and num_independent_weight != 20 and num_independent_weight != 400) 
            throw string("the number of independent weights should be 1, 20 or 400");

        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>("sigmoid_coupling", {n_residue}, [&](float* buffer) {
                for(int nr: range(n_residue)) {
                    wnumber[nr] = 0.0;
                    int ctype = coupling_types[nr];
                    for(int aa: range(n_restype)) 
                        wnumber[nr] += weights[ctype*n_restype+aa] * input.output(0, nr*n_restype+aa);
                    auto out = compact_sigmoid(wnumber[nr]-center[ctype], sharpness[ctype]);
                    buffer[nr] = out.x() * scale[ctype];
                }});
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("sigmoid_coupling");

        float pot = 0.f;
        for(int nr: range(n_residue)) {
            wnumber[nr] = 0.0;
            int ctype = coupling_types[nr];
            for(int aa: range(n_restype)) {
                wnumber[nr] += weights[ctype*n_restype+aa] * input.output(0, nr*n_restype+aa);
            }
            auto out = compact_sigmoid(wnumber[nr]-center[ctype], sharpness[ctype]);

            pot += scale[ctype] * out.x();
            for(int aa: range(n_restype)) 
                input.sens(0,nr*n_restype+aa) +=  weights[ctype*n_restype+aa] * scale[ctype] * out.y();
        }
        potential = pot;
    }

    virtual std::vector<float> get_param() const override {
        int wsize = n_restype*n_restype;
        vector<float> params(n_restype*3+wsize, 0.f);

        for(int ct: range(n_restype)) {
            params[0*n_restype+ct] = scale[ct];
            params[1*n_restype+ct] = center[ct];
            params[2*n_restype+ct] = sharpness[ct];
        }

        for(int ct: range(wsize)) 
            params[3*n_restype+ct] = weights[ct];

        return params;
    }

    virtual void set_param(const std::vector<float>& new_param) override {
        int csize = 3*n_restype;
        int wsize = n_restype*n_restype;
        if(new_param.size() != size_t(csize+wsize))
            throw string("the size of parameters should be 3*the number of amino acid type + (the number of amino acid type)^2");

        for(int ct: range(n_restype)) {
             scale[ct]     = new_param[0*n_restype+ct];
             center[ct]    = new_param[1*n_restype+ct];
             sharpness[ct] = new_param[2*n_restype+ct];
        }

        for(int ct: range(wsize)) 
             weights[ct]   = new_param[3*n_restype+ct];
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        int csize = 3*n_restype; // for scale. center and sharpness
        int wsize = n_restype*n_restype;

        vector<float> deriv(csize+wsize, 0.f);
        vector<float> sens(n_residue, 0.f);

        for(int nr: range(n_residue)) {
            wnumber[nr] = 0.0;
            int ctype = coupling_types[nr];

            for(int aa: range(n_restype)) 
                wnumber[nr] += weights[ctype*n_restype+aa] * input.output(0, nr*n_restype+aa);

            auto dist_coord = wnumber[nr]-center[ctype];
            auto out = compact_sigmoid(dist_coord, sharpness[ctype]);

            deriv[0*n_restype+ctype] +=  out.x();
            deriv[1*n_restype+ctype] += -scale[ctype] * out.y();
            deriv[2*n_restype+ctype] +=  scale[ctype] * dist_coord * out.y()/sharpness[ctype];

            if ( num_independent_weight > 1) 
                sens[nr]              =  out.y() * scale[ctype];
        }

        if ( num_independent_weight == 1 )
            return deriv;

        for(int nr: range(n_residue)) {
            int ctype = coupling_types[nr];
            for(int aa: range(n_restype)) {
                auto deriv_aa = sens[nr] * input.output(0, nr*n_restype+aa);
                if ( num_independent_weight == 400)
                    deriv[csize+ctype*n_restype+aa] += deriv_aa;
                else if ( num_independent_weight == 20 ) 
                    for (int aa2: range(n_restype)) 
                        deriv[csize+aa2*n_restype+aa] += deriv_aa;
            }
        }

        return deriv;
    }
#endif

};
static RegisterNodeType<SigmoidCoupling,1> sigmoid_coupling_node("sigmoid_coupling");

//struct BackboneWeightedBurialLevel : public CoordNode {
//
//    CoordNode& env_sc;
//
//    int n_restype;
//    int n_bb;
//    float scale;
//    float center;
//    float sharpness;
//    float hbond_weight;
//    vector<float> weights;
//
//    BackboneWeightedBurialLevel(hid_t grp, CoordNode& env_sc_):
//
//        PotentialNode(),
//        env_sc      (env_sc_),
//        n_restype   (read_attribute<int>(grp, ".", "num_aa_types")),
//        n_bb        (env_sc.n_elem/n_restype),
//        weights     (get_dset_size(1,grp, "weights")[0])
//    {
//        check_elem_width(env_sc, 1);  // could be generalized
//        traverse_dset<1,float>(grp,"weights",[&](size_t n, float x){weights[n]=x;});
//    }
//
//    virtual void compute_value(ComputeMode mode) override {
//        Timer timer("backbone_weighted_bl");
//
//        for(int nr: range(n_bb)) {
//            float bl = 0.0;
//            for(int aa: range(n_restype)) 
//                bl += weights[aa] * env_sc.output(0, nr*n_restype+aa);
//            bl += hbond_weight*env_hb.output(0, nr);
//
//            output(0, nr) = bl;
//        }
//    }
//
//    virtual void propagate_deriv() override {
//        Timer timer("d_backbone_weighted_bl");
//
//        for(int nr: range(n_bb)) {
//            for(int aa: range(n_restype)) 
//                env_sc.sens(0, nr*n_restype+aa) +=  weights[aa] * d_pot;
//        }
//    }
//};
//
//static RegisterNodeType<BackboneWeightedBurialLevel,1>_node("bb_sigmoid_coupling");

struct BackboneSigmoidCoupling : public PotentialNode {

    CoordNode& env_hb;
    CoordNode& env_sc;
    int n_restype;
    int n_bb;
    float scale;
    float center;
    float sharpness;
    float hbond_weight;
    vector<float> weights;

    BackboneSigmoidCoupling(hid_t grp, CoordNode& env_sc_, CoordNode& env_hb_):
        PotentialNode(),
        env_hb      (env_hb_),
        env_sc      (env_sc_),
        n_restype   (read_attribute<int>(grp, ".", "num_aa_types")),
        n_bb        (env_sc.n_elem/n_restype),
        scale       (read_attribute<float>(grp, ".", "scale")),
        center      (read_attribute<float>(grp, ".", "center")),
        sharpness   (read_attribute<float>(grp, ".", "sharpness")),
        hbond_weight(read_attribute<float>(grp, ".", "hbond_weight")),
        weights     (get_dset_size(1,grp, "weights")[0])
    {
        check_elem_width(env_sc, 1);  // could be generalized
        check_elem_width(env_hb, 1);  // could be generalized
        traverse_dset<1,float>(grp,"weights",[&](size_t n, float x){weights[n]=x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("backbone_sigmoid_coupling");
        float pot = 0.f;
        for(int nr: range(n_bb)) {
            float bl = 0.0;
            for(int aa: range(n_restype)) 
                bl += weights[aa] * env_sc.output(0, nr*n_restype+aa);
            bl += hbond_weight*env_hb.output(0, nr);
            auto out = compact_sigmoid(bl-center, sharpness);

            pot        += scale * out.x();
            auto d_pot  = scale * out.y();
            for(int aa: range(n_restype)) 
                env_sc.sens(0, nr*n_restype+aa) +=  weights[aa] * d_pot;
            env_hb.sens(0, nr) +=  hbond_weight* d_pot;
        }
        potential = pot;
    }

    virtual std::vector<float> get_param() const override {
        vector<float> params(4, 0.f);
        params[0] = scale;
        params[1] = center;
        params[2] = sharpness;
        params[3] = hbond_weight;
        return params;
    }
    virtual void set_param(const std::vector<float>& new_param) override {
        scale        = new_param[0];
        center       = new_param[1];
        sharpness    = new_param[2];
        hbond_weight = new_param[3];
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        vector<float> deriv(4, 0.f);
        for(int nr: range(n_bb)) {
            float bl = 0.0;
            for(int aa: range(n_restype)) 
                bl += weights[aa] * env_sc.output(0, nr*n_restype+aa);
            bl += hbond_weight*env_hb.output(0, nr);

            auto dist_coord = bl-center;
            auto out = compact_sigmoid(dist_coord, sharpness);
            deriv[0] +=  out.x();
            //deriv[1] += -scale * out.y();
            //deriv[2] +=  scale * dist_coord * out.y()/sharpness;
            //deriv[3] +=  scale * out.y() * env_hb.output(0, nr);
        }
        return deriv;
    }
#endif
};

static RegisterNodeType<BackboneSigmoidCoupling,2> bb_sigmoid_coupling_node("bb_sigmoid_coupling");

}
