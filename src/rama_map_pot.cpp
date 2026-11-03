#include "deriv_engine.h"
#include "timing.h"
#include "spline.h"
#include "state_logger.h"
#include "device_nodes.h"

using namespace h5;
using namespace std;

struct RamaMapParams {
    index_t residue;
    int       rama_map_id;
};


struct RamaMapPot : public PotentialNode
{
    int n_residue;
    CoordNode& rama;
    vector<RamaMapParams> params;
    LayeredPeriodicSpline2D<1> rama_map_data;
    vector<float> residue_potential;
    bool log_pot; // if false, never log potential
    vector<float> rama_pot; // the map values the spline is fitted to (the node's parameters)

    RamaMapPot(hid_t grp, CoordNode& rama_):
        PotentialNode(),
        n_residue(get_dset_size(1, grp, "residue_id")[0]), 
        rama(rama_), 
        params(n_residue),
        rama_map_data(
                get_dset_size(3, grp, "rama_pot")[0], 
                get_dset_size(3, grp, "rama_pot")[1], 
                get_dset_size(3, grp, "rama_pot")[2]),
        residue_potential(n_residue),
        log_pot(read_attribute<int>(grp,".","log_pot",1))
    {
        auto& r = rama_map_data;
        check_size(grp, "residue_id",  n_residue);
        check_size(grp, "rama_map_id", n_residue);
        check_size(grp, "rama_pot",    r.n_layer, r.nx, r.ny);

        if(r.nx != r.ny) throw string("must have same x and y grid spacing for Rama maps");
        check_elem_width(rama, 2);
        vector<double> raw_data(r.n_layer * r.nx * r.ny);

        traverse_dset<1,int>   (grp, "residue_id",  [&](size_t i, int x) {
                if(x<0 || x>=rama.n_elem) throw string("rama map residue_id outside its rama coordinate");
                params[i].residue = x;});
        traverse_dset<1,int>   (grp, "rama_map_id", [&](size_t i, int x) {
                if(x<0 || x>=r.n_layer) throw string("rama map rama_map_id outside rama_pot");
                params[i].rama_map_id = x;});
        traverse_dset<3,double>(grp, "rama_pot",    [&](size_t il, size_t ix, size_t iy, double x) {
                raw_data[(il*r.nx + ix)*r.ny + iy] = x;});
        set_param(vector<float>(raw_data.begin(), raw_data.end()));

        if(log_pot && logging(LOG_DETAILED))
            default_logger->add_logger<float>("rama_map_potential", {n_residue}, [&](float* buffer) {
                for(int nr: range(n_residue)) 
                    buffer[nr] = residue_potential[nr];
                    });
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("rama_map_pot"));

        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        VecArray ramac     = rama.output;
        VecArray rama_sens = rama.sens;
        if(pot) *pot = 0.f;

        // add a litte paranoia to make sure there are no rounding problems
        const float scale = rama_map_data.nx * (0.5f/M_PI_F - 1e-7f);
        const float shift = M_PI_F;

        if(pot) *pot = 0.f;
        for(int nr=0; nr<n_residue; ++nr) {
            const auto& p = params[nr];
            auto r = load_vec<2>(ramac, p.residue);

            float value,dx,dy;
            rama_map_data.evaluate_value_and_deriv(&value,&dx,&dy, p.rama_map_id, 
                    (r.v[0]+shift)*scale, (r.v[1]+shift)*scale);

            if(pot) {*pot += value; residue_potential[nr] = value;}
            rama_sens(0,p.residue) += dx * scale;
            rama_sens(1,p.residue) += dy * scale;
        }
    }

    virtual std::vector<float> get_param() const override {return rama_pot;}

    // the spline is fitted to the parameters as get_param returns them, so that set_param(get_param()) leaves
    // the node as it was
    virtual void set_param(const std::vector<float>& new_param) override {
        if(new_param.size() != size_t(rama_map_data.n_layer*rama_map_data.nx*rama_map_data.ny))
            throw string("wrong number of parameters");
        rama_pot = new_param;
        vector<double> raw_data(rama_pot.begin(), rama_pot.end());
        rama_map_data.fit_spline(raw_data.data());
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        throw string("the Ramachandran map parameter derivative is not implemented");
    }
#endif

    virtual DeviceNode* make_device(DeviceGraph& g, const vector<DerivComputation*>& systems) override {
        auto param = [&](int ns) -> const vector<RamaMapParams>& {return static_cast<RamaMapPot&>(*systems[ns]).params;};
        // each system's maps, those equal to an earlier system's (the replicas of one protein) uploaded once
        vector<const vector<float>*> distinct;
        vector<int> first_layer, map_id;
        vector<float> coefficients;
        int layer_size = rama_map_data.nx*rama_map_data.ny*16;
        for(int ns: range(g.n_system)) {
            auto& d = static_cast<RamaMapPot&>(*systems[ns]).rama_map_data;
            if(d.nx != rama_map_data.nx || d.ny != rama_map_data.ny)
                throw string("the device engine runs rama_map_pot's maps on one grid in every system");
            size_t k = 0;
            while(k<distinct.size() && *distinct[k] != d.coefficients) ++k;
            if(k == distinct.size()) {
                distinct.push_back(&d.coefficients);
                first_layer.push_back(coefficients.size()/layer_size);
                coefficients.insert(coefficients.end(), d.coefficients.begin(), d.coefficients.end());
            }
            for(auto& p: param(ns)) map_id.push_back(first_layer[k] + p.rama_map_id);
        }
        return device_rama_map_pot(g.potential_of(*this), g.coord(rama),
                g.offset_index(rama, [&](int ns) {vector<int> r; for(auto& p: param(ns)) r.push_back(p.residue); return r;}),
                map_id, coefficients, rama_map_data.nx, rama_map_data.ny,
                g.concat<int>([&](int ns) {return vector<int>(param(ns).size(), ns);}));
    }
};
static RegisterNodeType<RamaMapPot,1> rama_map_pot_node("rama_map_pot");

struct RamaMapPot2 : public PotentialNode
{
    int n_residue;
    CoordNode& rama;
    CoordNode& lambda;
    vector<RamaMapParams> params;
    LayeredPeriodicSpline2D<1> rama_map_data;
    vector<float> rama_pot; // the map values the spline is fitted to (the node's parameters)

    RamaMapPot2(hid_t grp, CoordNode& rama_, CoordNode& lambda_):
        PotentialNode(),
        n_residue(get_dset_size(1, grp, "residue_id")[0]), 
        rama(rama_), 
        lambda(lambda_), 
        params(n_residue),
        rama_map_data(
                get_dset_size(3, grp, "rama_pot")[0], 
                get_dset_size(3, grp, "rama_pot")[1], 
                get_dset_size(3, grp, "rama_pot")[2])
    {
        auto& r = rama_map_data;
        check_size(grp, "residue_id",  n_residue);
        check_size(grp, "rama_map_id", n_residue);
        check_size(grp, "rama_pot",    r.n_layer, r.nx, r.ny);

        if(r.nx != r.ny) throw string("must have same x and y grid spacing for Rama maps");
        check_elem_width(rama, 2);
        vector<double> raw_data(r.n_layer * r.nx * r.ny);

        traverse_dset<1,int>   (grp, "residue_id",  [&](size_t i, int x) {
                if(x<0 || x>=rama.n_elem) throw string("rama map residue_id outside its rama coordinate");
                params[i].residue = x;});
        traverse_dset<1,int>   (grp, "rama_map_id", [&](size_t i, int x) {
                if(x<0 || x>=r.n_layer) throw string("rama map rama_map_id outside rama_pot");
                params[i].rama_map_id = x;});
        traverse_dset<3,double>(grp, "rama_pot",    [&](size_t il, size_t ix, size_t iy, double x) {
                raw_data[(il*r.nx + ix)*r.ny + iy] = x;});
        set_param(vector<float>(raw_data.begin(), raw_data.end()));

        // row nr of lambda weights map row nr
        check_elem_width(lambda, 1);
        if(lambda.n_elem != n_residue) throw string("RamaMap2 needs one lambda row per residue_id row");
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("rama_map_pot"));

        float* pot = mode==PotentialAndDerivMode ? &potential : nullptr;
        VecArray ramac     = rama.output;
        VecArray rama_sens = rama.sens;
        VecArray lambdac     = lambda.output;
        VecArray lambda_sens = lambda.sens;
        if(pot) *pot = 0.f;

        // add a litte paranoia to make sure there are no rounding problems
        const float scale = rama_map_data.nx * (0.5f/M_PI_F - 1e-7f);
        const float shift = M_PI_F;

        if(pot) *pot = 0.f;
        for(int nr=0; nr<n_residue; ++nr) {
            const auto& p = params[nr];
            auto r = load_vec<2>(ramac, p.residue);

            float value,dx,dy;
            rama_map_data.evaluate_value_and_deriv(&value,&dx,&dy, p.rama_map_id, 
                    (r.v[0]+shift)*scale, (r.v[1]+shift)*scale);

	    float energy = value*lambdac(0, nr);
            if(pot) *pot += energy;
            rama_sens(0,p.residue) += dx * scale * lambdac(0, nr);
            rama_sens(1,p.residue) += dy * scale * lambdac(0, nr);
            lambda_sens(0,nr)      += value;
        }
    }

    virtual std::vector<float> get_param() const override {return rama_pot;}

    // the spline is fitted to the parameters as get_param returns them, so that set_param(get_param()) leaves
    // the node as it was
    virtual void set_param(const std::vector<float>& new_param) override {
        if(new_param.size() != size_t(rama_map_data.n_layer*rama_map_data.nx*rama_map_data.ny))
            throw string("wrong number of parameters");
        rama_pot = new_param;
        vector<double> raw_data(rama_pot.begin(), rama_pot.end());
        rama_map_data.fit_spline(raw_data.data());
    }

#ifdef PARAM_DERIV
    virtual std::vector<float> get_param_deriv() override {
        throw string("the Ramachandran map parameter derivative is not implemented");
    }
#endif
};
static RegisterNodeType<RamaMapPot2,2> rama_map_pot2_node("RamaMap2");

