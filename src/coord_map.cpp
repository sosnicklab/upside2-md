#include <math.h>
#include <string>
#include <algorithm>
#include <iostream>
#include "deriv_engine.h"
#include "timing.h"
#include "state_logger.h"
#include "vector_math.h"
#include <Eigen/Dense>

using namespace std;
using namespace h5;
using namespace Eigen;

struct PastePos : public CoordNode {
    CoordNode &pos1;
    CoordNode &pos2;

    int n_dim1;
    int n_dim2;

    struct Param {
        index_t index_pos1;
        index_t index_pos2;
    };
    vector<Param> params;

    PastePos(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(get_dset_size(1,grp,"index_pos1")[0], read_attribute<int>(grp, ".", "n_dim1")+read_attribute<int>(grp, ".", "n_dim2")),
        pos1(pos1_),
        pos2(pos2_),
        n_dim1(read_attribute<int>(grp, ".", "n_dim1")),
        n_dim2(read_attribute<int>(grp, ".", "n_dim2")),
        params(n_elem)
    {
        check_size(grp, "index_pos1", n_elem);
        check_size(grp, "index_pos2", n_elem);
        traverse_dset<1,int>(grp,"index_pos1" ,[&](size_t ne, int x){params[ne].index_pos1 = x;});
        traverse_dset<1,int>(grp,"index_pos2" ,[&](size_t ne, int x){params[ne].index_pos2 = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("paste_pos");

        for(int ne=0; ne<n_elem; ++ne) {
            auto p = params[ne];
            for(int nd=0; nd<n_dim1; ++nd) 
                output(nd, ne)    = pos1.output(nd, p.index_pos1);
            for(int nd=0; nd<n_dim2; ++nd)
                output(nd+n_dim1, ne) = pos2.output(nd, p.index_pos2);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_paste_pos");

        for(int ne=0; ne<n_elem; ++ne) {
            auto p = params[ne];
            for(int nd=0; nd<n_dim1; ++nd) 
                pos1.sens(nd,p.index_pos1) += sens(nd,ne);
            for(int nd=0; nd<n_dim2; ++nd)
                pos2.sens(nd,p.index_pos2) += sens(nd+n_dim1,ne);
        }
    }
};

static RegisterNodeType<PastePos,2> paste_pos_node("paste_pos");

struct PasteRama : public CoordNode {

    struct Param {
        index_t index_pos1;
        index_t index_pos2;
    };

    CoordNode &pos1;
    CoordNode &pos2;
    vector<Param> params;

    PasteRama(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(get_dset_size(1,grp,"index_pos1")[0], 2),
        pos1(pos1_),
        pos2(pos2_),
        params(n_elem)
    {
        check_size(grp, "index_pos1", n_elem);
        check_size(grp, "index_pos2", n_elem);
        traverse_dset<1,int>(grp,"index_pos1" ,[&](size_t ne, int x){params[ne].index_pos1 = x;});
        traverse_dset<1,int>(grp,"index_pos2" ,[&](size_t ne, int x){params[ne].index_pos2 = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("paste_rama");
        for(int ne=0; ne<n_elem; ++ne) {
            auto p = params[ne];
            output(0, ne) = pos1.output(0, p.index_pos1);
            output(1, ne) = pos2.output(0, p.index_pos2);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_paste_rama");
        for(int ne=0; ne<n_elem; ++ne) {
            auto p = params[ne];
            pos1.sens(0, p.index_pos1) += sens(0,ne);
            pos2.sens(0, p.index_pos2) += sens(1,ne);
        }
    }
};
static RegisterNodeType<PasteRama,2> paste_rama_node("rama_coord");

struct CatPos : public CoordNode {
    CoordNode &pos1;
    CoordNode &pos2;
    int pos1_size;
    int pos2_size;
    int n_dim;

    vector<int> index;

    CatPos(hid_t grp, CoordNode& pos1_, CoordNode& pos2_):
        CoordNode(get_dset_size(1,grp,"index_pos1")[0]+get_dset_size(1,grp,"index_pos2")[0], read_attribute<int>(grp, ".", "n_dim")),
        pos1(pos1_),
        pos2(pos2_),
        pos1_size(get_dset_size(1,grp,"index_pos1")[0]),
        pos2_size(get_dset_size(1,grp,"index_pos2")[0]),
        n_dim(read_attribute<int>(grp, ".", "n_dim")),
        index(n_elem)
    {
        traverse_dset<1,int>(grp,"index_pos1" ,[&](size_t ne, int x){index[ne] = x;});
        traverse_dset<1,int>(grp,"index_pos2" ,[&](size_t ne, int x){index[ne+pos1_size] = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("cat_pos");

        for(int ne=0; ne<pos1_size; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                output(nd, ne) = pos1.output(nd, index[ne]);
        }
        for(int ne=pos1_size; ne<n_elem; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                output(nd, ne) = pos2.output(nd, index[ne]);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_cat_pos");

        for(int ne=0; ne<pos1_size; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                pos1.sens(nd,index[ne]) += sens(nd,ne);
        }
        for(int ne=pos1_size; ne<n_elem; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                pos2.sens(nd,index[ne]) += sens(nd,ne);
        }
    }
};

static RegisterNodeType<CatPos,2> cat_pos_node("cat_pos");

struct Concat : public CoordNode
{
    vector<CoordNode*> coord_nodes;

    static int sum_n_elem(const vector<CoordNode*>& coord_nodes_) {
        int ne = 0;
        for(auto& cn: coord_nodes_) ne += cn->n_elem;
        return ne;
    }

    Concat(hid_t grp, const std::vector<CoordNode*> &coord_nodes_):
        CoordNode(sum_n_elem(coord_nodes_), coord_nodes[0]->elem_width),
        coord_nodes(coord_nodes_)
    {
        for(auto cn: coord_nodes)
            if(cn->n_elem != coord_nodes[0]->n_elem)
                throw string("Coord node n_elem mismatch");
    }

    virtual void compute_value(ComputeMode mode) override {
        int loc = 0;
        for(auto cn: coord_nodes) {
            int n_elem_cn = cn->n_elem;
            VecArray cn_output = cn->output;
            
            for(int ne=0; ne<n_elem_cn; ++ne){
                for(int nw=0; nw<elem_width; ++nw)
                    output(nw,loc) = cn_output(nw,ne);
                loc++;
            }
        }
        assert(loc==n_elem);
    } 

    virtual void propagate_deriv() override {
        int loc = 0;
        for(auto cn: coord_nodes) {
            int n_elem_cn = cn->n_elem;
            VecArray cn_sens = cn->sens;
            
            for(int ne=0; ne<n_elem_cn; ++ne){
                for(int nw=0; nw<elem_width; ++nw)
                    cn_sens(nw,loc) += sens(nw,ne);
                loc++;
            }
        }
        assert(loc==n_elem);
    }  
}; 
static RegisterNodeType<Concat,-1> concat_node("concat");

struct PosWeight : public CoordNode {
    CoordNode &pos;
    int n_dim;
    vector<int> index_pos;
    vector<float> weight;

    PosWeight(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1,grp,"index_pos")[0], read_attribute<int>(grp, ".", "n_dim")+1),
        pos(pos_),
        n_dim(read_attribute<int>(grp, ".", "n_dim")),
        index_pos(n_elem),
        weight(n_elem)
    {
        check_size(grp, "weight", n_elem);
        traverse_dset<1,int>(grp,"index_pos" ,[&](size_t ne, int x){index_pos[ne] = x;});
        traverse_dset<1,float>(grp,"weight" , [&](size_t ne, int x){weight[ne]    = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("weight_pos");
        for(int ne=0; ne<n_elem; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                output(nd, ne) = pos.output(nd, index_pos[ne]);
            output(n_dim, ne)  = weight[ne];
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_weight_pos");
        for(int ne=0; ne<n_elem; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                pos.sens(nd, index_pos[ne]) += sens(nd,ne);
        }
    }
};

static RegisterNodeType<PosWeight,1> pos_weight_node("weight_pos");

struct SelectPos : public CoordNode {
    CoordNode &pos;
    int n_dim;

    vector<int> index_pos;
    vector<int> index_dim;

    SelectPos(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1,grp,"index_pos")[0], read_attribute<int>(grp, ".", "n_dim")),
        pos(pos_),
        n_dim(read_attribute<int>(grp, ".", "n_dim")),
        index_pos(n_elem),
        index_dim(n_dim)
    {
        traverse_dset<1,int>(grp,"index_pos" ,[&](size_t ne, int x){index_pos[ne] = x;});
        traverse_dset<1,int>(grp,"index_dim" ,[&](size_t ne, int x){index_dim[ne] = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("select_pos");

        for(int ne=0; ne<n_elem; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                output(nd, ne) = pos.output(index_dim[nd], index_pos[ne]);
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_select_pos");

        for(int ne=0; ne<n_elem; ++ne) {
            for(int nd=0; nd<n_dim; ++nd)
                pos.sens(index_dim[nd], index_pos[ne]) += sens(nd, ne);
        }
    }
};

static RegisterNodeType<SelectPos,1> select_pos_node("select_pos");

struct Slice : public CoordNode
{
    int n_atom;
    vector<int> id;
    CoordNode& pos;

    Slice(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0], pos_.elem_width),
        n_atom(n_elem),
        id(n_atom),
        pos(pos_)
    {
        check_size(grp, "id", n_atom);
        traverse_dset<1,int> (grp, "id", [&](size_t i, int x) {id[i] = x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        for (int na = 0; na < n_atom; na++) {
            for (int d = 0; d < elem_width; d++) {
                output(d, na) = pos.output(d, id[na]);
            }
        }
    }

    virtual void propagate_deriv() override {
        for (int na = 0; na < n_atom; na++) {
            for (int d = 0; d < elem_width; d++) {
                pos.sens(d, id[na]) += sens(d, na);
            }
        }
    }
 };
static RegisterNodeType<Slice,1> slice_node("slice");


struct PosSum : public CoordNode {
    CoordNode &pos;
    int n_dim;

    vector<vector<int>> index_pos;
    vector<int> index_dim;

    PosSum(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(2,grp,"index_pos")[0], read_attribute<int>(grp, ".", "n_dim")),
        pos(pos_),
        n_dim(read_attribute<int>(grp, ".", "n_dim")),
        index_pos(n_elem),
        index_dim(n_dim)
    {
        traverse_dset<2, int>(grp, "index_pos", [&](size_t np, size_t na, int x) {index_pos[np].push_back(x);});
        traverse_dset<1, int>(grp, "index_dim" ,[&](size_t ne, int x){index_dim[ne] = x;});
        for(int ne=0; ne<n_elem; ++ne) {
            auto &ndx = index_pos[ne]; 
            vector<int>::iterator ip = std::unique(ndx.begin(), ndx.end()); 
            ndx.resize(std::distance(ndx.begin(), ip)); 
        }

        if(logging(LOG_DETAILED)) {
            default_logger->add_logger<float>("pos_sum", {n_elem, n_dim}, [&](float* buffer) {
                for(int ne=0; ne<n_elem; ++ne) {
                    for(int nd=0; nd<n_dim; ++nd) 
                        buffer[ne*n_dim + nd] = output(nd, ne);
                }
            });
        }
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("pos_sum");

        for(int ne=0; ne<n_elem; ++ne) {
            auto index = index_pos[ne];
            for(int nd=0; nd<n_dim; ++nd) {
                output(nd, ne) = 0.0;
                for( int na=0;na<(int)index.size(); ++na) 
                    output(nd, ne) += pos.output(index_dim[nd], index[na]);
            }
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_pos_sum");

        for(int ne=0; ne<n_elem; ++ne) {
            auto index = index_pos[ne];
            for(int nd=0; nd<n_dim; ++nd) {
                auto scale_sens = sens(nd, ne);
                for( int na=0;na<(int)index.size(); ++na) 
                    pos.sens(index_dim[nd], index[na]) += scale_sens;
            }
        }
    }
};

static RegisterNodeType<PosSum,1> pos_sum_node("sum_of_values");

struct ScaledSum: public PotentialNode
{
    CoordNode& input;
    float scale;

    ScaledSum(hid_t grp, CoordNode& input_):
        PotentialNode(),
        input(input_),
        scale(read_attribute<float>(grp, ".", "scale"))
    {

        if(input.elem_width != 1u) throw string("Sum only works on elem width 1");
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("scaled_sum")); 
        VecArray value = input.output;
        VecArray sens  = input.sens;
        int n_elem = input.n_elem;

        float pot = 0.f;
        for(int i=0; i<n_elem; ++i) pot += value(0,i);
        pot *= scale;
        potential = pot;

        for(int i=0; i<n_elem; ++i) sens(0,i) += scale;
    }
};
static RegisterNodeType<ScaledSum,1> scaled_sum_node("scaled_sum");


inline float relu(float x) {
    return x<0.f ? 0.f : x;
}

struct BackboneFeaturizer : public CoordNode
{
    CoordNode& rama;
    CoordNode& hbond;

    struct Param {
        int rama_idx;
        int donor_idx;
        int acceptor_idx;
    };
    vector<Param> params;

    BackboneFeaturizer(hid_t grp, CoordNode& rama_, CoordNode& hbond_):
        CoordNode(get_dset_size(1, grp, "rama_idx")[0], 6),
        rama(rama_),
        hbond(hbond_),
        params(n_elem)
    {
        check_size(grp, "rama_idx", n_elem);
        check_size(grp, "hbond_idx", n_elem, 2);

        traverse_dset<1,int>(grp, "rama_idx", [&](size_t ne, int x){params[ne].rama_idx=x;});

        // -1 is special code to indicate no hbond donor or acceptor
        traverse_dset<2,int>(grp, "hbond_idx", [&](size_t ne, size_t is_acceptor, int x){
                if(is_acceptor) params[ne].acceptor_idx=x;
                else            params[ne].donor_idx   =x;});
    }

    virtual void compute_value(ComputeMode mode) override {
        VecArray ramac = rama.output;
        VecArray hbondc = hbond.output;

        for(int ne=0; ne<n_elem; ++ne) {
            auto& p = params[ne];

            float phi = ramac(0,p.rama_idx);
            float psi = ramac(1,p.rama_idx);
            float don_hb = p.donor_idx   ==-1 ? 0.f : hbondc(6,p.donor_idx);
            float acc_hb = p.acceptor_idx==-1 ? 0.f : hbondc(6,p.acceptor_idx);

            output(0,ne) = sin(phi);
            output(1,ne) = cos(phi);
            output(2,ne) = sin(psi);
            output(3,ne) = cos(psi);
            output(4,ne) = don_hb;
            output(5,ne) = acc_hb;
        }
    }

    virtual void propagate_deriv() override {
        VecArray rama_s = rama.sens;
        VecArray hbond_s = hbond.sens;

        for(int ne=0; ne<n_elem; ++ne) {
            auto& p = params[ne];

            rama_s(0,p.rama_idx) += sens(1,ne)*(-output(0,ne)) + sens(0,ne)*output(1,ne);
            rama_s(1,p.rama_idx) += sens(3,ne)*(-output(2,ne)) + sens(2,ne)*output(3,ne);
            if(p.donor_idx   !=-1) hbond_s(6,p.donor_idx   ) += sens(4,ne);
            if(p.acceptor_idx!=-1) hbond_s(6,p.acceptor_idx) += sens(5,ne);
        }
    }
};
static RegisterNodeType<BackboneFeaturizer,2> backbone_featurizer_node("backbone_featurizer");


struct Conv1D : public CoordNode
{
    CoordNode& input;
    int n_elem_input;
    enum ActivationT {ReLU, Tanh, Identity};

    int conv_width, in_channels, out_channels;

    MatrixXf    weights;
    RowVectorXf bias;
    ActivationT activation;

    // temporaries
    MatrixXf input_conv_format;
    MatrixXf matmul_output;

    Conv1D(hid_t grp, CoordNode& input_):
        CoordNode   (input_.n_elem - get_dset_size(3, grp, "weights")[0] + 1,
                     get_dset_size(3, grp, "weights")[2]),
        input(input_),
        n_elem_input(input.n_elem),

        conv_width  (get_dset_size(3, grp, "weights")[0]),
        in_channels (get_dset_size(3, grp, "weights")[1]),
        out_channels(get_dset_size(3, grp, "weights")[2]),

        weights(conv_width*in_channels, out_channels),
        bias(out_channels),
        input_conv_format(n_elem, conv_width*in_channels),
        matmul_output    (n_elem, out_channels)
    {
        traverse_dset<3,float>(grp, "weights", [&](size_t nw, size_t in_c, size_t out_c, float x) {
                weights(nw*in_channels+in_c, out_c) = x;});

        traverse_dset<1,float>(grp, "bias", [&](size_t out_c, float x) {bias[out_c] = x;});

        // Read activation
        // This really should be a single string instead of a vector of them,
        // but I currently only have a reader for string vectors
        auto activation_str = read_attribute<vector<string>>(grp, ".", "activation");
        if(activation_str.size() != 1u) throw string("Invalid number of activations");

        if     (activation_str[0] == "ReLU") activation = ReLU;
        else if(activation_str[0] == "Tanh") activation = Tanh;
        else if(activation_str[0] == "Identity") activation = Identity;
        else throw string("Invalid activation name");
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer(string("conv1d")); 
        VecArray inputc = input.output;
        
        int n_elem_output = n_elem;

        // Copy input into proper format so that convolution becomes matmul
        for(int nr=0; nr<n_elem_output; ++nr)
            for(int nw=0; nw<conv_width; ++nw)
                for(int nc=0; nc<in_channels; ++nc)
                    input_conv_format(nr,nw*in_channels+nc) = inputc(nc,nr+nw);
       
        // Multiply by weights to make matmul_output
        matmul_output.noalias() = input_conv_format * weights;
        matmul_output.rowwise() += bias;

        // Apply activation and store in output
        // We also keep matmul_output around for storing information for backprop
        switch(activation) {
            case Identity:
                for(int nr=0; nr<n_elem_output; ++nr)
                    for(int nc=0; nc<out_channels; ++nc)
                        output(nc,nr) = matmul_output(nr,nc);
                break;
            case ReLU:
                for(int nr=0; nr<n_elem_output; ++nr)
                    for(int nc=0; nc<out_channels; ++nc)
                        output(nc,nr) = relu(matmul_output(nr,nc));
                break;
            case Tanh:
                for(int nr=0; nr<n_elem_output; ++nr)
                    for(int nc=0; nc<out_channels; ++nc)
                        output(nc,nr) = tanh(matmul_output(nr,nc));
                break;
        }
    }

    virtual void propagate_deriv() override {
        int n_elem_output = n_elem;

        // Backpropagate activations
        switch(activation) {
            case Identity:
                for(int nr=0; nr<n_elem_output; ++nr)
                    for(int nc=0; nc<out_channels; ++nc)
                        matmul_output(nr,nc) = sens(nc,nr);
                break;
            case ReLU:
                for(int nr=0; nr<n_elem_output; ++nr)
                    for(int nc=0; nc<out_channels; ++nc)
                        matmul_output(nr,nc) = output(nc,nr)>0.f ? sens(nc,nr) : 0.f;
                break;
            case Tanh:
                for(int nr=0; nr<n_elem_output; ++nr)
                    for(int nc=0; nc<out_channels; ++nc)
                        matmul_output(nr,nc) = sens(nc,nr) * (1.f-sqr(output(nc,nr)));
                break;
        }

        // Backpropagate convolution
        input_conv_format.noalias() = matmul_output * weights.transpose();

        // Backpropagte into sens
        VecArray inp_sens = input.sens;
        for(int nr=0; nr<n_elem_output; ++nr)
            for(int nw=0; nw<conv_width; ++nw)
                for(int nc=0; nc<in_channels; ++nc)
                    inp_sens(nc,nr+nw) += input_conv_format(nr,nw*in_channels+nc);
    }
};
static RegisterNodeType<Conv1D,1> conv1d_node("conv1d");

struct SigmoidCoord : public CoordNode {

    struct Params {
        index_t id;
        float c;
        float s;
	float sign;
    };

    int n_elem;
    CoordNode& pos;
    vector<Params> params;
    vector<float> deriv;
    int n_dim;
    int dim1;

    SigmoidCoord(hid_t grp, CoordNode& pos_):
        CoordNode(get_dset_size(1, grp, "id")[0],1),
        n_elem(get_dset_size(1, grp, "id")[0]), 
	pos(pos_), 
	params(n_elem),
	deriv(n_elem),
	n_dim(pos.output.row_width),
        dim1(read_attribute<int>(grp, ".", "dim1"))
    {
        check_size(grp, "c",    n_elem);
        check_size(grp, "s",    n_elem);
        check_size(grp, "sign", n_elem);

        auto& p = params;
        traverse_dset<1,int>  (grp, "id",   [&](size_t i, int   x) {p[i].id   = x;});
        traverse_dset<1,float>(grp, "c",    [&](size_t i, float x) {p[i].c    = x;});
        traverse_dset<1,float>(grp, "s",    [&](size_t i, float x) {p[i].s    = x;});
        traverse_dset<1,float>(grp, "sign", [&](size_t i, float x) {p[i].sign = x;});
	assert(dim1<n_dim);

        for(int ne=0; ne<n_elem; ++ne) {
            deriv[ne] = 0.f;
	}
    }

    virtual void compute_value(ComputeMode mode) override {
        Timer timer("sigmoid_coord");
        for(int ne=0; ne<n_elem; ++ne) {
            auto& p  = params[ne];
	    auto x   = pos.output(dim1, p.id);
	    auto sig = compact_sigmoid(x-p.c, p.s);
	    if (p.sign > 0) {
                output(0, ne) = sig.x();
		deriv[ne] = sig.y();
	    }
	    else {
                output(0, ne) = 1.f - sig.x();
		deriv[ne] = -sig.y();
            }
        }
    }

    virtual void propagate_deriv() override {
        Timer timer("d_sigmoid_coord");
        VecArray pos_sens = pos.sens;
        for(int ne=0; ne<n_elem; ++ne) 
            pos_sens(dim1, params[ne].id) += sens(0, ne)*deriv[ne];
    }
};

static RegisterNodeType<SigmoidCoord,1> sigmoid_coord_node("SigmoidCoord");

