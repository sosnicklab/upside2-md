// Device implementation of rotamer.cpp's RotamerSidechain: the side chains' belief propagation (BP) on the CPU's
// schedule (upside3.md 10.3).  Residues hold 1, 3 or 6 rotamers (classes 0, 1, 2).  BP iterates over the residue
// pairs of classes 1 and 2 (pairs 33, 36, 66), one thread block per system, each system stopping at its own test.
// As fill_holders does, a bead pair of two 1-rotamer residues adds to the energy, and one with a single 1-rotamer
// residue to the other residue's rotamer energies.  A system's residue pairs are found through a dense lookup of
// its pairs of residues, so they are numbered in (residue 1, residue 2) order, and each residue's messages are
// multiplied in the CPU's order: its pairs 33, then 36, then 66, each in that order.
#include "device_interaction_graph.h"
#include "bead_interaction.h"

using namespace std;

namespace {

constexpr int n_pair = 3;  // 33, 36, 66
UPSIDE_HD inline int class_n_rot(int c) {return c==0 ? 1 : c==1 ? 3 : 6;}
UPSIDE_HD inline int pair_c1(int p) {return p==2 ? 2 : 1;}
UPSIDE_HD inline int pair_c2(int p) {return p==0 ? 1 : 2;}

// The messages of an edge from the old node beliefs and messages (EdgeHolder::update_beliefs): into node 1, P
// applied to node 2's cavity; into node 2, P's transpose applied to node 1's
template <int NA, int NB>
UPSIDE_HD inline void edge_messages(const float* P, const float* b1, const float* b2, const float* m_old, float* m) {
    float v1[NA], v2[NB];
    for(int i=0; i<NA; ++i) v1[i] = b1[i] * rcp(1e-30f + m_old[i]);
    for(int j=0; j<NB; ++j) v2[j] = b2[j] * rcp(1e-30f + m_old[NA+j]);
    for(int i=0; i<NA; ++i) {float z = 0.f; for(int j=0; j<NB; ++j) z += P[i*NB+j]*v2[j]; m[i]    = z;}
    for(int j=0; j<NB; ++j) {float z = 0.f; for(int i=0; i<NA; ++i) z += v1[i]*P[i*NB+j]; m[NA+j] = z;}
}

// Everything the BP loops read and write.  Node arrays of BP class k (rotamer class k+1), pair arrays of pair p;
// beliefs and messages are double buffers, the second at stride elements after the first.
struct BPView {
    const float* prob[2];  float* belief[2];  float* marginal[2];  const int* start[2];  size_t belief_stride[2];
    const float* P[n_pair];  float* msg[n_pair];  float* edge_marginal[n_pair];  size_t msg_stride[n_pair];
    const int* node1[n_pair];  const int* node2[n_pair];  const int* n_edge[n_pair];  int cap[n_pair];
    const int* inc_start;  const int* inc_pair;  const int* inc_slot;  int n_node1;  // class 1's BP nodes come first
    float damping, tol;  int max_iter, chunk;
};

template <int k> struct NR {static constexpr int n = k==0 ? 3 : 6;};
template <int p> struct PairT {
    static constexpr int k1 = p==2 ? 1 : 0, k2 = p==0 ? 0 : 1;  // the BP classes of its nodes
    static constexpr int n1 = NR<k1>::n, n2 = NR<k2>::n;
};

// The items of BP's passes, each on its own data.  A pass runs the items of several kinds that touch disjoint data,
// each kind's items in one block_for_parts, so an update synchronises the block three times.

// f(part, i) for i in [0, n[part]) of every part, all parts' items divided among the block's threads in one pass
template <int N, typename F> UPSIDE_DEVICE inline void block_for_parts(const int (&n)[N], const F& f) {
    int total = 0;
    for(int part=0; part<N; ++part) total += n[part];
    block_for(total, [&](int i) {
        int part = 0;
        while(i >= n[part]) i -= n[part++];
        f(part, i);
    });
}

template <int k> UPSIDE_DEVICE inline int bp_n_node(const BPView& v, int ns) {return v.start[k][ns+1]-v.start[k][ns];}

// element i of class k's nodes' rotamers: the probability into buffer pn's beliefs
template <int k> UPSIDE_DEVICE inline void bp_copy_prob(const BPView& v, int ns, int pn, int i) {
    constexpr int nr = NR<k>::n;
    int a = v.start[k][ns];
    v.belief[k][pn*v.belief_stride[k] + a*nr+i] = v.prob[k][a*nr+i];
}

template <int p> UPSIDE_DEVICE inline void bp_messages(const BPView& v, int ns, int pn, int pe, int e) {
    typedef PairT<p> T;
    const float* b1 = v.belief[T::k1] + (1-pn)*v.belief_stride[T::k1];
    const float* b2 = v.belief[T::k2] + (1-pn)*v.belief_stride[T::k2];
    const float* m_old = v.msg[p] + (1-pe)*v.msg_stride[p];
    float*       m_cur = v.msg[p] +    pe *v.msg_stride[p];
    size_t s = size_t(ns)*v.cap[p] + e;
    edge_messages<T::n1,T::n2>(v.P[p] + s*T::n1*T::n2, b1 + v.node1[p][s]*T::n1, b2 + v.node2[p][s]*T::n2,
            m_old + s*(T::n1+T::n2), m_cur + s*(T::n1+T::n2));
}

// node i's belief: its probability times each of its messages in turn, normalised after each, as the CPU's edge
// loop forms it in place
template <int k> UPSIDE_DEVICE inline void bp_products(const BPView& v, int ns, int pn, int pe, int i) {
    constexpr int nr = NR<k>::n;
    float* cur = v.belief[k] + pn*v.belief_stride[k];
    int node = v.start[k][ns]+i;
    int f = (k==0 ? 0 : v.n_node1) + node;
    float b[nr];
    for(int r=0; r<nr; ++r) b[r] = cur[node*nr+r];
    for(int j=v.inc_start[f]; j<v.inc_start[f+1]; ++j) {
        int p = v.inc_pair[j]>>1, side = v.inc_pair[j]&1;
        int n1 = class_n_rot(pair_c1(p)), n2 = class_n_rot(pair_c2(p));
        const float* m = v.msg[p] + pe*v.msg_stride[p] + size_t(v.inc_slot[j])*(n1+n2) + (side ? n1 : 0);
        float s = 0.f;
        for(int r=0; r<nr; ++r) {b[r] = m[r]*b[r]; s += b[r];}
        float scale = rcp(s);
        for(int r=0; r<nr; ++r) b[r] *= scale;
    }
    for(int r=0; r<nr; ++r) cur[node*nr+r] = b[r];
}

// each half of a message to unit L1 norm (the CPU's scale is approximate; it cancels in the normalised beliefs)
template <int p> UPSIDE_DEVICE inline void bp_normalize(const BPView& v, int ns, int pe, int e) {
    typedef PairT<p> T;
    float* m = v.msg[p] + pe*v.msg_stride[p] + (size_t(ns)*v.cap[p]+e)*(T::n1+T::n2);
    float s1 = 0.f, s2 = 0.f;
    for(int i=0; i<T::n1; ++i) s1 += m[i];
    for(int j=0; j<T::n2; ++j) s2 += m[T::n1+j];
    float r1 = rcp(s1), r2 = rcp(s2);
    for(int i=0; i<T::n1; ++i) m[i] *= r1;
    for(int j=0; j<T::n2; ++j) m[T::n1+j] *= r2;
}

// NodeHolder::standardize_belief_update
template <int k> UPSIDE_DEVICE inline void bp_standardize(const BPView& v, int ns, int pn, float damping, int i) {
    constexpr int nr = NR<k>::n;
    int node = v.start[k][ns]+i;
    float* b = v.belief[k] + pn*v.belief_stride[k] + node*nr;
    const float* old = v.belief[k] + (1-pn)*v.belief_stride[k] + node*nr;
    float mx = b[0];
    for(int r=1; r<nr; ++r) mx = b[r] < mx ? mx : b[r];  // std::max(b[r], mx), which keeps a NaN
    float scale = rcp(mx);
    if(damping != 0.f) for(int r=0; r<nr; ++r) b[r] = (1.f-damping)*scale*b[r] + damping*old[r];
    else               for(int r=0; r<nr; ++r) b[r] = scale*b[r];  // zero damping keeps nothing, not even NaN
}

UPSIDE_DEVICE inline void bp_new_beliefs(const BPView& v, int ns, int& pn, int pe, float damping,
        bool swap_for_initial) {
    int n0 = bp_n_node<0>(v, ns), n1 = bp_n_node<1>(v, ns);
    int e0 = v.n_edge[0][ns], e1 = v.n_edge[1][ns], e2 = v.n_edge[2][ns];
    // the probabilities into the beliefs pn, and the messages pe from the old beliefs and messages
    block_for_parts<5>({n0*NR<0>::n, n1*NR<1>::n, e0, e1, e2}, [&](int part, int i) {
        if     (part==0) bp_copy_prob<0>(v, ns, pn, i);
        else if(part==1) bp_copy_prob<1>(v, ns, pn, i);
        else if(part==2) bp_messages<0>(v, ns, pn, pe, i);
        else if(part==3) bp_messages<1>(v, ns, pn, pe, i);
        else             bp_messages<2>(v, ns, pn, pe, i);
    });
    block_for_parts<2>({n0, n1}, [&](int part, int i) {
        if(part==0) bp_products<0>(v, ns, pn, pe, i);
        else        bp_products<1>(v, ns, pn, pe, i);
    });
    int pn_standardize = swap_for_initial ? 1-pn : pn;  // the first cavities then come from the probabilities
    block_for_parts<5>({e0, e1, e2, n0, n1}, [&](int part, int i) {
        if     (part==0) bp_normalize<0>(v, ns, pe, i);
        else if(part==1) bp_normalize<1>(v, ns, pe, i);
        else if(part==2) bp_normalize<2>(v, ns, pe, i);
        else if(part==3) bp_standardize<0>(v, ns, pn_standardize, damping, i);
        else             bp_standardize<1>(v, ns, pn_standardize, damping, i);
    });
    pn = pn_standardize;
}

// the largest change of any belief between the buffers
UPSIDE_DEVICE inline float bp_deviation(const BPView& v, int ns, int pn) {
    int n0 = bp_n_node<0>(v, ns)*NR<0>::n, n1 = bp_n_node<1>(v, ns)*NR<1>::n;
    int a0 = v.start[0][ns]*NR<0>::n, a1 = v.start[1][ns]*NR<1>::n;
    float d = block_max_keep_nan(n0+n1, [&](int i) {
        int k = i<n0 ? 0 : 1, j = (k ? a1 : a0) + (k ? i-n0 : i);
        const float* cur = v.belief[k] +    pn *v.belief_stride[k];
        const float* old = v.belief[k] + (1-pn)*v.belief_stride[k];
        return fabsf(cur[j]-old[j]);
    });
    return n0+n1 ? d : 0.f;
}

template <int k> UPSIDE_DEVICE inline void bp_node_marginals(const BPView& v, int ns, int pn, int i) {
    constexpr int nr = NR<k>::n;
    int node = v.start[k][ns]+i;
    const float* b = v.belief[k] + pn*v.belief_stride[k] + node*nr;
    float s = 0.f;
    for(int r=0; r<nr; ++r) s += b[r];
    float scale = rcp(s);
    for(int r=0; r<nr; ++r) v.marginal[k][node*nr+r] = b[r]*scale;
}

// EdgeHolder::calculate_marginals
template <int p> UPSIDE_DEVICE inline void bp_edge_marginals(const BPView& v, int ns, int pe, int e) {
    typedef PairT<p> T;
    size_t s = size_t(ns)*v.cap[p] + e;
    const float* b1 = v.marginal[T::k1] + v.node1[p][s]*T::n1;
    const float* b2 = v.marginal[T::k2] + v.node2[p][s]*T::n2;
    const float* m  = v.msg[p] + pe*v.msg_stride[p] + s*(T::n1+T::n2);
    const float* P  = v.P[p] + s*T::n1*T::n2;
    float* marg = v.edge_marginal[p] + s*T::n1*T::n2;
    float bc1[T::n1], bc2[T::n2];
    for(int i=0; i<T::n1; ++i) bc1[i] = b1[i] * rcp(1e-30f + m[i]);
    for(int j=0; j<T::n2; ++j) bc2[j] = b2[j] * rcp(1e-30f + m[T::n1+j]);
    float sum = 0.f;
    for(int i=0; i<T::n1; ++i)
        for(int j=0; j<T::n2; ++j) {marg[i*T::n2+j] = P[i*T::n2+j]*bc1[i]*bc2[j]; sum += marg[i*T::n2+j];}
    float scale = rcp(sum);
    for(int i=0; i<T::n1*T::n2; ++i) marg[i] *= scale;
}

// RotamerSidechain::solve_for_marginals for system ns
UPSIDE_DEVICE inline void bp_solve(const BPView& v, int ns) {
    int pn = 0, pe = 0;  // the buffers holding the current beliefs and messages
    int n0 = bp_n_node<0>(v, ns), n1 = bp_n_node<1>(v, ns);
    int e0 = v.n_edge[0][ns], e1 = v.n_edge[1][ns], e2 = v.n_edge[2][ns];

    // the old beliefs start as the probabilities, the old messages as 1
    block_for_parts<5>({n0*NR<0>::n, n1*NR<1>::n, e0*(3+3), e1*(3+6), e2*(6+6)}, [&](int part, int i) {
        if(part==0) bp_copy_prob<0>(v, ns, 1-pn, i);
        else if(part==1) bp_copy_prob<1>(v, ns, 1-pn, i);
        else {
            int p = part-2, w = class_n_rot(pair_c1(p)) + class_n_rot(pair_c2(p));
            v.msg[p][(1-pe)*v.msg_stride[p] + size_t(ns)*v.cap[p]*w + i] = 1.f;
        }
    });

    bp_new_beliefs(v, ns, pn, pe, 0.f, true);
    float max_deviation = 1e10f;
    for(int iter=0; max_deviation>v.tol && iter<v.max_iter; iter+=v.chunk) {
        for(int j=0; j<v.chunk; ++j) {
            pn = 1-pn;
            pe = 1-pe;
            bp_new_beliefs(v, ns, pn, pe, v.damping, false);
        }
        // the damped update moves each belief by (1-damping) times the residual of the undamped update
        max_deviation = bp_deviation(v, ns, pn) / (1.f-v.damping);
    }

    block_for_parts<2>({n0, n1}, [&](int part, int i) {
        if(part==0) bp_node_marginals<0>(v, ns, pn, i);
        else        bp_node_marginals<1>(v, ns, pn, i);
    });
    block_for_parts<3>({e0, e1, e2}, [&](int part, int i) {
        if     (part==0) bp_edge_marginals<0>(v, ns, pe, i);
        else if(part==1) bp_edge_marginals<1>(v, ns, pe, i);
        else             bp_edge_marginals<2>(v, ns, pe, i);
    });
}


struct DeviceRotamer : public DeviceNode {
    AccumView potential;
    int n_system;
    DeviceInteractionGraph<PosQuadSplineInteraction> igraph;
    int n_onebody;
    VecArray onebody_value[4];
    AccumView onebody_sens[4];
    int onebody_width[4];
    float damping, tol;
    int max_iter, chunk;

    DeviceArray<int> bead_class, bead_node, bead_rot;  // per bead; bead_node within its class
    vector<int> h_start[3];          // per class, each system's first node, then the class's count
    DeviceArray<int> start[3];
    int n_node[3];
    DeviceAccum node_energy[3];      // per node and rotamer
    DeviceArray<float> node_prob[3], node_offset[3], node_belief[3], node_marginal[3];
    DeviceArray<int> node_system[3];

    int cap[n_pair];                 // residue pairs per system
    DeviceArray<int> lookup_flag[n_pair], lookup_slot[n_pair], lookup_start[n_pair];  // a system's n1 x n2 lookup
    DeviceArray<int> n_edge[n_pair], node1[n_pair], node2[n_pair];
    DeviceAccum edge_energy[n_pair];
    DeviceArray<float> edge_prob[n_pair], edge_offset[n_pair], edge_msg[n_pair], edge_marg[n_pair];
    int overflow_code;

    DeviceArray<int> inc_count, inc_start, inc_pair, inc_slot;  // the BP nodes' edges (class 1's nodes first)
    DeviceArray<int> pair_kind, pair_target, pair_dim;  // per bead edge slot, for its sensitivity

    DeviceRotamer(DeviceGraph& graph, AccumView potential_, const DeviceInteractionData& ig,
            const vector<DeviceCoord*>& onebody, const vector<int>& bead_n_rot, const vector<int>& bead_residue,
            const vector<int>& bead_rot_, const vector<int>& n_residue, float damping_, int max_iter_, float tol_,
            int chunk_):
        potential(potential_), n_system(graph.n_system), igraph(graph, "rotamer", ig), n_onebody(onebody.size()),
        damping(damping_), tol(tol_), max_iter(max_iter_), chunk(chunk_),
        overflow_code(graph.add_check("rotamer's residue pairs exceed their capacity"))
    {
        if(n_onebody > 4) throw string("the device rotamer node takes at most 4 1-body energies");
        for(int i=0; i<n_onebody; ++i) {
            if(onebody[i]->offset != ig.pos1->offset)
                throw string("rotamer's 1-body energies need one element per bead in every system");
            onebody_value[i] = onebody[i]->value();
            onebody_sens[i]  = onebody[i]->sens.view();
            onebody_width[i] = onebody[i]->row_width;
        }

        // nodes: per class, each system's residues in turn
        for(int c=0; c<3; ++c) {
            h_start[c].assign(1, 0);
            vector<int> sys;
            for(int ns=0; ns<n_system; ++ns) {
                h_start[c].push_back(h_start[c].back() + n_residue[3*ns+c]);
                sys.insert(sys.end(), n_residue[3*ns+c], ns);
            }
            n_node[c] = h_start[c].back();
            int nr = class_n_rot(c);
            start[c] = DeviceArray<int>(h_start[c]);
            node_energy[c] = DeviceAccum(size_t(n_node[c])*nr);
            node_prob[c] = DeviceArray<float>(size_t(n_node[c])*nr);
            node_offset[c] = DeviceArray<float>(n_node[c]);
            node_belief[c] = DeviceArray<float>(2*size_t(n_node[c])*nr);
            node_marginal[c] = DeviceArray<float>(size_t(n_node[c])*nr);
            node_system[c] = DeviceArray<int>(sys);
        }

        // beads
        vector<int> cls, node;
        int bead = 0;
        for(int ns=0; ns<n_system; ++ns)
            for(int i=0; i<ig.n_elem1[ns]; ++i, ++bead) {
                int nr = bead_n_rot[bead];
                int c = nr==1 ? 0 : nr==3 ? 1 : nr==6 ? 2 : -1;
                if(c<0 || bead_residue[bead]<0 || bead_residue[bead]>=n_residue[3*ns+c] || bead_rot_[bead]>=nr)
                    throw string("invalid rotamer id for the device rotamer node");
                cls.push_back(c);
                node.push_back(h_start[c][ns] + bead_residue[bead]);
            }
        bead_class = DeviceArray<int>(cls);
        bead_node = DeviceArray<int>(node);
        bead_rot = DeviceArray<int>(bead_rot_);

        // residue pairs: a dense lookup per system, and as many pairs as the lookup or the bead pairs allow
        int n_inc = 0;
        for(int p=0; p<n_pair; ++p) {
            int c1 = pair_c1(p), c2 = pair_c2(p);
            vector<int> ls(1, 0);
            cap[p] = 0;
            for(int ns=0; ns<n_system; ++ns) {
                int n1 = n_residue[3*ns+c1], n2 = n_residue[3*ns+c2];
                ls.push_back(ls.back() + n1*n2);
                cap[p] = max(cap[p], min(n1*n2, ig.max_n_edge));
            }
            int w = class_n_rot(c1)*class_n_rot(c2);
            lookup_start[p] = DeviceArray<int>(ls);
            lookup_flag[p] = DeviceArray<int>(ls.back());
            lookup_slot[p] = DeviceArray<int>(ls.back());
            n_edge[p] = DeviceArray<int>(n_system);
            n_edge[p].zero();
            node1[p] = DeviceArray<int>(size_t(n_system)*cap[p]);
            node2[p] = DeviceArray<int>(size_t(n_system)*cap[p]);
            edge_energy[p] = DeviceAccum(size_t(n_system)*cap[p]*w);
            edge_prob[p] = DeviceArray<float>(size_t(n_system)*cap[p]*w);
            edge_offset[p] = DeviceArray<float>(size_t(n_system)*cap[p]);
            edge_msg[p] = DeviceArray<float>(2*size_t(n_system)*cap[p]*(class_n_rot(c1)+class_n_rot(c2)));
            edge_marg[p] = DeviceArray<float>(size_t(n_system)*cap[p]*w);
            n_inc += 2*n_system*cap[p];
        }
        inc_count = DeviceArray<int>(n_node[1]+n_node[2]);
        inc_start = DeviceArray<int>(n_node[1]+n_node[2]+1);
        inc_pair = DeviceArray<int>(n_inc);
        inc_slot = DeviceArray<int>(n_inc);
        pair_kind   = DeviceArray<int>(igraph.n_slot());
        pair_target = DeviceArray<int>(igraph.n_slot());
        pair_dim    = DeviceArray<int>(igraph.n_slot());
    }

    void compute_value(ComputeMode mode) override;
};

void DeviceRotamer::compute_value(ComputeMode mode) {
    bool want_pot = mode == PotentialAndDerivMode;
    auto pot = potential;
    int nsys = n_system;

    // fill_holders: bead pairs and the 1-body energies
    igraph.compute_edges();
    for(int c=0; c<3; ++c) node_energy[c].zero();
    for(int p=0; p<n_pair; ++p) {edge_energy[p].zero(); lookup_flag[p].zero();}

    const int *bc = bead_class.ptr, *bn = bead_node.ptr, *br = bead_rot.ptr;
    const int* loc = igraph.loc1.ptr;
    AccumView ne0 = node_energy[0].view(), ne1 = node_energy[1].view(), ne2 = node_energy[2].view();
    int nob = n_onebody;
    VecArray ob0 = onebody_value[0], ob1 = onebody_value[1], ob2 = onebody_value[2], ob3 = onebody_value[3];
    parallel_for(igraph.n_elem1, DEVICE_LAMBDA (int n) {
        float energy = 0.f;
        if(nob>0) energy += ob0(0,loc[n]);
        if(nob>1) energy += ob1(0,loc[n]);
        if(nob>2) energy += ob2(0,loc[n]);
        if(nob>3) energy += ob3(0,loc[n]);
        int nr = class_n_rot(bc[n]);
        (bc[n]==0 ? ne0 : bc[n]==1 ? ne1 : ne2).add(bn[n]*nr + br[n], energy);
    });

    // classify each bead pair: 11 to the energy, 1x to the other residue's rotamers, BP pairs marked in the lookup
    const int *e1 = igraph.pairlist.edge1.ptr, *e2 = igraph.pairlist.edge2.ptr, *ige = igraph.pairlist.n_edge.ptr;
    const float* val = igraph.edge_value.ptr;
    int icap = igraph.max_n_edge;
    const int *st1 = start[1].ptr, *st2 = start[2].ptr;
    int *lf0 = lookup_flag[0].ptr, *lf1 = lookup_flag[1].ptr, *lf2 = lookup_flag[2].ptr;
    const int *ls0 = lookup_start[0].ptr, *ls1 = lookup_start[1].ptr, *ls2 = lookup_start[2].ptr;
    int *pk = pair_kind.ptr, *pt = pair_target.ptr, *pd = pair_dim.ptr;
    slots_for(igraph.n_system, icap, ige, DEVICE_LAMBDA (int e) {
        int ns = e/icap;
        int a = e1[e], b = e2[e];
        if(class_n_rot(bc[a]) > class_n_rot(bc[b])) {int t = a; a = b; b = t;}
        int c1 = bc[a], c2 = bc[b], r1 = br[a], r2 = br[b];
        if(c1==0 && c2==0) {
            pk[e] = 0;
            if(want_pot) pot.add(ns, val[e]);
        } else if(c1==0) {
            pk[e] = 1; pt[e] = bn[b]; pd[e] = c2*8 + r2;
            (c2==1 ? ne1 : ne2).add(bn[b]*class_n_rot(c2) + r2, val[e]);
        } else {
            int p = (c1==1 && c2==1) ? 0 : (c1==1 ? 1 : 2);
            const int* sa = c1==1 ? st1 : st2;
            const int* sb = c2==1 ? st1 : st2;
            int l1 = bn[a]-sa[ns], l2 = bn[b]-sb[ns], n2 = sb[ns+1]-sb[ns];
            int idx = (p==0 ? ls0 : p==1 ? ls1 : ls2)[ns] + l1*n2 + l2;
            (p==0 ? lf0 : p==1 ? lf1 : lf2)[idx] = 1;
            pk[e] = 2+p; pt[e] = idx; pd[e] = r1*class_n_rot(c2) + r2;
        }
    });

    // number each system's residue pairs in lookup order
    int* err = igraph.pairlist.error;
    int code = overflow_code;
    for(int p=0; p<n_pair; ++p) {
        const int* sa = (pair_c1(p)==1 ? start[1] : start[2]).ptr;
        const int* sb = (pair_c2(p)==1 ? start[1] : start[2]).ptr;
        const int* ls = lookup_start[p].ptr;
        const int* lf = lookup_flag[p].ptr;
        int *lsl = lookup_slot[p].ptr, *nd1 = node1[p].ptr, *nd2 = node2[p].ptr, *ne = n_edge[p].ptr;
        int cp = cap[p];
        per_system(nsys, DEVICE_LAMBDA (int ns) {
            int n2 = sb[ns+1]-sb[ns];
            int o = ls[ns];
            int total = block_scan(ls[ns+1]-o, [=] (int i) {return lf[o+i];},
                    [=] (int i, int prefix, int v) {
                        lsl[o+i] = v ? prefix : -1;
                        if(v && prefix < cp) {
                            nd1[size_t(ns)*cp+prefix] = sa[ns] + i/n2;
                            nd2[size_t(ns)*cp+prefix] = sb[ns] + i%n2;
                        }});
            if(total > cp) *err = code;
            ne[ns] = total < cp ? total : cp;
        });
    }

    // the BP pairs' energies, by bead pair
    {
        const int *lsl0 = lookup_slot[0].ptr, *lsl1 = lookup_slot[1].ptr, *lsl2 = lookup_slot[2].ptr;
        AccumView ee0 = edge_energy[0].view(), ee1 = edge_energy[1].view(), ee2 = edge_energy[2].view();
        int cp0 = cap[0], cp1 = cap[1], cp2 = cap[2];
        slots_for(igraph.n_system, icap, ige, DEVICE_LAMBDA (int e) {
            int ns = e/icap;
            if(pk[e] < 2) return;
            int p = pk[e]-2;
            int k = (p==0 ? lsl0 : p==1 ? lsl1 : lsl2)[pt[e]];
            int cp = p==0 ? cp0 : p==1 ? cp1 : cp2;
            if(k >= cp) return;  // a capacity overflow, refused after this evaluation
            size_t slot = size_t(ns)*cp + k;
            int w = class_n_rot(pair_c1(p))*class_n_rot(pair_c2(p));
            (p==0 ? ee0 : p==1 ? ee1 : ee2).add(slot*w + pd[e], val[e]);
            pt[e] = int(slot);
        });
    }

    // energies to probabilities, each node's and pair's lowest energy taken out as its offset
    for(int c=0; c<3; ++c) {
        int nr = class_n_rot(c);
        auto en = node_energy[c].view();
        float *pr = node_prob[c].ptr, *off = node_offset[c].ptr;
        parallel_for(n_node[c], DEVICE_LAMBDA (int i) {
            float e_offset = en.get(i*nr);
            for(int d=1; d<nr; ++d) {float x = en.get(i*nr+d); e_offset = x < e_offset ? x : e_offset;}
            for(int d=0; d<nr; ++d) pr[i*nr+d] = expf(e_offset-en.get(i*nr+d));
            off[i] = e_offset;
        });
    }
    for(int p=0; p<n_pair; ++p) {
        int w = class_n_rot(pair_c1(p))*class_n_rot(pair_c2(p)), cp = cap[p];
        auto en = edge_energy[p].view();
        float *pr = edge_prob[p].ptr, *off = edge_offset[p].ptr;
        const int* ne = n_edge[p].ptr;
        slots_for(nsys, cp, ne, DEVICE_LAMBDA (int s) {
            float e_offset = en.get(size_t(s)*w);
            for(int d=1; d<w; ++d) {float x = en.get(size_t(s)*w+d); e_offset = x < e_offset ? x : e_offset;}
            for(int d=0; d<w; ++d) pr[size_t(s)*w+d] = expf(e_offset-en.get(size_t(s)*w+d));
            off[s] = e_offset;
        });
    }

    // each BP node's pairs in the CPU's order: 33, 36, 66, each in lookup (slot) order
    int nn1 = n_node[1], nn2 = n_node[2];
    const int *nsys1 = node_system[1].ptr, *nsys2 = node_system[2].ptr;
    const int *lsl0 = lookup_slot[0].ptr, *lsl1 = lookup_slot[1].ptr, *lsl2 = lookup_slot[2].ptr;
    int cp0 = cap[0], cp1 = cap[1], cp2 = cap[2];
    int *ic = inc_count.ptr, *is = inc_start.ptr, *ip = inc_pair.ptr, *isl = inc_slot.ptr;
    // node f's pairs in order, written as pair*2+side and slot from out_pair and out_slot (if given); returns their
    // count, on every lane
    auto for_pairs = [=] UPSIDE_HD (int f, const Lanes& lanes, int* out_pair, int* out_slot) {
        int k = 0;
        auto emit = [&](int p, int side, int k_slot, int ns, int cp) {
            int slot = lanes.place(k_slot >= 0, k);
            if(k_slot >= 0 && out_pair) {out_pair[slot] = 2*p+side; out_slot[slot] = ns*cp + k_slot;}
        };
        // the pairs of node l of pair p's lookup, in row order: as node 2 (column l) of rows a < l, as node 1 (all of
        // row l), as node 2 of rows a > l; a node of only one side has the whole column or only its row
        auto scan = [&](int p, int ns, int l, bool as1, bool as2, int n1, int n2, int o, const int* lsl, int cp) {
            auto column = [&](int a0, int a1) {
                lanes.chunks(a0, a1, [&](int a, bool in) {emit(p, 1, in ? lsl[o+a*n2+l] : -1, ns, cp);});
            };
            auto row = [&]() {
                lanes.chunks(0, n2, [&](int b, bool in) {emit(p, 0, in ? lsl[o+l*n2+b] : -1, ns, cp);});
            };
            if(as1 && as2) {column(0, l); row(); column(l+1, n1);}
            else if(as1)   row();
            else           column(0, n1);
        };
        if(f < nn1) {
            int ns = nsys1[f], l = f-st1[ns], n3 = st1[ns+1]-st1[ns], n6 = st2[ns+1]-st2[ns];
            scan(0, ns, l, true, true,  n3, n3, ls0[ns], lsl0, cp0);
            scan(1, ns, l, true, false, n3, n6, ls1[ns], lsl1, cp1);
        } else {
            int node = f-nn1;
            int ns = nsys2[node], l = node-st2[ns], n3 = st1[ns+1]-st1[ns], n6 = st2[ns+1]-st2[ns];
            scan(1, ns, l, false, true, n3, n6, ls1[ns], lsl1, cp1);
            scan(2, ns, l, true, true,  n6, n6, ls2[ns], lsl2, cp2);
        }
        return k;
    };
    lanes_for(nn1+nn2, DEVICE_LAMBDA (int f, Lanes lanes) {ic[f] = for_pairs(f, lanes, nullptr, nullptr);});
    per_system(1, DEVICE_LAMBDA (int ns) {
        int total = block_scan(nn1+nn2, [=] (int f) {return ic[f];}, [=] (int f, int prefix, int v) {is[f] = prefix;});
        is[nn1+nn2] = total;  // every thread writes the same total
    });
    lanes_for(nn1+nn2, DEVICE_LAMBDA (int f, Lanes lanes) {for_pairs(f, lanes, ip+is[f], isl+is[f]);});

    // BP
    BPView v;
    for(int k=0; k<2; ++k) {
        int c = k+1, nr = class_n_rot(c);
        v.prob[k] = node_prob[c].ptr; v.belief[k] = node_belief[c].ptr; v.marginal[k] = node_marginal[c].ptr;
        v.start[k] = start[c].ptr; v.belief_stride[k] = size_t(n_node[c])*nr;
    }
    for(int p=0; p<n_pair; ++p) {
        v.P[p] = edge_prob[p].ptr; v.msg[p] = edge_msg[p].ptr; v.edge_marginal[p] = edge_marg[p].ptr;
        v.msg_stride[p] = size_t(nsys)*cap[p]*(class_n_rot(pair_c1(p))+class_n_rot(pair_c2(p)));
        v.node1[p] = node1[p].ptr; v.node2[p] = node2[p].ptr; v.n_edge[p] = n_edge[p].ptr; v.cap[p] = cap[p];
    }
    v.inc_start = inc_start.ptr; v.inc_pair = inc_pair.ptr; v.inc_slot = inc_slot.ptr; v.n_node1 = nn1;
    v.damping = damping; v.tol = tol; v.max_iter = max_iter; v.chunk = chunk;
    per_system(nsys, DEVICE_LAMBDA (int ns) {bp_solve(v, ns);});

    // the energy from the marginals: the nodes' free energies, the 11 pairs' energies (added above) and the BP
    // pairs' free energies
    if(want_pot) {
        for(int c=0; c<3; ++c) {
            int nr = class_n_rot(c);
            const float *pr = node_prob[c].ptr, *off = node_offset[c].ptr, *mg = node_marginal[c].ptr;
            const int* sys = node_system[c].ptr;
            bool one = c==0;
            parallel_for(n_node[c], DEVICE_LAMBDA (int i) {
                float b[6], s = 0.f;
                for(int r=0; r<nr; ++r) {b[r] = one ? 1.f : mg[i*nr+r]; s += b[r];}
                float scale = rcp(s);
                float en = off[i];
                // free energy is average energy - entropy
                for(int r=0; r<nr; ++r) {
                    float br_ = b[r]*scale;
                    en += br_ * logf((1e-30f+br_)*rcp(1e-30f+pr[i*nr+r]));
                }
                pot.add(sys[i], en);
            });
        }
        for(int p=0; p<n_pair; ++p) {
            int c1 = pair_c1(p), c2 = pair_c2(p), n1 = class_n_rot(c1), n2 = class_n_rot(c2), cp = cap[p];
            const float *pr = edge_prob[p].ptr, *off = edge_offset[p].ptr, *mg = edge_marg[p].ptr;
            const float *m1 = node_marginal[c1].ptr, *m2 = node_marginal[c2].ptr;
            const int *nd1 = node1[p].ptr, *nd2 = node2[p].ptr, *ne = n_edge[p].ptr;
            slots_for(nsys, cp, ne, DEVICE_LAMBDA (int s) {
                int ns = s/cp;
                const float *b1 = m1 + nd1[s]*n1, *b2 = m2 + nd2[s]*n2;
                float en = off[s];
                for(int i=0; i<n1; ++i)
                    for(int j=0; j<n2; ++j) {
                        float q = mg[size_t(s)*n1*n2+i*n2+j];
                        // the average potential energy plus the mutual information; 1e-30f keeps exact zeros finite
                        en += q * logf((1e-30f+q) * rcp(1e-30f+pr[size_t(s)*n1*n2+i*n2+j]*b1[i]*b2[j]));
                    }
                pot.add(ns, en);
            });
        }
    }

    // derivatives by the envelope theorem: each bead pair weighted by its rotamer pair's marginal, each bead's 1-body
    // energy by its rotamer's
    const float *mg1 = node_marginal[1].ptr, *mg2 = node_marginal[2].ptr;
    const float *em0 = edge_marg[0].ptr, *em1 = edge_marg[1].ptr, *em2 = edge_marg[2].ptr;
    float* es = igraph.edge_sensitivity.ptr;
    slots_for(igraph.n_system, icap, ige, DEVICE_LAMBDA (int e) {
        int kind = pk[e];
        if(kind == 0)      es[e] = 1.f;
        else if(kind == 1) {int c2 = pd[e]/8, r2 = pd[e]%8; es[e] = (c2==1 ? mg1 : mg2)[pt[e]*class_n_rot(c2) + r2];}
        else {
            int p = kind-2, w = class_n_rot(pair_c1(p))*class_n_rot(pair_c2(p));
            es[e] = (p==0 ? em0 : p==1 ? em1 : em2)[size_t(pt[e])*w + pd[e]];
        }
    });
    igraph.propagate_derivatives();

    AccumView os0 = onebody_sens[0], os1 = onebody_sens[1], os2 = onebody_sens[2], os3 = onebody_sens[3];
    int ow0 = onebody_width[0], ow1 = onebody_width[1], ow2 = onebody_width[2], ow3 = onebody_width[3];
    parallel_for(igraph.n_elem1, DEVICE_LAMBDA (int n) {
        int c = bc[n], nr = class_n_rot(c);
        float b = c==0 ? 1.f : (c==1 ? mg1 : mg2)[bn[n]*nr + br[n]];
        if(nob>0) os0.add(loc[n]*ow0, b);
        if(nob>1) os1.add(loc[n]*ow1, b);
        if(nob>2) os2.add(loc[n]*ow2, b);
        if(nob>3) os3.add(loc[n]*ow3, b);
    });
}

}

DeviceNode* device_rotamer(DeviceGraph& graph, AccumView potential, const DeviceInteractionData& ig,
        const vector<DeviceCoord*>& onebody, const vector<int>& bead_n_rot, const vector<int>& bead_residue,
        const vector<int>& bead_rot, const vector<int>& n_residue, float damping, int max_iter, float tol,
        int iteration_chunk_size) {
    return new DeviceRotamer(graph, potential, ig, onebody, bead_n_rot, bead_residue, bead_rot, n_residue, damping,
            max_iter, tol, iteration_chunk_size);
}
