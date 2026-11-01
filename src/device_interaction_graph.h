#ifndef DEVICE_INTERACTION_GRAPH_H
#define DEVICE_INTERACTION_GRAPH_H

// The device counterparts of interaction_graph.h's PairlistComputation and InteractionGraph, for every system of a
// DeviceGraph; included by the device files only.  Elements pair only within their own system.  As on the CPU, a
// system's cache of pairs within cutoff + cache_buffer is rebuilt when one of its elements has moved more than
// half the buffer since the last rebuild, and each evaluation cuts the cache to the pairs within the cutoff.  The
// rebuild and the cut keep the pairs ordered (by first element, then second), so the result does not depend on the
// order the threads run in.  System s's pairs sit in the slots [s*max_n_edge, s*max_n_edge + n_edge[s]).

#include "device_nodes.h"

// A pair list over two sets of elements (one set paired with itself if symmetric).  Accept is a type with a
// device-callable static bool ok(int sequence_id1, int sequence_id2), the exclusion by sequence ids.
template <bool symmetric>
struct DevicePairlist {
    int n_system, n_elem1, n_elem2;
    int max_n_edge;  // per system
    float cutoff, cache_cutoff;
    DeviceArray<int> start1, start2;      // each system's first element in each set, then the set's size
    DeviceArray<int> system1, system2;    // each element's system
    DeviceArray<int> sequence_id1, sequence_id2;
    DeviceArray<float> cache_pos1, cache_pos2;   // the positions at the last rebuild, 3 per element
    DeviceArray<int> rebuild;             // per system, this evaluation
    DeviceArray<int> count, base;         // per element of set 1: its cache pairs, and their first slot
    DeviceArray<int> cache_n_edge, n_edge;
    DeviceArray<int> cache_edge1, cache_edge2, edge1, edge2;  // element indices, n_system*max_n_edge
    int* error;
    int overflow_code;
    bool first;

    DevicePairlist(DeviceGraph& graph, const std::string& name, const std::vector<int>& n_elem1_of_system,
            const std::vector<int>& n_elem2_of_system, const std::vector<int>& sequence_id1_,
            const std::vector<int>& sequence_id2_, float cutoff_, float cache_buffer, int max_n_edge_):
        n_system(graph.n_system), n_elem1(sequence_id1_.size()), n_elem2(sequence_id2_.size()),
        max_n_edge(max_n_edge_), cutoff(cutoff_), cache_cutoff(cutoff_+cache_buffer),
        start1(starts(n_elem1_of_system)), start2(starts(n_elem2_of_system)),
        system1(systems(n_elem1_of_system)), system2(systems(n_elem2_of_system)),
        sequence_id1(sequence_id1_), sequence_id2(sequence_id2_),
        cache_pos1(3*n_elem1), cache_pos2(3*n_elem2), rebuild(n_system),
        count(n_elem1), base(n_elem1), cache_n_edge(n_system), n_edge(n_system),
        cache_edge1(size_t(n_system)*max_n_edge), cache_edge2(size_t(n_system)*max_n_edge),
        edge1(size_t(n_system)*max_n_edge), edge2(size_t(n_system)*max_n_edge),
        error(graph.error.ptr),
        overflow_code(graph.add_check(name + "'s pair list needs more than its max_n_edge of " +
                    std::to_string(max_n_edge_))),
        first(true)
    {
        cache_n_edge.zero();
        n_edge.zero();
    }

    static std::vector<int> starts(const std::vector<int>& n) {
        std::vector<int> r(1, 0);
        for(int x: n) r.push_back(r.back()+x);
        return r;
    }
    static std::vector<int> systems(const std::vector<int>& n) {
        std::vector<int> r;
        for(int ns=0; ns<int(n.size()); ++ns) r.insert(r.end(), n[ns], ns);
        return r;
    }

    //! The pairs of elements of x1 and x2 (each element's position in its first three of stride floats)
    template <typename Accept>
    void find_edges(const float* x1, int stride1, const float* x2, int stride2) {
        if(symmetric) {x2 = x1; stride2 = stride1;}
        int* rb = rebuild.ptr;
        const int *s1 = system1.ptr, *s2 = system2.ptr;
        float *c1 = cache_pos1.ptr, *c2 = cache_pos2.ptr;

        // which systems must rebuild: every one at the first evaluation, then those with an element that moved
        // more than half the buffer
        if(first) {
            std::vector<int> all(n_system, 1);
            rebuild.upload(all);
            first = false;
        } else {
            rebuild.zero();
            float max_cache_dist2 = sqr(0.5f*(cache_cutoff - cutoff));
            parallel_for(n_elem1, DEVICE_LAMBDA (int i) {
                float d2 = sqr(x1[i*stride1+0]-c1[3*i+0]) + sqr(x1[i*stride1+1]-c1[3*i+1]) +
                           sqr(x1[i*stride1+2]-c1[3*i+2]);
                if(max_cache_dist2 < d2) rb[s1[i]] = 1;
            });
            if(!symmetric)
                parallel_for(n_elem2, DEVICE_LAMBDA (int i) {
                    float d2 = sqr(x2[i*stride2+0]-c2[3*i+0]) + sqr(x2[i*stride2+1]-c2[3*i+1]) +
                               sqr(x2[i*stride2+2]-c2[3*i+2]);
                    if(max_cache_dist2 < d2) rb[s2[i]] = 1;
                });
        }

        // rebuild: count each element's pairs within the cache cutoff, place them by a scan within the system,
        // and write them in order
        const int *st1 = start1.ptr, *st2 = start2.ptr;
        const int *q1 = sequence_id1.ptr, *q2 = symmetric ? sequence_id1.ptr : sequence_id2.ptr;
        int *cnt = count.ptr, *bs = base.ptr, *cne = cache_n_edge.ptr, *ce1 = cache_edge1.ptr, *ce2 = cache_edge2.ptr;
        float cache_cutoff2 = sqr(cache_cutoff);
        int cap = max_n_edge;
        int* err = error;
        int code = overflow_code;
        // the partners of element i in its system, within the cache cutoff and accepted, written in order to the
        // slots from k (if out1); returns k advanced past them, on every lane
        auto for_partners = [=] UPSIDE_HD (int i, const Lanes& l, int* out1, int* out2, int k) {
            int ns = s1[i];
            int j0 = symmetric ? i+1 : st2[ns];
            int j1 = symmetric ? st1[ns+1] : st2[ns+1];
            l.chunks(j0, j1, [&](int j, bool in) {
                bool keep = false;
                if(in) {
                    float d2 = sqr(x1[i*stride1+0]-x2[j*stride2+0]) + sqr(x1[i*stride1+1]-x2[j*stride2+1]) +
                               sqr(x1[i*stride1+2]-x2[j*stride2+2]);
                    keep = d2 < cache_cutoff2 && Accept::ok(q1[i], q2[j]);
                }
                int slot = l.place(keep, k);
                if(keep && out1 && slot<cap) {out1[slot] = i; out2[slot] = j;}
            });
            return k;
        };
        lanes_for(n_elem1, DEVICE_LAMBDA (int i, Lanes l) {
            if(rb[s1[i]]) cnt[i] = for_partners(i, l, nullptr, nullptr, 0);
        });
        per_system(n_system, DEVICE_LAMBDA (int ns) {
            if(!rb[ns]) return;
            int total = block_scan(st1[ns+1]-st1[ns], [=] (int i) {return cnt[st1[ns]+i];},
                    [=] (int i, int prefix, int v) {bs[st1[ns]+i] = prefix;});
            if(total > cap) *err = code;
            cne[ns] = total < cap ? total : cap;
        });
        lanes_for(n_elem1, DEVICE_LAMBDA (int i, Lanes l) {
            int ns = s1[i];
            if(rb[ns]) for_partners(i, l, ce1 + size_t(ns)*cap, ce2 + size_t(ns)*cap, bs[i]);
        });
        parallel_for(n_elem1, DEVICE_LAMBDA (int i) {
            if(rb[s1[i]]) for(int d=0; d<3; ++d) c1[3*i+d] = x1[i*stride1+d];
        });
        if(!symmetric)
            parallel_for(n_elem2, DEVICE_LAMBDA (int i) {
                if(rb[s2[i]]) for(int d=0; d<3; ++d) c2[3*i+d] = x2[i*stride2+d];
            });

        // the cache's pairs within the cutoff, in the cache's order
        int *ne = n_edge.ptr, *e1 = edge1.ptr, *e2 = edge2.ptr;
        float cutoff2 = sqr(cutoff);
        per_system(n_system, DEVICE_LAMBDA (int ns) {
            size_t o = size_t(ns)*cap;
            int total = block_scan(cne[ns],
                    [=] (int k) {
                        int i = ce1[o+k], j = ce2[o+k];
                        float d2 = sqr(x1[i*stride1+0]-x2[j*stride2+0]) + sqr(x1[i*stride1+1]-x2[j*stride2+1]) +
                                   sqr(x1[i*stride1+2]-x2[j*stride2+2]);
                        return int(d2 < cutoff2);},
                    [=] (int k, int prefix, int v) {
                        if(v) {e1[o+prefix] = ce1[o+k]; e2[o+prefix] = ce2[o+k];}});
            ne[ns] = total;
        });
    }
};


//! The device's InteractionGraph<IType>, built from a DeviceInteractionData (device_nodes.h): after compute_edges,
//! edge slot e of a system with k = e % max_n_edge < n_edge holds elements edge1[e] and edge2[e] (indices into the
//! graph's elements), its value and its derivatives by both elements' coordinates; propagate_derivatives adds each
//! edge's edge_sensitivity times them to the inputs' sensitivities.
template <typename IType>
struct DeviceInteractionGraph {
    constexpr static bool symmetric = IType::symmetric;
    constexpr static int  n_dim1 = IType::n_dim1, n_dim2 = IType::n_dim2, n_param = IType::n_param;
    constexpr static int  n_dim  = n_dim1 + n_dim2;

    struct Accept {UPSIDE_HD static bool ok(int a, int b) {return IType::template acceptable_id_pair<int>(a, b);}};

    int n_system;
    DeviceCoord &pos1, &pos2;
    int n_elem1, n_elem2, n_type2, max_n_edge;
    DeviceArray<int> loc1, loc2, type1, type2, id1, id2;
    DeviceArray<float> param;
    DeviceArray<float> x1, x2;   // the elements' coordinates, n_dim1 and n_dim2 per element
    DevicePairlist<symmetric> pairlist;
    DeviceArray<float> edge_value, edge_deriv, edge_sensitivity;  // per slot; n_dim per slot for edge_deriv

    DeviceInteractionGraph(DeviceGraph& graph, const std::string& name, const DeviceInteractionData& d):
        n_system(graph.n_system), pos1(*d.pos1), pos2(symmetric ? *d.pos1 : *d.pos2),
        n_elem1(d.loc1.size()), n_elem2(symmetric ? d.loc1.size() : d.loc2.size()),
        n_type2(d.n_type2), max_n_edge(d.max_n_edge),
        loc1(d.loc1), loc2(symmetric ? d.loc1 : d.loc2), type1(d.type1), type2(symmetric ? d.type1 : d.type2),
        id1(d.id1), id2(symmetric ? d.id1 : d.id2), param(d.interaction_param),
        x1(size_t(n_elem1)*n_dim1), x2(symmetric ? 0 : size_t(n_elem2)*n_dim2),
        pairlist(graph, name, d.n_elem1, symmetric ? d.n_elem1 : d.n_elem2, d.sequence_id1,
                symmetric ? d.sequence_id1 : d.sequence_id2, d.cutoff, d.cache_buffer, d.max_n_edge),
        edge_value(size_t(n_system)*max_n_edge), edge_deriv(size_t(n_system)*max_n_edge*n_dim),
        edge_sensitivity(size_t(n_system)*max_n_edge)
    {
        edge_sensitivity.zero();
    }

    int n_slot() const {return n_system*max_n_edge;}

    void compute_edges() {
        gather(pos1, loc1, x1, n_elem1, n_dim1);
        if(!symmetric) gather(pos2, loc2, x2, n_elem2, n_dim2);
        const float* c2 = symmetric ? x1.ptr : x2.ptr;
        pairlist.template find_edges<Accept>(x1.ptr, n_dim1, c2, symmetric ? n_dim1 : n_dim2);

        const float* c1 = x1.ptr;
        const int *e1 = pairlist.edge1.ptr, *e2 = pairlist.edge2.ptr, *ne = pairlist.n_edge.ptr;
        const int *t1 = type1.ptr, *t2 = type2.ptr;
        const float* p = param.ptr;
        float *val = edge_value.ptr, *der = edge_deriv.ptr;
        int cap = max_n_edge, nt2 = n_type2;
        slots_for(n_system, cap, ne, DEVICE_LAMBDA (int e) {
            int i1 = e1[e], i2 = e2[e];
            Vec<n_dim1> a = load_vec<n_dim1>(c1 + size_t(i1)*n_dim1);
            Vec<n_dim2> b = load_vec<n_dim2>(c2 + size_t(i2)*n_dim2);
            Vec<n_dim1> d1;
            Vec<n_dim2> d2;
            EdgeParams<float> ep = {p + (t1[i1]*nt2 + t2[i2])*n_param};
            val[e] = IType::template compute_edge<float>(d1, d2, ep, a, b);
            store_vec(der + size_t(e)*n_dim,          d1);
            store_vec(der + size_t(e)*n_dim + n_dim1, d2);
        });
    }

    void propagate_derivatives() {
        const int *e1 = pairlist.edge1.ptr, *e2 = pairlist.edge2.ptr, *ne = pairlist.n_edge.ptr;
        const int *l1 = loc1.ptr, *l2 = loc2.ptr;
        const float *der = edge_deriv.ptr, *sens = edge_sensitivity.ptr;
        auto s1 = pos1.sens.view(), s2 = pos2.sens.view();
        int w1 = pos1.row_width, w2 = pos2.row_width, cap = max_n_edge;
        slots_for(n_system, cap, ne, DEVICE_LAMBDA (int e) {
            float s = sens[e];
            for(int d=0; d<n_dim1; ++d) s1.add(l1[e1[e]]*w1+d, s*der[size_t(e)*n_dim+d]);
            for(int d=0; d<n_dim2; ++d) s2.add(l2[e2[e]]*w2+d, s*der[size_t(e)*n_dim+n_dim1+d]);
        });
    }

    // public because a CUDA device lambda cannot sit in a private member function
    static void gather(DeviceCoord& pos, DeviceArray<int>& loc, DeviceArray<float>& x, int n, int nd) {
        VecArray v = pos.value();
        const int* l = loc.ptr;
        float* xp = x.ptr;
        parallel_for(n, DEVICE_LAMBDA (int i) {for(int d=0; d<nd; ++d) xp[size_t(i)*nd+d] = v(d, l[i]);});
    }
};

#endif
