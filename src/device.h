#ifndef DEVICE_H
#define DEVICE_H

// The device backend: memory, kernel launches and force accumulation for the device engine (device_engine.h).
// Compiled by nvcc (an UPSIDE_CUDA build) the device files run on the GPU.  Compiled as C++ they run the same code
// serially on the host, which is how the device path is tested on a computer without a GPU.

#include "vector_math.h"
#include <vector>
#include <string>
#include <cstdint>

#ifdef __CUDACC__
#define UPSIDE_DEVICE __device__
#define DEVICE_LAMBDA [=] __device__
#else
#define UPSIDE_DEVICE
#define DEVICE_LAMBDA [=]
#endif

namespace device {
    // implemented in device.cu
    void* allocate(size_t n_bytes);
    void  release(void* ptr);
    void  copy_to_device(void* dst, const void* src, size_t n_bytes);
    void  copy_to_host  (void* dst, const void* src, size_t n_bytes);
    void  copy_on_device(void* dst, const void* src, size_t n_bytes);
    void  zero(void* ptr, size_t n_bytes);
    void  synchronize();
    void  check(const char* where);   // throws the message of the last failed launch
    std::string describe();           // the backend, and for a GPU its name, capability, driver and runtime
}

// An array in device memory, owned and move-only
template <typename T>
struct DeviceArray {
    T*     ptr;
    size_t n;

    DeviceArray(): ptr(nullptr), n(0) {}
    explicit DeviceArray(size_t n_): ptr(n_ ? static_cast<T*>(device::allocate(n_*sizeof(T))) : nullptr), n(n_) {}
    explicit DeviceArray(const std::vector<T>& host): DeviceArray(host.size()) {upload(host);}
    DeviceArray(const DeviceArray&) = delete;
    DeviceArray& operator=(const DeviceArray&) = delete;
    DeviceArray(DeviceArray&& o): ptr(o.ptr), n(o.n) {o.ptr = nullptr; o.n = 0;}
    DeviceArray& operator=(DeviceArray&& o) {
        if(this != &o) {device::release(ptr); ptr = o.ptr; n = o.n; o.ptr = nullptr; o.n = 0;}
        return *this;
    }
    ~DeviceArray() {device::release(ptr);}

    void upload(const std::vector<T>& host) {
        if(host.size() != n) throw std::string("device upload of the wrong size");
        device::copy_to_device(ptr, host.data(), n*sizeof(T));
    }
    std::vector<T> download() const {
        std::vector<T> host(n);
        device::copy_to_host(host.data(), ptr, n*sizeof(T));
        return host;
    }
    void zero() {device::zero(ptr, n*sizeof(T));}
};


// Deterministic accumulation.  Every sum whose terms arrive from several threads is kept in 64-bit fixed point
// with 32 fraction bits, so its value does not depend on the order of the terms and a GPU run repeats bit for bit.
// A term that is not finite or whose magnitude reaches 2^26 cannot be represented; it marks its element, which
// then reads back as NaN, so a blow-up reaches the coordinates as it does on the CPU.
typedef unsigned long long accum_t;
constexpr float  accum_scale = 4294967296.f;   // 2^32
constexpr double accum_unit  = 1./4294967296.; // 2^-32
constexpr float  accum_limit = 67108864.f;     // 2^26

//! a term in fixed point; false if it cannot be represented
UPSIDE_HD inline bool to_accum(accum_t& q, float v) {
    if(!(fabsf(v) < accum_limit)) return false;
    q = accum_t(llrintf(v*accum_scale));
    return true;
}

// A thread's partial sum of terms bound for one element.  The sum is of integers, so adding it to the element at
// once gives the element the bits that adding each term would.
struct AccumLocal {
    accum_t sum;
    int     bad;

    UPSIDE_HD AccumLocal(): sum(0), bad(0) {}
    UPSIDE_HD void add(float v) {
        accum_t q;
        if(to_accum(q, v)) sum += q;
        else bad = 1;
    }
};

struct AccumView {
    accum_t* sum;
    int*     bad;

    UPSIDE_HD void add_fixed(int i, accum_t q) const {
#ifdef __CUDA_ARCH__
        atomicAdd(sum+i, q);
#else
        sum[i] += q;
#endif
    }
    UPSIDE_HD void add(int i, float v) const {
        accum_t q;
        if(to_accum(q, v)) add_fixed(i, q);
        else bad[i] = 1;
    }
    UPSIDE_HD void add(int i, const AccumLocal& part) const {
        if(part.bad) bad[i] = 1;
        if(part.sum) add_fixed(i, part.sum);
    }
    UPSIDE_HD float get(int i) const {
        return bad[i] ? NAN : float(double((long long)(sum[i]))*accum_unit);
    }
};

struct DeviceAccum {
    DeviceArray<accum_t> sum;
    DeviceArray<int>     bad;

    DeviceAccum() {}
    explicit DeviceAccum(size_t n): sum(n), bad(n) {zero();}
    void zero() {sum.zero(); bad.zero();}
    AccumView view() const {return AccumView{sum.ptr, bad.ptr};}
    size_t size() const {return sum.n;}
    std::vector<float> download() const {
        auto s = sum.download();
        auto b = bad.download();
        std::vector<float> r(s.size());
        for(size_t i=0; i<s.size(); ++i) r[i] = b[i] ? NAN : float(double((long long)(s[i]))*accum_unit);
        return r;
    }
};


// One thread per index
#ifdef __CUDACC__
template <typename F>
__global__ void parallel_for_kernel(int n, F f) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i<n) f(i);
}

template <typename F>
void parallel_for(int n, const F& f) {
    if(n<=0) return;
    parallel_for_kernel<<<(n+127)/128, 128>>>(n, f);
    device::check("parallel_for");
}
#else
template <typename F>
void parallel_for(int n, const F& f) {
    for(int i=0; i<n; ++i) f(i);
}
#endif


// The used slots of a per-system slot array: slots_for(n_system, cap, n_used, f) runs f(e) for e in
// [s*cap, s*cap + n_used[s]) of each system s, n_used in device memory.  Each system has a row of a few blocks that
// stride over its used slots, so the launch is sized by them rather than by the capacity.
#ifdef __CUDACC__
template <typename F>
__global__ void slots_for_kernel(int cap, const int* n_used, F f) {
    int ns = blockIdx.y, n = n_used[ns];
    for(int k=blockIdx.x*blockDim.x + threadIdx.x; k<n; k+=gridDim.x*blockDim.x) f(ns*cap + k);
}

template <typename F>
void slots_for(int n_system, int cap, const int* n_used, const F& f) {
    if(n_system<=0 || cap<=0) return;
    constexpr int block = 128, max_blocks_per_system = 64;
    int blocks = (cap+block-1)/block;
    slots_for_kernel<<<dim3(blocks < max_blocks_per_system ? blocks : max_blocks_per_system, n_system), block>>>(
            cap, n_used, f);
    device::check("slots_for");
}
#else
template <typename F>
void slots_for(int n_system, int cap, const int* n_used, const F& f) {
    for(int ns=0; ns<n_system; ++ns)
        for(int k=0; k<n_used[ns]; ++k) f(ns*cap + k);
}
#endif


// One warp per index.  lanes_for(n, f) runs f(i, lanes) for i in [0,n) on the 32 threads of a warp, which divide the
// index's items among them through lanes; f's branches outside the items are uniform across the warp.  Each of
// lanes' combinations gives the serial loop's value: places keep the items' order, float sums add in it, and the
// fixed-point sums are of integers.  On the host backend one thread (width 1) runs each index, so the same code runs
// serially; the host bodies are that backend's.
#ifdef __CUDACC__
constexpr unsigned full_warp = 0xffffffffu;
#endif

struct Lanes {
    int lane;
#ifdef __CUDACC__
    static constexpr int width = 32;
#else
    static constexpr int width = 1;
#endif

    //! g(j, in) on every lane for each chunk of width items from j0 to j1, in order; in is false past j1
    template <typename G>
    UPSIDE_HD void chunks(int j0, int j1, const G& g) const {
        for(int base=j0; base<j1; base+=width) g(base+lane, base+lane<j1);
    }
    //! the slot of a kept item in an ordered output holding count items; count advances past the chunk's kept items
    UPSIDE_HD int place(bool keep, int& count) const {
#ifdef __CUDA_ARCH__
        unsigned kept = __ballot_sync(full_warp, keep);
        int slot = count + __popc(kept & ((1u<<lane)-1u));
        count += __popc(kept);
        return slot;
#else
        int slot = count;
        count += keep;
        return slot;
#endif
    }
    //! s plus every lane's v, added in lane order
    UPSIDE_HD float add_in_order(float s, float v) const {
#ifdef __CUDA_ARCH__
        for(int k=0; k<width; ++k) s += __shfl_sync(full_warp, v, k);
        return s;
#else
        return s + v;
#endif
    }
    //! every lane's partial sum, added to element i
    UPSIDE_HD void flush(const AccumLocal& part, const AccumView& a, int i) const {
        AccumLocal total = part;
#ifdef __CUDA_ARCH__
        for(int off=width/2; off>0; off>>=1) {
            total.sum += __shfl_down_sync(full_warp, total.sum, off);
            total.bad |= __shfl_down_sync(full_warp, total.bad, off);
        }
#endif
        if(!lane) a.add(i, total);
    }
};

#ifdef __CUDACC__
template <typename F>
__global__ void lanes_for_kernel(int n, F f) {
    int i = blockIdx.x*(blockDim.x/Lanes::width) + threadIdx.x/Lanes::width;
    if(i<n) f(i, Lanes{int(threadIdx.x%Lanes::width)});
}

template <typename F>
void lanes_for(int n, const F& f) {
    if(n<=0) return;
    constexpr int warps_per_block = 4;
    lanes_for_kernel<<<(n+warps_per_block-1)/warps_per_block, warps_per_block*Lanes::width>>>(n, f);
    device::check("lanes_for");
}
#else
template <typename F>
void lanes_for(int n, const F& f) {
    for(int i=0; i<n; ++i) f(i, Lanes{0});
}
#endif


// One thread block per system.  per_system(n, f) runs f(ns) on every thread of block ns; f is written for the
// whole block, with every branch uniform across it except inside block_for and block_scan, which divide their
// items among the block's threads and end with every thread synchronised.  On the host one "thread" runs each
// block, so the same code runs serially.
constexpr int device_block_size = 256;

#ifdef __CUDACC__
template <typename F>
__global__ void per_system_kernel(F f) {f(int(blockIdx.x));}

template <typename F>
void per_system(int n_system, const F& f) {
    if(n_system<=0) return;
    per_system_kernel<<<n_system, device_block_size>>>(f);
    device::check("per_system");
}

//! f(i) for i in [0,n), the items divided among the block's threads
template <typename F>
__device__ void block_for(int n, const F& f) {
    for(int i=threadIdx.x; i<n; i+=blockDim.x) f(i);
    __syncthreads();
}

//! write(i, prefix, value(i)) for i in [0,n), where prefix is the sum of value(j) over j<i, and the total to every
//! thread; value is called once per i.  Each chunk of the block's size is scanned within its warps by shuffles, then
//! across the warps' totals.
template <typename V, typename W>
__device__ int block_scan(int n, const V& value, const W& write) {
    constexpr int n_warp = device_block_size/32;
    __shared__ int warp_total[n_warp];
    int lane = threadIdx.x%32, warp = threadIdx.x/32;
    int carry = 0;  // the sum over the chunks before, the same on every thread
    for(int base=0; base<n; base+=blockDim.x) {
        int i = base + threadIdx.x;
        int v = i<n ? value(i) : 0;
        int x = v;  // inclusive scan within the warp
        for(int off=1; off<32; off<<=1) {
            int y = __shfl_up_sync(full_warp, x, off);
            if(lane>=off) x += y;
        }
        if(lane==31) warp_total[warp] = x;
        __syncthreads();
        if(warp==0) {
            int t = lane<n_warp ? warp_total[lane] : 0;
            for(int off=1; off<n_warp; off<<=1) {
                int y = __shfl_up_sync(full_warp, t, off);
                if(lane>=off) t += y;
            }
            if(lane<n_warp) warp_total[lane] = t;  // now inclusive over the warps
        }
        __syncthreads();
        if(i<n) write(i, carry + (warp ? warp_total[warp-1] : 0) + x - v, v);
        carry += warp_total[n_warp-1];
        __syncthreads();
    }
    return carry;
}

//! the largest value(i) over i in [0,n) to every thread, NaN if any value is NaN (as max_keep_nan); -inf for n=0
template <typename V>
__device__ float block_max_keep_nan(int n, const V& value) {
    __shared__ float buf[device_block_size];
    float m = -INFINITY;
    for(int i=threadIdx.x; i<n; i+=blockDim.x) {float v = value(i); m = (v!=v || m!=m) ? NAN : fmaxf(m,v);}
    buf[threadIdx.x] = m;
    __syncthreads();
    for(int off=blockDim.x/2; off>0; off>>=1) {
        if(threadIdx.x<off) {
            float a = buf[threadIdx.x], b = buf[threadIdx.x+off];
            buf[threadIdx.x] = (a!=a || b!=b) ? NAN : fmaxf(a,b);
        }
        __syncthreads();
    }
    float r = buf[0];
    __syncthreads();
    return r;
}
#else
template <typename F>
void per_system(int n_system, const F& f) {
    for(int ns=0; ns<n_system; ++ns) f(ns);
}

template <typename F>
void block_for(int n, const F& f) {
    for(int i=0; i<n; ++i) f(i);
}

template <typename V, typename W>
int block_scan(int n, const V& value, const W& write) {
    int sum = 0;
    for(int i=0; i<n; ++i) {int v = value(i); write(i, sum, v); sum += v;}
    return sum;
}

template <typename V>
float block_max_keep_nan(int n, const V& value) {
    float m = -INFINITY;
    for(int i=0; i<n; ++i) {float v = value(i); m = (v!=v || m!=m) ? NAN : fmaxf(m,v);}
    return m;
}
#endif

#endif
