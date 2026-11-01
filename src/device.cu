#include "device.h"
#include <cstring>
#include <cstdlib>

#ifdef __CUDACC__
#include <cuda_runtime.h>

namespace device {
    static void throw_on(cudaError_t err, const char* what) {
        if(err != cudaSuccess) throw std::string("CUDA ") + what + ": " + cudaGetErrorString(err);
    }

    void* allocate(size_t n_bytes) {
        void* ptr = nullptr;
        throw_on(cudaMalloc(&ptr, n_bytes), "allocation");
        return ptr;
    }
    void release(void* ptr) {if(ptr) cudaFree(ptr);}
    void copy_to_device(void* dst, const void* src, size_t n_bytes) {
        if(n_bytes) throw_on(cudaMemcpy(dst, src, n_bytes, cudaMemcpyHostToDevice), "copy to device");
    }
    void copy_to_host(void* dst, const void* src, size_t n_bytes) {
        if(n_bytes) throw_on(cudaMemcpy(dst, src, n_bytes, cudaMemcpyDeviceToHost), "copy to host");
    }
    void copy_on_device(void* dst, const void* src, size_t n_bytes) {
        if(n_bytes) throw_on(cudaMemcpy(dst, src, n_bytes, cudaMemcpyDeviceToDevice), "copy on device");
    }
    void zero(void* ptr, size_t n_bytes) {if(n_bytes) throw_on(cudaMemset(ptr, 0, n_bytes), "memset");}
    void synchronize() {throw_on(cudaDeviceSynchronize(), "synchronize");}
    void check(const char* where) {throw_on(cudaGetLastError(), where);}

    std::string describe() {
        int dev = 0, driver = 0, runtime = 0;
        throw_on(cudaGetDevice(&dev), "device query");
        cudaDeviceProp prop;
        throw_on(cudaGetDeviceProperties(&prop, dev), "device properties");
        cudaDriverGetVersion(&driver);
        cudaRuntimeGetVersion(&runtime);
        return std::string("CUDA device ") + std::to_string(dev) + " " + prop.name +
            ", compute capability " + std::to_string(prop.major) + "." + std::to_string(prop.minor) +
            ", " + std::to_string(prop.multiProcessorCount) + " multiprocessors, driver " + std::to_string(driver) +
            ", runtime " + std::to_string(runtime);
    }
}

#else

namespace device {
    // Memory starts as NaN, so that a read of memory no kernel wrote shows, as it would on a GPU, whose memory is
    // not cleared
    void* allocate(size_t n_bytes) {
        void* ptr = malloc(n_bytes);
        if(!ptr) throw std::string("host backend allocation failed");
        memset(ptr, 0xff, n_bytes);
        return ptr;
    }
    void release(void* ptr) {free(ptr);}
    void copy_to_device(void* dst, const void* src, size_t n_bytes) {if(n_bytes) memcpy(dst, src, n_bytes);}
    void copy_to_host  (void* dst, const void* src, size_t n_bytes) {if(n_bytes) memcpy(dst, src, n_bytes);}
    void copy_on_device(void* dst, const void* src, size_t n_bytes) {if(n_bytes) memcpy(dst, src, n_bytes);}
    void zero(void* ptr, size_t n_bytes) {if(n_bytes) memset(ptr, 0, n_bytes);}
    void synchronize() {}
    void check(const char* where) {}
    std::string describe() {return "host backend (the device code run serially on the CPU)";}
}

#endif
