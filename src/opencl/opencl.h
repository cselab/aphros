// Created by Petr Karnakov on 14.04.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#define CL_TARGET_OPENCL_VERSION 200
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>
#include <CL/cl_ext.h>

#include "geom/mesh.h"
#include "parse/vars.h"
#include "util/logger.h"

#define CLCALL(x)                                                     \
  do {                                                                \
    cl_int CLCALL_error;                                              \
    CLCALL_error = x;                                                 \
    if (CLCALL_error != CL_SUCCESS) {                                 \
      throw std::runtime_error(                                       \
          std::string() + __FILE__ + ":" + std::to_string(__LINE__) + \
          ": CL failed: " + std::to_string(CLCALL_error));            \
    }                                                                 \
  } while (0)

template <class M, class T>
void SharedToLocal(
    const FieldCell<T>& fc_shared, FieldCell<T>& fc_local, const M& m) {
  using MIdx = typename M::MIdx;
  auto& ms = m.GetShared();
  fc_local.Reinit(m);
  auto& ics = ms.GetIndexCells();
  for (auto c : m.CellsM()) {
    const auto cs = ics.GetIdx(MIdx(c));
    fc_local[c] = fc_shared[cs];
  }
}

template <class M, class T>
void LocalToShared(
    const FieldCell<T>& fc_local, FieldCell<T>& fc_shared, const M& m) {
  using MIdx = typename M::MIdx;
  auto& ms = m.GetShared();
  fc_shared.Reinit(ms);
  auto& ics = ms.GetIndexCells();
  for (auto c : m.CellsM()) {
    const auto cs = ics.GetIdx(MIdx(c));
    fc_shared[cs] = fc_local[c];
  }
}

template <class M>
struct OpenCL {
  static std::string GetErrorMessage(cl_int error);

  using Scal = typename M::Scal;
  using MIdx = typename M::MIdx;
  using MSize = generic::Vect<size_t, M::dim>;
  struct Device {
    Device() = default;
    Device(const Device&) = delete;
    Device(Device&&) = delete;
    Device& operator=(Device&) = delete;
    Device& operator=(Device&&) = delete;
    struct PlatformInfo {
      cl_platform_id id;
      std::string name;
      std::string vendor;
    };
    static std::vector<PlatformInfo> GetPlatformInfos();
    struct DeviceInfo {
      cl_device_id id;
      std::string name;
      std::string extensions;
    };
    static cl_device_id GetDevice(cl_platform_id platform);
    static DeviceInfo GetDeviceInfo(cl_platform_id platform);
    void Create(size_t pid);
    ~Device();
    operator cl_device_id() const {
      return handle;
    }
    cl_device_id handle = NULL;
    cl_platform_id platform = NULL;
  };

  struct Context {
    Context() = default;
    Context(const Context&) = delete;
    Context(Context&&) = delete;
    Context& operator=(Context&) = delete;
    Context& operator=(Context&&) = delete;
    void Create(const cl_device_id* device);
    ~Context();
    operator cl_context() const {
      return handle;
    }
    cl_context handle = NULL;
  };

  struct Queue {
    Queue() = default;
    Queue(const Queue&) = delete;
    Queue(Queue&&) = delete;
    Queue& operator=(Queue&) = delete;
    Queue& operator=(Queue&&) = delete;
    void Create(cl_context context, cl_device_id device);
    void Finish();
    ~Queue();
    operator cl_command_queue() const {
      return handle;
    }
    cl_command_queue handle = NULL;
  };

  struct Program {
    Program() = default;
    Program(const Program&) = delete;
    Program(Program&&) = delete;
    Program& operator=(Program&) = delete;
    Program& operator=(Program&&) = delete;
    void CreateFromString(
        std::string source, cl_context context, cl_device_id device);
    void CreateFromStream(
        std::istream& in, cl_context context, cl_device_id device);
    void CreateFromFile(
        std::string source_path, cl_context context, cl_device_id device);
    operator cl_program() const {
      return handle;
    }
    ~Program();
    cl_program handle = NULL;
  };
  template <class T>
  struct Buffer {
    Buffer() = default;
    Buffer(const Buffer&) = delete;
    Buffer(Buffer&&) = delete;
    Buffer& operator=(Buffer&) = delete;
    Buffer& operator=(Buffer&&) = delete;
    void Create(
        cl_context context, size_t size_,
        cl_mem_flags flags = CL_MEM_READ_WRITE) {
      size = size_;
      cl_int error;
      handle = clCreateBuffer(context, flags, sizeof(T) * size, NULL, &error);
      CLCALL(error);
    }
    ~Buffer() {
      if (handle) {
        clReleaseMemObject(handle);
      }
    }
    void EnqueueRead(cl_command_queue queue, T* buf) const {
      CLCALL(clEnqueueReadBuffer(
          queue, handle, CL_TRUE, 0, sizeof(T) * size, buf, 0, NULL, NULL));
    }
    void EnqueueRead(cl_command_queue queue, std::vector<T>& buf) const {
      EnqueueRead(queue, buf.data());
    }
    void EnqueueWrite(cl_command_queue queue, const T* buf) {
      CLCALL(clEnqueueWriteBuffer(
          queue, handle, CL_TRUE, 0, sizeof(T) * size, buf, 0, NULL, NULL));
    }
    void EnqueueWrite(cl_command_queue queue, const std::vector<T>& buf) {
      EnqueueWrite(queue, buf.data());
    }
    operator cl_mem() const {
      return handle;
    }
    void swap(Buffer& other) {
      std::swap(handle, other.handle);
    }
    cl_mem handle = NULL;
    size_t size;
  };

  template <class T>
  struct MirroredBuffer : public Buffer<T> {
    using Base = Buffer<T>;
    using Base::handle;
    using Base::size;
    void Create(
        cl_context context, size_t size_,
        cl_mem_flags flags = CL_MEM_READ_WRITE) {
      Base::Create(context, size_, flags);
      buf.resize(size);
    }
    void EnqueueRead(cl_command_queue queue) {
      Base::EnqueueRead(queue, buf);
    }
    void EnqueueWrite(cl_command_queue queue) {
      Base::EnqueueWrite(queue, buf);
    }
    const T& operator[](size_t i) const {
      return buf[i];
    }
    T& operator[](size_t i) {
      return buf[i];
    }

    std::vector<T> buf;
  };

  struct Kernel {
    Kernel() = default;
    Kernel(const Kernel&) = delete;
    Kernel(Kernel&&) = delete;
    Kernel& operator=(Kernel&) = delete;
    Kernel& operator=(Kernel&&) = delete;
    void Create(cl_program program, std::string name_);
    ~Kernel();
    template <class T>
    void SetArg(int pos, const T& value) {
      cl_int error = clSetKernelArg(handle, pos, sizeof(T), &value);
      fassert_equal(
          error, CL_SUCCESS,
          ". SetArg failed for kernel '" + name + "' at position " +
              std::to_string(pos));
    }
    template <class T>
    void SetArg(int pos, const Buffer<T>& value) {
      SetArg(pos, value.handle);
    }
    template <class T>
    void SetArg(int pos, const MirroredBuffer<T>& value) {
      SetArg(pos, value.handle);
    }
    void Enqueue(cl_command_queue queue, MSize global, MSize local) {
      CLCALL(clEnqueueNDRangeKernel(
          queue, handle, M::dim, NULL, global.data(), local.data(), 0, NULL,
          NULL));
    }
    void EnqueueWithArgs(
        int, cl_command_queue queue, MSize global, MSize local) {
      Enqueue(queue, global, local);
    }
    template <class T, class... Args>
    void EnqueueWithArgs(
        int pos, cl_command_queue queue, MSize global, MSize local,
        const T& value, const Args&... args) {
      SetArg(pos, value);
      EnqueueWithArgs(pos + 1, queue, global, local, args...);
    }
    template <class... Args>
    void EnqueueWithArgs(
        cl_command_queue queue, MSize global, MSize local,
        const Args&... args) {
      EnqueueWithArgs(0, queue, global, local, args...);
    }

    operator cl_kernel() const {
      return handle;
    }
    cl_kernel handle = NULL;
    std::string name;
  };

  struct HaloComm {
    HaloComm() = default;
    HaloComm(const HaloComm&) = delete;
    HaloComm(HaloComm&&) = delete;
    HaloComm& operator=(HaloComm&) = delete;
    HaloComm& operator=(HaloComm&&) = delete;
    void Create(cl_context context, const M& ms, const OpenCL& cl);
    void Comm(
        M& m, typename M::Sem& sem, Buffer<Scal>& d_field, Queue& queue,
        OpenCL& cl);

    FieldCell<Scal> fc_buf;
    MirroredBuffer<Scal> d_buf;
    std::vector<IdxCell> cells_inner;
    std::vector<IdxCell> cells_halo;
  };

  OpenCL(const M& ms, const Vars& var);
  Scal Sum(cl_mem d_u);
  Scal Dot(cl_mem d_u, cl_mem d_v);
  void Comm(M& m, typename M::Sem& sem, Buffer<Scal>& buf) {
    halocomm.Comm(m, sem, buf, queue, *this);
  }

  Context context;
  Device device;
  Program program;
  Queue queue;
  Kernel kernel_inner_to_buf;
  Kernel kernel_buf_to_halo;
  Kernel kernel_dot;
  Kernel kernel_sum;
  MirroredBuffer<Scal> d_buf_reduce;
  HaloComm halocomm;

  size_t size;
  MSize global_size;
  MSize local_size;
  size_t ngroups;
  int start; // offset of first element of inner cells
  int lead_x; // leading dimension in x
  MIdx msize; // number of inner cells
};
