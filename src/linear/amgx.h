// Created by Petr Karnakov on 15.11.2020
// Copyright 2020 ETH Zurich

#pragma once

#define AMGXCALL(amgxcall_expr)                                               \
  do {                                                                        \
    AMGX_RC err = amgxcall_expr;                                              \
    if (err != AMGX_RC_OK) {                                                  \
      const int kMaxLength = 4096;                                            \
      char buf[kMaxLength];                                                   \
      int buf_len = kMaxLength;                                               \
      AMGX_get_error_string(err, buf, buf_len);                               \
      throw std::runtime_error(                                               \
          FILELINE + ": AMGX failed with error code " + std::to_string(err) + \
          "\nError string: " + buf + "\nKnown error codes: " +             \
          Amgx::GetErrorCodes());                                             \
    }                                                                         \
  } while (0)

#define CUDACALL(cudacall_expr)                                    \
  do {                                                             \
    cudaError_t err = cudacall_expr;                               \
    if (err != cudaSuccess) {                                      \
      throw std::runtime_error(                                    \
          FILELINE + ": CUDA failed: " + cudaGetErrorString(err)); \
    }                                                              \
  } while (0)

namespace Amgx {

std::string GetErrorCodes() {
  std::stringstream res;
  res << "AMGX_RC_OK=" << AMGX_RC_OK;
  res << ", AMGX_RC_BAD_PARAMETERS=" << AMGX_RC_BAD_PARAMETERS;
  res << ", AMGX_RC_UNKNOWN=" << AMGX_RC_UNKNOWN;
  res << ", AMGX_RC_BAD_MODE=" << AMGX_RC_BAD_MODE;
  return res.str();
}

class Library {
 public:
  Library(std::ostream* log = &std::cout) {
    log_ = log;
    AMGXCALL(AMGX_register_print_callback(print_callback));
    AMGXCALL(AMGX_initialize());
    AMGXCALL(AMGX_initialize_plugins());
    AMGXCALL(AMGX_install_signal_handler());
  }
  ~Library() noexcept(false) {
    AMGXCALL(AMGX_finalize_plugins());
    AMGXCALL(AMGX_finalize());
  }
  auto GetVersion() const {
    struct {
      int major;
      int minor;
      std::string version;
      std::string date;
      std::string time;
    } res;
    AMGXCALL(AMGX_get_api_version(&res.major, &res.minor));
    char* version;
    char* date;
    char* time;
    AMGXCALL(AMGX_get_build_info_strings(&version, &date, &time));
    res.version = version;
    res.date = date;
    res.time = time;
    return res;
  }

 private:
  static std::ostream* log_;
  static void print_callback(const char* msg, int length) {
    static int cnt = 0;
    if (log_) {
      (*log_) << util::Format("==== {:} ====\n{}\n", cnt++, msg);
    }
  }
};

std::ostream* Library::log_ = nullptr;

class Config {
 public:
  Config(std::string path, std::string extra) {
    AMGXCALL(AMGX_config_create_from_file_and_string(
        &handle_, path.c_str(), extra.c_str()));
  }
  ~Config() noexcept(false) {
    AMGXCALL(AMGX_config_destroy(handle_));
  }
  operator AMGX_config_handle() const {
    return handle_;
  }
  int GetNumRings() const {
    int res;
    AMGXCALL(AMGX_config_get_default_number_of_rings(handle_, &res));
    return res;
  }
  void AddParameters(std::string add) {
    AMGXCALL(AMGX_config_add_parameters(&handle_, add.c_str()));
  }

 private:
  AMGX_config_handle handle_;
};

class Mode {
 public:
  Mode(AMGX_Mode mode) : mode_(mode) {}
  Mode(std::string s) : mode_(FromString(s)) {}
  operator AMGX_Mode() const {
    return mode_;
  }
  static AMGX_Mode FromString(std::string s) {
    if (s == "hDDI") return AMGX_mode_hDDI;
    if (s == "hDFI") return AMGX_mode_hDFI;
    if (s == "hFFI") return AMGX_mode_hFFI;
    if (s == "dDDI") return AMGX_mode_dDDI;
    if (s == "dDFI") return AMGX_mode_dDFI;
    if (s == "dFFI") return AMGX_mode_dFFI;
    fassert(false, "Unknown mode='" + s + "'");
  };
  bool IsVectorDouble() const {
    switch (mode_) {
      case AMGX_mode_hDDI:
      case AMGX_mode_dDDI:
        return true;
      default:
        return false;
    }
  }
  bool IsMatrixDouble() const {
    switch (mode_) {
      case AMGX_mode_hDDI:
      case AMGX_mode_hDFI:
      case AMGX_mode_dDDI:
      case AMGX_mode_dDFI:
        return true;
      default:
        return false;
    }
  }

 private:
  AMGX_Mode mode_;
};

class Vector {
 public:
  Vector(AMGX_resources_handle resources, AMGX_Mode mode) : mode_(mode) {
    AMGXCALL(AMGX_vector_create(&handle_, resources, mode));
  }
  Vector(const Vector&) = delete;
  Vector(Vector&&) = default;
  ~Vector() noexcept(false) {
    AMGXCALL(AMGX_vector_destroy(handle_));
  }
  Vector& operator=(const Vector&) = delete;
  Vector& operator=(Vector&&) = default;
  operator AMGX_vector_handle() const {
    return handle_;
  }
  struct Size {
    int n, block;
  };
  Size GetSize() const {
    Size res;
    AMGXCALL(AMGX_vector_get_size(handle_, &res.n, &res.block));
    return res;
  }
  template <class Scal>
  void Download(Scal* buf) const {
    fassert_equal(mode_.IsVectorDouble(), sizeof(Scal) == 8);
    AMGXCALL(AMGX_vector_download(handle_, (void*)buf));
  }
  template <class Scal>
  std::vector<Scal> Download() const {
    const auto size = GetSize();
    std::vector<Scal> buf(size.n * size.block);
    Download(buf.data());
    return buf;
  }
  void Upload(const double* buf, Size size) {
    fassert(mode_.IsVectorDouble());
    AMGXCALL(AMGX_vector_upload(handle_, size.n, size.block, (const void*)buf));
  }
  void Upload(const std::vector<double>& buf, Size size) {
    fassert_equal(size.n * size.block, buf.size());
    Upload(buf.data(), size);
  }
  void Bind(AMGX_matrix_handle matrix) {
    AMGXCALL(AMGX_vector_bind(handle_, matrix));
  }
  void SetZero(Size size) const {
    AMGXCALL(AMGX_vector_set_zero(handle_, size.n, size.block));
  }
  void SetZero() const {
    SetZero(GetSize());
  }
  template <class Scal>
  void PrintGlobal(std::ostream& out, MPI_Comm comm) const {
    const auto send = Download<Scal>();
    const bool isroot = MpiWrapper::GetCommRank(comm) == 0;
    const auto type = sizeof(Scal) == 8 ? MPI_DOUBLE : MPI_FLOAT;
    if (isroot) {
      const int commsize = MpiWrapper::GetCommSize(comm);
      std::vector<Scal> recv(send.size() * commsize);
      MPI_Gather(
          send.data(), send.size(), type, //
          recv.data(), send.size(), type, 0, comm);
      for (auto v : recv) {
        out << v << ' ';
      }
    } else {
      MPI_Gather(
          send.data(), send.size(), type, //
          nullptr, 0, MPI_DATATYPE_NULL, 0, comm);
    }
  }

 private:
  AMGX_vector_handle handle_;
  Mode mode_;
};

class Matrix {
 public:
  Matrix(AMGX_matrix_handle h) : handle_(h) {}
  Matrix(AMGX_resources_handle resources, AMGX_Mode mode) {
    AMGXCALL(AMGX_matrix_create(&handle_, resources, mode));
  }
  Matrix(const Matrix&) = delete;
  Matrix(Matrix&&) = default;
  ~Matrix() noexcept(false) {
    AMGXCALL(AMGX_matrix_destroy(handle_));
  }
  Matrix& operator=(const Matrix&) = delete;
  Matrix& operator=(Matrix&&) = default;
  operator AMGX_matrix_handle() const {
    return handle_;
  }
  struct Size {
    int n, block_x, block_y;
  };
  Size GetSize() const {
    Size res;
    AMGXCALL(AMGX_matrix_get_size(handle_, &res.n, &res.block_x, &res.block_y));
    return res;
  }

 private:
  AMGX_matrix_handle handle_;
};

class Solver {
 public:
  Solver(
      AMGX_resources_handle resources, AMGX_Mode mode,
      AMGX_config_handle config) {
    AMGXCALL(AMGX_solver_create(&handle_, resources, mode, config));
  }
  ~Solver() noexcept(false) {
    AMGXCALL(AMGX_solver_destroy(handle_));
  }
  operator AMGX_solver_handle() const {
    return handle_;
  }
  AMGX_SOLVE_STATUS GetStatus() const {
    AMGX_SOLVE_STATUS status;
    AMGXCALL(AMGX_solver_get_status(handle_, &status));
    return status;
  }
  int GetNumIters() const {
    int n;
    AMGXCALL(AMGX_solver_get_iterations_number(handle_, &n));
    return n;
  }
  int GetResidual(int iter, int idx_in_block) const {
    double res;
    AMGXCALL(
        AMGX_solver_get_iteration_residual(handle_, iter, idx_in_block, &res));
    return res;
  }
  int GetResidual() const {
    return GetResidual(GetNumIters() - 1, 0);
  }

 private:
  AMGX_solver_handle handle_;
};

class Resources {
 public:
  Resources(AMGX_config_handle config) {
    AMGXCALL(AMGX_resources_create_simple(&handle_, config));
  }
  Resources(AMGX_config_handle config, MPI_Comm comm, int device)
      : comm_(comm) {
    AMGXCALL(AMGX_resources_create(&handle_, config, &comm_, 1, &device));
  }
  ~Resources() noexcept(false) {
    AMGXCALL(AMGX_resources_destroy(handle_));
  }
  operator AMGX_resources_handle() const {
    return handle_;
  }

 private:
  AMGX_resources_handle handle_;
  MPI_Comm comm_; // copy of communicator for AMGX_resources_create()
                  // that does not accept a temporary
};

} // namespace Amgx

