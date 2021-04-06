// Created by Petr Karnakov on 20.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "util/logger.h"

namespace ReductionType {
struct Sum {};
struct Prod {};
struct Max {};
struct Min {};
struct MaxLoc {};
struct MinLoc {};
struct Concat {};
} // namespace ReductionType

namespace Reduction {
constexpr ReductionType::Sum sum;
constexpr ReductionType::Prod prod;
constexpr ReductionType::Max max;
constexpr ReductionType::Min min;
constexpr ReductionType::MaxLoc maxloc;
constexpr ReductionType::MinLoc minloc;
constexpr ReductionType::Concat concat;
} // namespace Reduction

template <class Scal>
class UReduce {
 public:
  // Reduction operation.
  class Op {
   public:
    virtual ~Op() {}
  };
  // Reduction on type T
  template <class T>
  class OpT : public Op {
   public:
    // v: buffer containing current value and used for result
    OpT(T* v) : v_(v) {}
    // Appends internal value to a
    virtual void Append(T& a) const {
      Append(a, *v_);
    }
    // Returns neutral value a such that Append(a, v) would set a=v
    virtual T Neutral() const = 0;
    virtual void Set(const T& v) {
      *v_ = v;
    }

   protected:
    // Appends v to a
    virtual void Append(T& a, const T& v) const = 0;

    T* v_;
  };

  // Reduction on Scal
  using OpS = OpT<Scal>;
  class OpSum : public OpS {
   public:
    using OpS::OpS;
    Scal Neutral() const override {
      return 0.;
    }

   protected:
    void Append(Scal& a, const Scal& v) const override {
      a += v;
    }
  };
  class OpProd : public OpS {
   public:
    using OpS::OpS;
    Scal Neutral() const override {
      return 1;
    }

   protected:
    void Append(Scal& a, const Scal& v) const override {
      a *= v;
    }
  };
  class OpMax : public OpS {
   public:
    using OpS::OpS;
    Scal Neutral() const override {
      return -std::numeric_limits<Scal>::max();
    }

   protected:
    void Append(Scal& a, const Scal& v) const override {
      a = std::max(a, v);
    }
  };
  class OpMin : public OpS {
   public:
    using OpS::OpS;
    Scal Neutral() const override {
      return std::numeric_limits<Scal>::max();
    }

   protected:
    void Append(Scal& a, const Scal& v) const override {
      a = std::min(a, v);
    }
  };

  // Reduction on std::pair<Scal, int>
  using OpSI = OpT<std::pair<Scal, int>>;
  class OpMinloc : public OpSI {
   public:
    using T = std::pair<Scal, int>;
    using OpSI::OpSI;
    T Neutral() const override {
      return std::make_pair(std::numeric_limits<Scal>::max(), int());
    }

   protected:
    void Append(T& a, const T& v) const override {
      if (v.first < a.first) {
        a = v;
      }
    }
  };
  class OpMaxloc : public OpSI {
   public:
    using T = std::pair<Scal, int>;
    using OpSI::OpSI;
    T Neutral() const override {
      return std::make_pair(-std::numeric_limits<Scal>::max(), int());
    }

   protected:
    void Append(T& a, const T& v) const override {
      if (v.first > a.first) {
        a = v;
      }
    }
  };

  // Reduction on std::vector<T>
  template <class T>
  using OpVT = OpT<std::vector<T>>;

  // Concatenation of std::vector<char>
  // Result only on root block.
  class OpCat : public OpVT<char> {
   public:
    using P = OpVT<char>;
    using Bytes = std::vector<char>;

    using OpVT<char>::OpVT;
    using P::Append;
    using P::Set;
    using P::v_;
    Bytes Neutral() const override {
      return {};
    }

   protected:
    void Append(Bytes& a, const Bytes& v) const override {
      a.insert(a.end(), v.begin(), v.end());
    }
  };

  // Raw copy of object.
  template <class T>
  class OpCatRaw : public OpCat {
   public:
    using Base = OpCat;
    using Bytes = std::vector<char>;
    OpCatRaw(T* elem) : OpCat(nullptr), elem_(elem) {}
    // Appends `buf` with serialized `elem_`
    void Append(Bytes& buf) const override {
      const Bytes elembuf = Serialize(*elem_);
      buf.insert(buf.end(), elembuf.begin(), elembuf.end());
    }
    // Sets `elem_` to deserialized `buf`
    void Set(const Bytes& buf) override {
      *elem_ = Deserialize(buf);
    }

   protected:
    using Base::Append;
    T* elem_;

   private:
    static Bytes Serialize(const T& elem) {
      Bytes res;
      const char* chars = (const char*)(const void*)&elem;
      for (size_t i = 0; i < sizeof(T); ++i) {
        res.push_back(chars[i]);
      }
      return res;
    }
    static T Deserialize(const Bytes& buf) {
      T res;
      char* chars = (char*)(void*)&res;
      for (size_t i = 0; i < sizeof(T); ++i) {
        chars[i] = buf[i];
      }
      return res;
    }
  };

  // Concatenation of std::vector<T>.
  // Result only on root block.
  template <class T>
  class OpCatT : public OpCat {
   public:
    using P = OpCat; // parent
    using Bytes = std::vector<char>;
    // vt: buffer containing current value and used for result
    // OpCat::v_ initialized with nullptr // TODO: revise
    OpCatT(std::vector<T>* vt) : OpCat(nullptr), vt_(vt) {}
    // Serializes vt_ and appends to a
    void Append(Bytes& a) const override {
      Bytes v = Serialize(*vt_);
      a.insert(a.end(), v.begin(), v.end());
    }
    // Deserializes v to vt_
    void Set(const Bytes& v) override {
      *vt_ = Deserialize(v);
    }

   protected:
    using P::Append;
    std::vector<T>* vt_;

   private:
    static Bytes Serialize(const std::vector<T>& vt) {
      Bytes r; // result
      for (const T& x : vt) {
        const char* c = (const char*)(const void*)&x;
        for (size_t i = 0; i < sizeof(T); ++i) {
          r.push_back(c[i]);
        }
      }
      return r;
    }
    static std::vector<T> Deserialize(const Bytes& v) {
      std::vector<T> r; // result
      assert(v.size() % sizeof(T) == 0);
      size_t k = 0;
      while (k < v.size()) {
        T x;
        char* c = (char*)(void*)&x;
        for (size_t i = 0; i < sizeof(T); ++i) {
          c[i] = v[k++];
        }
        r.push_back(x);
      }
      return r;
    }
  };

  // Concatenation of std::vector<std::vector<T>>.
  // Result only on root block.
  template <class T>
  class OpCatVT : public OpCat {
   public:
    using P = OpCat; // parent
    using Bytes = std::vector<char>;
    // vt: buffer containing current value and used for result
    // OpCat::v_ initialized with nullptr // TODO: revise
    OpCatVT(std::vector<std::vector<T>>* vt) : OpCat(nullptr), vt_(vt) {}
    void Append(Bytes& a) const override {
      Bytes v = Serialize(*vt_);
      a.insert(a.end(), v.begin(), v.end());
    }
    void Set(const Bytes& v) override {
      *vt_ = Deserialize(v);
    }

   protected:
    using P::Append;
    std::vector<std::vector<T>>* vt_;

   private:
    static Bytes Serialize(const std::vector<std::vector<T>>& vt) {
      Bytes res;
      for (auto& xx : vt) {
        // write size
        {
          size_t s = xx.size();
          const char* c = (const char*)(const void*)&s;
          for (size_t i = 0; i < sizeof(size_t); ++i) {
            res.push_back(c[i]);
          }
        }

        // write elems
        for (auto& x : xx) {
          const char* c = (const char*)(const void*)&x;
          for (size_t i = 0; i < sizeof(T); ++i) {
            res.push_back(c[i]);
          }
        }
      }
      return res;
    }
    static std::vector<std::vector<T>> Deserialize(const Bytes& buf) {
      std::vector<std::vector<T>> res;
      size_t k = 0; // index in buf
      while (k < buf.size()) {
        res.emplace_back();
        size_t size = 0;
        // read size
        {
          auto* bytes = reinterpret_cast<unsigned char*>(&size);
          for (size_t i = 0; i < sizeof(size_t); ++i) {
            bytes[i] = buf[k++];
          }
        }
        // read elems
        for (size_t i = 0; i < size; ++i) {
          assert(k < buf.size());
          T x;
          auto* bytes = reinterpret_cast<unsigned char*>(&x);
          for (size_t j = 0; j < sizeof(T); ++j) {
            bytes[j] = buf[k++];
          }
          res.back().push_back(x);
        }
      }
      return res;
    }
  };
  void Add(std::unique_ptr<Op>&& o) {
    reqs_.emplace_back(std::move(o));
  }
  // u: src and dst buffer
  // o: operation
  void Add(Scal* u, std::string o) {
    if (o == "sum") {
      reqs_.push_back(std::make_unique<OpSum>(u));
    } else if (o == "prod") {
      reqs_.push_back(std::make_unique<OpProd>(u));
    } else if (o == "max") {
      reqs_.push_back(std::make_unique<OpMax>(u));
    } else if (o == "min") {
      reqs_.push_back(std::make_unique<OpMin>(u));
    } else {
      fassert(false, "Reduce: unknown operation: '" + o);
    }
  }

  const std::vector<std::unique_ptr<Op>>& Get() const {
    return reqs_;
  }
  void Clear() {
    reqs_.clear();
  }

  template <class T>
  static std::unique_ptr<OpCatT<T>> Make(
      std::vector<T>* buf, ReductionType::Concat) {
    return std::make_unique<OpCatT<T>>(buf);
  }

 private:
  std::vector<std::unique_ptr<Op>> reqs_; // list of requests
};
