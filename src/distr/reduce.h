#include <vector>
#include <memory>
#include <utility>

// U: utility class

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
    virtual void Append(T& a) { Append(a, *v_); }
    // Returns neutral value a such that Append(a, v) would set a=v
    virtual T Neut() const = 0;
    virtual void Set(const T& v) { *v_ = v; }

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
    Scal Neut() const override { return 0.; }

   protected:
    void Append(Scal& a, const Scal& v) const override { a += v; }
  };
  class OpProd : public OpS {
   public:
    using OpS::OpS;
    Scal Neut() const override { return 1.; }
   
   protected:
    void Append(Scal& a, const Scal& v) const override { a *= v; }
  };
  class OpMax : public OpS {
   public:
    using OpS::OpS;
    Scal Neut() const override { return -std::numeric_limits<Scal>::max(); }

   protected:
    void Append(Scal& a, const Scal& v) const override { a = std::max(a, v); }
  };
  class OpMin : public OpS {
   public:
    using OpS::OpS;
    Scal Neut() const override { return std::numeric_limits<Scal>::max(); }

   protected:
    void Append(Scal& a, const Scal& v) const override { a = std::min(a, v); }
  };

  // Reduction on std::pair<Scal, int>
  using OpSI = OpT<std::pair<Scal, int>>;
  class OpMinloc : public OpSI {
   public:
    using T = std::pair<Scal, int>;
    using OpSI::OpSI;
    T Neut() const override { 
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
    T Neut() const override { 
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
    using V = std::vector<char>;

    using OpVT<char>::OpVT;
    using P::Append;
    using P::v_;
    using P::Set;
    V Neut() const override { return V(); }


   protected:
    void Append(V& a, const V& v) const override { 
      a.insert(a.end(), v.begin(), v.end());
    }
  };

  // Concatenation of std::vector<T>.
  // Result only on root block.
  template <class T>
  class OpCatT : public OpCat {
   public:
    using P = OpCat; // parent
    using V = std::vector<char>;
    using VT = std::vector<T>;
    // vt: buffer containing current value and used for result
    // OpCat::v_ initialized with nullptr // TODO: revise
    OpCatT(VT* vt) : OpCat(nullptr), vt_(vt) {}
    // Serializes vt_ and appends to a
    void Append(V& a) override { 
      V v = Ser(*vt_);
      a.insert(a.end(), v.begin(), v.end());
    }
    // Deserializes v to vt_
    void Set(const V& v) override { 
      *vt_ = Des(v); 
    }

   protected:
    using P::Append;
    VT* vt_;

   private:
    // Serialize vt 
    static V Ser(const VT& vt) {
      V r; // result
      for (const T& x : vt) {
        const char* c = (const char*)(const void*)&x;
        for (size_t i = 0; i < sizeof(T); ++i) {
          r.push_back(c[i]);
        }
      }
      return r;
    }
    // Deserialize v
    static VT Des(const V& v) {
      VT r; // result
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
    using V = std::vector<char>;
    using VT = std::vector<std::vector<T>>;
    // vt: buffer containing current value and used for result
    // OpCat::v_ initialized with nullptr // TODO: revise
    OpCatVT(VT* vt) : OpCat(nullptr), vt_(vt) {}
    // Serializes vt_ and appends to a
    void Append(V& a) override { 
      V v = Ser(*vt_);
      a.insert(a.end(), v.begin(), v.end());
    }
    // Deserializes v to vt_
    void Set(const V& v) override { 
      *vt_ = Des(v); 
    }

   protected:
    using P::Append;
    VT* vt_;

   private:
    // Serialize vt 
    static V Ser(const VT& vt) {
      V r; // result
      for (auto& xx : vt) {
        // write size
        {
          size_t s = xx.size();
          const char* c = (const char*)(const void*)&s;
          for (size_t i = 0; i < sizeof(size_t); ++i) {
            r.push_back(c[i]);
          }
        }

        // write elems
        for (auto& x : xx) {
          const char* c = (const char*)(const void*)&x;
          for (size_t i = 0; i < sizeof(T); ++i) {
            r.push_back(c[i]);
          }
        }
      }
      return r;
    }
    // Deserialize v
    static VT Des(const V& v) {
      VT r; // result
      //assert(v.size() % sizeof(T) == 0); // TODO: add assert
      size_t k = 0;
      while (k < v.size()) {
        r.emplace_back();
        size_t s;
        // read size
        {
          char* c = (char*)(void*)&s;
          for (size_t i = 0; i < sizeof(size_t); ++i) {
            c[i] = v[k++];
          }
        }
        // read elems
        for (size_t i = 0; i < s; ++i) {
          assert(k < v.size());
          T x;
          char* c = (char*)(void*)&x;
          for (size_t i = 0; i < sizeof(T); ++i) {
            c[i] = v[k++];
          }
          r.back().push_back(x);
        }
      }
      return r;
    }
  };


  void Add(const std::shared_ptr<Op>& o) {
    vrd_.push_back(o);
  }
  // u: src and dst buffer
  // o: operation 
  void Add(Scal* u, std::string o) {
    if (o == "sum") {
      vrd_.push_back(std::make_shared<OpSum>(u));
    } else if (o == "prod") {
      vrd_.push_back(std::make_shared<OpProd>(u));
    } else if (o == "max") {
      vrd_.push_back(std::make_shared<OpMax>(u));
    } else if (o == "min") {
      vrd_.push_back(std::make_shared<OpMin>(u));
    } else {
      throw std::runtime_error("Reduce: unknown operation: '" + o); 
    }
  }

  const std::vector<std::shared_ptr<Op>>& Get() const {
    return vrd_;
  }
  void Clear() {
    vrd_.clear();
  }

 private:
  std::vector<std::shared_ptr<Op>> vrd_; // list of requests
};

