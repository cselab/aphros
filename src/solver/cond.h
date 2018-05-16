namespace solver {

class CondFace {
 public:
  CondFace(size_t nci) : nci_(nci) {}
  virtual ~CondFace() {}
  // neighbour cell id
  virtual size_t GetNci() const {
    return nci_;
  }

 private:
  size_t nci_;
};

class CondFaceExtrap : public CondFace {
 public:
  CondFaceExtrap(size_t nci) : CondFace(nci) {}
};

template <class V>
class CondFaceVal : public CondFace {
 public:
  CondFaceVal(size_t nci) : CondFace(nci) {}
  virtual V GetValue() const = 0;
};

template <class Vect>
class CondFaceValComp :
    public CondFaceVal<typename Vect::value_type> {
 public:
  using Scal = typename Vect::value_type;
  using P = CondFaceVal<Scal>; // parent
  CondFaceValComp(CondFaceVal<Vect>* o, size_t d)
      : P(o->GetNci()), o_(o), d_(d) {}
  Scal GetValue() const override { return o_->GetValue()[d_]; }

 private:
  CondFaceVal<Vect>* o_;
  size_t d_;
};

template <class V>
class CondFaceValFixed : public CondFaceVal<V> {
 public:
  CondFaceValFixed(const V& v, size_t nci) : CondFaceVal<V>(nci), v_(v) {}
  V GetValue() const override { return v_; }
  void Set(const V& v) { v_ = v; }

 private:
  V v_;
};

template <class V>
class CondFaceGrad : public CondFace {
 public:
  CondFaceGrad(size_t nci) : CondFace(nci) {}
  virtual V GetGrad() const = 0;
};

template <class V>
class CondFaceGradFixed : public CondFaceGrad<V> {
 public:
  explicit CondFaceGradFixed(const V& v, size_t nci)
      : CondFaceGrad<V>(nci), v_(v) {}
  virtual V GetGrad() const override { return v_; }
  void Set(const V& v) { v = v; }

 private:
  V v_;
};

class CondCell {
 public:
  virtual ~CondCell() {}
};

template <class V>
class CondCellVal : public CondCell {
 public:
  virtual V GetValue() const = 0;
};

template <class V>
class CondCellValFixed : public CondCellVal<V> {
 public:
  explicit CondCellValFixed(const V& v) : v_(v) {}
  V GetValue() const override { return v_; }
  void Set(const V& v) { v_ = v; }

 private:
  V v_;
};

} // namespace solver
