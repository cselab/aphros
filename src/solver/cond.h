namespace solver {

class ConditionFace {
 public:
  ConditionFace(size_t nci) : nci_(nci) {}
  virtual ~ConditionFace() {}
  // neighbour cell id
  virtual size_t GetNci() const {
    return nci_;
  }
 private:
  size_t nci_;
};

class ConditionFaceExtrapolation : public ConditionFace {
 public:
  ConditionFaceExtrapolation(size_t nci) : ConditionFace(nci) {}
};

template <class V>
class ConditionFaceValue : public ConditionFace {
 public:
  ConditionFaceValue(size_t nci) : ConditionFace(nci) {}
  virtual V GetValue() const = 0;
};

template <class Vect>
class ConditionFaceValueExtractComponent :
    public ConditionFaceValue<typename Vect::value_type> {
 public:
  using Scal = typename Vect::value_type;
  using P = ConditionFaceValue<Scal>; // parent
  ConditionFaceValueExtractComponent(ConditionFaceValue<Vect>* o, size_t d)
      : P(o->GetNci()), o_(o), d_(d) {}
  Scal GetValue() const override { return o_->GetValue()[d_]; }

 private:
  ConditionFaceValue<Vect>* o_;
  size_t d_;
};

template <class V>
class ConditionFaceValueFixed : public ConditionFaceValue<V> {
 public:
  ConditionFaceValueFixed(const V& v, size_t nci)
      : ConditionFaceValue<V>(nci), v_(v) {}
  V GetValue() const override { return v_; }
  void Set(const V& v) { v_ = v; }

 private:
  V v_;
};

template <class V>
class ConditionFaceDerivative : public ConditionFace {
 public:
  ConditionFaceDerivative(size_t nci) : ConditionFace(nci) {}
  virtual V GetDerivative() const = 0;
};

template <class V>
class ConditionFaceDerivativeFixed : public ConditionFaceDerivative<V> {
 public:
  explicit ConditionFaceDerivativeFixed(const V& v, size_t nci)
      : ConditionFaceDerivative<V>(nci), v_(v) {}
  virtual V GetDerivative() const override { return v_; }
  void Set(const V& v) { v = v; }

 private:
  V v_;
};

class ConditionCell {
 public:
  virtual ~ConditionCell() {}
};

template <class V>
class ConditionCellValue : public ConditionCell {
 public:
  virtual V GetValue() const = 0;
};

template <class V>
class ConditionCellValueFixed : public ConditionCellValue<V> {
 public:
  explicit ConditionCellValueFixed(const V& v) : v_(v) {}
  V GetValue() const override { return v_; }
  void Set(const V& v) { v_ = v; }

 private:
  V v_;
};

} // namespace solver
