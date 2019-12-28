#pragma once

#include <functional>

#include "geom/field.h"
#include "parse/vars.h"

template <class M>
void InitVf(FieldCell<typename M::Scal>& fcu, const Vars& var, M& m);

template <class M>
std::function<void(
    FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&, const M&)>
CreateInitCl(const Vars& par, bool verb);

template <class M>
std::function<void(FieldCell<typename M::Scal>&, const M&)> CreateInitU(
    const Vars& par, bool verb);
