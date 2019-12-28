#pragma once

#include "geom/field.h"
#include "parse/vars.h"

template <class M>
void InitVf(FieldCell<typename M::Scal>& fcu, const Vars& var, M& m);
