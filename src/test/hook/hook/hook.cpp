// Created by Petr Karnakov on 03.05.2020
// Copyright 2020 ETH Zurich

#include <iostream>

#include <parse/config.h>
#include <util/posthook.h>

template <class M>
void InitVelHook(FieldCell<typename M::Vect>&, const Vars& var, const M&) {
  struct : public ConfigBase {
    VAR_INT(var0);
    VAR_DOUBLE(var1);
  } conf;
  conf.Read(var);

  std::cout << "[hook]";
  std::cout << " var0=" << conf.var0;
  std::cout << " var1=" << conf.var1;
  std::cout << std::endl;
}

using M = MeshCartesian<double, 3>;
template void InitVelHook(FieldCell<typename M::Vect>&, const Vars&, const M&);
