/*
 *  VectorOperator.cpp
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 03/18/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#include "VectorOperator.h"

const std::string Operator_Vort_4th::NAME = "Vort_4th";
const std::string Operator_IgradA2I_4th::NAME = "IgradA2I_4th";
const std::string Operator_divU_4th::NAME = "divU_4th";
const std::string Operator_SSDeviatoric_4th::NAME = "SSDeviatoric_4th";
const std::string Operator_Qcrit_4th::NAME = "Qcrit_4th";
const std::string Operator_gradU_4th::NAME = "gradU_4th";
template <> const std::string Operator_gradUij_4th<0>::NAME = "gradU11_4th";
template <> const std::string Operator_gradUij_4th<1>::NAME = "gradU12_4th";
template <> const std::string Operator_gradUij_4th<2>::NAME = "gradU13_4th";
template <> const std::string Operator_gradUij_4th<3>::NAME = "gradU21_4th";
template <> const std::string Operator_gradUij_4th<4>::NAME = "gradU22_4th";
template <> const std::string Operator_gradUij_4th<5>::NAME = "gradU23_4th";
template <> const std::string Operator_gradUij_4th<6>::NAME = "gradU31_4th";
template <> const std::string Operator_gradUij_4th<7>::NAME = "gradU32_4th";
template <> const std::string Operator_gradUij_4th<8>::NAME = "gradU33_4th";
const std::string Operator_Ucontraction_4th::NAME = "Ucontraction_4th";
const std::string Operator_PIC_4th::NAME = "PIC_4th";
