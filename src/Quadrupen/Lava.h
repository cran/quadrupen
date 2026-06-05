/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "SparseRegularizer.h"
#include "BaseLava.h"

using Rcpp::List;

template <typename matrix, SparseNorm norm>
class Lava : public BaseLava<matrix, SparseRegularizer<matrix,norm>> {

public:

  using Base = BaseLava<matrix, SparseRegularizer<matrix,norm>> ;

  Lava(RegressionData<matrix> orig_data, const List& regParam, const List& control)
    : Base(orig_data,
           Base::lava_preprocess(orig_data, Rcpp::as<double>(regParam["gamma"])),
           regParam, control) {}

} ;
