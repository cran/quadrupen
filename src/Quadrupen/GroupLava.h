/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "GroupSparseRegularizer.h"
#include "BaseLava.h"

using arma::uvec;
using Rcpp::List;

template <typename matrix, GroupSparseNorm norm>
class GroupLava : public BaseLava<matrix, GroupSparseRegularizer<matrix,norm>> {

public:

  using Base = BaseLava<matrix, GroupSparseRegularizer<matrix,norm>> ;

  GroupLava(RegressionData<matrix> orig_data, const uvec& group_ind,
            const List& regParam, const List& control)
    : Base(orig_data,
           Base::lava_preprocess(orig_data, Rcpp::as<double>(regParam["gamma"])),
           group_ind, regParam, control) {}

} ;
