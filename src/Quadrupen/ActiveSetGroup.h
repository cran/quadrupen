/*
 * Author: Julien CHIQUET
 */
#pragma once

using arma::vec;
using arma::mat;
using arma::uvec;
using arma::uword;
using arma::regspace;
using std::vector;

#include "RegressionData.h"
#include "ActiveSet.h"

template <typename matrix>
class ActiveSetGroup: public ActiveSet<matrix> {

public:
  using ActiveSet<matrix>::A_;
  using ActiveSet<matrix>::XATXA_;

  // VARIABLES FOR HANDLING THE ACTIVE SET
  uvec G_             ; // set of currently activated groups
  uvec is_grp_in_     ; // indicator of active groups (0/1)
  uvec grp_sizes_     ; // vector of group sizes
  vector<uvec> group_ ; // vector of parameters for each group
  vector<mat> V_      ; // eigen vectors of each active group
  vector<vec> D_      ; // eigen values of each active group
  bool use_evd_       ; // Whether to maintain the EVD cache

public:
  ActiveSetGroup() : ActiveSet<matrix>(), use_evd_(false) {};
  ActiveSetGroup(const RegressionData<matrix> &data, const uvec& group_ind, const bool use_evd);

  // ── Active set handling ──────────────────────────────────────────────────────────
  void add_group(uword, const RegressionData<matrix> &);
  void del_group(uword, vec&);
  void del_groups(uvec igrps_out, vec& beta);

  const uword size_grp() const { return G_.n_elem; }

  // ── EVD cache handling ────────────────────────────────────────────────────────────
  void update_EVD(uword grp_id);
  void downdate_EVD(uword igrp_out);
};

// ── Constructor ─────────────────────────────────────────────────────────────────────
template <typename matrix>
ActiveSetGroup<matrix>::ActiveSetGroup(const RegressionData<matrix>& data, const uvec& group_ind, const bool use_evd) :
  ActiveSet<matrix>(data, false), use_evd_(use_evd) {

    uvec grp = unique(group_ind);
    uword nb_grp = grp.n_elem;
    grp_sizes_.zeros(nb_grp);

    for (uword i = 0; i < group_ind.n_elem; ++i) {
      if (group_ind[i] == 0)
        Rcpp::stop("group indices must be 1-based; received 0 at position %d", (int)i);
      grp_sizes_[group_ind[i] - 1]++;
    }

    uword current_elt = 0;
    for (uword i = 0; i < grp_sizes_.n_elem; ++i) {
      group_.push_back(regspace<uvec>(current_elt, current_elt + grp_sizes_[i] - 1));
      current_elt += grp_sizes_[i];
    }

    is_grp_in_.zeros(grp_sizes_.n_elem);
  }

// ── Update EVD ────────────────────────────────────────
template <typename matrix>
void ActiveSetGroup<matrix>::update_EVD(uword grp_id) {
  uword sz = grp_sizes_(grp_id);
  // Le bloc se trouve toujours à la fin de XATXA_ suite à add_vars
  uword start = this->size() - sz;
  mat Hg = this->XATXA_.submat(start, start, this->size() - 1, this->size() - 1);

  vec d; mat v;
  eig_sym(d, v, Hg);

  D_.push_back(d);
  V_.push_back(v);
}

template <typename matrix>
void ActiveSetGroup<matrix>::downdate_EVD(uword igrp_out) {
  V_.erase(V_.begin() + igrp_out);
  D_.erase(D_.begin() + igrp_out);
}

// ──── ADD GROUP ────
template <typename matrix>
void ActiveSetGroup<matrix>::add_group(uword grp_in, const RegressionData<matrix>& data) {
  G_.resize(size_grp() + 1);
  G_.tail(1) = grp_in;
  is_grp_in_[grp_in] = 1;

  ActiveSet<matrix>::add_vars(group_[grp_in], data);

  if (use_evd_) update_EVD(grp_in);
}

// ──── DEL GROUP ────
template <typename matrix>
void ActiveSetGroup<matrix>::del_group(uword igrp_out, vec& beta) {
  uword ifirst_var_in = 0;
  for (uword i = 0; i < igrp_out; i++) ifirst_var_in += grp_sizes_(G_[i]);

  uword n_vars_to_del = grp_sizes_(G_[igrp_out]);
  // Ordre décroissant pour shedding dans del_vars
  uvec ivars_out = regspace<uvec>(ifirst_var_in + n_vars_to_del - 1, ifirst_var_in);

  this->del_vars(ivars_out, beta);

  is_grp_in_[G_[igrp_out]] = 0;
  if (use_evd_) downdate_EVD(igrp_out);

  G_.shed_row(igrp_out);
}

// ──── DEL GROUPS ────
template <typename matrix>
void ActiveSetGroup<matrix>::del_groups(uvec igrps_out, vec& beta) {
  if (igrps_out.is_empty()) return;
  uvec sorted_indices = sort(igrps_out, "descend");
  for (uword i = 0; i < sorted_indices.n_elem; ++i) {
    this->del_group(sorted_indices(i), beta);
  }
}
