/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "BoundedRegression.h"

using namespace Rcpp;
using namespace arma;

BoundedRegression::BoundedRegression(
  RegressionData<mat> data, const List& regParam, const List& control) :
  Regularizer<mat>::Regularizer(std::move(data), regParam) {
    
    // Set the penalty to L-infinity
    penalty_ = DensePenalty<DenseNorm::LINF>() ;
    get_lambda_seq(get_lambda_max(), regParam) ;

    // Set up the optimizer 
    solver_ = OptimizerLINF(penalty_, control) ;

    // Scale the structuring matrix according to the amount of l2 penalty 
    data_.scale_struct(gamma_) ;
    
    // Compute the Gram matrix (+ gamma * S)
    data_.precompute_XTX() ;

    // Initialize beta_ and the gradient
    beta_ = zeros<vec>(data_.p_) ; // vector of current parameters
    grad_ = -data_.XTy_          ; // vector of current gradient (smooth part)
  }

List BoundedRegression::to_list(const List& monitoring) {
  return List::create(
    Named("tuning_param") = List::create(
      Named("linf") = lambdas_,
      Named("l2")   = gamma_
    ),
    Named("coef")       = coef_,
    Named("active")     = unbounded_var(),
    Named("intercept")  = intercept_,
    Named("normx")      = data_.norm_X_,
    Named("df")         = df_,
    Named("monitoring") = monitoring
  ) ;
}

double BoundedRegression::get_df() {

  double df = data_.centered_ + unbounded_.size();

  if (gamma_ > 0 && !unbounded_.is_empty()) {
    mat C = inv_sympd(data_.XTX_(unbounded_,unbounded_));
    uword ku = unbounded_.size();
    // position of each unbounded variable, then a single pass over the non-zeros of S
    std::vector<long> pos(data_.p_, -1);
    for (uword i = 0; i < ku; i++) pos[unbounded_(i)] = i;
    mat SUU(ku, ku, fill::zeros);
    for (auto it = data_.S_.begin(); it != data_.S_.end(); ++it) {
      if (pos[it.row()] >= 0 && pos[it.col()] >= 0) SUU(pos[it.row()], pos[it.col()]) = *it;
    }
    df -= accu(SUU % C); // trace(SUU * C) for symmetric matrices
  }

  return(df);
}

double BoundedRegression::optimality_gap(const vec& grad, const double lambda) const {
  // KKT conditions of min 1/2 b'Hb - c'b + lambda max_i w_i |b_i|, with g = Hb - c:
  //  - b = 0 : sum_i |g_i| / w_i <= lambda
  //  - b != 0, B = {i : w_i |b_i| = max} : g_i = 0 outside B, sign(g_i) = -sign(b_i) on B
  //    and sum_{i in B} |g_i| / w_i = lambda
  const vec& w = lambda_factor_ ;
  const vec wb = abs(beta_) % w ;
  const double bound = wb.max() ;
  const vec gw = abs(grad) / w ;
  if (bound == 0.0) return std::max(0.0, accu(gw) - lambda) ;

  double gap = 0.0, dual_B = 0.0 ;
  for (uword i = 0; i < beta_.n_elem; ++i) {
    if (wb(i) >= bound * (1.0 - 1e-8)) {
      dual_B += gw(i) ;
      if (grad(i) * beta_(i) > 0.0) gap = std::max(gap, gw(i)) ;
    } else {
      gap = std::max(gap, gw(i)) ;
    }
  }
  return std::max(gap, std::abs(dual_B - lambda)) ;
}

List BoundedRegression::solution_path(const List& control) {

  // Parameters controlling the optimization
  const bool verbose(control["verbose"])      ; // verbosity level
  const double accuracy(control["threshold"]) ; // precision required
  const uword maxiter(control["maxiter"])     ; // max # of passes in the active set
  const uword maxfeat(control["maxfeat"])     ; // max # of variables activated

  SolverType algorithm = SolverType::QUADRA; // Optimizer (default to QUADRA)
  if (as<std::string>(control["method"]) == "FISTA") {
    algorithm = SolverType::FISTA;
  }

  // Variables monitoring the algorithm
  vector<double> gap, timing ; // timings and optimality measures
  vector<uword> status, iactive, ioptim ; // convergence and # of inner/outer iterates
  gap.reserve(lambdas_.size()); timing.reserve(lambdas_.size());
  status.reserve(lambdas_.size()); iactive.reserve(lambdas_.size()); ioptim.reserve(lambdas_.size());

  auto prox = [this](vec x, double l) {
    return(penalty_.proximal(x, l, this->lambda_factor_));
  } ;

  // LAMBDA LOOP
  coef_.set_size(data_.p_, lambdas_.size()) ;
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {
    if (verbose) {Rprintf("\n lambda_linf = %f",lambda_) ;}
    
    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    uword current_it = 0 ;
    double current_gap = datum::inf ;
    do {
      R_CheckUserInterrupt();
      current_it++;
      if (algorithm == SolverType::FISTA) {
        ioptim.push_back(
          solver_.fista(beta_, lambda_, data_.XTy_, data_.XTX_, prox, 1e-5, 10000)
        );
        grad_ = - data_.XTy_ + data_.XTX_ * beta_ ;
        current_gap = optimality_gap(grad_, lambda_) / std::max(1.0, lambda_) ;
        break;
      } else { // QUADRA solver
        try {
          if (current_it == maxiter) {
            throw std::runtime_error("Fail to converge...");
          } else {
            ioptim.push_back(
              solver_.quadratic_breg(beta_, grad_, lambda_, lambda_factor_, data_, unbounded_, accuracy, 10000)
            );
          }
        } catch (std::runtime_error& error) {
          if (verbose) {
            Rprintf("\nNumerical instability: switching to proximal algorithm (slower but safer).");
          }
          current_it = 0; // start this lambda all the way back, with FISTA algorithm
          algorithm = SolverType::FISTA ;
        }
      }

      // OPTIMALITY TESTING (relative to lambda, as the gradient scales with the data)
      grad_ = - data_.XTy_ + data_.XTX_ * beta_ ;
      current_gap = optimality_gap(grad_, lambda_) / std::max(1.0, lambda_) ;
    } while ((current_gap > accuracy) && (current_it <= maxiter));

    // Checking convergence status
    gap.push_back(current_gap) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_gap > accuracy) { status.back() = 1 ; }
    if ((unbounded_.n_elem > maxfeat) & 
        (algorithm == SolverType::QUADRA)) { status.back() = 2 ; }

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      coef_.col(df_.size()) = beta_/data_.norm_X_ ;
      intercept_.push_back(data_.y_bar_ - dot(beta_, data_.X_bar_));
      double max_abs = max(abs(beta_));
      bounded_.push_back(find(abs(beta_) >= max_abs * (1.0 - 1e-10))) ;
      df_.push_back(get_df()) ;
    }

    timing.push_back(timer.toc()) ;
  } // END OF THE LOOP OVER LAMBDA

  lambdas_.resize(df_.size()) ;
  coef_ = coef_.head_cols(df_.size()) ;
  
  return(
    List::create(
      Named("it_active")      = iactive,
      Named("it_optim")       = ioptim ,
      Named("max_grd")        = gap    ,
      Named("convergence")    = status ,
      Named("pensteps_timer") = timing 
    )
  );
}
