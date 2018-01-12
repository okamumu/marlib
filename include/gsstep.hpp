#pragma once

namespace marlib {

  /*
  !  SOR step for solving the following linear equation
  !
  !         alpha * trans(A - sigma I) * x = b
  !
  !         The step computes
  !
  !         fwd
  !         x := (D/omega + L)^(-1) (b/alpha - (U - D (1-omega)/omega - sigma I) * x)
  !
  !         bwd
  !         x := (D/omega + U)^(-1) (b/alpha - (L - D (1-omega)/omega - sigma I) * x)
  !
  !         A: square matrix
  !         x: vector (in; initial vector for the step, out; updated vector)
  !         b: constant vector
  */

  template <typename ValueT, typename RangeT>
  void gsstep_impl(const trans_t& trans, const ValueT& alpha,
    dense_matrix<ValueT,RangeT> A,
    const ValueT& sigma,
    const ValueT& omega,
    vector<ValueT,RangeT> b, vector<ValueT,RangeT> x);

  template <typename ValueT, typename RangeT>
  void gsstep_impl(const trans_t& trans, const ValueT& alpha,
    const csr_matrix<ValueT, RangeT>& A,
    const ValueT& sigma,
    const ValueT& omega,
    vector<ValueT,RangeT> b, vector<ValueT,RangeT> x);

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& gsstep(const trans_t& trans, const ValueT& alpha,
    const dense_matrix<ValueT,RangeT>& A, const ValueT& sigma, const ValueT& omega,
    const vector<ValueT,RangeT>& b, vector<ValueT,RangeT>& x) {
    gsstep_impl(trans, alpha, A, sigma, omega, b, x);
    return x;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& gsstep(const trans_t& trans, const ValueT& alpha,
    const csr_matrix<ValueT,RangeT>& A, const ValueT& sigma, const ValueT& omega,
    const vector<ValueT,RangeT>& b, vector<ValueT,RangeT>& x) {
    gsstep_impl(trans, alpha, A, sigma, omega, b, x);
    return x;
  }

}
