/*
! Description: matrix power
!
!        ME = MA^m
!
*/

#pragma once

namespace marlib {

  template <typename ValueT, typename RangeT, typename MatrixT>
  dense_matrix<ValueT,RangeT>& mexp_pade(const trans_t& trans,
    const MatrixT& MA, const ValueT& t, dense_matrix<ValueT,RangeT>& ME, const ValueT& eps);

  template <typename ValueT, typename RangeT, typename MatrixT, typename VectorT>
  VectorT& mexp_unif(const trans_t& trans, const MatrixT& P, const ValueT& qv,
  const poisson<ValueT,RangeT>& pois, const VectorT& x, VectorT& y,
  const ValueT& atol);

  template <typename ValueT, typename RangeT, typename MatrixT, typename VectorT>
  VectorT& mexpint_unif(const trans_t& trans, const MatrixT& P, const ValueT& qv,
  const poisson<ValueT,RangeT>& pois, const VectorT& x, VectorT& y, VectorT& cy);

  template <typename ValueT, typename RangeT, typename MatrixT, typename VectorT>
  VectorT& mexpconv_unif(const trans_t& transQ, const trans_t& transH,
  const MatrixT& P, const ValueT& qv, const poisson<ValueT,RangeT>& pois,
  const VectorT& x, const VectorT& y, VectorT& z, MatrixT& H);

}
