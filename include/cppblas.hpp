/*
  blas.hpp
*/

#pragma once

#include "types.h"


namespace marlib {

  // BLAS level 1

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& daxpy(const ValueT& alpha, const vector<ValueT,RangeT>& x, vector<ValueT,RangeT>& y);
  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& daxpy(const ValueT& alpha, const dense_matrix<ValueT,RangeT>& x, dense_matrix<ValueT,RangeT>& y);
  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& daxpy(const ValueT& alpha, const csr_matrix<ValueT,RangeT>& x, csr_matrix<ValueT,RangeT>& y);

  template <typename ValueT, typename RangeT>
  ValueT ddot(const vector<ValueT,RangeT>& x, const vector<ValueT,RangeT>& y);
  template <typename ValueT, typename RangeT>
  ValueT ddot(const dense_matrix<ValueT,RangeT>& x, const dense_matrix<ValueT,RangeT>& y);
  template <typename ValueT, typename RangeT>
  ValueT ddot(const csr_matrix<ValueT,RangeT>& x, const csr_matrix<ValueT,RangeT>& y);

  template <typename ValueT, typename RangeT>
  ValueT dnrm2(const vector<ValueT,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT dnrm2(const dense_matrix<ValueT,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT dnrm2(const csr_matrix<ValueT,RangeT>& x);

  template <typename ValueT, typename RangeT>
  ValueT dsum(const vector<ValueT,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT dsum(const dense_matrix<ValueT,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT dsum(const csr_matrix<ValueT,RangeT>& x);

  template <typename ValueT, typename RangeT>
  ValueT dasum(const vector<ValueT,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT dasum(const dense_matrix<ValueT,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT dasum(const csr_matrix<ValueT,RangeT>& x);

  template <typename ValueT, typename RangeT>
  RangeT iamax(const vector<ValueT,RangeT>& x);

  template <typename ValueT, typename RangeT>
  ValueT damax(const vector<ValueT,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT damax(const vector<ValueT*,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT damax(const dense_matrix<ValueT,RangeT>& x);
  template <typename ValueT, typename RangeT>
  ValueT damax(const csr_matrix<ValueT,RangeT>& x);

  // BLAS level 2

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& dgemv(const trans_t& trans, const ValueT& alpha,
    const dense_matrix<ValueT,RangeT>& A, const vector<ValueT,RangeT>& x,
    const ValueT& beta, vector<ValueT,RangeT>& y);
  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& dgemv(const trans_t& trans, const ValueT& alpha,
    const csr_matrix<ValueT,RangeT>& A, const vector<ValueT,RangeT>& x,
    const ValueT& beta, vector<ValueT,RangeT>& y);

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dger(const trans_t& trans,
    const ValueT& alpha, const vector<ValueT,RangeT>& x, const vector<ValueT,RangeT>& y,
    dense_matrix<ValueT,RangeT>& A);
  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& dger(const trans_t& trans,
    const ValueT& alpha, const vector<ValueT,RangeT>& x, const vector<ValueT,RangeT>& y,
    csr_matrix<ValueT,RangeT>& A);

  // BLAS level 3

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dgemm(
    const trans_t& transA,
    const trans_t& transB,
    const ValueT& alpha,
    const dense_matrix<ValueT,RangeT>& A,
    const dense_matrix<ValueT,RangeT>& B,
    const ValueT& beta,
    dense_matrix<ValueT,RangeT>& C);

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dgemm(
    const trans_t& transA,
    const trans_t& transB,
    const ValueT& alpha,
    const csr_matrix<ValueT,RangeT>& A,
    const dense_matrix<ValueT,RangeT>& B,
    const ValueT& beta,
    dense_matrix<ValueT,RangeT>& C);

  // dgemm for vector

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& dgemm(
    const trans_t& transA,
    const trans_t& transB,
    const ValueT& alpha,
    const dense_matrix<ValueT,RangeT>& A,
    const vector<ValueT,RangeT>& B,
    const ValueT& beta,
    vector<ValueT,RangeT>& C) {
    return dgemv<ValueT,RangeT>(transA, alpha, A, B, beta, C);
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& dgemm(
    const trans_t& transA,
    const trans_t& transB,
    const ValueT& alpha,
    const csr_matrix<ValueT,RangeT>& A,
    const vector<ValueT,RangeT>& B,
    const ValueT& beta,
    vector<ValueT,RangeT>& C) {
    return dgemv<ValueT,RangeT>(transA, alpha, A, B, beta, C);
  }

  // lapack

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dgesv(
    const dense_matrix<ValueT,RangeT>& A, dense_matrix<ValueT,RangeT>& B);

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& dgesv(
    const dense_matrix<ValueT,RangeT>& A, vector<ValueT,RangeT>& B);

}
