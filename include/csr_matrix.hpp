/*
  vector class
 */

#pragma once

#include <iostream>
#include <initializer_list>

#include "range.hpp"
#include "array.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT> class vector;
  template <typename ValueT, typename RangeT> class dense_matrix;

  template <typename ValueT, typename RangeT>
  class csr_matrix {
  public:
    using ValueType = ValueT;
    using RangeType = RangeT;

    csr_matrix(size_type row, size_type col, size_type nnz,
      RangeT* rowptr, RangeT* colind, ValueT* value);
    csr_matrix(const csr_matrix<ValueT,RangeT>& m);
    csr_matrix(const dense_matrix<ValueT,RangeT>& m);
    ~csr_matrix();

  private:
    csr_matrix(const range<RangeT>& row, const range<RangeT>& col,
      size_type nnz, const array<RangeT>& rowptr, const array<RangeT>& colind, const array<ValueT>& value);

    range<RangeT> m_row;
    range<RangeT> m_col;

    size_type m_nnz;
    array<RangeT> m_rowptr;
    array<RangeT> m_colind;
    array<ValueT> m_value;
    RangeT m_origin;

  public:
    ValueT* ptr();
    const ValueT* ptr() const;

    const RangeT& rbegin() const;
    const RangeT& rend() const;
    const RangeT& cbegin() const;
    const RangeT& cend() const;

    size_type nrow() const;
    size_type ncol() const;
    size_type nnz() const;

    csr_matrix<ValueT,RangeT> clone() const;

    const array<RangeT>& rowptr() const;
    const array<RangeT>& colind() const;
    const array<ValueT>& value() const;

    const RangeT& origin() const;
    size_type rowptr(const size_type& i) const;
    size_type colind(const size_type& i) const;
    ValueT& value(const size_type& i);
    const ValueT& value(const size_type& i) const;

    // equal
    csr_matrix<ValueT,RangeT>& operator=(const ValueT& v);
    csr_matrix<ValueT,RangeT>& operator=(const csr_matrix<ValueT,RangeT>& v);

    // arithmetic operators
    csr_matrix<ValueT,RangeT>& operator+=(const csr_matrix<ValueT,RangeT>& x);
    csr_matrix<ValueT,RangeT>& operator-=(const csr_matrix<ValueT,RangeT>& x);
    csr_matrix<ValueT,RangeT>& operator*=(const csr_matrix<ValueT,RangeT>& x);
    csr_matrix<ValueT,RangeT>& operator/=(const csr_matrix<ValueT,RangeT>& x);

    csr_matrix<ValueT,RangeT>& operator+=(const vector<ValueT,RangeT>& x);
    csr_matrix<ValueT,RangeT>& operator-=(const vector<ValueT,RangeT>& x);
    csr_matrix<ValueT,RangeT>& operator*=(const vector<ValueT,RangeT>& x);
    csr_matrix<ValueT,RangeT>& operator/=(const vector<ValueT,RangeT>& x);

    csr_matrix<ValueT,RangeT>& operator+=(const ValueT& v);
    csr_matrix<ValueT,RangeT>& operator-=(const ValueT& v);
    csr_matrix<ValueT,RangeT>& operator*=(const ValueT& v);
    csr_matrix<ValueT,RangeT>& operator/=(const ValueT& v);

    csr_matrix<ValueT,RangeT> operator+(const csr_matrix<ValueT,RangeT>& v) const;
    csr_matrix<ValueT,RangeT> operator-(const csr_matrix<ValueT,RangeT>& v) const;
    csr_matrix<ValueT,RangeT> operator*(const csr_matrix<ValueT,RangeT>& v) const;
    csr_matrix<ValueT,RangeT> operator/(const csr_matrix<ValueT,RangeT>& v) const;

    ////// print
    std::ostream& print(std::ostream& os) const;

    template <typename ValueTT, typename RangeTT>
    friend std::ostream& operator<< (std::ostream& os, const csr_matrix<ValueTT,RangeTT>& m);

    // utils

    const vector<ValueT*,RangeT> diag() const;
    vector<ValueT*,RangeT> diag();
    vector<ValueT*,RangeT>& diag(vector<ValueT*,RangeT>& x, int offset, const RangeT& xindex);
    csr_matrix<ValueT,RangeT>& eye();

    vector<ValueT,RangeT> row_sum() const;
    vector<ValueT,RangeT>& row_sum(vector<ValueT,RangeT>& x) const;

    csr_matrix<ValueT,RangeT>& copyfrom(const dense_matrix<ValueT,RangeT>& x);
    dense_matrix<ValueT,RangeT>& copyto(dense_matrix<ValueT,RangeT>& x) const;
  };

  template <typename ValueT, typename RangeT>
  std::ostream& operator<<(std::ostream& os, const csr_matrix<ValueT,RangeT>& m) {
    return m.print(os);
  }
}
