/*
  dense_matrix
*/

#pragma once

#include <iostream>
#include <initializer_list>

#include "range.hpp"
#include "array.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT> class vector;
  template <typename ValueT, typename RangeT> class csr_matrix;

  template <typename ValueT, typename RangeT>
  class dense_matrix {
  public:
    using ValueType = ValueT;
    using RangeType = RangeT;

    dense_matrix(const dense_matrix<ValueT,RangeT>& m);
    dense_matrix(const dense_matrix<ValueT,RangeT>& m, ValueT* p);
    dense_matrix(const csr_matrix<ValueT,RangeT>& m);

    dense_matrix(size_type nrow, size_type ncol);
    dense_matrix(size_type nrow, size_type ncol, ValueT* v, size_type ld);
    dense_matrix(size_type nrow, size_type ncol, std::initializer_list<ValueT> v);
    ~dense_matrix();

  private:
    dense_matrix(const range<RangeT>& row, const range<RangeT>& col, const array<ValueT>& a, size_type ld);
    range<RangeT> m_row;
    range<RangeT> m_col;
    array<ValueT> m_value;
    size_type m_ld;

  public:
    ValueT& operator()(const RangeT& i, const RangeT& j);
    const ValueT& operator()(const RangeT& i, const RangeT& j) const;

    ValueT* ptr();
    const ValueT* ptr() const;

    const RangeT& rbegin() const;
    const RangeT& rend() const;
    const RangeT& cbegin() const;
    const RangeT& cend() const;

    size_type nrow() const;
    size_type ncol() const;
    size_type nnz() const;

    size_type ld() const;

    dense_matrix<ValueT,RangeT> operator()(const range<RangeT>& row, const range<RangeT>& col);
    const dense_matrix<ValueT,RangeT> operator()(const range<RangeT>& row, const range<RangeT>& col) const;

    vector<ValueT,RangeT> operator()(const RangeT& row, const range<RangeT>& col);
    const vector<ValueT,RangeT> operator()(const RangeT& row, const range<RangeT>& col) const;
    vector<ValueT,RangeT> operator()(const range<RangeT>& row, const RangeT& col);
    const vector<ValueT,RangeT> operator()(const range<RangeT>& row, const RangeT& col) const;

    array<ValueT>& value();
    const array<ValueT>& value() const;

    void set_range(const range<RangeT>& row, const range<RangeT>& col);
    dense_matrix<ValueT,RangeT> clone() const;

    // equal
    dense_matrix<ValueT,RangeT>& operator=(const ValueT& v);
    dense_matrix<ValueT,RangeT>& operator=(const dense_matrix<ValueT,RangeT>& v);

    // arithmetic operators
    dense_matrix<ValueT,RangeT>& operator+=(const dense_matrix<ValueT,RangeT>& x);
    dense_matrix<ValueT,RangeT>& operator-=(const dense_matrix<ValueT,RangeT>& x);
    dense_matrix<ValueT,RangeT>& operator*=(const dense_matrix<ValueT,RangeT>& x);
    dense_matrix<ValueT,RangeT>& operator/=(const dense_matrix<ValueT,RangeT>& x);

    dense_matrix<ValueT,RangeT>& operator+=(const vector<ValueT,RangeT>& v);
    dense_matrix<ValueT,RangeT>& operator-=(const vector<ValueT,RangeT>& v);
    dense_matrix<ValueT,RangeT>& operator*=(const vector<ValueT,RangeT>& v);
    dense_matrix<ValueT,RangeT>& operator/=(const vector<ValueT,RangeT>& v);

    dense_matrix<ValueT,RangeT>& operator+=(const ValueT& v);
    dense_matrix<ValueT,RangeT>& operator-=(const ValueT& v);
    dense_matrix<ValueT,RangeT>& operator*=(const ValueT& v);
    dense_matrix<ValueT,RangeT>& operator/=(const ValueT& v);

    dense_matrix<ValueT,RangeT> operator+(const dense_matrix<ValueT,RangeT>& v) const;
    dense_matrix<ValueT,RangeT> operator-(const dense_matrix<ValueT,RangeT>& v) const;
    dense_matrix<ValueT,RangeT> operator*(const dense_matrix<ValueT,RangeT>& v) const;
    dense_matrix<ValueT,RangeT> operator/(const dense_matrix<ValueT,RangeT>& v) const;

    ////// print
    std::ostream& print(std::ostream& os) const;

    template <typename ValueTT, typename RangeTT>
    friend std::ostream& operator<< (std::ostream& os, const dense_matrix<ValueTT,RangeTT>& m);

    // utils

    static dense_matrix<ValueT,RangeT> eye(size_type n);
    static dense_matrix<ValueT,RangeT> zeros(size_type nrow, size_type ncol);

    const vector<ValueT,RangeT> diag() const;
    vector<ValueT,RangeT> diag();
    vector<ValueT*,RangeT>& diag(vector<ValueT*,RangeT>& x, int offset, const RangeT& xindex);
    dense_matrix<ValueT,RangeT>& eye();
    dense_matrix<ValueT,RangeT> tr() const;

    vector<ValueT,RangeT> row_sum() const;
    vector<ValueT,RangeT>& row_sum(vector<ValueT,RangeT>& x) const;

    dense_matrix<ValueT,RangeT>& copyfrom(const ValueT* x);
  };

  template <typename ValueT, typename RangeT>
  std::ostream& operator<<(std::ostream& os, const dense_matrix<ValueT,RangeT>& m) {
    return m.print(os);
  }

}
