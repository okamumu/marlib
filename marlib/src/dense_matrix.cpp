/*
  matrix class
 */

#include "marlib.hpp"
#include "dblas.h"
#include "dlapack.h"

namespace marlib {
  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>::dense_matrix(const range<RangeT>& row, const range<RangeT>& col, const array<ValueT>& a, size_type ld)
    : m_row(row), m_col(col), m_value(a), m_ld(ld) {}

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>::dense_matrix(const dense_matrix<ValueT,RangeT>& m)
  : m_row(m.nrow()), m_col(m.ncol()), m_value(m.m_value), m_ld(m.m_ld) {}

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>::dense_matrix(const dense_matrix<ValueT,RangeT>& m, ValueT* p)
  : m_row(m.nrow()), m_col(m.ncol()), m_value(m.m_value.size(),p), m_ld(m.m_ld) {}

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>::dense_matrix(const csr_matrix<ValueT,RangeT>& m)
  : m_row(m.nrow()), m_col(m.ncol()), m_value(m.nrow()*m.ncol()), m_ld(m.nrow()) {
    m.copyto(*this);
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>::dense_matrix(size_type nrow, size_type ncol)
  : m_row(nrow), m_col(ncol), m_value(nrow*ncol), m_ld(nrow) {}

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>::dense_matrix(size_type nrow, size_type ncol, ValueT* v, size_type ld)
  : m_row(nrow), m_col(ncol), m_value(ld*ncol,v), m_ld(ld) {}

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>::dense_matrix(size_type nrow, size_type ncol, std::initializer_list<ValueT> v)
  : m_row(nrow), m_col(ncol), m_value(v.size()), m_ld(nrow) {
    copyfrom(v.begin());
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>::~dense_matrix() { }

  // template <typename ValueT, typename RangeT>
  // ValueT& dense_matrix<ValueT,RangeT>::operator()(const RangeT i, const RangeT j) {
  //   const RangeT x = i - m_row.begin();
  //   const RangeT y = j - m_col.begin();
  //   return m_value[x + y * m_ld];
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const ValueT& dense_matrix<ValueT,RangeT>::operator()(const RangeT i, const RangeT j) const {
  //   const RangeT x = i - m_row.begin();
  //   const RangeT y = j - m_col.begin();
  //   return m_value[x + y * m_ld];
  // }
  //
  // template <typename ValueT, typename RangeT>
  // ValueT* dense_matrix<ValueT,RangeT>::ptr() {
  //   return &operator()(m_row.begin(), m_col.begin());
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const ValueT* dense_matrix<ValueT,RangeT>::ptr() const {
  //   return &operator()(m_row.begin(), m_col.begin());
  // }
  //
  // template <typename ValueT, typename RangeT>
  // size_type dense_matrix<ValueT,RangeT>::ld() const {
  //   return m_ld;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT dense_matrix<ValueT,RangeT>::rbegin() const {
  //   return m_row.begin();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT dense_matrix<ValueT,RangeT>::rend() const {
  //   return m_row.end();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT dense_matrix<ValueT,RangeT>::cbegin() const {
  //   return m_col.begin();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT dense_matrix<ValueT,RangeT>::cend() const {
  //   return m_col.end();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // size_type dense_matrix<ValueT,RangeT>::nrow() const {
  //   return m_row.size();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // size_type dense_matrix<ValueT,RangeT>::ncol() const {
  //   return m_col.size();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const vector<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator()(const range<RangeT>& row, const RangeT col) const {
  //   array<ValueT> tmp = m_value.subarray(row.begin() - m_row.begin() + (col - m_col.begin()) * m_ld);
  //   return vector<ValueT,RangeT>(row.size(), tmp, 1);
  // }
  //
  // template <typename ValueT, typename RangeT>
  // array<ValueT>& dense_matrix<ValueT,RangeT>::value() {
  //   return m_value;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const array<ValueT>& dense_matrix<ValueT,RangeT>::value() const {
  //   return m_value;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // void dense_matrix<ValueT,RangeT>::set_range(const range<RangeT>& row, const range<RangeT>& col) {
  //   assert(row.size() == m_row.size());
  //   assert(col.size() == m_col.size());
  //   m_row = row;
  //   m_col = col;
  // }

  template <typename ValueT, typename RangeT>
  size_type dense_matrix<ValueT,RangeT>::nnz() const {
    size_type n = 0;
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(); i<=rend(); i++) {
        if (operator()(i,j) != 0) {
          n += 1;
        }
      }
    }
    return n;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator()(const range<RangeT>& row, const range<RangeT>& col) {
    array<ValueT> tmp = m_value.subarray(row.begin() - m_row.begin() + (col.begin() - m_col.begin()) * m_ld);
    return dense_matrix<ValueT,RangeT>(range<RangeT>(row.size()), range<RangeT>(col.size()), tmp, m_ld);
  }

  template <typename ValueT, typename RangeT>
  const dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator()(const range<RangeT>& row, const range<RangeT>& col) const {
    array<ValueT> tmp = m_value.subarray(row.begin() - m_row.begin() + (col.begin() - m_col.begin()) * m_ld);
    return dense_matrix<ValueT,RangeT>(range<RangeT>(row.size()), range<RangeT>(col.size()), tmp, m_ld);
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator()(const RangeT row, const range<RangeT>& col) {
    array<ValueT> tmp = m_value.subarray(row - m_row.begin() + (col.begin() - m_col.begin()) * m_ld);
    return vector<ValueT,RangeT>(col.size(), tmp, m_ld);
  }

  template <typename ValueT, typename RangeT>
  const vector<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator()(const RangeT row, const range<RangeT>& col) const {
    array<ValueT> tmp = m_value.subarray(row - m_row.begin() + (col.begin() - m_col.begin()) * m_ld);
    return vector<ValueT,RangeT>(col.size(), tmp, m_ld);
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator()(const range<RangeT>& row, const RangeT col) {
    array<ValueT> tmp = m_value.subarray(row.begin() - m_row.begin() + (col - m_col.begin()) * m_ld);
    return vector<ValueT,RangeT>(row.size(), tmp, 1);
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::clone() const {
    dense_matrix<ValueT,RangeT> tmp(m_row, m_col, array<ValueT>(m_value.size()), m_ld);
    tmp = *this;
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::clone(ValueT* p) const {
    dense_matrix<ValueT,RangeT> tmp(m_row, m_col, array<ValueT>(m_value.size(), p), m_ld);
    tmp = *this;
    return tmp;
  }

  // equal
  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator=(const ValueT& v) {
    if (nrow() == ld()) {
      std::fill_n(ptr(), nrow()*ncol(), v);
    } else {
      for (RangeT j=cbegin(); j<=cend(); j++) {
        for (RangeT i=rbegin(); i<=rend(); i++) {
          operator()(i,j) = v;
        }
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator=(const dense_matrix<ValueT,RangeT>& x) {
    assert(nrow() == x.nrow());
    assert(ncol() == x.ncol());
    for (RangeT j=cbegin(), jx=x.cbegin(); j<=cend(); j++, jx++) {
      for (RangeT i=rbegin(), ix=x.rbegin(); i<=rend(); i++, ix++) {
        operator()(i,j) = x(ix,jx);
      }
    }
    return *this;
  }

  // arithmetic operators

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator+=(const dense_matrix<ValueT,RangeT>& x) {
    assert(nrow() == x.nrow());
    assert(ncol() == x.ncol());
    for (RangeT j=cbegin(), jx=x.cbegin(); j<=cend(); j++, jx++) {
      for (RangeT i=rbegin(), ix=x.rbegin(); i<=rend(); i++, ix++) {
        operator()(i,j) += x(ix,jx);
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator-=(const dense_matrix<ValueT,RangeT>& x) {
    assert(nrow() == x.nrow());
    assert(ncol() == x.ncol());
    for (RangeT j=cbegin(), jx=x.cbegin(); j<=cend(); j++, jx++) {
      for (RangeT i=rbegin(), ix=x.rbegin(); i<=rend(); i++, ix++) {
        operator()(i,j) -= x(ix,jx);
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator*=(const dense_matrix<ValueT,RangeT>& x) {
    assert(nrow() == x.nrow());
    assert(ncol() == x.ncol());
    for (RangeT j=cbegin(), jx=x.cbegin(); j<=cend(); j++, jx++) {
      for (RangeT i=rbegin(), ix=x.rbegin(); i<=rend(); i++, ix++) {
        operator()(i,j) *= x(ix,jx);
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator/=(const dense_matrix<ValueT,RangeT>& x) {
    assert(nrow() == x.nrow());
    assert(ncol() == x.ncol());
    for (RangeT j=cbegin(), jx=x.cbegin(); j<=cend(); j++, jx++) {
      for (RangeT i=rbegin(), ix=x.rbegin(); i<=rend(); i++, ix++) {
        operator()(i,j) /= x(ix,jx);
      }
    }
    return *this;
  }

  // arithmetic operators 2

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator+=(const vector<ValueT,RangeT>& x) {
    assert(nrow() == x.size());
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(), ix=x.begin(); i<=rend(); i++, ix++) {
        operator()(i,j) += x(ix);
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator-=(const vector<ValueT,RangeT>& x) {
    assert(nrow() == x.size());
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(), ix=x.begin(); i<=rend(); i++, ix++) {
        operator()(i,j) -= x(ix);
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator*=(const vector<ValueT,RangeT>& x) {
    assert(nrow() == x.size());
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(), ix=x.begin(); i<=rend(); i++, ix++) {
        operator()(i,j) *= x(ix);
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator/=(const vector<ValueT,RangeT>& x) {
    assert(nrow() == x.size());
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(), ix=x.begin(); i<=rend(); i++, ix++) {
        operator()(i,j) /= x(ix);
      }
    }
    return *this;
  }

  // arithmetic operators 3

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator+=(const ValueT& v) {
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(); i<=rend(); i++) {
        operator()(i,j) += v;
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator-=(const ValueT& v) {
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(); i<=rend(); i++) {
        operator()(i,j) -= v;
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator*=(const ValueT& v) {
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(); i<=rend(); i++) {
        operator()(i,j) *= v;
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::operator/=(const ValueT& v) {
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(); i<=rend(); i++) {
        operator()(i,j) /= v;
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator+(const dense_matrix<ValueT,RangeT>& v) const {
    dense_matrix<ValueT,RangeT> result = clone();
    result += v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator-(const dense_matrix<ValueT,RangeT>& v) const {
    dense_matrix<ValueT,RangeT> result = clone();
    result -= v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator*(const dense_matrix<ValueT,RangeT>& v) const {
    dense_matrix<ValueT,RangeT> result = clone();
    result *= v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::operator/(const dense_matrix<ValueT,RangeT>& v) const {
    dense_matrix<ValueT,RangeT> result = clone();
    result /= v;
    return result;
  }

  ////// print

  template <typename ValueT, typename RangeT>
  std::ostream& dense_matrix<ValueT,RangeT>::print(std::ostream& os) const {
    for (RangeT i=rbegin(); i<=rend(); i++) {
      for (RangeT j=cbegin(); j<=cend(); j++) {
        os << operator()(i,j) << " ";
      }
      os << std::endl;
    }
    return os;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::copyfrom(const ValueT* x) {
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(); i<=rend(); i++, x++) {
        operator()(i,j) = *x;
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  const vector<ValueT,RangeT> dense_matrix<ValueT,RangeT>::diag() const {
    assert(m_row == m_col);
    vector<ValueT,RangeT> tmp(nrow(), m_value, m_ld+1);
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> dense_matrix<ValueT,RangeT>::diag() {
    assert(m_row == m_col);
    vector<ValueT,RangeT> tmp(nrow(), m_value, m_ld+1);
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& dense_matrix<ValueT,RangeT>::diag(vector<ValueT*,RangeT>& x, int offset, const RangeT xindex) {
    RangeT ix=xindex, i, j;
    if (offset >= 0) {
      i=rbegin() + offset;
      j=cbegin();
    } else {
      i=rbegin();
      j=cbegin() - offset;
    }
    for (; ix<=x.end() && i<=rend() && j<=cend(); ix++, i++, j++) {
      x.ptr(ix) = &operator()(i,j);
    }
    return x;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::eye(size_type n) {
    dense_matrix<ValueT,RangeT> m(n,n);
    return m.eye();
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::zeros(size_type nrow, size_type ncol) {
    dense_matrix<ValueT,RangeT> m(nrow,ncol);
    m = 0;
    return m;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::eye() {
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(); i<=rend(); i++) {
        if (j == i) {
          operator()(i,j) = 1;
        } else {
          operator()(i,j) = 0;
        }
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT> dense_matrix<ValueT,RangeT>::tr() const {
    dense_matrix<ValueT,RangeT> tmp(ncol(), nrow());
    for (RangeT j=cbegin(), ii=tmp.rbegin(); j<=cend(); j++, ii++) {
      for (RangeT i=rbegin(), jj=tmp.cbegin(); i<=rend(); i++, jj++) {
        tmp(ii,jj) = operator()(i,j);
      }
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& dense_matrix<ValueT,RangeT>::row_sum(vector<ValueT,RangeT>& x) const {
    assert(x.size() == nrow());
    for (RangeT j=cbegin(); j<=cend(); j++) {
      for (RangeT i=rbegin(), ix=x.begin(); i<=rend(); i++, ix++) {
        x(ix) += operator()(i,j);
      }
    }
    return x;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> dense_matrix<ValueT,RangeT>::row_sum() const {
    vector<ValueT,RangeT> x(nrow());
    x = 0;
    return row_sum(x);
  }

#ifdef F77BLAS
  template <>
  dense_matrix<double,int>& dense_matrix<double,int>::operator+=(const dense_matrix<double,int>& v) {
    assert(nrow() == v.nrow());
    assert(ncol() == v.ncol());
    if (nrow() == ld() && v.nrow() == v.ld()) {
      dblas::daxpy(nrow()*ncol(), 1, v.ptr(), 1, ptr(), 1);
    } else {
      for (int j=cbegin(), jj=v.cbegin(); j<=cend(); j++, jj++) {
        dblas::daxpy(nrow(), 1, &v(v.rbegin(),jj), 1, &operator()(rbegin(),j), 1);
      }
    }
    return *this;
  }

  template <>
  dense_matrix<double,int>& dense_matrix<double,int>::operator-=(const dense_matrix<double,int>& v) {
    assert(nrow() == v.nrow());
    assert(ncol() == v.ncol());
    if (nrow() == ld() && v.nrow() == v.ld()) {
      dblas::daxpy(nrow()*ncol(), -1, v.ptr(), 1, ptr(), 1);
    } else {
      for (int j=cbegin(), jj=v.cbegin(); j<=cend(); j++, jj++) {
        dblas::daxpy(nrow(), -1, &v(v.rbegin(),jj), 1, &operator()(rbegin(),j), 1);
      }
    }
    return *this;
  }

  template <>
  dense_matrix<double,int>& dense_matrix<double,int>::operator*=(const double& v) {
    if (nrow() == ld()) {
      dblas::dscal(nrow()*ncol(), v, ptr(), 1);
    } else {
      for (int j=cbegin(); j<=cend(); j++) {
        dblas::dscal(nrow(), v, &operator()(rbegin(),j), 1);
      }
    }
    return *this;
  }

  template <>
  dense_matrix<double,int>& dense_matrix<double,int>::operator/=(const double& v) {
    if (nrow() == ld()) {
      dblas::dscal(nrow()*ncol(), 1/v, ptr(), 1);
    } else {
      for (int j=cbegin(); j<=cend(); j++) {
        dblas::dscal(nrow(), 1/v, &operator()(rbegin(),j), 1);
      }
    }
    return *this;
  }

  template <>
  dense_matrix<double,int>& dense_matrix<double,int>::operator=(const dense_matrix<double,int>& v) {
    assert(nrow() == v.nrow());
    assert(ncol() == v.ncol());
    if (nrow() == ld() && v.nrow() == v.ld()) {
      dblas::dcopy(nrow()*ncol(), v.ptr(), 1, ptr(), 1);
    } else {
      for (int j=cbegin(), jj=v.cbegin(); j<=cend(); j++, jj++) {
        dblas::dcopy(nrow(), &v(v.rbegin(),jj), 1, &operator()(rbegin(),j), 1);
      }
    }
    return *this;
  }
#endif

  // instance
  template class dense_matrix<double,int>;
}
