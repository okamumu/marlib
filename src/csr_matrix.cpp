/*
  matrix class
 */

#include <cassert>

#include "marlib.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT>
  const RangeT csr_matrix<ValueT,RangeT>::default_origin = range<RangeT>::default_origin;

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>::csr_matrix(const range<RangeT>& row, const range<RangeT>& col,
    size_type nnz, const array<RangeT>& rowptr, const array<RangeT>& colind, const array<ValueT>& value,
    const RangeT origin)
  : m_row(row), m_col(col), m_nnz(nnz), m_rowptr(rowptr), m_colind(colind), m_value(value),
  m_origin(origin) {}

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>::csr_matrix(const dense_matrix<ValueT,RangeT>& m, const RangeT origin)
  : m_row(m.nrow()), m_col(m.ncol()), m_nnz(m.nnz()), m_rowptr(m.nrow()+1), m_colind(m_nnz), m_value(m_nnz),
  m_origin(origin) {
    copyfrom(m);
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>::csr_matrix(const csr_matrix<ValueT,RangeT>& m)
  : m_row(m.nrow()), m_col(m.ncol()), m_nnz(m.m_nnz), m_rowptr(m.m_rowptr), m_colind(m.m_colind), m_value(m.m_value),
  m_origin(m.m_origin) {}

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>::csr_matrix(size_type row, size_type col, size_type nnz,
    RangeT* rowptr, RangeT* colind, ValueT* value, const RangeT origin)
  : m_row(row), m_col(col), m_nnz(nnz), m_rowptr(row+1,rowptr), m_colind(nnz,colind), m_value(nnz,value), m_origin(origin) { }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>::~csr_matrix() { }

  // template <typename ValueT, typename RangeT>
  // ValueT* csr_matrix<ValueT,RangeT>::ptr() {
  //   return &m_value[0];
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const ValueT* csr_matrix<ValueT,RangeT>::ptr() const {
  //   return &m_value[0];
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT csr_matrix<ValueT,RangeT>::rbegin() const {
  //   return m_row.begin();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT csr_matrix<ValueT,RangeT>::rend() const {
  //   return m_row.end();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT csr_matrix<ValueT,RangeT>::cbegin() const {
  //   return m_col.begin();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT csr_matrix<ValueT,RangeT>::cend() const {
  //   return m_col.end();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // size_type csr_matrix<ValueT,RangeT>::nrow() const {
  //   return m_row.size();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // size_type csr_matrix<ValueT,RangeT>::ncol() const {
  //   return m_col.size();
  // }
  //
  // template <typename ValueT, typename RangeT>
  // size_type csr_matrix<ValueT,RangeT>::nnz() const {
  //   return m_nnz;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const array<RangeT>& csr_matrix<ValueT,RangeT>::rowptr() const {
  //   return m_rowptr;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const array<RangeT>& csr_matrix<ValueT,RangeT>::colind() const {
  //   return m_colind;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const array<ValueT>& csr_matrix<ValueT,RangeT>::value() const {
  //   return m_value;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const RangeT csr_matrix<ValueT,RangeT>::origin() const {
  //   return m_origin;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // size_type csr_matrix<ValueT,RangeT>::rowptr(const size_type i) const {
  //   return m_rowptr[i] - m_origin;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // size_type csr_matrix<ValueT,RangeT>::colind(const size_type i) const {
  //   return m_colind[i] - m_origin;
  // }
  //
  // template <typename ValueT, typename RangeT>
  // ValueT csr_matrix<ValueT,RangeT>::value(const size_type i) {
  //   return m_value[i];
  // }
  //
  // template <typename ValueT, typename RangeT>
  // const ValueT csr_matrix<ValueT,RangeT>::value(const size_type i) const {
  //   return m_value[i];
  // }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT> csr_matrix<ValueT,RangeT>::clone() const {
    csr_matrix<ValueT,RangeT> tmp(m_row, m_col, m_nnz, m_rowptr, m_colind, array<ValueT>(m_value.size()), m_origin);
    tmp = *this;
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT> csr_matrix<ValueT,RangeT>::clone(ValueT* p) const {
    csr_matrix<ValueT,RangeT> tmp(m_row, m_col, m_nnz, m_rowptr, m_colind, array<ValueT>(m_value.size(), p), m_origin);
    tmp = *this;
    return tmp;
  }

  // equal
  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator=(const ValueT& v) {
    for (size_type i=0; i<m_value.size(); i++) {
      m_value[i] = v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator=(const csr_matrix<ValueT,RangeT>& x) {
    assert(nnz() == x.nnz());
    for (size_type i=0; i<m_value.size(); i++) {
      m_value[i] = x.m_value[i];
    }
    return *this;
  }

  // arithmetic operators 1

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator+=(const csr_matrix<ValueT,RangeT>& x) {
    assert(nnz() == x.nnz());
    for (size_type i=0; i<nnz(); i++) {
      m_value[i] += x.m_value[i];
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator-=(const csr_matrix<ValueT,RangeT>& x) {
    assert(nnz() == x.nnz());
    for (size_type i=0; i<nnz(); i++) {
      m_value[i] -= x.m_value[i];
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator*=(const csr_matrix<ValueT,RangeT>& x) {
    assert(nnz() == x.nnz());
    for (size_type i=0; i<nnz(); i++) {
      m_value[i] *= x.m_value[i];
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator/=(const csr_matrix<ValueT,RangeT>& x) {
    assert(nnz() == x.nnz());
    for (size_type i=0; i<nnz(); i++) {
      m_value[i] /= x.m_value[i];
    }
    return *this;
  }

  // arithmetic operators 2

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator+=(const vector<ValueT,RangeT>& x) {
    assert(nrow() == x.size());
    for (size_type i=0; i<nrow(); i++) {
      for (RangeT z=rowptr(i); z<rowptr(i+1); z++) {
        value(z) += x(i + x.begin());
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator-=(const vector<ValueT,RangeT>& x) {
    assert(nrow() == x.size());
    for (size_type i=0; i<nrow(); i++) {
      for (RangeT z=rowptr(i); z<rowptr(i+1); z++) {
        value(z) -= x(i + x.begin());
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator*=(const vector<ValueT,RangeT>& x) {
    assert(nrow() == x.size());
    for (size_type i=0; i<nrow(); i++) {
      for (RangeT z=rowptr(i); z<rowptr(i+1); z++) {
        value(z) *= x(i + x.begin());
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator/=(const vector<ValueT,RangeT>& x) {
    assert(nrow() == x.size());
    for (size_type i=0; i<nrow(); i++) {
      for (RangeT z=rowptr(i); z<rowptr(i+1); z++) {
        value(z) /= x(i + x.begin());
      }
    }
    return *this;
  }

  // arithmetic operators 3

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator+=(const ValueT& v) {
    for (size_type i=0; i<nnz(); i++) {
      m_value[i] += v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator-=(const ValueT& v) {
    for (size_type i=0; i<nnz(); i++) {
      m_value[i] -= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator*=(const ValueT& v) {
    for (size_type i=0; i<nnz(); i++) {
      m_value[i] *= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::operator/=(const ValueT& v) {
    for (size_type i=0; i<nnz(); i++) {
      m_value[i] /= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT> csr_matrix<ValueT,RangeT>::operator+(const csr_matrix<ValueT,RangeT>& v) const {
    csr_matrix<ValueT,RangeT> result = clone();
    result += v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT> csr_matrix<ValueT,RangeT>::operator-(const csr_matrix<ValueT,RangeT>& v) const {
    csr_matrix<ValueT,RangeT> result = clone();
    result -= v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT> csr_matrix<ValueT,RangeT>::operator*(const csr_matrix<ValueT,RangeT>& v) const {
    csr_matrix<ValueT,RangeT> result = clone();
    result *= v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT> csr_matrix<ValueT,RangeT>::operator/(const csr_matrix<ValueT,RangeT>& v) const {
    csr_matrix<ValueT,RangeT> result = clone();
    result /= v;
    return result;
  }

  ////// print
  template <typename ValueT, typename RangeT>
  std::ostream& csr_matrix<ValueT,RangeT>::print(std::ostream& os) const {
    for (size_type i=0; i<nrow(); i++) {
      for (size_type z=rowptr(i); z<rowptr(i+1); z++) {
        os << "(" << i + origin() << "," << colind(z) + origin() << ")=" << m_value[z] << std::endl;
      }
    }
    return os;
  }

  // utils

  template <typename ValueT, typename RangeT>
  const vector<ValueT*,RangeT> csr_matrix<ValueT,RangeT>::diag() const {
    assert(m_row == m_col);
    vector<ValueT*,RangeT> tmp(nrow());
    RangeT v = tmp.begin();
    for (size_type i=0; i<nrow(); i++) {
      bool find_diag = false;
      for (size_type z=rowptr(i); z<rowptr(i+1); z++) {
        size_type j = colind(z);
        if (i == j) {
          tmp.ptr(v) = const_cast<ValueT*>(&m_value[z]);
          v++;
          find_diag = true;
          break;
        }
      }
      assert(find_diag == true);
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT> csr_matrix<ValueT,RangeT>::diag() {
    assert(m_row == m_col);
    vector<ValueT*,RangeT> tmp(nrow());
    RangeT v = tmp.begin();
    for (size_type i=0; i<nrow(); i++) {
      bool find_diag = false;
      for (size_type z=rowptr(i); z<rowptr(i+1); z++) {
        size_type j = colind(z);
        if (i == j) {
          tmp.ptr(v) = const_cast<ValueT*>(&m_value[z]);
          v++;
          find_diag = true;
          break;
        }
      }
      assert(find_diag == true);
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& csr_matrix<ValueT,RangeT>::diag(vector<ValueT*,RangeT>& x, int offset, const RangeT xindex) {
    size_type i;
    if (offset >= 0) {
      i = offset;
    } else {
      i = 0;
    }
    RangeT v = xindex;
    for (; i<nrow() && i<nrow()+offset; i++) {
      bool find_diag = false;
      for (size_type z=rowptr(i); z<rowptr(i+1); z++) {
        size_type j = colind(z);
        if (i - offset == j) {
          x.ptr(v) = &m_value[z];
          v++;
          find_diag = true;
          break;
        }
      }
      assert(find_diag == true);
    }
    return x;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::eye() {
    assert(m_row == m_col);
    for (size_type i=0; i<nrow(); i++) {
      bool find_diag = false;
      for (size_type z=rowptr(i); z<rowptr(i+1); z++) {
        size_type j = colind(z);
        if (i == j) {
          m_value[z] = 1;
          find_diag = true;
        } else {
          m_value[z] = 0;
        }
      }
      assert(find_diag == true);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::row_sum(vector<ValueT,RangeT>& x) const {
    assert(nrow() == x.size());
    for (size_type i=0; i<nrow(); i++) {
      for (size_type z=rowptr(i); z<rowptr(i+1); z++) {
        x(i + x.begin()) += m_value[z];
      }
    }
    return x;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> csr_matrix<ValueT,RangeT>::row_sum() const {
    vector<ValueT,RangeT> x(nrow());
    x = 0;
    return row_sum(x);
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::copyfrom(const dense_matrix<ValueT,RangeT>& m) {
    assert(nrow() == m.nrow());
    assert(ncol() == m.ncol());
    size_type z = 0;
    size_type zi = 0;
    m_rowptr[zi] = z + m_origin;
    for (RangeT i=m.rbegin(); i<=m.rend(); i++) {
      for (RangeT j=m.cbegin(); j<=m.cend(); j++) {
        if (m(i,j) != 0) {
          m_value[z] = m(i,j);
          m_colind[z] = (j-m.cbegin()) + m_origin;
          z++;
        }
      }
      zi++;
      m_rowptr[zi] = z + m_origin;
    }
    assert(nnz() == z);
    return *this;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& csr_matrix<ValueT,RangeT>::copyto(dense_matrix<ValueT,RangeT>& x) const {
    x = 0;
    for (size_type i=0; i<nrow(); i++) {
      for (size_type z=rowptr(i); z<rowptr(i+1); z++) {
        size_type j = colind(z);
        x(i + x.rbegin(), j + x.cbegin()) = m_value[z];
      }
    }
    return x;
  }

  // instance
  template class csr_matrix<double,int>;
}
