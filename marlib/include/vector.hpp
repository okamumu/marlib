/*
  vector class
 */

namespace marlib {

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>::vector(size_type size, const array<ValueT>& a, size_type inc)
  : m_range(size), m_value(a), m_inc(inc) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>::vector(const range<RangeT>& r, const array<ValueT>& a, size_type inc)
  : m_range(r), m_value(a), m_inc(inc) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>::vector(const vector<ValueT,RangeT>& v)
  : m_range(v.m_range), m_value(v.m_value), m_inc(v.m_inc) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>::vector(const vector<ValueT,RangeT>& v, ValueT* p)
  : m_range(v.m_range), m_value(v.m_value.size(),p), m_inc(v.m_inc) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>::vector(size_type size)
  : m_range(size), m_value(size), m_inc(1) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>::vector(size_type size, ValueT* v, size_type inc)
  : m_range(size), m_value(size*inc,v), m_inc(inc) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>::vector(std::initializer_list<ValueT> v)
  : m_range(v.size()), m_value(v.size()), m_inc(1) {
    copyfrom(v.begin());
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>::~vector() { }

  template <typename ValueT, typename RangeT>
  inline ValueT& vector<ValueT,RangeT>::operator()(const RangeT i) {
    return m_value[(i - m_range.begin()) * m_inc];
  }

  template <typename ValueT, typename RangeT>
  inline const ValueT& vector<ValueT,RangeT>::operator()(const RangeT i) const {
    return m_value[(i - m_range.begin()) * m_inc];
  }

  template <typename ValueT, typename RangeT>
  inline ValueT* vector<ValueT,RangeT>::ptr() {
    return &operator()(m_range.begin());
  }

  template <typename ValueT, typename RangeT>
  inline const ValueT* vector<ValueT,RangeT>::ptr() const {
    return &operator()(m_range.begin());
  }

  template <typename ValueT, typename RangeT>
  inline size_type vector<ValueT,RangeT>::inc() const {
    return m_inc;
  }

  template <typename ValueT, typename RangeT>
  inline const RangeT vector<ValueT,RangeT>::begin() const {
    return m_range.begin();
  }

  template <typename ValueT, typename RangeT>
  inline const RangeT vector<ValueT,RangeT>::end() const {
    return m_range.end();
  }

  template <typename ValueT, typename RangeT>
  inline size_type vector<ValueT,RangeT>::size() const {
    return m_range.size();
  }

  template <typename ValueT, typename RangeT>
  inline const array<ValueT>& vector<ValueT,RangeT>::value() const {
    return m_value;
  }

  template <typename ValueT, typename RangeT>
  inline void vector<ValueT,RangeT>::set_range(const range<RangeT>& r) {
    assert(r.size() == m_range.size());
    m_range = r;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT,RangeT>::operator()(const range<RangeT>& r) {
    array<ValueT> tmp = m_value.subarray((r.begin() - m_range.begin()) * m_inc);
    return vector<ValueT,RangeT>(range<RangeT>(r.size()), tmp, m_inc);
  }

  template <typename ValueT, typename RangeT>
  const vector<ValueT,RangeT> vector<ValueT,RangeT>::operator()(const range<RangeT>& r) const {
    array<ValueT> tmp = m_value.subarray((r.begin() - m_range.begin()) * m_inc);
    return vector<ValueT,RangeT>(range<RangeT>(r.size()), tmp, m_inc);
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT,RangeT>::clone() const {
    vector<ValueT,RangeT> tmp(m_range, array<ValueT>(m_value.size()), m_inc);
    tmp = *this;
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT,RangeT>::clone(ValueT* p) const {
    vector<ValueT,RangeT> tmp(m_range, array<ValueT>(m_value.size(), p), m_inc);
    tmp = *this;
    return tmp;
  }

  // equal
  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator=(const ValueT& v) {
    if (inc() == 1) {
      std::fill_n(ptr(), size(), v);
    } else {
      for (RangeT i=begin(); i<=end(); i++) {
        operator()(i) = v;
      }
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator=(const vector<ValueT,RangeT>& x) {
    assert(size() == x.size());
    for (RangeT i=begin(), j=x.begin(); i<=end(); i++, j++) {
      operator()(i) = x(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator=(const vector<ValueT*,RangeT>& x) {
    assert(size() == x.size());
    for (RangeT i=begin(), j=x.begin(); i<=end(); i++, j++) {
      operator()(i) = x(j);
    }
    return *this;
  }

  // arithmetic operators

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator+=(const vector<ValueT,RangeT>& v) {
    assert(size() == v.size());
    for (RangeT i=begin(), j=v.begin(); i<=end(); i++, j++) {
      operator()(i) += v(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator-=(const vector<ValueT,RangeT>& v) {
    assert(size() == v.size());
    for (RangeT i=begin(), j=v.begin(); i<=end(); i++, j++) {
      operator()(i) -= v(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator*=(const vector<ValueT,RangeT>& v) {
    assert(size() == v.size());
    for (RangeT i=begin(), j=v.begin(); i<=end(); i++, j++) {
      operator()(i) *= v(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator/=(const vector<ValueT,RangeT>& v) {
    assert(size() == v.size());
    for (RangeT i=begin(), j=v.begin(); i<=end(); i++, j++) {
      operator()(i) /= v(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator+=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) += v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator-=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) -= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator*=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) *= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::operator/=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) /= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT,RangeT>::operator+(const vector<ValueT,RangeT>& v) const {
    assert(size() == v.size());
    vector<ValueT,RangeT> result = clone();
    result += v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT,RangeT>::operator-(const vector<ValueT,RangeT>& v) const {
    assert(size() == v.size());
    vector<ValueT,RangeT> result = clone();
    result -= v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT,RangeT>::operator*(const vector<ValueT,RangeT>& v) const {
    assert(size() == v.size());
    vector<ValueT,RangeT> result = clone();
    result *= v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT,RangeT>::operator/(const vector<ValueT,RangeT>& v) const {
    assert(size() == v.size());
    vector<ValueT,RangeT> result = clone();
    result /= v;
    return result;
  }

  ////// print

  template <typename ValueT, typename RangeT>
  std::ostream& vector<ValueT,RangeT>::print(std::ostream& os) const {
    for (RangeT i=begin(); i<=end(); i++) {
      os << operator()(i) << " ";
    }
    return os;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& vector<ValueT,RangeT>::copyfrom(const ValueT* x) {
    for (RangeT i=begin(); i<=end(); i++, x++) {
      operator()(i) = *x;
    }
    return *this;
  }

  // specialized

#ifdef F77BLAS

  template <>
  inline
  vector<double,int>& vector<double,int>::operator+=(const vector<double,int>& v) {
    assert(size() == v.size());
    dblas::daxpy(size(), 1, v.ptr(), v.inc(), ptr(), inc());
    return *this;
  }

  template <>
  inline
  vector<double,int>& vector<double,int>::operator-=(const vector<double,int>& v) {
    assert(size() == v.size());
    dblas::daxpy(size(), -1, v.ptr(), v.inc(), ptr(), inc());
    return *this;
  }

  template <>
  inline
  vector<double,int>& vector<double,int>::operator*=(const double& v) {
    dblas::dscal(size(), v, ptr(), inc());
    return *this;
  }

  template <>
  inline
  vector<double,int>& vector<double,int>::operator/=(const double& v) {
    dblas::dscal(size(), 1/v, ptr(), inc());
    return *this;
  }

  template <>
  inline
  vector<double,int>& vector<double,int>::operator=(const vector<double,int>& x) {
    assert(size() == x.size());
    dblas::dcopy(size(), x.ptr(), x.inc(), ptr(), inc());
    return *this;
  }

  // template <>
  // vector<double,int>& vector<double,int>::operator=(const double* x) {
  //   dblas::dcopy(size(), x, 1, ptr(), inc());
  //   return *this;
  // }

#endif

}
