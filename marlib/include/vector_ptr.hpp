/*
  vector class
 */

#include "marlib.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>::vector(size_type size, const array<ValueT*>& a, size_type inc)
  : m_range(size), m_value(a), m_inc(inc) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>::vector(const range<RangeT>& r, const array<ValueT*>& a, size_type inc)
  : m_range(r), m_value(a), m_inc(inc) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>::vector(const vector<ValueT*,RangeT>& v)
  : m_range(v.m_range), m_value(v.m_value), m_inc(v.m_inc) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>::vector(size_type size)
  : m_range(size), m_value(size), m_inc(1) {}

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>::~vector() { }

  template <typename ValueT, typename RangeT>
  inline ValueT*& vector<ValueT*,RangeT>::ptr(const RangeT i) {
    return m_value.ptr((i - m_range.begin()) * m_inc);
  }

  template <typename ValueT, typename RangeT>
  inline ValueT& vector<ValueT*,RangeT>::operator()(const RangeT i) {
    return m_value[(i - m_range.begin()) * m_inc];
  }

  template <typename ValueT, typename RangeT>
  inline const ValueT& vector<ValueT*,RangeT>::operator()(const RangeT i) const {
    return m_value[(i - m_range.begin()) * m_inc];
  }

  template <typename ValueT, typename RangeT>
  inline const RangeT vector<ValueT*,RangeT>::begin() const {
    return m_range.begin();
  }

  template <typename ValueT, typename RangeT>
  inline const RangeT vector<ValueT*,RangeT>::end() const {
    return m_range.end();
  }

  template <typename ValueT, typename RangeT>
  inline size_type vector<ValueT*,RangeT>::size() const {
    return m_range.size();
  }

  template <typename ValueT, typename RangeT>
  inline void vector<ValueT*,RangeT>::set_range(const range<RangeT>& r) {
    assert(r.size() == m_range.size());
    m_range = r;
  }

  // equal
  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) = v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator=(const vector<ValueT,RangeT>& x) {
    assert(size() == x.size());
    for (RangeT i=begin(), j=x.begin(); i<=end(); i++, j++) {
      operator()(i) = x(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator=(const vector<ValueT*,RangeT>& x) {
    assert(size() == x.size());
    for (RangeT i=begin(), j=x.begin(); i<=end(); i++, j++) {
      operator()(i) = x(j);
    }
    return *this;
  }

  // arithmetic operators

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator+=(const vector<ValueT,RangeT>& v) {
    assert(size() == v.size());
    for (RangeT i=begin(), j=v.begin(); i<=end(); i++, j++) {
      operator()(i) += v(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator-=(const vector<ValueT,RangeT>& v) {
    assert(size() == v.size());
    for (RangeT i=begin(), j=v.begin(); i<=end(); i++, j++) {
      operator()(i) -= v(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator*=(const vector<ValueT,RangeT>& v) {
    assert(size() == v.size());
    for (RangeT i=begin(), j=v.begin(); i<=end(); i++, j++) {
      operator()(i) *= v(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator/=(const vector<ValueT,RangeT>& v) {
    assert(size() == v.size());
    for (RangeT i=begin(), j=v.begin(); i<=end(); i++, j++) {
      operator()(i) /= v(j);
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator+=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) += v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator-=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) -= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator*=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) *= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT*,RangeT>& vector<ValueT*,RangeT>::operator/=(const ValueT& v) {
    for (RangeT i=begin(); i<=end(); i++) {
      operator()(i) /= v;
    }
    return *this;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT*,RangeT>::operator+(const vector<ValueT,RangeT>& v) const {
    assert(size() == v.size());
    vector<ValueT,RangeT> result(size());
    result = *this;
    result += v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT*,RangeT>::operator-(const vector<ValueT,RangeT>& v) const {
    assert(size() == v.size());
    vector<ValueT,RangeT> result(size());
    result = *this;
    result -= v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT*,RangeT>::operator*(const vector<ValueT,RangeT>& v) const {
    assert(size() == v.size());
    vector<ValueT,RangeT> result(size());
    result = *this;
    result *= v;
    return result;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT> vector<ValueT*,RangeT>::operator/(const vector<ValueT,RangeT>& v) const {
    assert(size() == v.size());
    vector<ValueT,RangeT> result(size());
    result = *this;
    result /= v;
    return result;
  }

  ////// print

  template <typename ValueT, typename RangeT>
  std::ostream& vector<ValueT*,RangeT>::print(std::ostream& os) const {
    for (RangeT i=begin(); i<=end(); i++) {
      os << operator()(i) << " ";
    }
    return os;
  }

  template class vector<int*,int>;
  template class vector<size_type*,int>;
  template class vector<double*,int>;
}