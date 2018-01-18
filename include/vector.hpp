/*
  vector class
 */

#pragma once

#include <iostream>
#include <initializer_list>

#include "range.hpp"
#include "array.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT>
  class vector {
  public:
    using ValueType = ValueT;
    using RangeType = RangeT;

    vector(size_type size, const array<ValueT>& a, size_type inc);
    vector(const vector<ValueT,RangeT>& v);
    vector(const vector<ValueT,RangeT>& v, ValueT* p);
    vector(size_type size);
    vector(size_type size, ValueT* v, size_type inc);
    vector(std::initializer_list<ValueT> v);
    ~vector();

  private:
    vector(const range<RangeT>& r, const array<ValueT>& a, size_type inc);
    range<RangeT> m_range;
    array<ValueT> m_value;
    size_type m_inc;

  public:
    ValueT& operator()(const RangeT i);
    const ValueT& operator()(const RangeT i) const;

    const ValueT* ptr() const;
    ValueT* ptr();

    size_type inc() const;
    const RangeT begin() const;
    const RangeT end() const;
    size_type size() const;
    vector<ValueT,RangeT> operator()(const range<RangeT>& r);
    const vector<ValueT,RangeT> operator()(const range<RangeT>& r) const;
    const array<ValueT>& value() const;
    void set_range(const range<RangeT>& r);
    vector<ValueT,RangeT> clone() const;
    vector<ValueT,RangeT> clone(ValueT* p) const;

    // equal
    vector<ValueT,RangeT>& operator=(const ValueT& v);
    vector<ValueT,RangeT>& operator=(const vector<ValueT,RangeT>& v);
    vector<ValueT,RangeT>& operator=(const vector<ValueT*,RangeT>& v);

    // arithmetic operators
    vector<ValueT,RangeT>& operator+=(const vector<ValueT,RangeT>& v);
    vector<ValueT,RangeT>& operator-=(const vector<ValueT,RangeT>& v);
    vector<ValueT,RangeT>& operator*=(const vector<ValueT,RangeT>& v);
    vector<ValueT,RangeT>& operator/=(const vector<ValueT,RangeT>& v);

    vector<ValueT,RangeT>& operator+=(const ValueT& v);
    vector<ValueT,RangeT>& operator-=(const ValueT& v);
    vector<ValueT,RangeT>& operator*=(const ValueT& v);
    vector<ValueT,RangeT>& operator/=(const ValueT& v);

    vector<ValueT,RangeT> operator+(const vector<ValueT,RangeT>& v) const;
    vector<ValueT,RangeT> operator-(const vector<ValueT,RangeT>& v) const;
    vector<ValueT,RangeT> operator*(const vector<ValueT,RangeT>& v) const;
    vector<ValueT,RangeT> operator/(const vector<ValueT,RangeT>& v) const;

    ////// print
    std::ostream& print(std::ostream& os) const;

    template <typename ValueTT, typename RangeTT>
    friend std::ostream& operator<< (std::ostream& os, const vector<ValueTT,RangeTT>& v);

    vector<ValueT,RangeT>& copyfrom(const ValueT* x);
  };

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
  std::ostream& operator<<(std::ostream& os, const vector<ValueT,RangeT>& v) {
    return v.print(os);
  }

  //////////////////////////

  template <typename ValueT, typename RangeT>
  class vector<ValueT*,RangeT> {
  public:
    vector(size_type size, const array<ValueT*>& a, size_type inc);
    vector(const vector<ValueT*,RangeT>& v);
    vector(size_type size);
    ~vector();

  private:
    vector(const range<RangeT>& r, const array<ValueT*>& a, size_type inc);
    range<RangeT> m_range;
    array<ValueT*> m_value;
    size_type m_inc;

  public:
    ValueT*& ptr(const RangeT i);
    ValueT& operator()(const RangeT i);
    const ValueT& operator()(const RangeT i) const;

    const RangeT begin() const;
    const RangeT end() const;
    size_type size() const;

    void set_range(const range<RangeT>& r);

    // equal
    vector<ValueT*,RangeT>& operator=(const ValueT& v);
    vector<ValueT*,RangeT>& operator=(const vector<ValueT,RangeT>& v);
    vector<ValueT*,RangeT>& operator=(const vector<ValueT*,RangeT>& v);

    // arithmetic operators
    vector<ValueT*,RangeT>& operator+=(const vector<ValueT,RangeT>& v);
    vector<ValueT*,RangeT>& operator-=(const vector<ValueT,RangeT>& v);
    vector<ValueT*,RangeT>& operator*=(const vector<ValueT,RangeT>& v);
    vector<ValueT*,RangeT>& operator/=(const vector<ValueT,RangeT>& v);

    vector<ValueT*,RangeT>& operator+=(const ValueT& v);
    vector<ValueT*,RangeT>& operator-=(const ValueT& v);
    vector<ValueT*,RangeT>& operator*=(const ValueT& v);
    vector<ValueT*,RangeT>& operator/=(const ValueT& v);

    vector<ValueT,RangeT> operator+(const vector<ValueT,RangeT>& v) const;
    vector<ValueT,RangeT> operator-(const vector<ValueT,RangeT>& v) const;
    vector<ValueT,RangeT> operator*(const vector<ValueT,RangeT>& v) const;
    vector<ValueT,RangeT> operator/(const vector<ValueT,RangeT>& v) const;

    ////// print
    std::ostream& print(std::ostream& os) const;

    template <typename ValueTT, typename RangeTT>
    friend std::ostream& operator<< (std::ostream& os, const vector<ValueTT,RangeTT>& v);
  };

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
}
