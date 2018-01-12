/*
  array class
*/

#pragma once

#include <iostream>
#include <memory>
#include <cassert>

#include "types.h"

namespace marlib {

  template <typename T>
  class array_base {
  public:
    array_base(size_type size)
    : m_delete_ptr(true), m_ptr(new T[size]) {}

    array_base(const array_base<T>& a)
    : m_delete_ptr(a.m_delete_ptr), m_ptr(a.m_ptr) { }

    array_base(T* v)
    : m_delete_ptr(false), m_ptr(v) { }

    ~array_base() {
      // std::cout << "destructor" << std::endl;
      if (m_delete_ptr) {
        delete [] m_ptr;
      }
    }

  private:
    bool m_delete_ptr;
    T* m_ptr;

  public:
    T& operator[](size_type i) {
      return m_ptr[i];
    }

    const T& operator[](size_type i) const {
      return m_ptr[i];
    }
  };

  template <typename T>
  class array {
  public:
    array(size_type size) : m_ptr(new array_base<T>(size)), m_size(size) {
      m_base = &(*m_ptr)[0];
    }

    array(const array<T>& a) : m_ptr(a.m_ptr), m_size(a.m_size), m_base(a.m_base) { }

    array(size_type size, T* v) : m_ptr(new array_base<T>(v)), m_size(size), m_base(v) { }

    ~array() {}

  private:
    std::shared_ptr< array_base<T> > m_ptr;
    size_type m_size;
    T* m_base;

    // create sub-array
    array(size_type size, const std::shared_ptr<array_base<T>>& ptr, T* base) : m_ptr(ptr), m_size(size), m_base(base) {}

  public:
    T& operator[](size_type i) {
      assert(i >= 0 && i < m_size);
      return m_base[i];
    }

    const T& operator[](size_type i) const {
      assert(i >= 0 && i < m_size);
      return m_base[i];
    }

    size_type size() const {
      return m_size;
    }

    array<T> subarray(size_type i) {
      return array(m_size - i, m_ptr, m_base + i);
    }

    const array<T> subarray(size_type i) const {
      return array(m_size - i, m_ptr, m_base + i);
    }

    // print
    std::ostream& print(std::ostream& os) const {
      for (size_type i=0; i<m_size; i++) {
        os << m_base[i] << " ";
      }
      return os;
    }

    template <typename TT>
    friend std::ostream& operator<< (std::ostream& os, const array<TT>& v);

  };

  template <typename T>
  std::ostream& operator<<(std::ostream& os, const array<T>& v) {
    return v.print(os);
  }

  template <typename T>
  class array<T*> {
  public:
    array(size_type size)
    : m_ptr(new array_base<T*>(size)), m_size(size), m_base(&((*m_ptr)[0])) { }

    array(const array<T>& a)
    : m_ptr(a.m_ptr), m_size(a.m_size), m_base(a.m_base) { }

    array(size_type size, T** v)
    : m_ptr(new array_base<T*>(v)), m_size(size), m_base(v) { }

    ~array() {}

  private:
    std::shared_ptr<array_base<T*>> m_ptr;
    size_type m_size;
    T** m_base;

    // create sub-array
    array(size_type size, const std::shared_ptr<array_base<T*>>& ptr, T** base) : m_ptr(ptr), m_size(size), m_base(base) {}

  public:
    T*& ptr(size_type i) {
      assert(i >= 0 && i < m_size);
      return m_base[i];
    }

    T& operator[](size_type i) {
      assert(i >= 0 && i < m_size);
      return *m_base[i];
    }

    const T& operator[](size_type i) const {
      assert(i >= 0 && i < m_size);
      return *m_base[i];
    }

    size_type size() const {
      return m_size;
    }

    array<T*> subarray(size_type i) {
      return array(m_size - i, m_ptr, m_base + i);
    }

    const array<T*> subarray(size_type i) const {
      return array(m_size - i, m_ptr, m_base + i);
    }

    // print
    std::ostream& print(std::ostream& os) const {
      for (size_type i=0; i<m_size; i++) {
        os << *m_base[i] << std::endl;
      }
      return os;
    }

    template <typename TT>
    friend std::ostream& operator<< (std::ostream& os, const array<TT>& v);

  };
}
