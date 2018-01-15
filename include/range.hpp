////////
///  range
///////

#pragma once

#include "types.h"

namespace marlib {

  template <typename T>
  class range {
  public:
    static const T default_origin = 1;

    range(const T& size) : m_begin(default_origin), m_end(size-1+default_origin) {}
    range(const T& begin, const T& end) : m_begin(begin), m_end(end) {}
    range(const range<T>& r) : m_begin(r.m_begin), m_end(r.m_end) {}
    ~range() {}

  private:
    T m_begin;
    T m_end;

  public:
    size_type size() const {
      return static_cast<size_type>(m_end - m_begin + 1);
    }

    const T& begin() const {
      return m_begin;
    }

    const T& end() const {
      return m_end;
    }

    range<T>& operator=(const range<T>& v) {
      m_begin = v.m_begin;
      m_end = v.m_end;
      return *this;
    }

    bool operator==(const range<T>& v) const {
      return m_begin == v.m_begin && m_end == v.m_end;
    }
  };

}
