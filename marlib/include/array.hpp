/*
  array class
*/

namespace marlib {

  template <typename T>
  class array_base {
  public:
    inline array_base(size_type size)
    : m_delete_ptr(true), m_ptr(new T[size]) {}

    inline array_base(const array_base<T>& a)
    : m_delete_ptr(a.m_delete_ptr), m_ptr(a.m_ptr) { }

    inline array_base(T* v)
    : m_delete_ptr(false), m_ptr(v) { }

    inline ~array_base() {
      // std::cout << "destructor" << std::endl;
      if (m_delete_ptr) {
        delete [] m_ptr;
      }
    }

  private:
    bool m_delete_ptr;
    T* m_ptr;

  public:
    inline T& operator[](size_type i) {
      return m_ptr[i];
    }

    inline const T& operator[](size_type i) const {
      return m_ptr[i];
    }
  };

  template <typename T>
  class array {
  public:
    inline array(size_type size) : m_ptr(new array_base<T>(size)), m_size(size) {
      m_base = &(*m_ptr)[0];
    }

    inline array(const array<T>& a) : m_ptr(a.m_ptr), m_size(a.m_size), m_base(a.m_base) { }

    inline array(size_type size, T* v) : m_ptr(new array_base<T>(v)), m_size(size), m_base(v) { }

    inline ~array() {}

  private:
    std::shared_ptr< array_base<T> > m_ptr;
    size_type m_size;
    T* m_base;

    // create sub-array
    inline array(size_type size, const std::shared_ptr<array_base<T>>& ptr, T* base) : m_ptr(ptr), m_size(size), m_base(base) {}

  public:
    inline T& operator[](size_type i) {
      assert(i >= 0 && i < m_size);
      return m_base[i];
    }

    inline const T& operator[](size_type i) const {
      assert(i >= 0 && i < m_size);
      return m_base[i];
    }

    inline size_type size() const {
      return m_size;
    }

    inline array<T> subarray(size_type i) {
      return array(m_size - i, m_ptr, m_base + i);
    }

    inline const array<T> subarray(size_type i) const {
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
    inline array(size_type size)
    : m_ptr(new array_base<T*>(size)), m_size(size), m_base(&((*m_ptr)[0])) { }

    inline array(const array<T>& a)
    : m_ptr(a.m_ptr), m_size(a.m_size), m_base(a.m_base) { }

    inline array(size_type size, T** v)
    : m_ptr(new array_base<T*>(v)), m_size(size), m_base(v) { }

    inline ~array() {}

  private:
    std::shared_ptr<array_base<T*>> m_ptr;
    size_type m_size;
    T** m_base;

    // create sub-array
    inline array(size_type size, const std::shared_ptr<array_base<T*>>& ptr, T** base) : m_ptr(ptr), m_size(size), m_base(base) {}

  public:
    inline T*& ptr(size_type i) {
      assert(i >= 0 && i < m_size);
      return m_base[i];
    }

    inline T& operator[](size_type i) {
      assert(i >= 0 && i < m_size);
      return *m_base[i];
    }

    inline const T& operator[](size_type i) const {
      assert(i >= 0 && i < m_size);
      return *m_base[i];
    }

    inline size_type size() const {
      return m_size;
    }

    inline array<T*> subarray(size_type i) {
      return array(m_size - i, m_ptr, m_base + i);
    }

    inline const array<T*> subarray(size_type i) const {
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