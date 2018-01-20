/*
  vector class
 */

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
  std::ostream& operator<<(std::ostream& os, const vector<ValueT,RangeT>& v) {
    return v.print(os);
  }
}
