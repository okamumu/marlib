/*
  vector class
 */

namespace marlib {

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

}
