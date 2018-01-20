
namespace marlib {

/*
!   Description: Gauss quadrature for the following integral

!      | b
!      |  f(x) dx
!      | a

!   gaussinte_w: make points and weights for n discrete points
!     n (in): the number of points. This is the size of both x and w.
!     x (out): x-axis points in the interval [-1, 1].
!     w (out): weights for the points.
!     eps (in): tolerance error.

!   gaussinte_fx: make x points for the interval [a,b]
!     n (in): the number of points.
!     x (in): the x points for interval [-1, 1].
!     a, b (in): lower and upper points for the integral
!     fx (out): x point for the interval [a, b]
!     return value: (b-a)/2

!   gaussinte_fv: compute the integral
!     n (in): the number of points.
!     w (in): weights for the x points in [-1,1]
!     c (in): (b-a)/2 ?
!     fv (in): function values at x points derived by gauss_inte_fx
!     return value: the interal value
*/

  template <typename ValueT, typename RangeT>
  class gaussinte {
  public:
    gaussinte(size_type n, const ValueT eps = 1.0e-8);

    ValueT comp_fx(const ValueT a, const ValueT b, vector<ValueT,RangeT>& fx);
    ValueT comp_value(const ValueT c, const vector<ValueT,RangeT>& fv);
    const vector<ValueT,RangeT>& get_w() const;

    template <typename Func>
    ValueT comp(const ValueT a, const ValueT b, Func&& f) {
      vector<ValueT,RangeT> fx(m_n);
      ValueT c = comp_fx(a, b, fx);
      for (RangeT i=fx.begin(); i<=fx.end(); i++) {
        fx(i) = f(fx(i));
      }
      return comp_value(c, fx);
    }

  private:
    size_type m_n;
    vector<ValueT,RangeT> m_x;
    vector<ValueT,RangeT> m_w;

    void comp_w(const ValueT eps);
  };

}
