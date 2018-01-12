
#include <cmath>
#include "marlib.hpp"

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
  gaussinte<ValueT,RangeT>::gaussinte(size_type n, const ValueT& eps)
  : m_n(n), m_x(n), m_w(n) {
      comp_w(eps);
  }

  template <typename ValueT, typename RangeT>
  void gaussinte<ValueT,RangeT>::comp_w(const ValueT& eps) {
    constexpr ValueT pai = 3.14159265358979324;
    switch (m_n) {
      case 1:
        m_x(1) = 0.0;
        m_w(1) = 2.0;
        return;
      case 2:
        m_x(1) = -std::sqrt(1.0/3.0);
        m_w(1) = 1.0;
        m_x(2) = -m_x(2);
        m_w(2) = m_w(1);
        return;
      case 3:
        m_x(1) = -std::sqrt(0.6);
        m_w(1) = 5.0/9.0;
        m_x(2) = 0.0;
        m_w(2) = 8.0/9.0;
        m_x(3) = -m_x(1);
        m_w(3) = m_w(1);
        return;
    }
    RangeT m = m_n / 2;
    ValueT npai = pai / (m_n + 0.5);
    ValueT p0, p1, p2, q0, q1, q2;
    for (RangeT i=1; i<=m; i++) {
      ValueT tmp = cos((i - 0.25) * npai);
      ValueT dt = tmp;
      while (std::abs(dt) > std::abs(tmp) * eps) {
        p1 = tmp;
        p2 = (3.0 * tmp * tmp - 1.0) * 0.5;
        q1 = 1.0;
        q2 = 3.0 * tmp;
        for (size_type l=3; l<=m_n; l++) {
          p0 = p1;
          p1 = p2;
          p2 = ((l + l - 1) * tmp * p1 - (l-1) * p0) / l;
          q0 = q1;
          q1 = q2;
          q2 = ((l + l - 1) * (tmp * q1 + p1) - (l-1) * q0) / l;
        }
        dt = p2 / q2;
        tmp = tmp - dt;
      }
      m_x(i) = tmp;
      m_w(i) = 2.0 / (m_n * p1 * q2);
    }

    if (m_n % 2 == 1) {
      ValueT tmp = m_n;
      for (RangeT i=1; i<=m; i++) {
        tmp = tmp * (0.5 - i) / i;
      }
      m_x(m+1) = 0.0;
      m_w(m+1) = 2.0 / (tmp * tmp);
    }

    for (RangeT i=1; i<=m; i++) {
      m_x(m_n+1-i) = m_x(i);
      m_x(i) = -m_x(i); // reverse order
      m_w(m_n+1-i) = m_w(i);
    }
  }

  template <typename ValueT, typename RangeT>
  const vector<ValueT,RangeT>& gaussinte<ValueT,RangeT>::get_w() const {
      return m_w;
  }

  template <typename ValueT, typename RangeT>
  ValueT gaussinte<ValueT,RangeT>::comp_fx(const ValueT& a, const ValueT& b, vector<ValueT,RangeT>& fx) {
    ValueT t1 = (b - a)/2;
    ValueT t2 = (b + a)/2;
    for (RangeT i=m_x.begin(), ifx=fx.begin(); i<=m_x.end(); i++, ifx++) {
      fx(ifx) = t1 * m_x(i) + t2;
    }
    return t1;
  }

  template <typename ValueT, typename RangeT>
  ValueT gaussinte<ValueT,RangeT>::comp_value(const ValueT& c, const vector<ValueT,RangeT>& fv) {
    ValueT s = 0;
    for (RangeT iw=m_w.begin(), iv=fv.begin(); iw<=m_w.end(); iw++, iv++) {
      s += m_w(iw) * fv(iv);
    }
    s *= c;
    return s;
  }

  template class gaussinte<double,int>;

}
