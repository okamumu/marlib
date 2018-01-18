/*
  Gamma Functions
*/

#include <cmath>

namespace marlib {

  template <typename ValueT, typename IntegerT>
  struct gamma {
    static constexpr IntegerT FACTMAX = 20;
    static constexpr ValueT PI = 3.14159265358979324; // pi
    static constexpr ValueT LOG_2PI = 1.83787706640934548;
    static constexpr ValueT LOG_PI = 1.14472988584940017; // log(pi)

    static constexpr IntegerT N = 8;
    static constexpr ValueT B0 = 1;
    static constexpr ValueT B1 = (-1.0 / 2.0);
    static constexpr ValueT B2 = ( 1.0 / 6.0);
    static constexpr ValueT B4 = (-1.0 / 30.0);
    static constexpr ValueT B6 = ( 1.0 / 42.0);
    static constexpr ValueT B8 = (-1.0 / 30.0);
    static constexpr ValueT B10 = ( 5.0 / 66.0);
    static constexpr ValueT B12 = (-691.0 / 2730.0);
    static constexpr ValueT B14 = ( 7.0 / 6.0);
    static constexpr ValueT B16 = (-3617.0 / 510.0);

    static const ValueT nfact[];
    static const ValueT lognfact[];

  };

  template <typename ValueT, typename IntegerT>
  const ValueT gamma<ValueT,IntegerT>::nfact[] = {
    1.0,                        // 0
    1.0,                        // 1
    2.0,                        // 2
    6.0,                        // 3
    24.0,                       // 4
    120.0,                      // 5
    720.0,                      // 6
    5040.0,                     // 7
    40320.0,                    // 8
    362880.0,                   // 9
    3628800.0,                  // 10
    39916800.0,                 // 11
    479001600.0,                // 12
    6227020800.0,               // 13
    87178291200.0,              // 14
    1307674368000.0,            // 15
    20922789888000.0,           // 16
    355687428096000.0,          // 17
    6402373705728000.0,         // 18
    121645100408832000.0,       // 19
    2432902008176640000.0       // 20
  };

  template <typename ValueT, typename IntegerT>
  const ValueT gamma<ValueT,IntegerT>::lognfact[] = {
    0.0,
    0.0,
    0.6931471805599453,
    1.791759469228055,
    3.1780538303479458,
    4.787491742782046,
    6.579251212010101,
    8.525161361065415,
    10.60460290274525,
    12.801827480081469,
    15.104412573075516,
    17.502307845873887,
    19.987214495661885,
    22.552163853123425,
    25.19122118273868,
    27.89927138384089,
    30.671860106080672,
    33.50507345013689,
    36.39544520803305,
    39.339884187199495,
    42.335616460753485
  };

  template <typename ValueT, typename IntegerT>
  ValueT lgamma(const ValueT xinit) {
    ValueT x = xinit;
    ValueT v = 1;
    while (x < gamma<ValueT,IntegerT>::N) { v *=x; x++; }
    ValueT w = 1 / (x * x);
    return ((((((((gamma<ValueT,IntegerT>::B16 / (16 * 15)) * w + (gamma<ValueT,IntegerT>::B14 / (14 * 13))) * w
      + (gamma<ValueT,IntegerT>::B12 / (12 * 11))) * w + (gamma<ValueT,IntegerT>::B10 / (10 * 9))) * w
    + (gamma<ValueT,IntegerT>::B8 / (8 * 7))) * w + (gamma<ValueT,IntegerT>::B6 / (6 * 5))) * w
    + (gamma<ValueT,IntegerT>::B4 / (4 * 3))) * w + (gamma<ValueT,IntegerT>::B2 / (2 * 1))) / x
    + 0.5 * gamma<ValueT,IntegerT>::LOG_2PI - log(v) - x + (x - 0.5) * log(x);
  }

  template <typename ValueT, typename IntegerT>
  ValueT tgamma(const ValueT x) {
    if (x < 0) {
      return gamma<ValueT,IntegerT>::PI / (sin(gamma<ValueT,IntegerT>::PI * x) * exp(lgamma<ValueT,IntegerT>(1-x)));
    }
    return exp(lgamma<ValueT,IntegerT>(x));
  }

  template <typename ValueT, typename IntegerT>
  ValueT psi(const ValueT xinit) {
    ValueT x = xinit;
    ValueT v = 0;
    while (x < gamma<ValueT,IntegerT>::N) { v += 1 / x; x++; }
    ValueT w = 1 / (x * x);
    v += ((((((((gamma<ValueT,IntegerT>::B16 / 16) * w + (gamma<ValueT,IntegerT>::B14 /14)) * w
      + (gamma<ValueT,IntegerT>::B12 / 12)) * w + (gamma<ValueT,IntegerT>::B10 / 10)) * w
    + (gamma<ValueT,IntegerT>::B8 / 8)) * w + (gamma<ValueT,IntegerT>::B6 / 6)) * w
    + (gamma<ValueT,IntegerT>::B4 / 4)) * w + (gamma<ValueT,IntegerT>::B2 / 2)) * w + 0.5 / x;
    return log(x) - v;
  }

  template <typename ValueT, typename IntegerT>
  ValueT polygamma(const IntegerT n, const ValueT xinit) {
    ValueT x = xinit;
    ValueT u = 1;
    for(IntegerT k=1-n; k<0; k++) u *= k;
    ValueT v = 0;
    while (x < gamma<ValueT,IntegerT>::N) { v +=1 / pow(x, n+1); x++; }
    ValueT w = x * x;
    ValueT t = (((((((gamma<ValueT,IntegerT>::B16
      * (n + 15.0) * (n + 14) / (16 * 15 * w) + gamma<ValueT,IntegerT>::B14)
    * (n + 13.0) * (n + 12) / (14 * 13 * w) + gamma<ValueT,IntegerT>::B12)
    * (n + 11.0) * (n + 10) / (12 * 11 * w) + gamma<ValueT,IntegerT>::B10)
    * (n + 9.0) * (n + 8) / (10 * 9 * w) + gamma<ValueT,IntegerT>::B8)
    * (n + 7.0) * (n + 6) / (8 * 7 * w) + gamma<ValueT,IntegerT>::B6)
    * (n + 5.0) * (n + 4) / (6 * 5 * w) + gamma<ValueT,IntegerT>::B4)
    * (n + 3.0) * (n + 2) / (4 * 3 * w) + gamma<ValueT,IntegerT>::B2)
    * (n + 1.0) * n / (2 * 1 * w)
    + 0.5 * n / x + 1;
    return u * (t / pow(x, n) + n * v);
  }

  template <typename ValueT, typename IntegerT>
  ValueT tfact(const IntegerT s) {
    if (s <= gamma<ValueT,IntegerT>::FACTMAX) {
      return gamma<ValueT,IntegerT>::nfact[s];
    } else {
      return exp(lgamma<ValueT,IntegerT>(1.0 + s));
    }
  }

  template <typename ValueT, typename IntegerT>
  ValueT lfact(const IntegerT s) {
    if (s <= gamma<ValueT,IntegerT>::FACTMAX) {
      return gamma<ValueT,IntegerT>::lognfact[s];
    } else {
      return lgamma<ValueT,IntegerT>(1.0 + s);
    }
  }

  template struct gamma<double,int>;
  template double lgamma<double,int>(const double x);
  template double tgamma<double,int>(const double x);
  template double psi<double,int>(const double x);
  template double polygamma<double,int>(const int n, const double x);
  template double tfact<double,int>(const int s);
  template double lfact<double,int>(const int s);

}
