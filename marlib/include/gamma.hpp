/*
  Gamma Functions
*/

namespace marlib {

  template <typename ValueT, typename IntegerT>
  ValueT lgamma(const ValueT x);

  template <typename ValueT, typename IntegerT>
  ValueT tgamma(const ValueT x);

  template <typename ValueT, typename IntegerT>
  ValueT psi(const ValueT x);

  template <typename ValueT, typename IntegerT>
  ValueT polygamma(const IntegerT n, const ValueT x);

  template <typename ValueT, typename IntegerT>
  ValueT tfact(const IntegerT s);

  template <typename ValueT, typename IntegerT>
  ValueT lfact(const IntegerT s);

}
