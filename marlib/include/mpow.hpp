
namespace marlib {

  template <typename MatrixT, typename IntegerT>
  dense_matrix<typename MatrixT::ValueType, typename MatrixT::RangeType>& mpow(const trans_t& trans,
    const MatrixT& MA, IntegerT m, dense_matrix<typename MatrixT::ValueType, typename MatrixT::RangeType>& ME);

}
