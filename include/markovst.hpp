#pragma once

namespace marlib {

  template <typename VectorT, typename MatrixT>
  void gth_impl(const MatrixT& Q, VectorT x);

  template <typename VectorT, typename MatrixT>
  VectorT& gth(const MatrixT& Q, VectorT& x) {
    gth_impl(Q, x);
    return x;
  }

  template <typename ValueT, typename VectorT>
  struct gsconf {
    size_type iter;
    size_type maxiter;
    size_type steps;
    ValueT atol;
    ValueT rtol;
    ValueT aerror;
    ValueT rerror;
    int info;

    gsconf() :
    iter(0),
    maxiter(1000),
    steps(5),
    atol(0),
    rtol(1.0e-8),
    aerror(0),
    rerror(0),
    info(0) { }

    void callback(const VectorT& x) {
      std::cout << "iter=" << iter << " rerror=" << rerror;
      std::cout << " " << x << std::endl;
    }
  };

  template <typename ValueT, typename RangeT, typename MatrixT, typename ConfigT>
  vector<ValueT,RangeT>& stgs(const MatrixT& Q, vector<ValueT,RangeT>& x, ConfigT& conf);

}
