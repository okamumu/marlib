# marlib

C++ library for Markov chains

This is a library for computing Markov chaings. The package consists of header files on templates. The package requires f77blas.

## How to use

Include the header 
```c
#include <marlib.h>
```

Write the matrix traits for your own environment. For example, the traits for Rcpp are
```
#ifndef MARLIB_RCPP_H
#define MARLIB_RCPP_H

#include <Rcpp.h>
#include <marlib.h>

namespace marlib {

template <>
struct vector_traits<Rcpp::NumericVector> {
  static int size(const Rcpp::NumericVector& v) { return v.size(); }
  static const double* value(const Rcpp::NumericVector& v) { return &v[0]; }
  static double* value(Rcpp::NumericVector& v) { return &v[0]; }
  static int inc(const Rcpp::NumericVector& v) { return 1; }
};

template <>
struct vector_traits<Rcpp::NumericMatrix> {
  static int size(const Rcpp::NumericMatrix& v) { return v.size(); }
  static const double* value(const Rcpp::NumericMatrix& v) { return &v[0]; }
  static double* value(Rcpp::NumericMatrix& v) { return &v[0]; }
  static int inc(const Rcpp::NumericMatrix& v) { return 1; }
};

template <>
struct dense_matrix_traits<Rcpp::NumericMatrix> {
  static int nrow(const Rcpp::NumericMatrix& m) { return m.nrow(); }
  static int ncol(const Rcpp::NumericMatrix& m) { return m.ncol(); }
  static const double* value(const Rcpp::NumericMatrix& m) { return &m[0]; }
  static double* value(Rcpp::NumericMatrix& m) { return &m[0]; }
  static int ld(const Rcpp::NumericMatrix& m) { return nrow(m); }
};

template <>
struct vector_traits<Rcpp::S4> {
  static int size(const Rcpp::S4& v) { return Rcpp::as<Rcpp::NumericVector>(v.slot("x")).length(); }
  static const double* value(const Rcpp::S4& v) { return &Rcpp::as<Rcpp::NumericVector>(v.slot("x"))[0]; }
  static double* value(Rcpp::S4& v) { return &Rcpp::as<Rcpp::NumericVector>(v.slot("x"))[0]; }
  static int inc(const Rcpp::S4& v) { return 1; }
};

template <>
struct dense_matrix_traits<Rcpp::S4> {
  static int nrow(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
  static int ncol(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[1]; }
  static const double* value(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static double* value(Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static int ld(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
};

template <>
struct csr_matrix_traits<Rcpp::S4> {
  static int nrow(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
  static int ncol(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[1]; }
  static int nnz(const Rcpp::S4& m) { return Rcpp::as<Rcpp::NumericVector>(m.slot("x")).length(); }
  static int base(const Rcpp::S4& m) { return 0; }
  static double* value(Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const double* value(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const int* rowptr(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("p"))[0]; }
  static const int* colind(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("j"))[0]; }
};

template <>
struct csc_matrix_traits<Rcpp::S4> {
  static int nrow(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
  static int ncol(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[1]; }
  static int nnz(const Rcpp::S4& m) { return Rcpp::as<Rcpp::NumericVector>(m.slot("x")).length(); }
  static int base(const Rcpp::S4& m) { return 0; }
  static double* value(Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const double* value(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const int* colptr(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("p"))[0]; }
  static const int* rowind(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("i"))[0]; }
};

template <>
struct coo_matrix_traits<Rcpp::S4> {
  static int nrow(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[0]; }
  static int ncol(const Rcpp::S4& m) { return Rcpp::as<Rcpp::IntegerVector>(m.slot("Dim"))[1]; }
  static int nnz(const Rcpp::S4& m) { return Rcpp::as<Rcpp::NumericVector>(m.slot("x")).length(); }
  static int base(const Rcpp::S4& m) { return 0; }
  static double* value(Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const double* value(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::NumericVector>(m.slot("x"))[0]; }
  static const int* rowind(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("i"))[0]; }
  static const int* colind(const Rcpp::S4& m) { return &Rcpp::as<Rcpp::IntegerVector>(m.slot("j"))[0]; }
};

}

#endif
```

## Functions

```c
void ctmc_st_gth(const T1& Q, T2& x, MatT);
```

This function computes the steady-state vector of CTMC with GTH algorithm.

- Q (in): An infinitesimal generator.
- x (out): The steady-state vector.
- MatT: Type for matrix. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.

```c
void ctmc_st_power(T1& Q, T2& x, Params& params, const Func& callback, MatT);
```

This function computes the steady-state vector of CTMC with power method.

- Q (inout): An infinitesimal generator. This is changed to the probability matrix after computation.
- x (out): The steady-state vector.
- params (inout): The parameters for computation. unif, rtol, steps and maxiter are required. info, iter and rerror are stoted after the computation.
- callback: A callback function. The argument is Params. This is called after every `steps` times iteration.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.

```c
void ctmc_st_gs(const T1& Q, T2& x, Params& params, const Func& callback, MatT);
```

This function computes the steady-state vector of CTMC with Gauss-Seidal method.

- Q (in): An infinitesimal generator.
- x (out): The steady-state vector.
- params (inout): The parameters for computation. rtol, steps and maxiter are required. info, iter and rerror are stoted after the computation.
- callback: A callback function. The argument is Params. This is called after every `steps` times iteration.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.

```c
void ctmc_stsen_gs(const T1& Q, T2& x, const T3& b, const T4& pis, Params& params, const Func& callback, MatT);
```

This function computes the sensitivity of steady-state vector of CTMC with Gauss-Seidal method.

- Q (in): An infinitesimal generator.
- x (out): The sensitivity vector.
- b (in): The vector corresponding to the terms except for xQ.
- pis (in): The steady-state vector.
- params (inout): The parameters for computation. rtol, steps and maxiter are required. info, iter and rerror are stoted after the computation.
- callback: A callback function. The argument is Params. This is called after every `steps` times iteration.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.

```c
double unif(T1& Q, double ufact, MatT);
```

This function provides the unformized DTMC from an infinitesimal generator.

- Q (inout): An infinitesimal generator of CTMC. The output is the uniformized generator of DTMC.
- ufact (in): The uniformization factor. This is set as 1.01 as usual.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.

```c
void mexp_func(TR, T1& Q, const T2& x, T3& y, double t, Params& params, const Func& callback, MatT, VecT):
```

This function computes the matrix exponential:
```
y = exp(Qt) x
```
where x and y are vectors or (dense) matrices.

- TR: Type indicating to use a transpose of Q. An instance of marlib::NoTrans or marlib::Trans.
- Q (inout): An infinitesimal generator of CTMC. The output is the generator of DTMC.
- x (in): A vector or matrix which is multiplied by exp(Qt)
- y (out): A vector or matrix.
- t (in): A value as a scale.
- params (inout): The parameters for computation. ufact, eps and rmax are required. r is given as the upper bound of Poisson prob after the computation.
- callback: A callback function. The argument is Params. This is called when the uppoer bound of Poisson prob exceeds rmax of params.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.
- VecT: Type for x. An instance of marlib::ArrayT ro marlib::DenseMatrixT.

```c
void mexpint_func(TR, T1& Q, const T2& x, T3& y, T4& cy, double t, Params& params, const Func& callback, MatT, VecT);
```

This function computes the matrix exponential with cumulative value:
```
y = exp(Qt) x
cy = \int_0^t exp(Qs) ds x
```
where x, y, cy are vectors or (dense) matrices.

- TR: Type indicating to use a transpose of Q. An instance of marlib::NoTrans or marlib::Trans.
- Q (inout): An infinitesimal generator of CTMC. The output is the generator of DTMC.
- x (in): A vector or matrix which is multiplied by exp(Qt)
- y (out): A vector or matrix.
- cy (out): A vector or matrix as a cumulative value.
- t (in): A value as a scale.
- params (inout): The parameters for computation. ufact, eps and rmax are required. r is given as the upper bound of Poisson prob after the computation.
- callback: A callback function. The argument is Params. This is called when the uppoer bound of Poisson prob exceeds rmax of params.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.
- VecT: Type for x. An instance of marlib::ArrayT ro marlib::DenseMatrixT.

```c
void mexp_mix(TR, T1& Q, T2& x, T2& y, const T3& w, const T4& t, Params& params, const Func1& callback1, const Func2& callback2, MatT, VecT);
```

This function computes the integral of matrix exponential:
```
y = \int_0^\infty exp(Qt) f(t) dt x
```
where x and y are vectors or (dense) matrices and f is a probability density function.

- TR: Type indicating to use a transpose of Q. An instance of marlib::NoTrans or marlib::Trans.
- Q (inout): An infinitesimal generator of CTMC. The output is the generator of DTMC.
- x (inout): A vector or matrix which is multiplied by exp(Qt)
- y (out): A vector or matrix.
- w (in): The weight vector of f(t)
- t (in): The time instance for f(t)
- params (inout): The parameters for computation. ufact, eps and rmax are required. r is given as the upper bound of Poisson prob after the computation.
- callback1: A callback function. The argument is Params. This is called when the uppoer bound of Poisson prob exceeds rmax of params.
- callback2: A callback function. The argument is Params. This is called for each time instance.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.
- VecT: Type for x. An instance of marlib::ArrayT ro marlib::DenseMatrixT.

```c
void mexpint_mix(TR, T1& A, T2& x, T2& cx, T2& y, T2& cy, const T3& w, const T4& t, Params& params, const Func1& callback1, const Func2& callback2, MatT, VecT);
```

This function computes the integral of matrix exponential:
```
y = \int_0^\infty exp(Qt) f(t) dt x
cy = cx + \int_0^\infty \int_0^t exp(Qs) ds f(t) dt x
```
where x, cx, y and cy are vectors or (dense) matrices and f is a probability density function.

- TR: Type indicating to use a transpose of Q. An instance of marlib::NoTrans or marlib::Trans.
- Q (inout): An infinitesimal generator of CTMC. The output is the generator of DTMC.
- x (inout): A vector or matrix which is multiplied by exp(Qt). The output is not guaranteed.
- cx (inout): A vector or matrix.
- y (out): A vector or matrix.
- cy (out): A vector or matrix.
- w (in): The weight vector of f(t)
- t (in): The time instance for f(t)
- params (inout): The parameters for computation. ufact, eps and rmax are required. r is given as the upper bound of Poisson prob after the computation.
- callback1: A callback function. The argument is Params. This is called when the uppoer bound of Poisson prob exceeds rmax of params.
- callback2: A callback function. The argument is Params. This is called for each time instance.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.
- VecT: Type for x. An instance of marlib::ArrayT ro marlib::DenseMatrixT.

```c
void ctmc_tran(TR, T1& Q, T2& x, T3& cx, const T4& t, T5& res_x, T6& res_cx, Params& params, const Func1& callback1, const Func2& callback2, MatT, VecT);
```

This function provides a transient solution of CTMC, i.e., compute the following vector or matrices for each time point t:
```
res_x(t) = exp(Qt) x
res_cx(t) = cx + \int_0^t exp(Qs) ds x
```
where x and cx are vectors or (dense) matrices.

- TR: Type indicating to use a transpose of Q. An instance of marlib::NoTrans or marlib::Trans.
- Q (inout): An infinitesimal generator of CTMC. The output is the generator of DTMC.
- x (inout): An initial vector or matrix at time 0. The output is the probability vector at the last time of t.
- cx (inout): A vector or matrix. The output is the cumulative value at the last time of t.
- t (in): The time instance for f(t)
- res_x (out): A result for the probability vector. The result is given by a vector. If the result is a matrix, the result is stored as a column-major matrix format.
- res_cx (out): A result for the cumulative value. The result is given by a vector. If the result is a matrix, the result is stored as a column-major matrix format.
- params (inout): The parameters for computation. ufact, eps and rmax are required. r is given as the upper bound of Poisson prob after the computation.
- callback1: A callback function. The argument is Params. This is called when the uppoer bound of Poisson prob exceeds rmax of params.
- callback2: A callback function. The argument is Params. This is called for each time instance.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.
- VecT: Type for x. An instance of marlib::ArrayT ro marlib::DenseMatrixT.

```c
void ctmc_tran_rwd(TR, T1& Q, T2& x, T3& cx, const T4& rwd, const T5& t, T6& res_irwd, T7& res_crwd, Params& params, const Func1& callback1, const Func2& callback2, MatT, VecT1, VecT2);
```

This function provides a transient solution of CTMC, i.e., compute the following vector or matrices for each time point t:
```
res_irwd(t) = rwd^T exp(Qt) x
res_crwd(t) = rwd^T (cx + \int_0^t exp(Qs) ds) x
```
where x and cx are vectors or (dense) matrices.

- TR: Type indicating to use a transpose of Q. An instance of marlib::NoTrans or marlib::Trans.
- Q (inout): An infinitesimal generator of CTMC. The output is the generator of DTMC.
- x (inout): An initial vector or matrix at time 0. The output is the probability vector at the last time of t.
- cx (inout): A vector or matrix. The output is the cumulative value at the last time of t.
- rwd (in): The reward vector or matrix.
- t (in): The time instance for f(t)
- res_irwd (out): A result for the probability vector. The result is given by a vector. If the result is a matrix, the result is stored as a column-major matrix format.
- res_crwd (out): A result for the cumulative value. The result is given by a vector. If the result is a matrix, the result is stored as a column-major matrix format.
- params (inout): The parameters for computation. ufact, eps and rmax are required. r is given as the upper bound of Poisson prob after the computation.
- callback1: A callback function. The argument is Params. This is called when the uppoer bound of Poisson prob exceeds rmax of params.
- callback2: A callback function. The argument is Params. This is called for each time instance.
- MatT: Type for matrix Q. An instance of marlib::DenseMatrixT, marlib::CSRMatrixT, marlib::CSCMatrixT and marlib::COOMatrixT.
- VecT1: Type for x. An instance of marlib::ArrayT ro marlib::DenseMatrixT.
- VecT2: Type for rwd. An instance of marlib::ArrayT ro marlib::DenseMatrixT. Note that VecT2 is marlib::DenseMatrixT if VecT1 is marlib::DenseMatrixT.