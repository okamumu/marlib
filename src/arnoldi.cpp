#include <cmath>

#include "marlib.hpp"

namespace marlib {

  /*
  !     Description:

  !     Arnoldi process. The matrix A is approximated by h on the subspace

  !       h = v^T A v on Krylov subspace  {x, A*x, A^2*x, ... A^(m-1)*x}

  !     Parameters:

  !     A (in): a n-by-n square sparse matrix
  !     x (in): a vector with n elements
  !     h (out): m-by-m square Hessenburg matrix
  !     v (out): n-by-m matrix of orthogonal vectors on Krylov subspace {x, A*x, A^2*x, ... A^(m-1)*x}

  !     tol (in): tolerance error for happy breakdown
  !     ite (in): the number of times of iteration for Arnoldi process to reduce roundoff errors
  !               ite=1 is enough in most cases
  !     rnorm (out): 2-norm of vector x (and is used for checking for happy breakdown)
  !     info (out): the size of effective dimensions of h when a happy breakdown occurs

  !     work: working area. The required size is n.
  */

  template <typename ValueT, typename RangeT, typename MatrixT>
  size_type arnoldi(const trans_t& trans, const MatrixT& A,
    const vector<ValueT,RangeT>& x, dense_matrix<ValueT,RangeT> H, dense_matrix<ValueT,RangeT> V,
    ValueT& beta, ValueT& rnorm, const ValueT& tol, const size_type& ite) {
    assert(H.ncol() == V.ncol());
    assert(A.nrow() == A.ncol());
    size_type n = A.nrow();
    size_type m = H.ncol();
    H.set_range(range<RangeT>(1,m),range<RangeT>(1,m));
    V.set_range(range<RangeT>(1,n),range<RangeT>(1,m));
    H = 0;
    V = 0;

    rnorm = dnrm2(A);
    beta = dnrm2(x);

    range<RangeT> _(1,n);
    vector<ValueT,RangeT> tmp(n);

    vector<ValueT,RangeT> V1 = V(_,1);
    daxpy(1/beta, x, V1);
    for (RangeT j=1; j<=m; j++) {
      dgemv(trans, ValueT(1), A, V(_,j), ValueT(0), tmp);
      for (RangeT i=1; i<=j; i++) {
        vector<ValueT,RangeT> Vi = V(_,i);
        for (size_type l=1; l<=ite; l++) {
          ValueT r = ddot(Vi, tmp);
          H(i,j) += r;
          daxpy(-r, Vi, tmp);
        }
      }
      if (j != m) {
        H(j+1,j) = dnrm2(tmp);
        if (H(j+1,j) < rnorm * tol) {
          return j;
        }
        vector<ValueT,RangeT> Vj1 = V(_,j+1);
        daxpy(1/H(j+1,j), tmp, Vj1);
      }
    }
    return m;
  }

  template size_type arnoldi(const trans_t& trans, const dense_matrix<double,int>& A,
    const vector<double,int>& x, dense_matrix<double,int> H, dense_matrix<double,int> V,
    double& beta, double& rnorm, const double& tol, const size_type& ite);

  template size_type arnoldi(const trans_t& trans, const csr_matrix<double,int>& A,
    const vector<double,int>& x, dense_matrix<double,int> H, dense_matrix<double,int> V,
    double& beta, double& rnorm, const double& tol, const size_type& ite);
}
