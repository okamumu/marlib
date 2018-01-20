
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
    ValueT& beta, ValueT& rnorm, ValueT tol, size_type ite);

}
