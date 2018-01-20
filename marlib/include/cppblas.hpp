/*
  blas.hpp
*/

namespace marlib {

  // BLAS level 1

  // daxpy

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& daxpy(const ValueT alpha, const vector<ValueT,RangeT>& x, vector<ValueT,RangeT>& y) {
    assert(x.size() == y.size());
    for (RangeT i=y.begin(), j=x.begin(); i<=y.end(); i++, j++) {
      y(i) += alpha * x(j);
    }
    return y;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& daxpy(const ValueT alpha, const dense_matrix<ValueT,RangeT>& x,
    dense_matrix<ValueT,RangeT>& y) {
    assert(y.nrow() == x.nrow());
    assert(y.ncol() == x.ncol());
    for (RangeT j=y.cbegin(), jx=x.cbegin(); j<=y.cend(); j++, jx++) {
      for (RangeT i=y.rbegin(), ix=x.rbegin(); i<=y.rend(); i++, ix++) {
        y(i,j) += alpha * x(ix,jx);
      }
    }
    return y;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& daxpy(const ValueT alpha, const csr_matrix<ValueT,RangeT>& x,
    csr_matrix<ValueT,RangeT>& y) {
    assert(y.nnz() == x.nnz());
    for (size_type i=0; i<x.nnz(); i++) {
      y.value(i) += alpha * x.value(i);
    }
    return y;
  }

#ifdef F77BLAS
  template <>
  inline
  vector<double,int>& daxpy(const double alpha, const vector<double,int>& x, vector<double,int>& y) {
    assert(x.size() == y.size());
    dblas::daxpy(x.size(), alpha, x.ptr(), x.inc(), y.ptr(), y.inc());
    return y;
  }

  template <>
  inline
  dense_matrix<double,int>& daxpy(const double alpha, const dense_matrix<double,int>& x,
    dense_matrix<double,int>& y) {
    assert(x.nrow() == y.nrow());
    assert(x.ncol() == y.ncol());
    if (x.nrow() == x.ld()) {
      dblas::daxpy(x.nrow()*x.ncol(), alpha, x.ptr(), 1, y.ptr(), 1);
    } else {
      for (int j=y.cbegin(), jx=x.cbegin(); j<=y.cend(); j++, jx++) {
        dblas::daxpy(x.nrow(), alpha, &x(x.rbegin(),jx), 1, &y(y.rbegin(),j), 1);
      }
    }
    return y;
  }

  template <>
  inline
  csr_matrix<double,int>& daxpy(const double alpha, const csr_matrix<double,int>& x,
    csr_matrix<double,int>& y) {
    assert(y.nnz() == x.nnz());
    dblas::daxpy(x.nnz(), alpha, &x.value(0), 1, &y.value(0), 1);
    return y;
  }
#endif

  // ddot

  template <typename ValueT, typename RangeT>
  ValueT ddot(const vector<ValueT,RangeT>& x, const vector<ValueT,RangeT>& y) {
    assert(x.size() == y.size());
    ValueT tmp = 0;
    for (RangeT i=y.begin(), j=x.begin(); i<=y.end(); i++, j++) {
      tmp += y(i) * x(j);
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  ValueT ddot(const dense_matrix<ValueT,RangeT>& x, const dense_matrix<ValueT,RangeT>& y) {
    assert(y.nrow() == x.nrow());
    assert(y.ncol() == x.ncol());
    ValueT tmp = 0;
    for (RangeT j=y.cbegin(), jx=x.cbegin(); j<=y.cend(); j++, jx++) {
      for (RangeT i=y.rbegin(), ix=x.rbegin(); i<=y.rend(); i++, ix++) {
        tmp += y(i,j) * x(ix,jx);
      }
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  ValueT ddot(const csr_matrix<ValueT,RangeT>& x, const csr_matrix<ValueT,RangeT>& y) {
    assert(y.nnz() == x.nnz());
    ValueT tmp = 0;
    for (size_type i=0; i<x.nnz(); i++) {
      tmp += y.value(i) * x.value(i);
    }
    return tmp;
  }

#ifdef F77BLAS
  template <>
  inline
  double ddot<double,int>(const vector<double,int>& x, const vector<double,int>& y) {
    assert(x.size() == y.size());
    return dblas::ddot(x.size(), x.ptr(), x.inc(), y.ptr(), y.inc());
  }

  template <>
  inline
  double ddot<double,int>(const dense_matrix<double,int>& x, const dense_matrix<double,int>& y) {
    assert(x.nrow() == y.nrow());
    assert(x.ncol() == y.ncol());
    if (x.nrow() == x.ld()) {
      return dblas::ddot(x.nrow()*x.ncol(), x.ptr(), 1, y.ptr(), 1);
    } else {
      double tmp = 0;
      for (int j=y.cbegin(), jx=x.cbegin(); j<=y.cend(); j++, jx++) {
        tmp += dblas::ddot(x.nrow(), &x(x.rbegin(),jx), 1, &y(y.rbegin(),j), 1);
      }
      return tmp;
    }
  }

  template <>
  inline
  double ddot<double,int>(const csr_matrix<double,int>& x, const csr_matrix<double,int>& y) {
    assert(y.nnz() == x.nnz());
    return dblas::ddot(x.nnz(), &x.value(0), 1, &y.value(0), 1);
  }
#endif

  // dnrm2

  template <typename ValueT, typename RangeT>
  ValueT dnrm2(const vector<ValueT,RangeT>& x) {
    return std::sqrt(ddot(x,x));
  }

  template <typename ValueT, typename RangeT>
  ValueT dnrm2(const dense_matrix<ValueT,RangeT>& x) {
    return std::sqrt(ddot(x,x));
  }

  template <typename ValueT, typename RangeT>
  ValueT dnrm2(const csr_matrix<ValueT,RangeT>& x) {
    return std::sqrt(ddot(x,x));
  }

  // dsum

  template <typename ValueT, typename RangeT>
  ValueT dsum(const vector<ValueT,RangeT>& x) {
    ValueT tmp = 0;
    for (RangeT i=x.begin(); i<=x.end(); i++) {
      tmp += x(i);
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  ValueT dsum(const dense_matrix<ValueT,RangeT>& x) {
    ValueT tmp = 0;
    for (RangeT j=x.cbegin(); j<=x.cend(); j++) {
      for (RangeT i=x.rbegin(); i<=x.rend(); i++) {
        tmp += x(i,j);
      }
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  ValueT dsum(const csr_matrix<ValueT,RangeT>& x) {
    ValueT tmp = 0;
    for (size_type i=0; i<x.nnz(); i++) {
      tmp += x.value(i);
    }
    return tmp;
  }

  // dasum

  template <typename ValueT, typename RangeT>
  ValueT dasum(const vector<ValueT,RangeT>& x) {
    ValueT tmp = 0;
    for (RangeT i=x.begin(); i<=x.end(); i++) {
      tmp += std::abs(x(i));
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  ValueT dasum(const dense_matrix<ValueT,RangeT>& x) {
    ValueT tmp = 0;
    for (RangeT j=x.cbegin(); j<=x.cend(); j++) {
      for (RangeT i=x.rbegin(); i<=x.rend(); i++) {
        tmp += std::abs(x(i,j));
      }
    }
    return tmp;
  }

  template <typename ValueT, typename RangeT>
  ValueT dasum(const csr_matrix<ValueT,RangeT>& x) {
    ValueT tmp = 0;
    for (size_type i=0; i<x.nnz(); i++) {
      tmp += std::abs(x.value(i));
    }
    return tmp;
  }

  // iamax

  template <typename ValueT, typename RangeT>
  RangeT iamax(const vector<ValueT,RangeT>& x) {
    RangeT maxi = x.begin();
    ValueT maxv = std::abs(x(maxi));
    for (RangeT i=x.begin(); i<=x.end(); i++) {
      ValueT tmp = std::abs(x(i));
      if (tmp > maxv) {
        maxi = i;
        maxv = tmp;
      }
    }
    return maxi;
  }

  // damax

  template <typename ValueT, typename RangeT>
  ValueT damax(const vector<ValueT,RangeT>& x) {
    ValueT maxv = 0;
    for (RangeT i=x.begin(); i<=x.end(); i++) {
      ValueT tmp = std::abs(x(i));
      if (tmp > maxv) {
        maxv = tmp;
      }
    }
    return maxv;
  }

  template <typename ValueT, typename RangeT>
  ValueT damax(const vector<ValueT*,RangeT>& x) {
    ValueT maxv = 0;
    for (RangeT i=x.begin(); i<=x.end(); i++) {
      ValueT tmp = std::abs(x(i));
      if (tmp > maxv) {
        maxv = tmp;
      }
    }
    return maxv;
  }

  template <typename ValueT, typename RangeT>
  ValueT damax(const dense_matrix<ValueT,RangeT>& x) {
    ValueT mx = 0;
    for (RangeT j=x.cbegin(); j<=x.cend(); j++) {
      for (RangeT i=x.rbegin(); i<=x.rend(); i++) {
        ValueT tmp = std::abs(x(i,j));
        if (mx < tmp) {
          mx = tmp;
        }
      }
    }
    return mx;
  }

  template <typename ValueT, typename RangeT>
  ValueT damax(const csr_matrix<ValueT,RangeT>& x) {
    ValueT mx = 0;
    for (size_type i=0; i<x.nnz(); i++) {
      ValueT tmp = std::abs(x.value(i));
      if (mx < tmp) {
        mx = tmp;
      }
    }
    return mx;
  }

  // BLAS level 2

  // dgemv

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& dgemv(
    const trans_t& trans,
    const ValueT alpha,
    const dense_matrix<ValueT,RangeT>& A,
    const vector<ValueT,RangeT>& x,
    const ValueT beta,
    vector<ValueT,RangeT>& y) {

    switch(trans) {
    case NoTrans:
      assert(A.ncol() == x.size());
      assert(A.nrow() == y.size());
      y *= beta;
      for (RangeT j=y.begin(), ia=A.rbegin(); j<=y.end(); j++, ia++) {
        for (RangeT i=x.begin(), ja=A.cbegin(); i<=x.end(); i++, ja++) {
          y(j) += alpha * A(ia,ja) * x(i);
        }
      }
      break;
    case Trans:
      assert(A.nrow() == x.size());
      assert(A.ncol() == y.size());
      y *= beta;
      for (RangeT j=y.begin(), ja=A.cbegin(); j<=y.end(); j++, ja++) {
        for (RangeT i=x.begin(), ia=A.rbegin(); i<=x.end(); i++, ia++) {
          y(j) = alpha * A(ia,ja) * x(i);
        }
      }
      break;
    }
    return y;
  }

  template <typename ValueT, typename RangeT>
  vector<ValueT,RangeT>& dgemv(const trans_t& trans,
    const ValueT alpha, const csr_matrix<ValueT,RangeT>& A,
    const vector<ValueT,RangeT>& x, const ValueT beta, vector<ValueT,RangeT>& y) {

    switch (trans) {
      case NoTrans:
      y *= beta;
      for (size_type i=0; i<A.nrow(); i++) {
        for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
          size_type j = A.colind(z);
          y(i + y.begin()) += alpha * A.value(z) * x(j + x.begin());
        }
      }
      break;
      case Trans:
      y *= beta;
      for (size_type i=0; i<A.nrow(); i++) {
        for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
          size_type j = A.colind(z);
          y(j + y.begin()) += alpha * A.value(z) * x(i + x.begin());
        }
      }
      break;
    }
    return y;
  }

#ifdef F77BLAS
  template <>
  inline
  vector<double,int>& dgemv(
    const trans_t& trans, const double alpha,
    const dense_matrix<double,int>& A,
    const vector<double,int>& x, const double beta, vector<double,int>& y) {
    dblas::dgemv(trans, A.nrow(), A.ncol(), alpha, A.ptr(), A.ld(),
      x.ptr(), x.inc(), beta, y.ptr(), y.inc());
    return y;
  }
#endif

  // dger

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dger(const trans_t& trans,
    const ValueT alpha, const vector<ValueT,RangeT>& x, const vector<ValueT,RangeT>& y,
    dense_matrix<ValueT,RangeT>& A) {

    switch (trans) {
      case NoTrans:
      assert(A.nrow() == x.size());
      assert(A.ncol() == y.size());
      for (RangeT ja=A.cbegin(), j=y.begin(); j<=y.end(); ja++, j++) {
        for (RangeT ia=A.rbegin(), i=x.begin(); i<=x.end(); ia++, i++) {
          A(ia,ja) += alpha * x(i) * y(j);
        }
      }
      break;
      case Trans:
      assert(A.ncol() == x.size());
      assert(A.nrow() == y.size());
      for (RangeT ja=A.cbegin(), j=x.begin(); j<=x.end(); ja++, j++) {
        for (RangeT ia=A.rbegin(), i=y.begin(); i<=y.end(); ia++, i++) {
          A(ia,ja) += alpha * x(j) * y(i);
        }
      }
      break;
    }
    return A;
  }

  template <typename ValueT, typename RangeT>
  csr_matrix<ValueT,RangeT>& dger(const trans_t& trans,
    const ValueT alpha, const vector<ValueT,RangeT>& x, const vector<ValueT,RangeT>& y,
    csr_matrix<ValueT,RangeT>& A) {

    switch (trans) {
      case NoTrans:
      for (size_type i=0; i<A.nrow(); i++) {
        for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
          size_type j = A.colind(z);
          ValueT tmp = alpha;
          tmp *= x(i + x.begin());
          tmp *= y(j + y.begin());
          A.value(z) += tmp;
        }
      }
      break;
      case Trans:
      for (size_type i=0; i<A.nrow(); i++) {
        for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
          size_type j = A.colind(z);
          ValueT tmp = alpha;
          tmp *= x(j + x.begin());
          tmp *= y(i + y.begin());
          A.value(z) += tmp;
        }
      }
      break;
    }
    return A;
  }

#ifdef F77BLAS
  template <>
  inline
  dense_matrix<double,int>& dger(
    const trans_t& trans,
    const double alpha, const vector<double,int>& x, const vector<double,int>& y,
    dense_matrix<double,int>& A) {
    switch (trans) {
      case NoTrans:
        dblas::dger(A.nrow(), A.ncol(), alpha, x.ptr(), x.inc(), y.ptr(), y.inc(), A.ptr(), A.ld());
        break;
      case Trans:
        dblas::dger(A.nrow(), A.ncol(), alpha, y.ptr(), y.inc(), x.ptr(), x.inc(), A.ptr(), A.ld());
        break;
    }
    return A;
  }
#endif

  // BLAS level 3

  // dgemm

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dgemm(
    const trans_t& transA,
    const trans_t& transB,
    const ValueT alpha,
    const dense_matrix<ValueT,RangeT>& A,
    const dense_matrix<ValueT,RangeT>& B,
    const ValueT beta,
    dense_matrix<ValueT,RangeT>& C) {

    switch (transA) {
      case NoTrans:
      switch (transB) {
        case NoTrans:
        assert(C.nrow() == A.nrow() && C.ncol() == B.ncol() && A.ncol() == B.nrow());
        for (RangeT jc=C.cbegin(), jb=B.cbegin(); jc<=C.cend(); jc++, jb++) {
          for (RangeT ic=C.rbegin(), ia=A.rbegin(); ic<=C.rend(); ic++, ia++) {
            ValueT tmp = 0;
            for (RangeT ja=A.cbegin(), ib=B.rbegin(); ja<=A.cend(); ja++, ib++) {
              tmp += A(ia,ja) * B(ib,jb);
            }
            tmp *= alpha;
            C(ic,jc) *= beta;
            C(ic,jc) += tmp;
          }
        }
        break;
        case Trans:
        assert(C.nrow() == A.nrow() && C.ncol() == B.nrow() && A.ncol() == B.ncol());
        for (RangeT jc=C.cbegin(), ib=B.rbegin(); jc<=C.cend(); jc++, ib++) {
          for (RangeT ic=C.rbegin(), ia=A.rbegin(); ic<=C.rend(); ic++, ia++) {
            ValueT tmp = 0;
            for (RangeT ja=A.cbegin(), jb=B.cbegin(); ja<=A.cend(); ja++, jb++) {
              tmp += A(ia,ja) * B(ib,jb);
            }
            tmp *= alpha;
            C(ic,jc) *= beta;
            C(ic,jc) += tmp;
          }
        }
        break;
      }
      break;
      case Trans:
      switch (transB) {
        case NoTrans:
        assert(C.nrow() == A.ncol() && C.ncol() == B.ncol() && A.nrow() == B.nrow());
        for (RangeT jc=C.cbegin(), jb=B.cbegin(); jc<=C.cend(); jc++, jb++) {
          for (RangeT ic=C.rbegin(), ja=A.cbegin(); ic<=C.rend(); ic++, ja++) {
            ValueT tmp = 0;
            for (RangeT ia=A.rbegin(), ib=B.rbegin(); ia<=A.rend(); ia++, ib++) {
              tmp += A(ia,ja) * B(ib,jb);
            }
            tmp *= alpha;
            C(ic,jc) *= beta;
            C(ic,jc) += tmp;
          }
        }
        break;
        case Trans:
        assert(C.nrow() == A.ncol() && C.ncol() == B.nrow() && A.nrow() == B.ncol());
        for (RangeT jc=C.cbegin(), ib=B.rbegin(); jc<=C.cend(); jc++, ib++) {
          for (RangeT ic=C.rbegin(), ja=A.cbegin(); ic<=C.rend(); ic++, ja++) {
            ValueT tmp = 0;
            for (RangeT ia=A.rbegin(), jb=B.cbegin(); ia<=A.rend(); ia++, jb++) {
              tmp += A(ia,ja) * B(ib,jb);
            }
            tmp *= alpha;
            C(ic,jc) *= beta;
            C(ic,jc) += tmp;
          }
        }
        break;
      }
      break;
    }
    return C;
  }

  template <typename ValueT, typename RangeT>
  dense_matrix<ValueT,RangeT>& dgemm(
    const trans_t& transA,
    const trans_t& transB,
    const ValueT alpha,
    const csr_matrix<ValueT,RangeT>& A,
    const dense_matrix<ValueT,RangeT>& B,
    const ValueT beta,
    dense_matrix<ValueT,RangeT>& C) {

    switch (transB) {
      case NoTrans:
      switch (transA) {
        case NoTrans:
        C *= beta;
        for (size_type i=0; i<A.nrow(); i++) {
          for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            for (RangeT jx=B.cbegin(), jy=C.cbegin(); jx<=B.cend(); jx++, jy++) {
              C(i + C.rbegin(), jy) += alpha * A.value(z) * B(j + B.rbegin(), jx);
            }
          }
        }
        break;
        case Trans:
        C *= beta;
        for (size_type i=0; i<A.nrow(); i++) {
          for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            for (RangeT jx=B.cbegin(), jy=C.cbegin(); jx<=B.cend(); jx++, jy++) {
              C(j + C.rbegin(), jy) += alpha * A.value(z) * B(i + B.rbegin(), jx);
            }
          }
        }
        break;
      }
      break;
      case Trans:
      switch (transA) {
        case NoTrans:
        C *= beta;
        for (size_type i=0; i<A.nrow(); i++) {
          for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            for (RangeT jx=B.rbegin(), jy=C.rbegin(); jx<=B.rend(); jx++, jy++) {
              C(jy, i + C.cbegin()) += alpha * A.value(z) * B(jx, j + B.cbegin());
            }
          }
        }
        break;
        case Trans:
        C *= beta;
        for (size_type i=0; i<A.nrow(); i++) {
          for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            for (RangeT jx=B.rbegin(), jy=C.rbegin(); jx<=B.rend(); jx++, jy++) {
              C(jy, j + C.cbegin()) += alpha * A.value(z) * C(jx, i + C.cbegin());
            }
          }
        }
        break;
      }
      break;
    }
    return C;
  }

#ifdef F77BLAS
  template <>
  inline
  dense_matrix<double,int>& dgemm(
    const trans_t& transA,
    const trans_t& transB,
    const double alpha,
    const dense_matrix<double,int>& A,
    const dense_matrix<double,int>& B,
    const double beta,
    dense_matrix<double,int>& C) {
      switch (transA) {
        case NoTrans:
        switch (transB) {
          case NoTrans:
          assert(C.nrow() == A.nrow() && C.ncol() == B.ncol() && A.ncol() == B.nrow());
          dblas::dgemm(transA, transB, C.nrow(), C.ncol(), A.ncol(), alpha, A.ptr(), A.ld(), B.ptr(), B.ld(), beta, C.ptr(), C.ld());
          break;
          case Trans:
          assert(C.nrow() == A.nrow() && C.ncol() == B.nrow() && A.ncol() == B.ncol());
          dblas::dgemm(transA, transB, C.nrow(), C.ncol(), A.ncol(), alpha, A.ptr(), A.ld(), B.ptr(), B.ld(), beta, C.ptr(), C.ld());
          break;
        }
        break;
        case Trans:
        switch (transB) {
          case NoTrans:
          assert(C.nrow() == A.ncol() && C.ncol() == B.ncol() && A.nrow() == B.nrow());
          dblas::dgemm(transA, transB, C.nrow(), C.ncol(), A.nrow(), alpha, A.ptr(), A.ld(), B.ptr(), B.ld(), beta, C.ptr(), C.ld());
          break;
          case Trans:
          assert(C.nrow() == A.ncol() && C.ncol() == B.nrow() && A.nrow() == B.ncol());
          dblas::dgemm(transA, transB, C.nrow(), C.ncol(), A.nrow(), alpha, A.ptr(), A.ld(), B.ptr(), B.ld(), beta, C.ptr(), C.ld());
          break;
        }
        break;
      }
    return C;
  }
#endif

  // lapack

  // dgesv: solve A X = B
  // The solution is assgined to B

#ifdef F77BLAS
  template <>
  inline
  dense_matrix<double,int>& dgesv(
    const dense_matrix<double,int>& A, dense_matrix<double,int>& B) {
    assert(A.nrow() == A.ncol());
    dense_matrix<double,int> MA = A.clone();
    array<int> ipiv(A.nrow());
    int info = dblas::dgesv(A.nrow(), B.ncol(), MA.ptr(), MA.ld(), &ipiv[0], B.ptr(), B.ld());
    assert(info == 0);
    return B;
  }

  template <>
  inline
  vector<double,int>& dgesv(
    const dense_matrix<double,int>& A, vector<double,int>& B) {
    assert(A.nrow() == A.ncol());
    dense_matrix<double,int> MA = A.clone();
    array<int> ipiv(A.nrow());
    int info = dblas::dgesv(A.nrow(), 1, MA.ptr(), MA.ld(), &ipiv[0], B.ptr(), B.size());
    assert(info == 0);
    return B;
  }
#endif

}
