
#ifndef _DLAPACK_H_
#define _DLAPACK_H_

namespace dblas {

  int dgesv(const int& n, const int& nrhs, double *A, const int& lda,
    int *ipiv, double *B, const int& ldb);

}

#endif
