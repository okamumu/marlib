namespace dblas {

  int dgesv(const int& n, const int& nrhs, double *A, const int& lda,
    int *ipiv, double *B, const int& ldb);

#ifdef F77LAPACK

extern "C" {
	void dgesv_(const int* n, const int* nrhs, double* A, const int* lda,
		int* ipiv, double* B, const int* ldb, int* info);
}

#define __DGESV__ dgesv_
#endif

inline int dgesv(const int& n, const int& nrhs, double *A, const int& lda,
	int *ipiv, double *B, const int& ldb) {
	int info;
	__DGESV__(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
	return info;
}

}
