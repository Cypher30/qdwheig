extern void dgemv_(char *TRANS, int *M, int *N, double *ALPHA, double *A, int *LDA, double *X, int *INCX, double *BETA, double *Y, int *INCY);
extern double dnrm2_(int *n, double *x, int *incx);
extern void dscal_(int *N, double *DA, double *DX, int *INCX);
extern double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
extern void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM, double *RCOND, double *WORK, int *IWORK, int *INFO);
extern void dgeqrf_(int *M, int *N, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
extern void dorgqr_(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
extern void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA, double *C, int *LDC);
extern void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
