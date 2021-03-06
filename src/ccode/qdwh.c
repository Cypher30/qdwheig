#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<complex.h>
#include<float.h>
#include "../include/reflapack.h"


#define DEBUG 0


//Estimate the 2-norm of a symmetric matrix A (fill upper triangular part)
double dsynormest(int M, const double *A, int lda, double tol, double *WORK)
{
	//Initialization
	srand(time(NULL));
	double *y = WORK;
	double *x = WORK + M;
	
	//Construct initial vector
	for (int i = 0; i < M; i++)
	{
		y[i] = (double) rand() / (double) RAND_MAX;
	}
	double nest = 0.0;
	double n0, normx;
	int INCX = 1, INCY = 1;
	double alpha = 1.0, beta = 0.0;

	while(1)
	{
		n0 = nest;
		dsymv_("U", &M, &alpha, A, &lda, y, &INCY, &beta, x, &INCX);
		normx = dnrm2_(&M, x, &INCX);
		if (normx == 0.0)
		{
			for (int i = 0; i < M; i++)
			{
				x[i] = rand();
			}
		}
		else
		{
			double quo = 1.0 / normx;
			dscal_(&M, &quo, x, &INCX);
		}
		dsymv_("U", &M, &alpha, A, &lda, x, &INCX, &beta, y, &INCY);
		nest = dnrm2_(&M, y, &INCY);
		if (abs(n0 - nest) <= tol * nest) break;
	}
	return nest;
}


//QDWH algorithm to compute polar factor U of a symmetric matrix A
//INPUT: 
//M: Dimension of A
//A: M * M symmetric matrix with leading dimension lda, column major, fill upper triangular part and store polar factor (upper triangular part)
//lda: Leading dimenstion of A
//U: Work space dimension M * M
//B: Working space, double, dimenstion 2M * M
//IPIV: Working space, int, dimension M
//QWORK: Working space, double, dimension M
//WORK: Working space, dimension LWORK
//LWORK: Size of WORK, if LWORK = -1, it returns the size of WORK in WORK[0]
//Output:
//U: polar factor of matrix A
void dsyqdwh(int M, double *A, int lda, double *U, double *B, int *IPIV, double *QWORK, int *IWORK, double *WORK, int LWORK)
{
	if (LWORK == -1)
	{
		int BM = 2 * M;
		int BN = M;
		int INFO;
		dgeqrf_(&BM, &BN, B, &BM, QWORK, WORK, &LWORK, &INFO);
		double temp = *WORK;
		dsytrf_("U", &M, A, &lda, IPIV, WORK, &LWORK, &INFO);
		*WORK = (*WORK <= temp) ? temp : *WORK;
		return;
	}

	//Construct initial matrix U = A / dsynormest(A, 3e-1)
	//and choose initial L = 0.9 / condest(A)
	double *tempA = A; 
	double *tempU = U;
	double quo = 1.0 / dsynormest(M, A, lda, 3e-1, WORK);

	for (int j = 0; j < M; j++)
	{
		memcpy(tempU, tempA, (j + 1) * sizeof(double));
		tempA += lda;
		tempU += M;
	}
	
	double Anorm = dlansy_("1", "U", &M, A, &lda, WORK);
	double RCOND;
	int INFO;
	dsytrf_("U", &M, A, &lda, IPIV, WORK, &LWORK, &INFO);
	dsycon_("U", &M, A, &lda, IPIV, &Anorm, &RCOND, WORK, IWORK, &INFO);
	double L = 0.9 * RCOND;
	//printf("RCOND = %f\n", RCOND);
	
	/*
	tempU = U;
	tempA = A;
	for (int j = 0; j < M; j++)
	{
		memcpy(tempU, tempA, (j + 1) * sizeof(double));
		tempA += lda;
		tempU += ldu;
	}
	*/

	//Fill the lower triangular part of U
	tempU = U;
	for (int j = 0; j < M - 1; j++)
	{
		for (int i = j + 1; i < M; i++)
		{
			tempU[i] = U[j + i * M];
		}
		tempU += M;
	}

	//Scaling
	tempU = U;
	int INCU = 1;
	for (int j = 0; j < M; j++)
	{
		dscal_(&M, &quo, tempU, &INCU);
		tempU += M;
	}

	
	//Start iteration
	double tol1 = 10 * DBL_EPSILON / 2;
	double tol2 = pow(tol1, 1.0 / 3.0);
	int it = 0;
	double delta = 100.0;
	double *tempB;

	//printf("tol1 = %e, tol2 = %e\n", tol1, tol2);
	
	
	while(it == 0 || delta > tol2 || abs(1 - L) > tol1)
	{
		it = it + 1;
		
		tempA = A;
		tempU = U;
		for (int j = 0; j < M; j++)
		{
			memcpy(tempA, tempU, (j + 1) * sizeof(double));
			tempA += lda;
			tempU += M;
		}

		double L2 = L * L;
		double dd = pow(4.0 * (1.0 - L2) / (L2 * L2), 1.0 / 3.0);
		double sqd = sqrt(1 + dd);
		double a = sqd + sqrt(8.0 - 4.0 * dd + 8.0 * (2.0 - L2) / (L2 * sqd)) / 2.0;
		double b = pow(a - 1.0, 2) / 4.0; 
		double c = a + b - 1.0;

		//update L
		L = L * (a + b * L2) / (1.0 + c * L2);
		L = (L > 1) ? 1.0 : L;
		
		if (c > 100)
		{
			tempB = B;
			tempU = U;
			for (int j = 0; j < M; j++)
			{
				memcpy(tempB, tempU, M * sizeof(double));
				double DA = sqrt(c);
				int INCX = 1;
				dscal_(&M, &DA, tempB, &INCX);
				tempB += 2 * M;
				tempU += M;
			}
			tempB = B;
			for (int j = 0; j < M; j++)
			{
				memset(tempB + M, 0, M * sizeof(double));
				tempB[M + j] = 1.0;
				tempB += 2 * M;
			}

			int BM = 2 * M;
			int BN = M;
			int LDB = BM;
			int INFO;
			dgeqrf_(&BM, &BN, B, &LDB, QWORK, WORK, &LWORK, &INFO);
			dorgqr_(&BM, &BN, &BN, B, &LDB, QWORK, WORK, &LWORK, &INFO);

			double ALPHA = (a - b / c) / sqrt(c);
			double BETA = b / c;
			dgemm_("N", "T", &M, &M, &BN, &ALPHA, B, &LDB, B + M, &LDB, &BETA, U, &M);
		}
		else
		{
			int LDB = 2 * M;
			tempB = B;
			tempU = U;
			for (int j = 0; j < M; j++)
			{
				memcpy(tempB, tempU, M * sizeof(double));
				tempB += 2 * M;
				tempU += M;
			}
			tempB = B;
			for (int j = 0; j < M; j++)
			{
				memset(tempB + M, 0, M * sizeof(double));
				tempB[M + j] = 1.0;
				tempB += 2 * M;
			}
#if DEBUG
			for (int i = 0; i < 2 * M; i++)
			{
				for (int j = 0; j < M; j++)
				{
					printf("%.6f ", B[i + j * LDB]);
				}
				printf("\n");
			}
			printf("\n");
#endif
			double BETA = 1.0;
			int INFO;
			dsyrk_("U", "T", &M, &M, &c, B, &LDB, &BETA, B + M, &LDB);
#if DEBUG
			printf("c = %.6f\n", c);
			for (int i = 0; i < 2 * M; i++)
			{
				for (int j = 0; j < M; j++)
				{
					printf("%.6f ", B[i + j * LDB]);
				}
				printf("\n");
			}
			printf("\n");
#endif
			dpotrf_("U", &M, B + M, &LDB, &INFO);
#if DEBUG
			for (int i = 0; i < 2 * M; i++)
			{
				for (int j = 0; j < M; j++)
				{
					printf("%.6f ", B[i + j * LDB]);
				}
				printf("\n");
			}
			printf("\n");
#endif
			for (int i = 0; i < M; i++)
			{
				IPIV[i] = i + 1;
			}
			dgetrs_("T", &M, &M, B + M, &LDB, IPIV, B, &LDB, &INFO);
			dgetrs_("N", &M, &M, B + M, &LDB, IPIV, B, &LDB, &INFO);
#if DEBUG
			for (int i = 0; i < 2 * M; i++)
			{
				for (int j = 0; j < M; j++)
				{
					printf("%.6f ", B[i + j * LDB]);
				}
				printf("\n");
			}
			printf("\n");
#endif
			double ALPHA = a - b / c;
			BETA = b / c;
			int INCX = 1, INCY = 1;
			tempB = B;
			tempU = U;
			for (int i = 0; i < M; i++)
			{
				dscal_(&M, &BETA, tempU, &INCY);
				daxpy_(&M, &ALPHA, tempB, &INCX, tempU, &INCY);
				tempB += 2 * M;
				tempU += M;
			}
		}
		tempU = U;
		tempA = A;
		for (int j = 0; j < M; j++)
		{
			for (int i = 0; i < j + 1; i++)
			{
				tempA[i] -= tempU[i];
			}
			tempU += M;
			tempA += lda;
		}
		delta = dlansy_("F", "U", &M, A, &lda, WORK);
		//printf("delta = %e, a = %e, b = %e, c = %e, L = %e, dd = %e, iter = %d\n", delta, a, b, c, L, dd, it);
	}
	
	
	//printf("%d\n", it);
	tempU = U;
	tempA = A;
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < j + 1; i++)
		{
			tempA[i] = (tempU[i] + U[j + i * M]) / 2;
		}
		tempU += M;
		tempA += lda;
	}
}
