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
//A: M * M symmetric matrix with leading dimension lda, column major, fill upper triangular part
//lda: Leading dimenstion of A
//U: Space to store polar factor
//ldu: Leading dimension of U
//B: Working space, double, dimenstion 2M * M
//Uprev: Working space, double, dimension M * M
//IPIV: Working space, int, dimension M
//QWORK: Working space, double, dimension M
//WORK: Working space, dimension LWORK
//LWORK: Size of WORK, if LWORK = -1, it returns the size of WORK in WORK[0]
//Output:
//U: polar factor of matrix A
void dsyqdwh(int M, double *A, int lda, double *U, int ldu, double *B, double *Uprev, int *IPIV, double *QWORK, int *IWORK, double *WORK, int LWORK)
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
		tempU += ldu;
	}

#if DEBUG
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < M; i++)
		{
			printf("%6f", A[i + j * lda]);
		}
		printf("\n");
	}
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i <= j; i++)
		{
			printf("%6f", U[i + j * ldu]);
		}
		printf("\n");
	}
#endif
	
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
			tempU[i] = U[j + i * ldu];
		}
		tempU += ldu;
	}

	//Scaling
	tempU = U;
	int INCU = 1;
	for (int j = 0; j < M; j++)
	{
		dscal_(&M, &quo, tempU, &INCU);
		tempU += M;
	}

#if DEBUG
	printf("quo = %6f\n", quo);
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < M; i++)
		{
			printf("%6f", U[i + j * ldu]);
		}
		printf("\n");
	}
#endif
	
	//Start iteration
	double tol1 = 10 * DBL_EPSILON / 2;
	double tol2 = pow(tol1, 1.0 / 3.0);
	int it = 0;
	double delta = 100.0;
	double *tempUprev;
	double *tempB;

	//printf("tol1 = %e, tol2 = %e\n", tol1, tol2);
	
	
	while(it == 0 || delta > tol2 || abs(1 - L) > tol1)
	{
		it = it + 1;
		
		tempUprev = Uprev;
		tempU = U;
		for (int j = 0; j < M; j++)
		{
			memcpy(tempUprev, tempU, (j + 1) * sizeof(double));
			tempUprev += M;
			tempU += ldu;
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
		
		tempB = B;
		tempU = U;
		for (int j = 0; j < M; j++)
		{
			memcpy(tempB, tempU, M * sizeof(double));
			double DA = sqrt(c);
			int INCX = 1;
			dscal_(&M, &DA, tempB, &INCX);
			tempB += 2 * M;
			tempU += ldu;
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
		dgemm_("N", "T", &M, &M, &BN, &ALPHA, B, &LDB, B + M, &LDB, &BETA, U, &ldu);

		tempU = U;
		tempUprev = Uprev;
		for (int j = 0; j < M; j++)
		{
			for (int i = 0; i < j + 1; i++)
			{
				tempUprev[i] = tempUprev[i] - tempU[i];
			}
			tempU += ldu;
			tempUprev += M;
		}
		delta = dlansy_("F", "U", &M, Uprev, &M, WORK);
		//printf("delta = %e, a = %e, b = %e, c = %e, L = %e, dd = %e, iter = %d\n", delta, a, b, c, L, dd, it);
	}
	//printf("%d\n", it);
	tempU = U;
	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < M; i++)
		{
			if (i > j)
			{
				tempU[i] = (tempU[i] + U[j + i * ldu]) / 2;
			}
			else
			{
				tempU[i] = U[j + i * ldu];
			}
		}
		tempU += ldu;
	}
}
