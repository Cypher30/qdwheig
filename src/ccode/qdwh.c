#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<complex.h>
#include<float.h>
#include "../include/reflapack.h"




//Estimate the 2-norm of a symmetric matrix A (fill upper triangular part)
double dsynormest(int M, const double *A, int lda, double tol)
{
	//Initialization
	srand(time(NULL));
	double *y = (double *)malloc(M * sizeof(double));
	double *x = (double *)malloc(M * sizeof(double));
	
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
//Output:
//U: polar factor of matrix A
void dsyqdwh(int M, const double *A, int lda, double *U, int ldu)
{
	//Construct initial matrix U = A / dsynormest(A, 3e-1)
	//and choose initial L = 0.9 / condest(A)
	const double *tempA = A; 
	double *tempU = U;

	for (int j = 0; j < M; j++)
	{
		memcpy(tempU, tempA, (j + 1) * sizeof(double));
		tempA += lda;
		tempU += ldu;
	}
	
	double *WORK = (double *)malloc(M * sizeof(double));
	int *IWORK = (int *)malloc(M * sizeof(int));
	double Anorm = dlansy_("1", "U", &M, U, &ldu, WORK);
	int *IPIV = (int *)malloc(M * sizeof(int));
	double RCOND;
	int INFO;
	int LWORK = -1;
	dsytrf_("U", &M, U, &ldu, IPIV, WORK, &LWORK, &INFO);
	LWORK = *WORK;
	free(WORK);
	WORK = (double *)malloc(LWORK * sizeof(double));
	dsytrf_("U", &M, U, &ldu, IPIV, WORK, &LWORK, &INFO);
	free(WORK);
	WORK = (double *)malloc(2 * M * sizeof(double));
	dsycon_("U", &M, U, &ldu, IPIV, &Anorm, &RCOND, WORK, IWORK, &INFO);
	double L = 0.9 * RCOND;
	//printf("RCOND = %f\n", RCOND);
	
	tempU = U;
	tempA = A;
	for (int j = 0; j < M; j++)
	{
		memcpy(tempU, tempA, (j + 1) * sizeof(double));
		tempA += lda;
		tempU += ldu;
	}

	//Fill the lower triangular part of U
	tempU = U;
	for (int j = 0; j < M - 1; j++)
	{
		for (int i = j + 1; i < M; i++)
		{
			tempU[i] = U[j + i * ldu];
		}
	}

	//Scaling
	tempU = U;
	int INCU = 1;
	double quo = 1.0 / dsynormest(M, A, lda, 3e-1);
	for (int j = 0; j < M; j++)
	{
		dscal_(&M, &quo, tempU, &INCU);
		tempU += M;
	}
	
	
	//Start iteration
	double *Uprev = (double *)malloc(M * M * sizeof(double));
	double *B = (double *)malloc(2 * M * M * sizeof(double));
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
			memcpy(tempUprev, tempU, M * sizeof(double));
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
		
		memset(B, 0, 2 * M * M * sizeof(double));
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
		for (int j = 0; j < M; j++)
		{
			B[M + j + j * 2 * M] = 1.0;
		}

		int BM = 2 * M;
		int BN = M;
		int LDB = BM;
		double *TAU = (double *)malloc(BN * sizeof(double));
		int LWORK = -1;
		double *WORK = (double *)malloc(1 * sizeof(double));
		int INFO;
		dgeqrf_(&BM, &BN, B, &LDB, TAU, WORK, &LWORK, &INFO);
		LWORK = *WORK;
		WORK = (double *)malloc(LWORK * sizeof(double));
		dgeqrf_(&BM, &BN, B, &LDB, TAU, WORK, &LWORK, &INFO);
		dorgqr_(&BM, &BN, &BN, B, &LDB, TAU, WORK, &LWORK, &INFO);

		double ALPHA = (a - b / c) / sqrt(c);
		double BETA = b / c;
		dgemm_("N", "T", &M, &M, &BN, &ALPHA, B, &LDB, B + M, &LDB, &BETA, U, &ldu);

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

		tempU = U;
		tempUprev = Uprev;
		for (int j = 0; j < M; j++)
		{
			for (int i = 0; i < M; i++)
			{
				tempUprev[i] = tempUprev[i] - tempU[i];
			}
			tempU += ldu;
			tempUprev += M;
		}
		delta = dlange_("F", &M, &M, Uprev, &M, WORK);
		//printf("delta = %e, iter = %d\n", delta, it);
	}
}
