#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include "../include/reflapack.h"
#include "../include/qdwheig.h"


#define NORMEST_TEST 1
#define QR_TEST 0
#define DSYQDWH_TEST 0 
#define DGECON_TEST 0


int main(int argc, char** argv)
{
	//Initialization
	if (argc != 4)
	{
		printf("Test function usage: ./Main M N lda\n");
		exit(-1);
	}

	printf("Welcome to branch v1.1!\n");
	int M = atoi(argv[1]);
	int N = atoi(argv[2]);
	int lda = atoi(argv[3]);
	double *A = (double *)malloc(N * lda * sizeof(double));
	srand(time(NULL));
	FILE *fp;
	
	//Testing function normest
#if NORMEST_TEST
	fp = fopen("TEST.m", "w");
	fprintf(fp, "function [norm0, norm1] = TEST\n");
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			if (i >= j)
			{
				A[i + j * lda] = (double) rand() / (double) RAND_MAX;
				fprintf(fp, "A(%d, %d) = %.12f;\n", i + 1, j + 1, A[i + j * lda]);
			}
			else
			{
				A[i + j * lda] = A[j + i * lda];
				fprintf(fp, "A(%d, %d) = %.12f;\n", i + 1, j + 1, A[i + j * lda]);
			}
		}
	}

	fprintf(fp, "norm0 = %.12f;\n", dsynormest(M, A, lda, 3e-1));
	fprintf(fp, "norm1 = normest(A, 3e-1);\n");
	fclose(fp);
#endif

	//QR factorization test
#if QR_TEST
	fp = fopen("TEST.m", "w");
	fprintf(fp, "function [Q0, Q1, error] = TEST\nA = zeros(%d, %d);\nQ1 = zeros(%d, %d);\n", M, N, M, N);	
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			A[i + j * lda] = (double) rand() / (double) RAND_MAX;
			fprintf(fp, "A(%d, %d) = %.12f;\n", i + 1, j + 1, A[i + j * lda]);
		}
	}
	int TAUsize = (M < N) ? M : N;
	double *TAU = (double *)malloc(TAUsize * sizeof(double));
	double *WORK = (double *)malloc(1 * sizeof(double));
	int LWORK = -1;
	int INFO;
	dgeqrf_(&M, &N, A, &lda, TAU, WORK, &LWORK, &INFO);
	LWORK = WORK[0];
	WORK = (double *)malloc(LWORK * sizeof(double));
	dgeqrf_(&M, &N, A, &lda, TAU, WORK, &LWORK, &INFO);
	//int ldc = M;
	//double *C = (double *)malloc(ldc * N * sizeof(double));
	//dormqr_("L", "N", &M, &N, &N, A, &lda, TAU, C, &ldc, WORK, &LWORK, &INFO);
	dorgqr_(&M, &N, &N, A, &lda, TAU, WORK, &LWORK, &INFO);
	
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			fprintf(fp, "Q1(%d, %d) = %.12f;\n", i + 1, j + 1, A[i + j * lda]);
		}
	}
	fprintf(fp, "[Q0, ~] = qr(A);\n");
	fprintf(fp, "error = norm(Q0 - Q1);\n");
	fclose(fp);
#endif

	//dsyqdwh test
#if DSYQDWH_TEST		
	fp = fopen("TEST.m", "w");
	fprintf(fp, "function [A, U0, U1, error] = TEST\n");
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			if (i >= j)
			{
				A[i + j * lda] = (double) rand() / (double) RAND_MAX;
				fprintf(fp, "A(%d, %d) = %.12f;\n", i + 1, j + 1, A[i + j * lda]);
			}
			else
			{
				A[i + j * lda] = A[j + i * lda];
				fprintf(fp, "A(%d, %d) = %.12f;\n", i + 1, j + 1, A[i + j * lda]);
			}
		}
	}
	double *U = (double *)malloc(N * N * sizeof(double));
	dsyqdwh(N, A, lda, U, N);
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < N; i++)
		{
			fprintf(fp, "U1(%d, %d) = %.12f;\n", i + 1, j + 1, U[i + j * N]);
		}
	}

	fprintf(fp, "U0 = qdwh(A, normest(A, 3e-1), 0.9 / condest(A));\n");
	fprintf(fp, "error = norm(U0 - U1, 'fro');\n");
	fprintf(fp, "end");
	fclose(fp);
	free(U);
#endif

	//dlange & dgecon test
#if DGECON_TEST 
	fp = fopen("TEST.m", "w");
	fprintf(fp, "function [A, error1, error2] = TEST\n");
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < N; i++)
		{
			if (i >= j)
			{
				A[i + j * lda] = (double) rand() / (double) RAND_MAX;
				fprintf(fp, "A(%d, %d) = %.12f;\n", i + 1, j + 1, A[i + j * lda]);
			}
			else
			{
				A[i + j * lda] = A[j + i * lda];
				fprintf(fp, "A(%d, %d) = %.12f;\n", i + 1, j + 1, A[i + j * lda]);
			}
		}
	}
			
	double *WORK = (double *)malloc(4 * N * sizeof(double));
	double *IWORK = (double *)malloc(N * sizeof(double));
	int INFO;
	double RCOND;
	double ANORM = dlange_("1", &N, &N, A, &lda, WORK);
	int *IPIV = (int *)malloc(N * sizeof(int));
	dgetrf_(&N, &N, A, &lda, IPIV, &INFO);
	dgecon_("1", &N, A, &lda, &ANORM, &RCOND, WORK, IWORK, &INFO);
	fprintf(fp, "norm1 = %.16f;\n", ANORM);
	fprintf(fp, "norm0 = norm(A, 1);\n");
	fprintf(fp, "error1 = abs(norm1 - norm0) / norm0;\n");
	fprintf(fp, "condest1 = %.16f;\n", 1.0 / RCOND);
	fprintf(fp, "condest0 = condest(A);\n");
	fprintf(fp, "error2 = abs(condest1 - condest0) / condest0;\n");
	fprintf(fp, "end");
	fclose(fp);
#endif

	
	
	free(A);
	return 0;
}
