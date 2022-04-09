//the purpose of this C++ file is to replicate the NAG function E02CBF
#include <vector>

#include "K_e02cbf.h"

//function that fills vector T of length n+1 with chebyshev values of the first kind { T0(x), T1(x), ..., T_n(x) }
//use recurrence relation from https://en.wikipedia.org/wiki/Chebyshev_polynomials#Definition
void K_e02cbf::chebyshevSeries(double x, int n, std::vector<double>& T)
{
	if (n > 0)
		T[0] = 1.0;
	if (n > 1)
		T[1] = x;
	for (int i = 2; i <= n; i++)
		T[i] = 2.0 * x * T[i - 1] - T[i - 2];
}


//for m = MFIRST, MFIRST+1, ..., MLAST, function sets FF[m] = sum of A[i,j] * Ti(X[m]) * Tj(Y) over i=0...K and j=0...L
//Ti means Chebyshev polynomial of the first kind of degree i
//There are MLAST-MFIRST+1 values for X and a single value for Y
void K_e02cbf::sumPolynomials(int MFIRST, int MLAST, int K, int L, std::vector<double>& X, double Y, double* FF,
	const double* A)
{
	//assume that A[i,j] is stored in A[i * (L+1) + j]
	//the function description is correct that indices start at 0, so A contains (0 to K) times L+1 (0 to L) values

	//create vector of length L+1 containing chebyshev polynomials T0(Y)...Tn(Y)
	std::vector<double> TY(L+1);
	K_e02cbf::chebyshevSeries(Y, L, TY);

	//create vector of length K+1 such that ATY[i] = sum of A[i,j] * Tj(Y) over j=0...L
	std::vector<double> ATY(K + 1, 0.0);
	double increment;
	for (int i = 0; i <= K; i++)
	{
		for (int j = 0; j <= L; j++)
		{
			increment = A[i*(L + 1) + j] * TY[j]; // increment by A[i, j] * Tj(Y)
			if (j == 0)
				increment *= 0.5;
			ATY[i] += increment; // increment by A[i, j] * Tj(Y)
		}
	}

	//prepare vector of length K+1 that will be filled with chebyshev polynomials of X[m] in each loop below
	std::vector<double> TXm(K + 1);
	//fill FF
	for (int m = MFIRST; m <= MLAST; m++)
	{
		FF[m] = 0;
		//fill TXm with chebyshev polynomials of X[m]
		K_e02cbf::chebyshevSeries(X[m], K, TXm);
		for (int i = 0; i <= K; i++)
		{
			increment = TXm[i] * ATY[i];
			if (i == 0)
				increment *= 0.5;
			FF[m] += increment;
		}
	}
	return;
}

//first scales the values in the array X, and also scales the single Y value. Then does sumPolynomials.
void K_e02cbf::e02cbf(int MFIRST, int MLAST, int K, int L, const double* X, double XMIN, double XMAX, double Y, double YMIN,
	double YMAX, double* FF, const double* A, int NA, const double* /*WORK*/, int NWORK, int* IFAIL)
{
	*IFAIL = 0;

	//make *IFAIL nonzero if something wrong
	if (MFIRST > MLAST || K < 0 || L < 0 || NA < (K + 1) * (L + 1) || NWORK < K + 1)
	{
		*IFAIL = 1;
		return; //no point in doing function
	}
	if (YMIN >= YMAX || Y<YMIN || Y>YMAX)
		*IFAIL = 2;
	for (int m = MFIRST; m <= MLAST; m++)
	{
		if (X[m]<XMIN || X[m]>XMAX)
			*IFAIL = 3;
	}
	if (XMIN >= XMAX)
		*IFAIL = 3;

	//scale x and y
	double newY = (2.0 * Y - (YMAX + YMIN)) / (YMAX - YMIN);
	int numX = MLAST - MFIRST + 1;
	std::vector<double> newX(numX);
	for (int m = 0; m < numX; m++)
		newX[m] = (2 * X[m] - (XMAX + XMIN)) / (XMAX - XMIN);

	//call sumPolynomials
	K_e02cbf::sumPolynomials(0, numX - 1, K, L, newX, newY, FF, A);
}