#pragma once

#include <vector>

class K_e02cbf
{
	// member functions
public:

	//Mimics NAG function e02cbf. First scales the values in the array X, and also scales the single Y value. Then does
	//sumPolynomials.
	static void e02cbf(int MFIRST, int MLAST, int K, int L, const double* X, double XMIN, double XMAX, double Y,
		double YMIN, double YMAX, double* FF, const double* A, int NA, const double* WORK, int NWORK, int* IFAIL);

private:

	//function that fills vector T of length n+1 with chebyshev values of the first kind { T0(x), T1(x), ..., T_n(x) }
	//use recurrence relation from https://en.wikipedia.org/wiki/Chebyshev_polynomials#Definition
	static void chebyshevSeries(double x, int n, std::vector<double>& T);

	//for m = MFIRST, MFIRST+1, ..., MLAST, function sets FF[m] = sum of A[i,j] * Ti(X[m]) * Tj(Y) over i=0...K and j=0...L
	//Ti means Chebyshev polynomial of the first kind of degree i
	//There are MLAST-MFIRST+1 values for X and a single value for Y
	static void sumPolynomials(int MFIRST, int MLAST, int K, int L, std::vector<double>& X, double Y, double* FF,
		const double* A);
};
