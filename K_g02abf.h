#pragma once
//This is the header file for K_g02abf.cpp
#include <vector>

class K_g02abf
{
public:
	static void g02abf(double* G, const int* N, char OPT, const double* W, const double* ERRTOL, const int* MAXITS, const int* MAXIT, double* X, const int* ITER, int* IFAIL);

private:
	//The following three functions have been copied from https://nickhigham.wordpress.com/2013/02/13/the-nearest-correlation-matrix/
	static std::vector<double> nearcorr(std::vector<double>& A, int n, double tol, std::vector<double>& nWeights, bool useWeights, int maxits, int* iter);
	static void proj_spd(std::vector<double>& A, std::vector<double>& result, int n); //used by nearcorr()
	static void proj_unitdiag(std::vector<double>& A, int n); //used by nearcorr()
};