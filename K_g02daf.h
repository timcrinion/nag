//This is the header file for K_g02daf.cpp
#pragma once
#include <vector>

class K_g02daf
{
private:

public:

	static void multLinReg(int n, std::vector<double> X, int p, std::vector<double> Y, std::vector<double>& B, double* se_array, double tol, int* ifail);
	static void g02daf(bool mean, int n, const double* X_array, int m, const int* isx, int ip, double* Y_array, double* B_array, double* se_array, double tol, int* ifail);
	static double getVariance(std::vector<double> Y, std::vector<double> X, std::vector<double> B, int n, int p, int rank, int* ifail);

};