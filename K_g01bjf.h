#pragma once

class K_g01bjf
{
public:
	static void pd_g01bjf(const int* N, const double* P, const int* K, double* PLEK, double* PGTK, double* PEQK,
	                      const int* IFAIL);

private:
	static double BinomialCDF(int n, int k, double p);
	static double BinomialPMF(int n, int k, double p);
};

