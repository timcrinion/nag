#include "K_g01bjf.h"
#include <cmath> //so we can use pow(,)


//Binomial Probability Mass Function
//Returns probability of exactly k successes, given a total of n tries, and probability p of success per try
double K_g01bjf::BinomialPMF(int n, int k, double p)
{
	const double q = 1 - p;
	if (k == 0)
	{
		return std::pow(q, n);
	}
	if (k == n)
	{
		return std::pow(p, n);
	}
	double result = 1; //will multiply this result k times to get n(n-1)...(n-k+1) / k! p^k (1-p)^(n-k)
	const double constMult = std::pow(q, static_cast<double>(n - k) / static_cast<double>(k) );
	int topMult = n;
	int bottomMult = k;
	for (int i = 0; i < k; i++) // k times
	{
		result = result * topMult * p * constMult / bottomMult;
		topMult--;
		bottomMult--;
	}
	return result;
	//could return Choose(k, n) * std::pow(p, k) * std::pow(1 - p, n - k); but is inaccurate for large values of n
}

//Binomial Cumulative Density Function
//Returns probability of less than or equal to k successes, given a total of n tries, and probability p of success per try
//This function could be made far more efficient using arrays to store the nth row of Pascal's triangle, and powers of p and 1-p
double K_g01bjf::BinomialCDF(int n, int k, double p)
{
	//initialise these values to be updated in next loop
	double result = 0;
	for (int i = 0; i <= k; i++)
	{
		result += BinomialPMF(n, i, p);
	}
	return result;
}

/*
The purpose of this function is to mimic the NAG function G01BJF from https://www.nag.co.uk/numeric/fl/nagdoc_fl24/html/g01/g01bjf.html
All inputs and outputs are pointers of these variables:
In:
		N, total number of attempts
		P, probability of success
		K, an integer
Out:
		PLEK, probability that total number of successful attempts is less than or equal to K
		PGTK, probability that total number of successful attempts is more than K
		PEQK, probability that total number of successful attempts is exactly K		*/
void K_g01bjf::pd_g01bjf(const int* N, const double* P, const int* K, double* PLEK, double* PGTK, double* PEQK, const int* IFAIL)
{
		*PLEK = BinomialCDF(*N, *K, *P);
		*PEQK = BinomialPMF(*N, *K, *P);
		*PGTK = 1 - *PLEK;
		IFAIL = nullptr;
}