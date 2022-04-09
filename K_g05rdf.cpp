//The purpose of this file is to mimic the NAG function g05rdf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05rdf.html

#include "K_g01eaf.h"
#include "K_g05rdf.h"
#include "K_g05rzf.h"
#include <vector>


//Function that returns a random vector (length m) with values between 0 and 1
//Inputs include the mean vector (m) and covariance matrix (mxm) of X
//Another input is state array to generate random numbers
std::vector<double> K_g05rdf::copulaNormal(int* state, std::vector<double>& meanX, std::vector<double>& covX, int m)
{
	std::vector<double> X = K_g05rzf::qq_g05rzf(state, meanX, covX, m);
	//convert to probabilities between 0 and 1
	double Xi;
	int ifail = 0;
	for (int i = 0; i < m; i++)
	{
		//This seems to be what happens by looking at the examples here
		//https://www.nag.co.uk/numeric/fl/nagdoc_fl26.0/pdf/g05/g05rdf.pdf
		//https://www.nag.com/numeric/fl/nagdoc_fl22/pdf/g05/g05rzf.pdf
		const double Cii = covX[i * m + i]; //variance of Xi
		Xi = (X[i] - meanX[i]) / sqrt(Cii);
		X[i] = K_g01eaf::g01eaf('L', &Xi, &ifail); //probability that standard normal < Xi
	}
	return X;
}


//This function mimics the NAG function g05rdf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05rdf.html assuming mode=2 and ldc=ldx=m
void K_g05rdf::qq_g05rdf(const int* mode [[maybe_unused]], const int* n, const int* m, const double* C, const int* ldc [[maybe_unused]], const double* r [[maybe_unused]],
	const int* ldr [[maybe_unused]], int* state, double* X, const int* ldx [[maybe_unused]], const int* ifail [[maybe_unused]])//if mode=2 then r irrelevant
{
	//create vector version of C
	std::vector<double> Cvec((*m) * (*m), 0.0);
	for (int i = 0; i < (*m) * (*m); i++)
		Cvec[i] = C[i];
	//Let the mean of X be the zero vector
	std::vector<double> meanX(*m, 0.0);
	//fill X (nxm)
	for (int i = 0; i < *n; i++)
	{
		std::vector<double> Xi;
		Xi = copulaNormal(state, meanX, Cvec, *m);
		for (int j = 0; j < *m; j++)
			X[i * (*m) + j] = Xi[j];
	}
}
