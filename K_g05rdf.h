//This is the header file for K_g05rdf.cpp
#pragma once
#include <vector>

class K_g05rdf
{
public:
	//This function mimics the NAG function g05rdf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05rdf.html assuming mode=2 and ldc=ldx=m
	static void qq_g05rdf(const int* mode [[maybe_unused]] , const int* n, const int* m, const double* C, const int* ldc [[maybe_unused]] , const double* r [[maybe_unused]] ,
		const int* ldr [[maybe_unused]] , int* state, double* X, const int* ldx [[maybe_unused]] , const int* ifail [[maybe_unused]] );//if mode=2 then r irrelevant5/g05rdf.html
private:
	static std::vector<double> copulaNormal(int* state, std::vector<double>& meanX, std::vector<double>& covX, int m);
};

