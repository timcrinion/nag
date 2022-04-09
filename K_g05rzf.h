//This is the header file for K_g05rzf.cpp
#pragma once
#include <vector>

class K_g05rzf
{
private:
public:
	static std::vector<double> CholeskyDecomposition_vec(std::vector<double>& C, int n, bool* fail);
	//Multiplies matrix M by vector v
	static std::vector<double> MatrixByVector(std::vector<double> M, std::vector<double> v, int numRowsM, int numColsM);
	//Add two vectors of length n
	static std::vector<double> addVectors(std::vector<double>& vec1, std::vector<double>& vec2, int n);
	//mimic NAG function g05rzf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05rzf.html with n=1
	//meanX is mean of random vector X of length m
	//covX is mxm covariance matrix of X
	//indRand is a vector of m independent random standard normal variables N(0,1)^m
	static std::vector<double> qq_g05rzf(int* state, std::vector<double>& xMean, std::vector<double>& covX, int m);
	//Puts the elements from a double array into a vector
	static std::vector<double> doubleArrayToVector(const double* inArray, int length);
	//Check a n by n matrix M is symmetric. Returns true if symmetric, false otherwise.
	static bool isSymmetric(std::vector<double>& M, int n, double tolerance);
};
