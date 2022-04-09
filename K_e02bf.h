//This is the header file for K_e02bf.cpp
#pragma once
#include <vector>

class K_e02bf
{
	/*
	Given a sequence of knots t[0] < t[1] < ... < t[n]
	define B-splines iteratively as follows:

	w(i,k,x) = (x-t[i]) / (t[i+k]-t[i])       that is, the line passing through (0,t[i]) and (1,t[i+k])

	B(i,1,x) = 1 if t[i] <= x < t[i+1]
			   0 otherwise

	for k>0:
	B(i,k+1,x) = w(i,k,x)B(i,k,x) + (1-w(i+1,k,x))B(i+1,k,x)
	*/

private:
	//define w as above, with knots t
	static double w(int i, int k, double x, std::vector<double>& t);

public:
	//define B as above, with knots t
	static double B(int i, int k, double x, std::vector<double>& t);

	static std::vector<double> leastSquares(int m, std::vector<double>& x, std::vector<double>& y, std::vector<double>& w, int q, std::vector<double>& t, bool& possible);

	static void e02baf (const int m, const int ncap7, const double x[], const double y[], const double w[], double lamda[],
		const double work1[], const double work2[], double c[], const double* ss, int* ifail);

	static std::vector<double> e02baf_inner(int m, std::vector<double>& x, std::vector<double>& y, std::vector<double>& w, int ncap7, std::vector<double>& t, int& ifail);

	//mimics nag function e02bbf https://www.nag.com/numeric/fl/nagdoc_fl25/pdf/e02/e02bbf.pdf
	static void e02bbf_vec(int ncap7, std::vector<double>& lambda, std::vector<double>& C, double x, double& s, int& ifail);

	//mimics nag function e02bbf https://www.nag.com/numeric/fl/nagdoc_fl25/pdf/e02/e02bbf.pdf
	static void e02bbf(int ncap7, double* lambda, double* C, double x, double& s, int& ifail);
};