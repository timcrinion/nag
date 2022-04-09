//This is the header file for K_e02df.cpp
#pragma once
#include <vector>

class K_e02daf
{
private:

public:

	static std::vector<double> leastSquares(int m, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, std::vector<double>& w,
		int px, std::vector<double>& lambda, int py, std::vector<double>& mu, bool& possible, double* DL, double eps);
		static double s(double x, double y, int px, std::vector<double>& lambda, int py, std::vector<double>& mu, std::vector<double>& c);
		static void daf(int m, int px, int py, double* x, double* y, double* f, double* w, double* lambda, double* mu, const int* point [[maybe_unused]] , int npoint [[maybe_unused]] , double* dl,
			double* c, int nc [[maybe_unused]] , const double* ws [[maybe_unused]] , int nws [[maybe_unused]] , double eps, const double* sigma [[maybe_unused]] , const int* rank [[maybe_unused]] , int* ifail);

	static void RREF(std::vector<double>& M, int rowCount, int columnCount);
	static void reduce(std::vector<double>& M, int numRows, int numCols, double eps);

};