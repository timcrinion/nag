//This is the header file for K_borsdorf.cpp
#pragma once

class K_borsdorf
{
private:


public:

	//following functions copied from section A.2 of http://eprints.maths.manchester.ac.uk/1085/1/thesis.pdf

	static std::vector<double> diagEnlarge(const std::vector<double>& v, int n); //v length n, output length n^2
	static std::vector<double> diagReduce(const std::vector<double>& M, int n); //M length n^2, output length n
	static std::vector<double> cornewton(std::vector<double>& G, int n, double error_tol, int diaones, int maxitmeth, int maxit, bool* ifail); //G length n^2, output is nearest correlation matrix to G, length n^2
	static std::vector<double> computeX(std::vector<double>& P, std::vector<double>& lambda, int n, int diaones); //P length n^2, lambda length n, output length n^2
	static void eigdecomp(std::vector<double>& P, std::vector<double>& lambda, const std::vector<double>& G, std::vector<double>& y, int n); //y and lambda length n, P and G length n^2
	static void gradient(double* f, std::vector<double>& userFy, std::vector<double>& y, std::vector<double>& lambda, std::vector<double>& P, std::vector<double>& b0, int n);//inputs length n except P of length n^2
	static std::vector<double> omega_matrix(std::vector<double>& lambda, int n); //input length n, output length n^2
	static double maxTwoDoubles(double a, double b);
	static double minTwoDoubles(double a, double b);
	static double minVec(std::vector<double> vec, int vecLen);
	static std::vector<double> precond_matrix(std::vector<double>& Omega, std::vector<double>& P, int n); //inputs length n^2, output length n
	static std::vector<double> Jacobian_matrix(std::vector<double>& x, std::vector<double>& Omega, std::vector<double>& P, int n); //Omega and P length n^2, x and output length n

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Following two functions are copied and pasted from minres.cpp
	//Same as originals except A*x is replaced by Jacobean_matrix(x)
	//Also, new variables bool flag and int iterk introduced
	//Therefore K_borsdorf::gmres() will modify x such that b=Jacobean_matrix(x,Omega,P,n)
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	static int gmres(std::vector<double>& Omega, std::vector<double>& P, std::vector<double>& b, std::vector<double>& x, int n, int max_iterations, double threshold);
	static void arnoldi(std::vector<double>& H, std::vector<double>& Omega, std::vector<double>& P, std::vector<double>& Q, int k, int n, int m);

};
