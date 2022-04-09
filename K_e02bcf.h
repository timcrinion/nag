#pragma once

#include <vector>

class K_e02bcf
{
	// member functions
public:

	//Calculate derivation of a cublic B-spline
	//use http://mat.fsv.cvut.cz/gcg/sbornik/prochazkova.pdf
	static void e02bcf(int NCAP7, const double* LAMBDA, const double* C, double X, int LEFT, double* S, int* IFAIL);

private:

	//function identical to K_e02baf, but derivative d is added on end
	//returns result of B(i,k,x,lambda) differentiated d times w.r.t x
	static double dB(int i, int k, double x, std::vector<double>& t, int d);
};
