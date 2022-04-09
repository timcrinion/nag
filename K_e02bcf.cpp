#include "K_e02bcf.h"
#include "K_e02bf.h"

//returns result of the spline B(i,k,x,lambda) differentiated d times w.r.t x
//inputs identical to those of K_e02baf::B, but derivative d is added on end
//copied from http://mat.fsv.cvut.cz/gcg/sbornik/prochazkova.pdf
double K_e02bcf::dB(int i, int k, double x, std::vector<double>& t, int d)
{
	if (x < t[i])
		return 0.0;
	if (x > t[i + 4])
		return 0.0;
	//from this point can assume t[i] < x < t[i+4]
	if (d == 0)
		return K_e02bf::B(i, k, x, t);
	//from this point can assume d>0
	if (k == 1)
		return 0.0;
	//from this point can assume k>1
	double div1 = t[i + k - 1] - t[i]; //want to divide by this number unless it is zero
	double mult1 = 0.0;
	if (div1 != 0.0)
		mult1 = 1.0 / div1;
	double div2 = t[i + k] - t[i + 1]; //want to divide by this number unless it is zero
	double mult2 = 0.0;
	if (div2 != 0.0)
		mult2 = 1.0 / div2;
	return (k - 1) *mult1 * dB(i, k - 1, x, t, d - 1) - (k - 1) * mult2 *dB(i + 1, k - 1, x, t, d - 1);
}

//Calculate derivation of a cublic B-spline
//use http://mat.fsv.cvut.cz/gcg/sbornik/prochazkova.pdf
void K_e02bcf::e02bcf(int NCAP7, const double* LAMBDA, const double* C, double X, int LEFT [[maybe_unused]], double* S,
	int* IFAIL)
{
	//spline s(x) = sum from i=0 to q-1 of C[i] * Ni(x)
	//   where Ni(x) defined as the cubic spline B(i, 4, x, LAMBDA), only positive when LAMBDA[i] <= x < LAMBDA[i+4]
	//LAMBDA has length NCAP7 = q + 4
	//C has length NCAP7
	//S has length 4

	if (NCAP7 < 8)
		*IFAIL = 1;
	if (LAMBDA[3] >= LAMBDA[NCAP7 - 4] || X < LAMBDA[3] || X>LAMBDA[NCAP7 - 4])
		*IFAIL = 2;

	int q = NCAP7 - 4;

	//lambda as a vector
	std::vector<double> lambdaVec(NCAP7);
	for (int i = 0; i < NCAP7; i++)
		lambdaVec[i] = LAMBDA[i];

	//fill S
	for (int d = 0; d < 4; d++)
	{
		S[d] = 0.0;
		for (int i = 0; i < q; i++)
			S[d] += C[i] * K_e02bcf::dB(i, 4, X, lambdaVec, d);
	}
}