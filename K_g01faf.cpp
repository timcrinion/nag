#include "K_g01faf.h"
#include <algorithm>

//Mimic NAG function g01faf https://www.nag.co.uk/numeric/fl/nagdoc_fl24/html/g01/g01faf.html
//return x such that N (standard normal), P and x satisfy one of the cases below
double K_g01faf::ndd_g01faf(char TAIL, const double* P, int* IFAIL)
{
	*IFAIL = 0;
	if ((*P <= 0) || (*P >= 1))
		*IFAIL = 2;
	else if (TAIL == 'L') //P = Prob(N<x)
		return NormalDeviate(*P);
	else if (TAIL == 'U') //P = Prob(x<N)
	{
		const double q = 1 - *P; // q = Prob(N<x)
		return NormalDeviate(q);
	}
	else if (TAIL == 'S') //P = Prob(|x|<N) + Prob(N<-|x|)
	{
		const double q = 1 - 0.5 * (*P); //q = Prob(N<|x|)
		return NormalDeviate(q); //always returns positive value
	}
	else if (TAIL == 'C') //P = Prob(N<|x|) - Prob(N<-|x|)
	{
		const double q = 0.5 * (*P + 1); //q = Prob(N<|x|)
		return NormalDeviate(q); //always returns positive value
	}
	else
		*IFAIL = 1;
	return 0;
}

//Wichura's function PPND16, described in http://csg.sph.umich.edu/abecasis/gas_power_calculator/algorithm-as-241-the-percentage-points-of-the-normal-distribution.pdf
//Returns x such that p is the probability that Normal(mean 0, variance 1) < x
//should be accurate to sixteen figures for 10^-316 < p < 1-10^-316
double K_g01faf::NormalDeviate(double p)
{
	const double q = p - 0.5;
	double r;

	const double split1 = 0.425;
	const double split2 = 5;
	const double const2 = 1.6;
	if (abs(q) <= split1) // p within 0.425 of 0.5
	{
		r = split1 * split1 - q * q;
		return q * polynomialA(r) / polynomialB(r);
	}
		// p futher than 0.425 from 0.5
	r = sqrt(-log(std::min(p, 1 - p)));
	if (r <= split2)    // p further than exp(-25) from 0 and 1
	{
		r = r - const2;
		return sign(q)* polynomialC(r) / polynomialD(r);
	}
	//p within exp(-25) of 0 or 1
	r = r - split2;
	return sign(q)* polynomialE(r) / polynomialF(r);
	// end if
	// end if
}

int K_g01faf::sign(double input)
{
	int output = 1;
	if (input < 0) {
		output = -1;
	}//end if
	return output;
}