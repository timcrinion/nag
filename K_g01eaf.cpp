#include "K_g01eaf.h"
#include <iostream>
#include <cmath> //so we can use pow(,)
#include "K_s15aef.h"
#include "K_s15abf.h"


//Probability that the standard normal distribution (mean 0, variance 1) is less than abs(x) and more than -abs(x)
//https://www.nag.co.uk/numeric/fl/nagdoc_fl24/html/g01/g01eaf.html
double K_g01eaf::NormalSandwich(double x, int* IFAIL)
{
	return K_s15aef::s15aef_erfx(abs(x) / sqrt(2), IFAIL);
} // end function

/*
The purpose of this function is to mimic the NAG function G01EAF https://www.nag.co.uk/numeric/fl/nagdoc_fl24/html/g01/g01eaf.html
*/
double K_g01eaf::g01eaf(char TAIL, const double* X, int* IFAIL)
{
	if (TAIL == 'L') //lower tail probability P(Normal<x) is returned
	{
		return K_s15abf::s15abf_cdf(X, IFAIL);
	}
	if (TAIL == 'U') //upper tail probability P(x<Normal) is returned
	{
		return 1.0 - K_s15abf::s15abf_cdf(X, IFAIL);
	}
	if (TAIL == 'S') //two tail significance probability P(N>|x|) + P(N<-|x|)
	{
		return 1.0 - NormalSandwich(*X, IFAIL);
	}
	if (TAIL == 'C') //two tail confidence probability P(N<|x|) - P(N<-|x|)
	{
		return NormalSandwich(*X, IFAIL);
	}
	*IFAIL = 1;
	return 0;
}
