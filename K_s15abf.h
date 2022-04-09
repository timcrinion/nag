#pragma once
#include <cmath>
#include "K_s15adf.h"
class K_s15abf
{
	//Normal Cumulative Distribution Function
	//Returns probability that the standard Normal distribution (mean 0, variance 1) is less than x
	//https://www.nag.co.uk/numeric/fl/nagdoc_fl26/html/s/s15abf.html
public:
	static double s15abf_cdf(const double* X, int* IFAIL) { return 0.5* K_s15adf::s15adf_erfcx(-(*X) / sqrt(2), IFAIL); };
};