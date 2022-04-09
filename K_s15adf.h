#pragma once

class K_s15adf
{
public:
	//Complimentary error function s15adf. Equal to 1-erf(x). Equivalently, equal to integral from x to infinity of 2*exp(-t^2)/sqrt(pi) dt
	//https://www.nag.co.uk/numeric/fl/nagdoc_fl24/html/s/s15adf.html
	static double s15adf_erfcx(double x, int* IFAIL);

	//Partial version of error function erf(): Only works when 0 < x < 0.5
	//Returns integral from 0 to x of 2*exp(-t^2)/sqrt(pi) dt
	//https://www.nag.co.uk/numeric/fl/nagdoc_fl24/html/s/s15aef.html
	static double erfBetweenZeroAndHalf(double x);
};

