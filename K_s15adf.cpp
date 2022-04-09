#include "K_s15adf.h"
#include <cmath>
#include "KUtils.h"

//Partial version of error function erf(): Only works when 0 < x < 0.5
//Returns integral from 0 to x of 2*exp(-t^2)/sqrt(pi) dt
//https://www.nag.co.uk/numeric/fl/nagdoc_fl24/html/s/s15aef.html
double K_s15adf::erfBetweenZeroAndHalf(double x)
{
	//numbers used for the calculation
	const double p0 = 3.209377589138469472562 * pow(10, 3);
	const double p1 = 3.774852376853020208137 * pow(10, 2);
	const double p2 = 1.138641541510501556495 * pow(10, 2);
	const double p3 = 3.161123743870565596947 * pow(10, 0);
	const double p4 = 1.857777061846031526730 * pow(10, -1);
	const double q0 = 2.844236833439170622273 * pow(10, 3);
	const double q1 = 1.282616526077372275645 * pow(10, 3);
	const double q2 = 2.440246379344441733056 * pow(10, 2);
	const double q3 = 2.360129095234412093499 * pow(10, 1);
	const double q4 = 1.000000000000000000000 * pow(10, 0);
	//powers of x
	const double x2 = x * x;
	const double x4 = x2 * x2;
	const double x6 = x4 * x2;
	const double x8 = x6 * x2;
	const double top = p0 + p1 * x2 + p2 * x4 + p3 * x6 + p4 * x8;
	const double bottom = q0 + q1 * x2 + q2 * x4 + q3 * x6 + q4 * x8;
	return x * top / bottom;
}

//Complimentary error function erfc(). Equal to 1-erf(x). Equivalently, equal to integral from x to infinity of 2*exp(-t^2)/sqrt(pi) dt
//https://www.nag.co.uk/numeric/fl/nagdoc_fl24/html/s/s15adf.html
double K_s15adf::s15adf_erfcx(double x, int* IFAIL)
{
	const double split1 = 0.46865;
	const double split2 = 4;
	//Function acts differently depending on which of these cases x falls into:
	//Case 1:          x < 0
	//Case 2:      0 <=x <= split1
	//Case 3: split1 < x <= split2
	//Case 4: split2 < x
	if (x < 0) //Case 1
	{
		return 2 - s15adf_erfcx(-x, IFAIL);  //2-erfc(-x)
	}
	if (x <= split1) //Case 2
	{
		return 1 - erfBetweenZeroAndHalf(x);
	}
	if (x <= split2) //Case 3
	{
		//numbers used for the calculation
		const double p0 = 1.23033935479799725272 * pow(10, 3);
		const double p1 = 2.05107837782607146532 * pow(10, 3);
		const double p2 = 1.71204761263407058314 * pow(10, 3);
		const double p3 = 8.81952221241769090411 * pow(10, 2);
		const double p4 = 2.98635138197400131132 * pow(10, 2);
		const double p5 = 6.61191906371416294775 * pow(10, 1);
		const double p6 = 8.88314979438837594118 * pow(10, 0);
		const double p7 = 5.64188496988670089180 * pow(10, -1);
		const double p8 = 2.15311535474403846343 * pow(10, -8);
		const double q0 = 1.23033935480374942043 * pow(10, 3);
		const double q1 = 3.43936767414372163696 * pow(10, 3);
		const double q2 = 4.36261909014324715820 * pow(10, 3);
		const double q3 = 3.29079923573345962678 * pow(10, 3);
		const double q4 = 1.62138957456669018874 * pow(10, 3);
		const double q5 = 5.37181101862009857509 * pow(10, 2);
		const double q6 = 1.17693950891312499305 * pow(10, 2);
		const double q7 = 1.57449261107098347253 * pow(10, 1);
		const double q8 = 1.00000000000000000000 * pow(10, 0);
		//powers of x
		const double x2 = x * x;
		const double x3 = x2 * x;
		const double x4 = x3 * x;
		const double x5 = x4 * x;
		const double x6 = x5 * x;
		const double x7 = x6 * x;
		const double x8 = x7 * x;
		const double top = p0 + p1 * x + p2 * x2 + p3 * x3 + p4 * x4 + p5 * x5 + p6 * x6 + p7 * x7 + p8 * x8;
		const double bottom = q0 + q1 * x + q2 * x2 + q3 * x3 + q4 * x4 + q5 * x5 + q6 * x6 + q7 * x7 + q8 * x8;
		return exp(-x2) * top / bottom;
	}
		//Case 4
	//numbers used for the calculation
	const double p0 = -6.58749161529837803157 * pow(10, -4);
	const double p1 = -1.60837851487422766278 * pow(10, -2);
	const double p2 = -1.25781726111229246204 * pow(10, -1);
	const double p3 = -3.60344899949804439429 * pow(10, -1);
	const double p4 = -3.05326634961232344035 * pow(10, -1);
	const double p5 = -1.63153871373020978498 * pow(10, -2);
	const double q0 = 2.33520497626869185443 * pow(10, -3);
	const double q1 = 6.05183413124413191178 * pow(10, -2);
	const double q2 = 5.27905102951428412248 * pow(10, -1);
	const double q3 = 1.87295284992346047209 * pow(10, 0);
	const double q4 = 2.56852019228982242072 * pow(10, 0);
	const double q5 = 1.00000000000000000000 * pow(10, 0);
	//negative powers of x
	const double x1 = 1 / x;
	const double x2 = x1 * x1;
	const double x4 = x2 * x2;
	const double x6 = x4 * x2;
	const double x8 = x6 * x2;
	const double x10 = x8 * x2;
	const double top = p0 + p1 * x2 + p2 * x4 + p3 * x6 + p4 * x8 + p5 * x10;
	const double bottom = q0 + q1 * x2 + q2 * x4 + q3 * x6 + q4 * x8 + q5 * x10;
	return exp(-x * x) * x1 * (pow(KUtils::Pi(), -0.5) + x2 * top / bottom); //check that KUtils::Pi() is accurate enough
	// end if
}


