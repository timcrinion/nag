#pragma once
#include <cmath>

class K_g01faf
{
public:
	static double ndd_g01faf(char TAIL, const double* P, int* IFAIL);

private:
	//Coefficients used by polynomials:
//near to 0.5
	static double A0() { return 3.3871328727963666080* pow(10, 0); }
	static double A1() { return 1.3314166789178437745* pow(10, 2); }
	static double A2() { return 1.9715909503065514427* pow(10, 3); }
	static double A3() { return 1.3731693765509461125* pow(10, 4); }
	static double A4() { return 4.5921953931549871457* pow(10, 4); }
	static double A5() { return 6.7265770927008700853* pow(10, 4); }
	static double A6() { return 3.3430575583588128105* pow(10, 4); }
	static double A7() { return 2.5090809287301226727* pow(10, 3); }
	static double B0() { return 1; }
	static double B1() { return 4.2313330701600911252* pow(10, 1); }
	static double B2() { return 6.8718700749205790830* pow(10, 2); }
	static double B3() { return 5.3941960214247511077* pow(10, 3); }
	static double B4() { return 2.1213794301586595867* pow(10, 4); }
	static double B5() { return 3.9307895800092710610* pow(10, 4); }
	static double B6() { return 2.8729085735721942674* pow(10, 4); }
	static double B7() { return 5.2264952788528545610* pow(10, 3); }
	//far from 0.5, and not extremely close to 0 or 1
	static double C0() { return 1.42343711074968357734* pow(10, 0); }
	static double C1() { return 4.63033784615654529590* pow(10, 0); }
	static double C2() { return 5.76949722146069140550* pow(10, 0); }
	static double C3() { return 3.64784832476320460504* pow(10, 0); }
	static double C4() { return 1.27045825245236838258* pow(10, 0); }
	static double C5() { return 2.41780725177450611770* pow(10, -1); }
	static double C6() { return 2.27238449892691845833* pow(10, -2); }
	static double C7() { return 7.74545014278341407640* pow(10, -4); }
	static double D0() { return 1; }
	static double D1() { return 2.05319162663775882187* pow(10, 0); }
	static double D2() { return 1.67638483018380384940* pow(10, 0); }
	static double D3() { return 6.89767334985100004550* pow(10, -1); }
	static double D4() { return 1.48103976427480074590* pow(10, -1); }
	static double D5() { return 1.51986665636164571966* pow(10, -2); }
	static double D6() { return 5.47593808499534494600* pow(10, -4); }
	static double D7() { return 1.05075007164441684324* pow(10, -9); }
	//extremely close to 0 or 1
	static double E0() { return 6.65790464350110377720* pow(10, 0); }
	static double E1() { return 5.46378491116411436990* pow(10, 0); }
	static double E2() { return 1.78482653991729133580* pow(10, 0); }
	static double E3() { return 2.96560571828504891230* pow(10, -1); }
	static double E4() { return 2.65321895265761230930* pow(10, -2); }
	static double E5() { return 1.24266094738807843860* pow(10, -3); }
	static double E6() { return 2.71155556874348757815* pow(10, -5); }
	static double E7() { return 2.01033439929228813265* pow(10, -7); }
	static double F0() { return 1; }
	static double F1() { return 5.99832206555887937690* pow(10, -1); }
	static double F2() { return 1.36929880922735805310* pow(10, -1); }
	static double F3() { return 1.48753612908506148525* pow(10, -2); }
	static double F4() { return 7.86869131145613259100* pow(10, -4); }
	static double F5() { return 1.84631831751005468180* pow(10, -5); }
	static double F6() { return 1.42151175831644588870* pow(10, -7); }
	static double F7() { return 2.04426310338993978564* pow(10, -15); }

	//polynomial functions
	static double polynomialA(double x) { return ((((((A7() * x + A6()) * x + A5()) * x + A4()) * x + A3()) * x + A2()) * x + A1()) * x + A0(); }
	static double polynomialB(double x) { return ((((((B7() * x + B6()) * x + B5()) * x + B4()) * x + B3()) * x + B2()) * x + B1()) * x + B0(); }
	static double polynomialC(double x) { return ((((((C7() * x + C6()) * x + C5()) * x + C4()) * x + C3()) * x + C2()) * x + C1()) * x + C0(); }
	static double polynomialD(double x) { return ((((((D7() * x + D6()) * x + D5()) * x + D4()) * x + D3()) * x + D2()) * x + D1()) * x + D0(); }
	static double polynomialE(double x) { return ((((((E7() * x + E6()) * x + E5()) * x + E4()) * x + E3()) * x + E2()) * x + E1()) * x + E0(); }
	static double polynomialF(double x) { return ((((((F7() * x + F6()) * x + F5()) * x + F4()) * x + F3()) * x + F2()) * x + F1()) * x + F0(); }

	static double NormalDeviate(double p);
	static int sign(double input);
};
