//The purpose of this C++ file is to imitate the NAG function g01edf() from https://www.nag.com/numeric/fl/nagdoc_latest/html/g01/g01edf.html
//and g01ebf() from https://www.nag.com/numeric/fl/nagdoc_latest/html/g01/g01ebf.html

#include "K_g01edf.h"
#include "K_s15abf.h"
#include "KUtils.h"
#include <vector>

//lower tail probability of the F or variance-ratio distribution, with df1 and df2 degrees of freedom
//	0 < df1 < 100000
//	0 < df2 < 100000
//	0 <= f
//returns Prob(F <= f)
double K_g01edf::G01edf(char tail, double f, double df1, double df2, int* ifail)
{
	//set ifail to 0 unless input values are unacceptable
	*ifail = 0;
	if (f < 0.0)
	{
		*ifail = 2;
		return -1.0;
	}
	if (df1 <= 0.0 || df2 <= 0.0)
	{
		*ifail = 3;
		return -1.0;
	}

	//use beta distribution
	const double x = df1 * f / (df1 * f + df2); //0 <= x < 1
	const double a = df1 * 0.5; // 0 < a < 50000
	const double b = df2 * 0.5; // 0 < b < 50000
	double result = K_g01edf::BetaCumulativeProb(x, a, b);
	if (tail == 'L') // return prob(F<f)
		return result;
	else if (tail == 'U') //return prob(F>f)
		return 1.0 - result;
	else
	{
		*ifail = 1;
		return -1.0;
	}

}


//mimics NAG function g01ebf
double K_g01edf::Ebf(char tail, double t, double df, int* ifail)
{
	*ifail = 0;
	if (df < 1.0)
	{
		*ifail = 2;
		return -1.0;
	}

	//calculate probability
	double prob = -1.0;
	if (tail == 'U') //user wants upper tail probability P(T>t)
		prob = 1.0 - K_g01edf::StudentProb(t, df);
	else if (tail == 'S') //user wants two tail (significance level) probability P(T<-|t|) + P(T>|t|)
		prob = K_g01edf::StudentProb(-abs(t), df) + 1.0 - K_g01edf::StudentProb(abs(t), df);
	else if (tail == 'C') //user wants two tail (confidence interval) probability P(T<|t|) - P(T<-|t|)
		prob = K_g01edf::StudentProb(abs(t), df) - K_g01edf::StudentProb(-abs(t), df);
	else if (tail == 'L') //user wants lower tail probability P(T<t)
		prob = K_g01edf::StudentProb(t, df);
	else
		*ifail = 1;
	return prob;
}

//returns cumulative probability P(T<t:df) of Student's t-distribution with df degrees of freedom
//must have df >= 1
double K_g01edf::StudentProb(double t, double df)
{
	if (t < 0.0)
		return 1.0 - K_g01edf::StudentProb(-t, df);
	if (df < 20.0)
		return 1.0 - 0.5 * K_g01edf::BetaCumulativeProb(df / (df + t * t), 0.5 * df, 0.5);
	else // df>=20 and t>=0
	{//use Cornish-Fisher expansion from https://dl.acm.org/doi/10.1145/355598.362775
		const double a = df - 0.5;
		const double b = 48.0 * a * a;
		const double z2 = a * log(1.0 + t * t / df);
		const double z = sqrt(z2);
		const double z3 = z * z2;
		const double z4 = z2 * z2;
		const double z5 = z2 * z3;
		const double z7 = z2 * z5;
		const double z9 = z2 * z7;
		const double z11 = z2 * z9;
		const double thirdDivisor = 10.0 * b * (b + 0.8 * z4 + 100.0);
		double X = z +
			(z3 + 3.0 * z) / b -
			(4.0 * z7 + 33.0 * z5 + 240.0 * z3 + 855.0 * z) / thirdDivisor +
			(64.0 * z11 + 788.0 * z9 + 9801.0 * z7 + 89775.0 * z5 + 543375.0 * z3 + 1788885.0 * z) / (210.0 * b * b * b);

		int ifail;
		return K_s15abf::s15abf_cdf(&X, &ifail);
	}
}

/*
Cumulative probability distribution for beta function.
Identical to regular incomplete beta function.
Conditions:
	0 <=x < 1
	0 < a < 50000
	0 < b < 50000
*/
double K_g01edf::BetaCumulativeProb(double x, double a, double b)
{//https://en.wikipedia.org/wiki/Beta_distribution#Cumulative_distribution_function
	if (x >= 1.0)
		return 1.0;
	if (x <= 0.0)
		return 0.0;
	return K_g01edf::RegIncBeta(x, a, b);
}



//Beta function
//returns x! * y! / (x+y)!
double K_g01edf::Beta(double x, double y)
{
	return K_g01edf::LanczosGamma(x)* K_g01edf::LanczosGamma(y) / K_g01edf::LanczosGamma(x + y);
}


//Lanczos approximation for Gamma function of a real number z
//see https://en.wikipedia.org/wiki/Lanczos_approximation#Simple_implementation
//should be correct to 15 decimal places
//Gamma function = integral from zero to infinity of x^(z-1) e^(-x) dx
double K_g01edf::LanczosGamma(double z)
{
	//from cmath import sin, sqrt, pi, exp

	std::vector<double> coeff = { 676.5203681218851
		,-1259.1392167224028
		,771.32342877765313
		,-176.61502916214059
		,12.507343278686905
		,-0.13857109526572012
		,0.0000099843695780195716
		,0.00000015056327351493116 };
	int numCoeff = 8;

	if (z < 0.5)
		return KUtils::Pi() / (sin(KUtils::Pi() * z) * K_g01edf::LanczosGamma(1.0 - z));  //Reflection formula
	else
	{
		double z1 = z - 1.0;
		double x = 0.99999999999980993;
		for (int i = 0; i < numCoeff; i++)
			x += coeff[i] / (z1 + static_cast<double>(i) + 1.0);
		double t = z1 + static_cast<double>(numCoeff) - 0.5;
		return sqrt(2 * KUtils::Pi()) * pow(t, z1 + 0.5) * exp(-t) * x;
	}
}//another approximation is sqrt(2 pi x) * (x/e)^x (x sinh(1/x))^(x/2) * e^(7 / (324 (35x^2 + 33) x^3)) from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5840229/



/*Regularised incomplete beta function
Conditions:
	0 <=x < 1
	0 < a < 50000
	0 < b < 50000
function requires memory O(b-a) and time O(b) */
double K_g01edf::RegIncBeta(double x, double a, double b)
{
	if (x > 0.5)
		return 1.0 - K_g01edf::RegIncBeta(1.0 - x, b, a); //using https://dlmf.nist.gov/8.17#i
	if (b >= a && b >= 2)
		return K_g01edf::RegIncBeta(x, a + 1.0, b - 1.0) + pow(x, a) * pow((1.0 - x), b - 1.0) / (a * K_g01edf::Beta(a, b)); //using https://dlmf.nist.gov/8.17#iv
	//at this point, 0 <= x <= 0.5 and either 0 < b < a < 50000 or b < 2
	//incomplete beta function
	double Bxab = (pow(x, a) * pow(1.0 - x, b) / a) * K_g01edf::Hypergeometric_for_regIncBeta(a + b, a + 1, x); // x^a * (1-x)^b / a * F(a+b, 1; a+1; x) using https://dlmf.nist.gov/8.17#ii
	//regular incomplete beta function
	return Bxab / K_g01edf::Beta(a, b);
}



/*
Special hypergoemetric function

This function is intended only to be called from regIncBeta()
since that guarantees 0 <= x <= 0.5 and 0 < p < 2q+1

Function returns F(p,1;q;x) =
1 +
px / q +
p(p+1) x^2 / (q(q+1)) +
p(p+1)(p+2) x^3 / (q(q+1)(q+2)) + etc

where
	1 < q
	0 < p < 2q+1
	0 <= x <= 0.5

Given above conditions, function accurate to at least 8 decimal places
*/
double K_g01edf::Hypergeometric_for_regIncBeta(double p, double q, double x)
{
	double result = 1.0;
	double mult;
	double increment = 1.0;
	int consecutiveDrops = 0;
	for (int i = 0; i < 999; i++)
	{//each time, multiply increment and add it to result
		mult = (p + i) * x / (q + i);
		increment *= mult;
		result += increment;
		//if we have had mult<0.6 fifty times in a row, break
		if (mult < 0.6)
			consecutiveDrops++;
		else
			consecutiveDrops = 0;
		if (consecutiveDrops == 50)
			break;
	}
	return result;
}



//integral from 0 to x of t^(s-1)*e^(-t) dt
double K_g01edf::LowerIncompleteGamma(double s, double x)
{
	double result = 0.0;

	//use https://en.wikipedia.org/wiki/Incomplete_gamma_function#Evaluation_formulae
	//loop below adds increments to result. Want number of terms to be such that increment decreases rapidly
	double increment = 1.0 / s;
	int numIncrements = static_cast<int>(2.0 * x - s) + 20; //at this point the remaining terms will sum to less than 0.00000001 times the maximum term
	for (int i = 0; i < numIncrements; i++)
	{
		result += increment;
		increment *= x / (s + static_cast<double>(i) + 1.0);
	}

	result *= pow(x, s) * exp(-x);
	return result;
}