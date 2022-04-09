//This is the header file for K_g01edf.cpp
#pragma once
#include <vector>

class K_g01edf
{
private:
	static double Beta(double x, double y);
	static double BetaCumulativeProb(double x, double a, double b);
	static double Hypergeometric_for_regIncBeta(double p, double q, double x);
	static double LanczosGamma(double z);
	static double LowerIncompleteGamma(double s, double x);
	static double RegIncBeta(double x, double a, double b);
	static double StudentProb(double t, double df);

public:

	static double Ebf(char tail, double t, double df, int* ifail);
	static double G01edf(char tail, double f, double df1, double df2, int* ifail);

};