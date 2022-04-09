//This is the header file for K_g05saf.cpp
#pragma once
#include <vector>

class K_g05saf
{
private:
	static void updateState(int* state);
	static double randFromState(int* state);
public:
	//Purpose of this function is to mimic NAG function g05saf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05saf.html
	//Fill array x with n random numbers between 0 and 1 using state array
	//We assume that the Wichmann-Hill-II generator is used
	static void g05saf(const int* n, int* state, double* x, const int* ifail [[maybe_unused]] );
	//Function that generates a random standard normal variable (mean 0, variance 1) using state array
	static double rndStdNormal(int* state);
	//Function that generates a vector of n random independent standard normal variables (mean 0, variance 1) using state array
	static std::vector<double> idtStdNormals(int* state, int n);
};

