#pragma once
class K_g05kff
{
private:
	static void WichmannHillSecondInitializer(const int* seedArray, int* stateArray);
public:
	//Function that populates STATE (array of doubles) with repeatable random numbers between 0 and 1
	//This function mimics https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05kff.html
	static void g05kff(const int & GENID, const int* SEED, const int & LSEED, int* STATE, const int & LSTATE, int* IFAIL);
};
