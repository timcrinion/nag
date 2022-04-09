#pragma once
class K_g05kgf
{
private:

public:
	//Function that populates STATE (array of doubles) with non-repeatable random numbers between 0 and 1
	//mimics NAG function g05kgf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05kgf.html
	static void g05kgf(int* GENID, int* SUBID, int* STATE, int* LSTATE, int* IFAIL);
};
