//The purpose of this C++ file is to mimic the NAG function g05sqf

#include "K_g05saf.h"
#include "K_g05sqf.h"
#include "KUtils.h"


//fill array x with n random numbers between a and b
void K_g05sqf::sqf(int* n, double* a, double* b, int* state, double* x, int* ifail)
{
	*ifail = 0;

	if (n < 0)
		* ifail = 1;

	if (b < a)
		* ifail = 3;

	K_g05saf::g05saf(n, state, x, ifail);
	double diff = *b - *a;
	for (int i = 0; i < *n; i++)
		x[i] = *a + diff * x[i];

}
