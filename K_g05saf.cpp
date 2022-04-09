//The purpose of this file is to mimic the NAG function g05saf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05saf.html
//This is a function that generates pseudorandom numbers between 0 and 1. It does this using a state array in the computer's memory:
//Every time the computer requires a random number, it calculates that random number from the state array
//Immediately after calculating the random number, the state array is changed so that the next random number will be different.


#include <cmath>
#include "K_g01faf.h"
#include "K_g05saf.h"
#include <vector>


//Function that updates state array exactly once
//We assume that the Wichmann-Hill-II generator is used
void K_g05saf::updateState(int* state)
{
	//Pull values from state array
	int w = state[4]; // was 123456 when state created
	int x = state[5]; // was 102938 when state created
	int y = state[6]; // was 394875 when state created
	int z = state[7]; // was 765432 when state created
	const int a1 = state[8]; // 11600	a1
	const int wMult = state[9]; // 10379
	const int a2 = state[10]; // 47003
	const int xMult = state[11]; // 10479
	const int a3 = state[12]; // 23000
	const int yMult = state[13]; // 19423
	const int a4 = state[14]; // 33000
	const int zMult = state[15]; // 8123
	const int wMod = state[16]; // 185127
	const int xMod = state[17]; // 45688
	const int yMod = state[18]; // 93368
	const int zMod = state[19]; // 65075
	const int m1 = state[20]; // 2147483579
	const int m2 = state[21];	// 2147483543
	const int m3 = state[22]; // 2147483423
	const int m4 = state[23];	// 2147483123

	//Calculate next w,x,y,z
	w = a1 * (w % wMod) - wMult * (w / wMod);
	if (w < 0)
		w += m1; //same effect as w = (a1*w) mod m1 but safer
	x = a2 * (x % xMod) - xMult * (x / xMod);
	if (x < 0)
		x += m2; //same effect as x = (a2*x) mod m2 but safer
	y = a3 * (y % yMod) - yMult * (y / yMod);
	if (y < 0)
		y += m3; //same effect as y = (a3*y) mod m3 but safer
	z = a4 * (z % zMod) - zMult * (z / zMod);
	if (z < 0)
		z += m4; //same effect as z = (a4*z) mod m4 but safer

	//Update state
	state[4] = w;
	state[5] = x;
	state[6] = y;
	state[7] = z;
}

//Pulls values from state array and outputs pseudorandom number between 0 and 1
//Also updates state array so next random number will be different
//We assume that the Wichmann-Hill-II generator is used
double K_g05saf::randFromState(int* state)
{
	//Update state array
	updateState(state);

	//Pull values from state array
	const int w = state[4]; // was 123456 when state created
	const int x = state[5]; // was 102938 when state created
	const int y = state[6]; // was 394875 when state created
	const int z = state[7]; // was 765432 when state created
	const int m1 = state[20]; // 2147483579
	const int m2 = state[21]; // 2147483543
	const int m3 = state[22]; // 2147483423
	const int m4 = state[23]; // 2147483123

	//Return random number
	return std::fmod(w / static_cast<double>(m1) + x / static_cast<double>(m2) + y / static_cast<double>(m3) + z / static_cast<double>(m4), 1.0);
}

//Purpose of this function is to mimic NAG function g05saf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05saf.html
//Fill array x with n independent random numbers between 0 and 1
void K_g05saf::g05saf(const int* n, int* state, double* x, const int* ifail [[maybe_unused]] )
{
	//Fill x one element at a time
	for (int i = 0; i < *n; i++)
		x[i] = randFromState(state);
}


//Function that generates a random standard normal variable (mean 0, variance 1) using state array
double K_g05saf::rndStdNormal(int* state)
{
	double prob = randFromState(state);
	int ifail = 0;
	return K_g01faf::ndd_g01faf('L', &prob, &ifail);
}

//Function that generates a vector of n random independent standard normal variables (mean 0, variance 1) using state array
std::vector<double> K_g05saf::idtStdNormals(int*state, int n)
{
	//initialise sample
	std::vector<double> sample(n, 0.0);
	//fill sample
	for (int i = 0; i < n; i++)
		sample[i] = rndStdNormal(state);
	return sample;
}