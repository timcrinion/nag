//The purpose of this C++ file is to mimic the NAG function g05kgf from https://www.nag.co.uk/numeric/fl/nagdoc_latest/html/g05/g05kgf.html

#include <iostream>
#include <cmath>
#include "K_g05kff.h"
#include "K_g05kgf.h"
#include <ctime>
#include <chrono>

//Function that populates STATE (array of doubles) with non-repeatable random numbers between 0 and 1
void K_g05kgf::g05kgf(int* GENID, int* SUBID [[maybe_unused]], int* STATE, int* LSTATE, int* IFAIL)
{

	//want SIX independent nonrepeatable random numbers

	//first nonrepeatable random number
	int rand1 = static_cast<int>(time(0)); //time in seconds since 1st Jan 1970

	//second nonrepeatable random number
	auto now = std::chrono::high_resolution_clock::now();
	typedef std::chrono::high_resolution_clock::period period_t;
	auto dur = now.time_since_epoch();
	int rand2 = static_cast<int>(dur.count());

	//third nonrepeatable random number
	int rand3 = static_cast<int>(clock());

	//seed
	int seed[6] = { rand1, rand2, rand3, 1234, 1234, 1234 };  //look for a way of getting time of higher precision
	
	//length of seed
	int LSeed;

	if (*GENID == 1) //Basic NAG generator
	{
		LSeed = 1;
	}
	else if (*GENID == 2)
	{
		LSeed = 3;
	}
	else if (*GENID == 3)
	{
		//Mersenne twister generator
		LSeed = 1;
	}
	else if (*GENID == 4)
	{
		LSeed = 4;
	}
	else if (*GENID == 5)
	{
		//ACORN generator
		LSeed = 1;
	}
	else if (*GENID == 6)
	{
		LSeed = 6;
	}

	//update STATE
	K_g05kff::g05kff(*GENID, seed, LSeed, STATE, *LSTATE, IFAIL); //gareth draper wants int not double
}
