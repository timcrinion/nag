//The purpose of this C++ file is to mimic the NAG function g05kff from https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05kff.html
//This function creates a state array that is stored in the computer's memory, and is edited by g05saf whenever it needs to calculate a random number.

#include "K_g05kff.h"

/****************************************************************************************************************************************************/

//Initializer for Wichmann-Hill II Generator from https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05intro.html#iigenerator
void K_g05kff::WichmannHillSecondInitializer(const int* seedArray, int* stateArray)
{
	//Fill array with same values as original NAG function so that SAM can't tell the difference. Discovered these values using debug points and studying STATE in Immediate Window:
	stateArray[0] = 29; //length of array?
	stateArray[1] = 4826; //?
	stateArray[2] = 4; //number of seeds, or representing Wichmann-Hill-II?
	stateArray[3] = 0; //?
	stateArray[4] = seedArray[0]; //seed for w, usually 123456, or current state?
	stateArray[5] = seedArray[1]; //seed for x, usually 102938, or current state?
	stateArray[6] = seedArray[2]; //seed for y, usually 394875, or current state?
	stateArray[7] = seedArray[3]; //seed for z, usually 765432, or current state?
	stateArray[8] = 11600; //w multiplier
	stateArray[9] = 10379; //w small multiplier
	stateArray[10] = 47003; //x multiplier
	stateArray[11] = 10479; //x samll multiplier
	stateArray[12] = 23000; //y multiplier
	stateArray[13] = 19423; //y small multiplier
	stateArray[14] = 33000; //z multiplier
	stateArray[15] = 8123; //z small multiplier
	stateArray[16] = 185127; //w small mod
	stateArray[17] = 45688; //x small mod
	stateArray[18] = 93368; //y small mod
	stateArray[19] = 65075; //z small mod
	stateArray[20] = 2147483579; //w mod
	stateArray[21] = 2147483543; //x mod
	stateArray[22] = 2147483423; //y mod
	stateArray[23] = 2147483123; //z mod
	stateArray[24] = 0; //?
	stateArray[25] = 4826; //?
	stateArray[26] = 1; //?
	stateArray[27] = 1; //?
	stateArray[28] = 4826; //?
}

/****************************************************************************************************************************************************/

//Function that populates STATE array with initial values so that it can generate pseudorandom numbers later on
//This function mimics https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05kff.html
void K_g05kff::g05kff(const int & GENID, const int* SEED, const int & LSEED, int* STATE, const int & LSTATE, int* IFAIL)
{
	if (GENID == 1) //Basic NAG generator
	{
		//SAM does not use initializer for Basic NAG generator, so this is not done here. https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05intro.html#nagbasicgenerator
	}
	else if (GENID == 2) //Wichmann-Hill-I generator
	{
		//SAM does not use initializer for Wichmann-Hill-I generator, so this is not done here. https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05intro.html#igenerator
	}
	else if (GENID == 3) //Mersenne Twister generator
	{
		//SAM does not use initializer for Mersenne Twister generator, so this is not done here. https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05intro.html#merseenetwister
	}
	else if (GENID == 4) //Wichmann-Hill-II generator
	{
		WichmannHillSecondInitializer(SEED, STATE); //https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05intro.html#iigenerator
		if (LSEED != 4 || LSTATE != 29)
		{
			*IFAIL = 1;
		}
	}
	else if (GENID == 5) //ACORN generator
	{
		//SAM does not use initializer for ACORN generator, so this is not done here. https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05intro.html#acorn
	}
	else if (GENID == 6) //L'Ecuyer generator
	{
		//SAM does not use initializer for L'Ecuyer generator, so this is not done here. https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g05/g05intro.html#lecuyer
	}
	else
	{
		*IFAIL = 1;
		//should also set *LSTATE to 1 but the value being passed to the function is constant
	}
}
