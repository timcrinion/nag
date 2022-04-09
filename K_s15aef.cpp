#include "K_s15adf.h"
#include "K_s15aef.h"

double K_s15aef::s15aef_erfx(double x, int* IFAIL)
{
	//Function acts differently depending on which of these cases x falls into:
	//Case 1:        x <  0
	//Case 2:   0 <= x <= 0.5
	//Case 3: 0.5 <  x
	if (x < 0.0) //Case 1
	{
		return -s15aef_erfx(-x, IFAIL); //set to -erf(-x)
	}
	if (x <= 0.5) //Case 2
	{
		return K_s15adf::erfBetweenZeroAndHalf(x);
	}
	//Case 3
	return 1.0 - s15aef_erfx(x, IFAIL); //1-efrc(x)
	// end if
}