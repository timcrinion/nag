//Purpose of C++ file is to mimic the NAG function G01AAF https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g01/g01aaf.html
#include <iostream>
#include <cmath> //so we can use pow(,)
#include "K_g01aaf.h"

/**********************************************************************************************************************************************************************************************/

//Function that returns the sum of an array arr = { arr[0]...arr[length-1] }
double K_g01aaf::sumArray(const double* arr, int length)
{
	double sum = 0;
	for (int i = 0; i < length; i++)
	{
		sum += arr[i];
	}//end for
	return sum;
} //end function

/**********************************************************************************************************************************************************************************************/

//Function that returns the minimum value of an array arr = { arr[0]...arr[length-1] }
double K_g01aaf::minArray(const double* arr, int length)
{
	double minValue = arr[0];
	for (int i = 1; i < length; i++)
	{
		if (arr[i] < minValue)
		{
			minValue = arr[i];
		} //end if
	} //end for
	return minValue;
} //end function

/**********************************************************************************************************************************************************************************************/

//Function that returns the maximum value of an array arr = { arr[0]...arr[length-1] }
double K_g01aaf::maxArray(const double* arr, int length)
{
	double maxValue = arr[0];
	for (int i = 1; i < length; i++)
	{
		if (arr[i] > maxValue)
		{
			maxValue = arr[i];
		} //end if
	} //end for
	return maxValue;
} //end function

/**********************************************************************************************************************************************************************************************/

//Function that accepts two arrays A and B of length n
//Modifies array AB. For example A={0,2,4} and B={1,10,100} will cause AB={0,20,400}
void K_g01aaf::mult2Arrays(int n, const double* A, const double* B, double* AB)
{
	//AB will be updated as it's external to function.
	for (int i = 0; i < n; i++)
		AB[i] = A[i] * B[i];
} // end function (helped by https://stackoverflow.com/questions/6422993/returning-an-array-from-a-function)

/**********************************************************************************************************************************************************************************************/

//Dot product of two arrays of equal length. Eg {1,10,100} and {1,2,3} return 321
double K_g01aaf::dotProductArrays(const double* array1, const double* array2, int n)
{
	double result = 0;
	for (int i = 0; i < n; i++)
		result += array1[i] * array2[i];
	return result;
}

/**********************************************************************************************************************************************************************************************/

/*
Purpose of function is to mimic G01AAF https://www.nag.co.uk/numeric/fl/nagdoc_fl25/html/g01/g01aaf.html
All inputs and outputs are pointers of these variables
In:
		n (length of X and WT)
		X = array of length n containing sample values
		WT = array of length n containing weights of the values in X
		IWT = 0 if WT not provided, 1 if provided
		IFAIL = -1, 0, or 1. Not sure what this means
Out:
		XMEAN = mean of X
		S2 = standard deviation of X
		S3 = coefficient of skewness of X
		S4 = coefficient of kurtosis of X
		XMIN = smallest value of X
		XMAX = largest value of X
		WTSUM = sum of WT
		IFAIL =	1 if N<1
				2 if number of valid cases, m, is 1. I don't understand this
				3 if either the number of valid cases is 0, or a weight is negative
				-99 if an unexpected error has been triggered by this function
				-399 if your licence key does not work
				-999 if dynamic memory allocation failed
*/
void K_g01aaf::getStats(int n, const double* X, int IWT, const double* WT, double* XMEAN, double* S2, double* S3, double* S4, double* XMIN, double* XMAX, double* WTSUM, int* IFAIL)
{
	*IFAIL = 0;
	if (n < 1) //if N<1 return fail
	{
		*IFAIL = 1;
		return;
	}
	//Set weights to either WT or {1,1,1...} depending on IWT
	if (IWT == 1) //user provides weights inside WT
	{
		//sum the weights
		*WTSUM = sumArray(WT, n);
		//get mean of weighted X
		*XMEAN = dotProductArrays(WT, X, n) / *WTSUM;
		//Redefine X as X-XMEAN and declare following sums of products so we can calculate S2, S3 and S4
		double sumWX2 = 0, sumWX3 = 0, sumWX4 = 0, sumW2 = 0;
		for (int i = 0; i < n; i++)
		{
			const double x = X[i] - *XMEAN;
			const double wx2 = WT[i] * x * x;
			const double wx3 = wx2 * x;
			const double wx4 = wx3 * x;
			sumWX2 += wx2;
			sumWX3 += wx3;
			sumWX4 += wx4;
			sumW2 += WT[i] * WT[i];
			if (WT[i] < 0)
			{
				*IFAIL = 3;
			}
		}
		double denominator = *WTSUM - sumW2 / *WTSUM;
		if (abs(denominator) <= 0.000000000001)//to avoid division by zero
		{
			*IFAIL = -99;
			return;
		}
		const double varX = sumWX2 / denominator; //danger of division by zero
		//Standard deviation
		*S2 = sqrt(varX);
		//Coefficient of skewness
		denominator *= (*S2) * varX;
		*S3 = sumWX3 / denominator;
		//Coefficient of kurtosis
		denominator *= (*S2);
		*S4 = sumWX4 / denominator;
	}
	else //user provides no weights inside WT, use all ones instead
	{
		//sum the weights
		*WTSUM = n;
		//get mean of weighted X
		*XMEAN = sumArray(X, n) / *WTSUM;
		//Redefine X as X-XMEAN and declare following sums of powers of x so we can calculate S2, S3 and S4
		double sumX2 = 0, sumX3 = 0, sumX4 = 0;
		for (int i = 0; i < n; i++)
		{
			const double x = X[i] - *XMEAN;
			const double x2 = x * x;
			const double x3 = x2 * x;
			const double x4 = x3 * x;
			sumX2 += x2;
			sumX3 += x3;
			sumX4 += x4;
		}
		if (n == 1)//to avoid division by zero
		{
			*IFAIL = 2;
			return;
		}
		double denominator = n - 1;
		const double varX = sumX2 / denominator; //danger of division by zero
		//Standard deviation
		*S2 = sqrt(varX);
		//Coefficient of skewness
		denominator *= (*S2) * varX;
		*S3 = sumX3 / denominator;
		//Coefficient of kurtosis
		denominator *= (*S2);
		*S4 = sumX4 / denominator;
	}
	//end if

	//get min and max of X
	*XMIN = minArray(X, n);
	*XMAX = maxArray(X, n);
}// end function
