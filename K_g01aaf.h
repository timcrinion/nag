#pragma once
//This is the header file for g01aaf.cpp

class K_g01aaf
{
public:

	//Function that returns the sum of an array arr = { arr[0]...arr[length-1] }
	static double sumArray(const double* arr, int length);

	//Function that returns the minimum value of an array arr = { arr[0]...arr[length-1] }
	static double minArray(const double* arr, int length);

	//Function that returns the maximum value of an array arr = { arr[0]...arr[length-1] }
	static double maxArray(const double* arr, int length);

	//Function that accepts two arrays A and B of length n
	//Modifies array AB. For example A={0,2,4} and B={1,10,100} would cause AB = {0,20,400}
	static void mult2Arrays(int n, const double* A, const double* B, double* AB);

	//Dot product of two arrays of equal length. Eg {1,10,1000} and {1,2,3} return 3021
	static double dotProductArrays(const double* array1, const double* array2, int n);


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
	static void getStats(int n, const double* X, int IWT, const double* WT, double* XMEAN, double* S2, double* S3, double* S4, double* XMIN, double* XMAX, double* WTSUM, int* IFAIL);

private:

};