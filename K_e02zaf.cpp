#include "K_e02zaf.h"

/* Sorts m coordinates into panels.
Panels defined by lambda and mu, coordinates by x and y
For example:

 y
 ^
 |
 |
  ----> x

	   3 |  7  | 11
mu[6] ---|-----|---
	   2 |  6  | 10
mu[5] ---|-----|---
	   1 |  5  | 9
mu[4] ---|-----|---
	   0 |  4  | 8
   lambda[4] lambda[5] (in this example px=10 and py=11)

lambda, length px. Only lambda[4], lambda[5] ... lambda[px-5] exist at start of function, nondecreasing order
mu, length py. Only mu[4], mu[5] ... mu[py-5] exist at start of function, nondecreasing order
x, length m
y, length m
*/
void K_e02zaf::e02zaf(int px, int py, const double* lambda, const double* mu, int m, const double* x, const double* y, int* point, int npoint, const int* adres [[maybe_unused]], int nadres, int* ifail)
{
	//Panels will be numbered 0 to numPanels-1
	const int numPanels = (px - 7) * (py - 7); //px-7 panels horizontally, py-7 verticallly
	/*point has length m+numPanels and is filled using a confusing mix of zero-based and one-based numbering.
	point[0]...point[m-1] are used differently from point[m]...point[m + numPanels - 1]
	point[m+i] = i1 means the first point in panel i is the (i1)th (defined by x[i1-1] y[i1-1])
	point[i1-1] = i2 means the next point in panel i is the (i2)th (defined by x[i2-1] y[i2-1])
	point[i2-1] = i3 means the next point in panel i is the (i3)th (defined by x[i3-1] y[i3-1])
	...
	point[in-1] = 0 means the (in)th point (x[in-1]y[in-1]) is the last in panel i
	*/

	*ifail = 0;

	//check lambda nondecreasing
	for (int i = 4; i < px - 5; i++)
	{
		if (lambda[i] > lambda[i + 1])
			*ifail = 1;
	}
	//check mu nondecreasing
	for (int i = 4; i < py - 5; i++)
	{
		if (mu[i] > mu[i + 1])
			*ifail = 1;
	}
	//check lengths of arrays are ok
	if ((px < 8) || (py < 8) || (m <= 0) || (nadres != numPanels) || (npoint < m + numPanels))
		*ifail = 2;

	if (*ifail != 0)
		return;

	//fill point with zeros
	for (int i = 0; i < m + numPanels; i++)
		point[i] = 0;

	//fill point
	int i = -1;
	int j = -1;
	int panel;
	for (int r = m - 1; r >= 0; r--) // x[m-1]y[m-1] to x[0]y[0]
	{
		//There are px-7 panels horizontally
		//First is          x < lambda[4]
		//then lambda[4] <= x < lambda[5]
		//...
		//then lambda[px-6] <= x < lambda[px-5]
		//then lambda[px-5] <= x
		//find i such that lambda[i-1] <= x < lambda[i]
		for (i = 4; i < px - 4; i++) //check from lambda[4] to lambda[px-5]
		{
			if (x[r] < lambda[i])
				break;
		}//now i has px-7 options, 4 <= i <= px-4

		//There are py-7 panels vertically
		//First is      y < mu[4]
		//then mu[4] <= y < mu[5]
		//...
		//then mu[py-6] <= y < mu[py-5]
		//then mu[py-5] <= y
		//find j such that mu[j-1] <= y < mu[j]
		for (j = 4; j < py - 4; j++) //check from mu[4] to mu[py-5]
		{
			if (y[r] < mu[j])
				break;
		}//now j has py-7 options, 4 <= j <= py-4

		//get panel number
		i -= 4; // 0 <= i <= px-8
		j -= 4; // 0 <= j <= py-8
		panel = i * (py - 7) + j;
		//update point
		point[r] = point[m + panel];
		point[m + panel] = r + 1;
	}
}
