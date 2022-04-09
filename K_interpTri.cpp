//This C++ file imitates the NAG functions e01saf and e01sbf
//Details in https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770025881.pdf
#include <algorithm>
#include "K_Cartesian2D.h"
#include "K_Thiessen.h"
#include "K_interpTri.h"
#include "K_minres.h"
#include "K_qshep2.h"
#include <vector>


/*a weighted least squares method
Inputs are three vectors w,x,y,z of length n and double f.
Outputs are coefficients a,b,c,d,e such that z approximated by ax^2 + bxy + cy^2 +dx + ey + f
To find these, define residue R = sum over i of wi * (zi - (axi^2 + bxiyi + cyi^2 +dxi + eyi + f))^2 for i=0 to n-1
Want to minimise R by adjusting coefficients. Therefore set its derivatives to zero:
dR/da = 2 sum_i( wi (zi - (axi^2 + bxiyi + cyi^2 +dxi + eyi + f)) * -xi^2 ) = 0
dR/db = 2 sum_i( wi (zi - (axi^2 + bxiyi + cyi^2 +dxi + eyi + f)) * -xiyi ) = 0
dR/dc = 2 sum_i( wi (zi - (axi^2 + bxiyi + cyi^2 +dxi + eyi + f)) * -yi^2 ) = 0
dR/dd = 2 sum_i( wi (zi - (axi^2 + bxiyi + cyi^2 +dxi + eyi + f)) * -xi   ) = 0
dR/de = 2 sum_i( wi (zi - (axi^2 + bxiyi + cyi^2 +dxi + eyi + f)) * -yi   ) = 0
dR/df = 2 sum_i( wi (zi - (axi^2 + bxiyi + cyi^2 +dxi + eyi + f)) * -1    ) = 0
This is equivalent to solving the following matrix equation for a,b,c,d,e:

|sum(w x^4) sum(w x^3 y)   sum(w x^2 y^2) sum(w x^3)   sum(w x^2 y)|  |a|   |sum(w x^2 (z-f))|
|           sum(w x^2 y^2) sum(w x y^3)   sum(w x^2 y) sum(w x y^2)|  |b|   |sum(w x y (z-f))|
|                          sum(w y^4)     sum(w x y^2) sum(w y^3)  |  |c|   |sum(w y^2 (z-f))|
|        (symmetric)                      sum(w x^2)   sum(w x y)  |  |d| = |sum(w x (z-f))  |
|                                                      sum(w y^2)  |  |e|   |sum(w y (z-f))  |

Therefore a,b,c,d,e can be found by multiplying both sides by the inverse of the 5-by-5 symmetric matrix on the LHS
Returns true if function successful, false otherwise
*/
bool K_interpTri::WeightedLeastSquares(int n, std::vector<double>& W, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z, double& a, double& b, double& c, double& d, double& e, double f)
{
	//Initiate sums:
	//order 3
	double sumWXX = 0.0; //sum over i of W[i]*X[i]*X[i]
	double sumWXY = 0.0;
	double sumWYY = 0.0;
	double sumWXZ = 0.0;
	double sumWYZ = 0.0;
	//order 4
	double sumWXXX = 0.0;
	double sumWXXY = 0.0;
	double sumWXYY = 0.0;
	double sumWYYY = 0.0;
	double sumWXXZ = 0.0;
	double sumWYYZ = 0.0;
	double sumWXYZ = 0.0;
	//order 5
	double sumWXXXX = 0.0;
	double sumWXXXY = 0.0;
	double sumWXXYY = 0.0;
	double sumWXYYY = 0.0;
	double sumWYYYY = 0.0;

	//increment sums:
	for (int i = 0; i < n; i++)
	{
		//order 1
		double w = W[i];
		double x = X[i];
		double y = Y[i];
		double z = Z[i] - f;
		//order 3
		double wxx = w * x * x;
		double wyy = w * y * y;
		double wxy = w * x * y;
		double wxz = w * x * z;
		double wyz = w * y * z;
		sumWXX += wxx;
		sumWYY += wyy;
		sumWXY += wxy;
		sumWXZ += wxz;
		sumWYZ += wyz;
		//order 4
		double wxxx = wxx * x;
		double wxxy = wxx * y;
		double wyyy = wyy * y;
		sumWXXX += wxxx;
		sumWXXY += wxxy;
		sumWXYY += wxy * y;
		sumWYYY += wyyy;
		sumWXXZ += wxx * z;
		sumWYYZ += wyy * z;
		sumWXYZ += wxy * z;
		//order 5
		sumWXXXX += wxxx * x;
		sumWXXXY += wxxx * y;
		sumWXXYY += wxxy * y;
		sumWXYYY += wyyy * x;
		sumWYYYY += wyyy * y;
	}

	//Set M as symmetric 6 by 6 matrix on the LHS, and v as vector on the RHS
	std::vector<double> M = {
	sumWXXXX, sumWXXXY, sumWXXYY, sumWXXX, sumWXXY,
	sumWXXXY, sumWXXYY, sumWXYYY, sumWXXY, sumWXYY,
	sumWXXYY, sumWXYYY, sumWYYYY, sumWXYY, sumWYYY,
	sumWXXX,  sumWXXY,  sumWXYY,  sumWXX,  sumWXY,
	sumWXXY,  sumWXYY,  sumWYYY,  sumWXY,  sumWYY }; //length 25
	std::vector<double> v = { sumWXXZ, sumWXYZ, sumWYYZ, sumWXZ, sumWYZ }; //length 5

	//want abcde = M^-1 * v
	std::vector<double> abcde = v;
	bool ok = K_minres::gmres(M, v, abcde, 5, 5, 0.000000000001); //or this if convenient
	a = abcde[0];
	b = abcde[1];
	c = abcde[2];
	d = abcde[3];
	e = abcde[4];
	return ok;
}



//following function fills DX and DY with estimated partial differentials df/dx and df/dy
//returns true if everything is fine, false if at least one problem with least squares method
bool K_interpTri::GetPartials(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& DX, std::vector<double>& DY,
	std::vector<int>& randTri, std::vector<int>& TriAdj, std::vector<int>& TriPts)
{
//set up grid for finding nearest points
	auto gridDim = static_cast<int>(sqrt(n / 3));
	if (gridDim < 1)
		gridDim = 1;
	std::vector<int> grid(gridDim * gridDim, -1);
	std::vector<int> next(n, -1);
	double dx, dy, xMin, yMin;
	K_Cartesian2D::SetupGrid(n, gridDim, x, y, grid, next, dx, dy, xMin, yMin);

	//declare vectors of length n now to prevent loop below having order O(n)
	std::vector<int> pNeighbours(n, -1);
	std::vector<int> nInt(n, -1);
	std::vector<bool> nBool(n, false);

	//for each node p
	bool all_ok = true;
	for (int p = 0; p < n; p++)
	{
		//set k to be 2+number of nodes that affect partials (+2 because answer will include p, and we also want to set radius as distance from nearest node not affecting partials)
		int k = std::min(n, 10); //10
		bool leastSquaresFailed = true;
		//repeat next loop until least squares method succeeds and DX[p] and DY[p] filled
		while (leastSquaresFailed)
		{
			//fill k_near_p with k nodes nearest to p. This will include p itself, plus k-1 others
			std::vector<int> kNearP = K_Cartesian2D::Nearest_k(x[p], y[p], k, n, gridDim, x, y, grid, next, dx, dy, xMin, yMin);

			//fill pNeighbours with all nodes 1 or 2 triangles away from p, including p itself
			int len_pNeighbours = K_interpTri::Neighbours2(n, p, nBool, pNeighbours, nInt, randTri, TriAdj, TriPts, 2); //nBool and nInt only used inside this function

			//Merge these into nearP, the complete collection of nearby points we wish to use (it will include p) of length numNear
			int numNear = 0;
			std::vector<int> nearP = K_interpTri::MergeIntVectors(kNearP, k, pNeighbours, len_pNeighbours, numNear);

			//sort by distance from p, this will cause nearP[0]=p
			K_Cartesian2D::OrderPointsByDistance(x[p], y[p], x, y, nearP, numNear);

			//exclude p from nearP
			for (int i = 0; i < numNear - 1; i++)
				nearP[i] = nearP[i + 1];
			numNear--; //because we excluded p

			//set radius as a bit further than the furthest in nearP
			double radius = sqrt(K_Cartesian2D::DistSquared(x[p], y[p], x[nearP[numNear]], y[nearP[numNear]]));
			numNear--; //because we want to ignore furthest point from now on

			//normally would only take k nearest, and ignore neighbours, but do this because some tina tests have hundreds of points close together in lines that are far apart

			//Apply least squares method to the following points to calculate formula that predicts z from x and y, using weighting w
			std::vector<double> xNear(numNear);
			std::vector<double> yNear(numNear);
			std::vector<double> zNear(numNear);
			std::vector<double> wNear(numNear);
			for (int i = 0; i < numNear; i++)//for each nearby point
			{
				xNear[i] = x[nearP[i]] - x[p]; //relative to p
				yNear[i] = y[nearP[i]] - y[p]; //relative to p
				zNear[i] = z[nearP[i]];
				wNear[i] = sqrt(xNear[i] * xNear[i] + yNear[i] * yNear[i]); //distance from p
				wNear[i] = (radius - wNear[i]) / (radius * wNear[i]);
			}
			//get a,b,c,d,e such that z approximated by a(x-xp)^2 + b(x-xp)(y-yp) + c(y-yp)^2 +d(x-xp) + e(y-yp) + zp
			double a = 0.0;
			double b = 0.0;
			double c = 0.0;
			double d = 0.0;
			double e = 0.0;
			leastSquaresFailed = !K_interpTri::WeightedLeastSquares(numNear, wNear, xNear, yNear, zNear, a, b, c, d, e, z[p]);
			//according to this approximation, dz/dx=d and dz/dy=e at the point p
			DX[p] = d;
			DY[p] = e;
			//prepare for next loop
			k += static_cast<int>(sqrt(n) / 3) + 1;
			if (k >= n - 3)
			{
				all_ok = false;
				break;
			}
		}
	}
	return all_ok;
}

//function that outputs n+1 mod 3
int add1(int n)
{
	if (n == 2)
		return 0;
	return n + 1;
}

//function that outputs n-1 mod 3
int sub1(int n)
{
	if (n == 0)
		return 2;
	return n - 1;
}

//interpolate f at point P=(xx,yy) within a triangle. Algorithm copied from:
//https://ntrs.nasa.gov/search.jsp?R=19760016832
//https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19760016832.pdf
//"C1-Compatible Interpolation Over A Triangle" by C. L. Lawson 1976
double K_interpTri::Interpolate(double xx, double yy, std::vector<double> & x, std::vector<double> & y, std::vector<double> & f, std::vector<double> & fx, std::vector<double> & fy, bool linear)
{
	//phase 1

	std::vector<double> u(3);
	std::vector<double> v(3);
	std::vector<double> L2(3);
	for (int i = 0; i < 3; i++)
	{
		u[i] = x[sub1(i)] - x[add1(i)];
		v[i] = y[sub1(i)] - y[add1(i)];
		L2[i] = u[i] * u[i] + v[i] * v[i];
	}
	double delta = u[0] * v[1] - u[1] * v[0]; //assume this is not zero since that would mean 3 points colinear

	//barycentric coordinates from chapter 4
	std::vector<double> r(3);
	for (int i = 0; i < 3; i++)
		r[i] = (u[i] * (yy - y[add1(i)]) - v[i] * (xx - x[add1(i)])) / delta; //eq 22

	if (linear)
	{
		double w = 0.0;
		for (int i = 0; i < 3; i++)
			w += f[i] * r[i];
		return w;
	}

	//phi from chapter 7
	std::vector<double> phi(3);
	for (int i = 0; i < 3; i++)
		phi[i] = r[add1(i)] * r[sub1(i)]; //eq 40

	//phase 2

	//m from chapter 8
	int m = 0;
	if (r[1] < r[m])
		m = 1;
	if (r[2] < r[m])
		m = 2;

	//rho from chapter 8
	std::vector<double> rho(3);
	double a = 0.5 * r[m] * r[m];
	double b = r[m] / 3;
	rho[m] = r[m] * (phi[m] + 5 * a / 3) - a;
	rho[add1(m)] = a * (r[sub1(m)] - b);
	rho[sub1(m)] = a * (r[add1(m)] - b);

	//phase 3 version 1

	std::vector<double> h(3); //h tilde
	std::vector<double> k(3); //k tilde
	std::vector<double> g(3); //g tilde
	for (int i = 0; i < 3; i++)
	{
		h[i] = u[i] * fx[add1(i)] + v[i] * fy[add1(i)];
		k[i] = u[i] * fx[sub1(i)] + v[i] * fy[sub1(i)];
		g[i] = (r[add1(i)] - r[sub1(i)]) * phi[i] + 3 * (L2[add1(i)] - L2[sub1(i)]) * rho[i] / L2[i] - rho[add1(i)] + rho[sub1(i)];
	}
	double w = 0.0;
	for (int i = 0; i < 3; i++)
		w += f[i] * r[i] + 0.5 * (h[i] - k[i]) * phi[i] + (0.5 * (h[i] + k[i]) - f[sub1(i)] + f[add1(i)]) * g[i];
	return w;
}

//Imitates NAG function e01saf
//inputs:
//   nodes represented by vectors x,y,f of length n (must be no duplicate coordinates (x,y))
//outputs:
//   triangulated grid on x,y in TriAdj,TriPts,randTri
//   estimated lists of df/dx and df/dy at each node (if linear=false. User may want derivatives depending on how they plan to use e01sbf)
//   grid (returned) of length gridDim^2 and values dx,dy,xMin,xMax corressponding to grid
void K_interpTri::E01saf(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, std::vector<int> & TriAdj, std::vector<int> & TriPts,
	std::vector<double> & dfdx, std::vector<double> & dfdy, std::vector<int>& grid, int gridDim, double& dx, double& dy, double& xMin, double& yMin, int& ifail, bool linear)
{
	//must have at least 3 nodes
	ifail = 0;
	if (n < 3)
	{
		ifail = 1;
		return;
	}

	//build triangular grid, store in TriAdj and TriPts (and randTri)
	for (int i = 0; i < 6 * n; i++)
	{
		TriAdj[i] = -1;
		TriPts[i] = -1;
	}
	std::vector<int> randTri(n, -1);
	K_Thiessen::TRIGRD(n, x, y, randTri, TriAdj, TriPts);

	//get partial derivatives if desired by user
	if (!linear)
	{
		if (K_interpTri::GetPartials(n, x, y, f, dfdx, dfdy, randTri, TriAdj, TriPts))
			ifail = 0;
		else
			ifail = 2;
	}

	//construct grid such that grid[i,j] = some triangle t inside cell [i,j]
	//however, if mid-point of grid[i,j] outside triangular grid, it will be set to a triangle facing t
	std::vector<int> next(n, -1);
	//set up grid such that a point (x,y) is in grid[i,j] = grid[i*gridDim + j] if
	//	i <= (x - xMin) / dx < i + 1 and
	//	j <= (y - yMin) / dy < j + 1
	K_Cartesian2D::SetupGrid(n, gridDim, x, y, grid, next, dx, dy, xMin, yMin);
	//fill each cell [i,j] in grid
	int t; //triangle
	int tTopRow = 0; //triangle last time i was zero
	int duffInt;
	for (int i = 0; i < gridDim; i++)
	{
		t = tTopRow;
		for (int j = 0; j < gridDim; j++)
		{
			//replace t with a triangle facing or containing mid-point of cell [i,j]
			K_Thiessen::TRFIND(duffInt, t, xMin + (i + 0.5) * dx, yMin + (j + 0.5) * dy, x, y, TriAdj, TriPts);
			//store t in grid[i,j]
			grid[i * gridDim + j] = t;
			//if this is the top row, set tTopRow
			if (j == 0)
				tTopRow = t;
		}
	}
}

//estimates f at given point P=(Px,Py)
//Input triangular grid defined by TriAdj and TriPts
//Input partial derivatives given by dfdx and dfdy
//ifail set to 0 if a triangle contains P, 3 otherwise
//user can set linear=true if they want f on plane of triangle (partial derivatives ignored), false if they want to use partial derivatives. True results in "jaggedy" surface, false in smooth surface
double K_interpTri::E01sbf(double Px, double Py, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, std::vector<int>& TriAdj, std::vector<int>& TriPts,
	std::vector<double>& dfdx, std::vector<double>& dfdy, std::vector<int>& grid, int gridDim, double dx, double dy, double xMin, double yMin, int& ifail, bool linear)
{
	//Locate triangle triP containing P

	//find cell [i,j] of grid in which p lies
	auto i = static_cast<int>((Px - xMin) / dx);
	auto j = static_cast<int>((Py - yMin) / dy);

	//find a triangle in/near that cell
	i = std::max(i, 0);
	i = std::min(i, gridDim - 1);
	j = std::max(j, 0);
	j = std::min(j, gridDim - 1);
	int triP = grid[i * gridDim + j];

	//triP is currently a triangle near P
	//Set it to either the triangle containing P (if P inside triangular grid) or a triangle with edge facing P (if P outside triangular grid)
	int edge;
	bool P_in_triP = K_Thiessen::TRFIND(edge, triP, Px, Py, x, y, TriAdj, TriPts);

	//set ifail to 3 if P not in any triangle, 0 otherwise
	ifail = 3;
	if (P_in_triP)
		ifail = 0;

	//get 3 points of triangle a,b,c
	int a = TriPts[3 * triP];
	int b = TriPts[3 * triP + 1];
	int c = TriPts[3 * triP + 2];
	//interpolate inside triangle
	std::vector<double> xTri = { x[a], x[b], x[c] };
	std::vector<double> yTri = { y[a], y[b], y[c] };
	std::vector<double> fTri = { f[a], f[b], f[c] };
	std::vector<double> dxTri = { dfdx[a], dfdx[b], dfdx[c] };
	std::vector<double> dyTri = { dfdy[a], dfdy[b], dfdy[c] };
	double Pf = K_interpTri::Interpolate(Px, Py, xTri, yTri, fTri, dxTri, dyTri, linear);
	return Pf;
}



//sets input vector to all -1s
//assume on input that list consists only of nonnegatives until a certain point, after which it is only -1s. Eg {2,4,0,12,-1,-1,-1}
//This function efficient, <O(n), because it only bothers changing first few entries, not all n of them
void resetListOfNeighbours(int n, std::vector<int>& list)
{
	for (int i = 0; i < n; i++)
	{
		if (list[i] == -1)
			return;
		else
			list[i] = -1;
	}
}




//returns number of neighbours of p, and inserts neighbours into vector neighboursOut of length n
//Tp is any triangle with node p
//When function called, neighboursOut should never contain a -1 followed by a non-negative number
int K_interpTri::Neighbours(int p, int Tp, std::vector<int> & TriAdj, std::vector<int> & TriPts, std::vector<int> & neighboursOut, int n)
{
	//reset neighboursOut to all -1s
	resetListOfNeighbours(n, neighboursOut);

	int numNeighbours = 0;
	int i, nextEdge;

	//Starting with Tp, record nodes clockwise around p until either we hit boundary or arrive back at Tp
	int t = Tp;
	while (t > -1) //while t is a triangle and not -1
	{
		//find next edge to cross, travelling clockwise around p
		nextEdge = K_Thiessen::PtIndex(p, t, TriPts); //index of node p in t
		//get index i in t of neighbour joining p along nextEdge
		i = nextEdge + 1;
		if (i > 2)
			i = 0;
		//record neighbour
		neighboursOut[numNeighbours++] = TriPts[3 * t + i]; //TriPts[t,i]
		//move to next triangle
		t = TriAdj[3 * t + nextEdge]; //next triangle = TriAdj[t,nextEdge]
		//if next triangle is Tp, we have gone full circle and all neighbours collected
		if (t == Tp)
			return numNeighbours;
	}

	//If this point in the function is reached, p is a point on the boundary

	//Starting with Tp, record nodes *anticlockwise* around p until we hit boundary
	t = Tp;
	while (t > -1) //while t is a triangle, not -1
	{
		//find next edge to cross, travelling anticlockwise around p
		nextEdge = K_Thiessen::PtIndex(p, t, TriPts); //index of node p in t
		nextEdge--;
		if (nextEdge < 0)
			nextEdge = 2;
		//record neighouring node, who has index nextEdge
		neighboursOut[numNeighbours++] = TriPts[3 * t + nextEdge]; //TriPts[t,nextEdge]
		//move to next triangle
		t = TriAdj[3 * t + nextEdge]; //next triangle = TriAdj[t,nextEdge]
	}
	return numNeighbours;
}


//output is ordered vector C of length lenC containing A[0], A[1] ... A[lenA-1] and B[0], B[1] ... B[lenB-1] with no repetitions
//assume lenA+lenB > 0
std::vector<int> K_interpTri::MergeIntVectors(std::vector<int>& A, int lenA, std::vector<int>& B, int lenB, int& lenC)
{
	//Let D be A followed by B
	std::vector<int> D(lenA + lenB);
	for (int i = 0; i < lenA; i++)
		D[i] = A[i];
	for (int i = 0; i < lenB; i++)
		D[i + lenA] = B[i];

	//sort D
	std::sort(std::begin(D), std::end(D), [&](int i1, int i2) { return i1 < i2; });

	//count repetitions
	int numRep = 0;
	for (int i = 1; i < lenA + lenB; i++)
	{
		if (D[i - 1] == D[i])
			numRep++;
	}

	//output C of length lenA + lenB - numRep
	std::vector<int> C(lenA + lenB - numRep);
	C[0] = D[0];
	lenC = 1;
	for (int i = 1; i < lenA + lenB; i++)
	{
		if (D[i] != C[lenC - 1])
			C[lenC++] = D[i];
	}
	return C;
}


//function that fills output (length n) with all neighbours, and neighbours of neighbours, of node p
//Returns size of output (number of elements in output that are not -1)
//isOutput[i] = true if node i is in output. Starts all false. Length n
//qNeighbours is a working vector of length n and consists only of nonnegatives followed by -1s. It is only ever used inside this function
//   but is declared outside this function to stop the efficiency having order O(n).
//Except qNeighbours, only output enters and exits function changed. output should consist only of nonnegatives followed by -1s.
//isOutput should be all false on entry and exit. Only used in fucntion, but declared outside to avoid efficiency order O(n)
//may include p in answer
//user can set discSize to 0 to output empty vector, 1 to output immediate neighbours, 2 to output immediate neighbours and their neighbours
int K_interpTri::Neighbours2(int n, int p, std::vector<bool>& isOutput, std::vector<int>& output, std::vector<int>& qNeighbours,
	std::vector<int>& randTri, std::vector<int>& TriAdj, std::vector<int>& TriPts, int discSize)
{
	if (discSize == 0) //user wants nothing done
		return 0;

	//put immediate neighbours of p into output
	int pNumNeighbours = K_interpTri::Neighbours(p, randTri[p], TriAdj, TriPts, output, n);

	if (discSize == 1) //user only wants immediate neighbours
		return pNumNeighbours;

	//mark each neighbour
	for (int i = 0; i < pNumNeighbours; i++)
		isOutput[output[i]] = true;

	//add neighbours of each neighbour q to output
	int sizeOutput = pNumNeighbours;
	int q, qj;
	for (int i = 0; i < pNumNeighbours; i++)
	{
		q = output[i];
		if (q < 0)
			break; //checked all neighbours q of p
		//get neighbours of q
		int qNumNeighbours = K_interpTri::Neighbours(q, randTri[q], TriAdj, TriPts, qNeighbours, n);
		//place each neighbour of q into output if it is not already in there
		for (int j = 0; j < qNumNeighbours; j++)
		{
			qj = qNeighbours[j]; //jth neighbour of q
			if (!isOutput[qj]) //if qj not in output, add to output
			{
				//add qj to current disc
				output[sizeOutput++] = qj;
				isOutput[qj] = true;
			}
		}
	}

	//before exiting function, reset isOutput to its original state, all false
	for (int i = 0; i < n; i++)
	{
		if (output[i] < 0)
			break;
		isOutput[output[i]] = false;
	}

	return sizeOutput;
}