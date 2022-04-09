//This C++ file is used for constructing a triangular grid, used in function e01saf() and e01sbf()
//details in https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770025881.pdf

#include "K_Thiessen.h"
#include "K_Cartesian2D.h"
#include "KUtils.h"
#include <vector>


/*algorithm from section 2.2 of https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770025881.pdf
inputs:
n, number of nodes (this means there can't be more than 2n triangles)
x, vector of length n such that x[i] is x coordinate of ith node (x horizontal, increases from left to right)
y, vector of length n such that y[i] is y coordinate of ith node (y vertical, increases from down to up)
outputs:
TriAdj, a 2n by 3 vector such that TriAdj[i,:] = indices of triangles adjacent to ith triangle (counterclockwise, -1 for none, as per Fig 2.1)
TriPts, a 2n by 3 vector such that TriPts[i,:] = indices of nodes, counterclockwise. TriPts[i,j] joins TriAdj[i,j] and TriAdj[i,j+1] as shown

		 2
		/ \
	  2/   \1
	  /     \
	 0_______1
		 0
randTri, vector of length n such that rantTri[i] = some triangle connected to node i
We write tables as vectors such that TriAdj[i,j] = TriAdj[3*i + j] and TriPts[i,j] = TriPts[3*i + j]
Function makes the following assumptions:
   n>=3
   There are no duplicate points
   It is not the case that all points are colinear
   All three output vectors are full of -1s when function starts
*/
void K_Thiessen::TRIGRD(int n, std::vector<double>& x, std::vector<double>& y, std::vector<int>& randTri, std::vector<int>& TriAdj, std::vector<int>& TriPts)
{
	//Find leftmost node q0
	int q0 = 0;
	for (int i = 1; i < n; i++)
	{
		if (x[i] < x[q0]) //if ith node to the left of q0
			q0 = i;
		if (KUtils::AreDoublesEqual(x[i], x[q0]) && y[i] < y[q0]) //if ith node has same x as q0, but smaller y
			q0 = i;
	}


	//Create vector q of length n. q[0] will be the leftmost point q0 and q[1:n - 1] will be the points ordered by their distance from q0
	std::vector<int> q(n);
	std::vector<double> dLeftmost(n); //dLeftmost[i] will be distance^2 of ith node from leftmost node q0
	for (int i = 0; i < n; i++)
		q[i] = i; //q={0,1,2...n}
	K_Cartesian2D::OrderPointsByDistance(x[q0], y[q0], x, y, q, n);


	//First triangle
	//find next point q[k] not colinear with q[0] and q[1]
	int k = 2; //q[k] will be node that completes the first triangle
	while (K_Cartesian2D::Left(x[q[0]], y[q[0]], x[q[1]], y[q[1]], x[q[k]], y[q[k]]) == 0) //while q[0], q[1], q[k] colinear
		k++; //if k reaches n here then all nodes are colinear. We assume this does not happen!
	//relabel q[2] to q[k]
	int qk = q[k];
	for (int i = k; i > 2; i--) //i=k...3
		q[i] = q[i - 1];
	q[2] = qk;
	//points of first triangle
	TriPts[0] = q[0];
	if (K_Cartesian2D::Left(x[q[0]], y[q[0]], x[q[1]], y[q[1]], x[q[2]], y[q[2]]) == 1) //if q2 to the left of an observer at q0 facing q1
	{
		TriPts[1] = q[1]; //anticlockwise order is 012
		TriPts[2] = q[2];
	}
	else //if q2 to the right of an observer at q0 facing q1, then anticlockwise order is 021
	{
		TriPts[1] = q[2]; //anticlockwise order is 012
		TriPts[2] = q[1];
	}

	int numTriangles = 1;


	//define c = (cx,cy) as centroid of first triangle q[0]q[1]q[2]. All angles measured around c.
	double cx = (x[q[0]] + x[q[1]] + x[q[2]]) / 3;
	double cy = (y[q[0]] + y[q[1]] + y[q[2]]) / 3;
	//all angles will be measured anticlockwise around c from the line from c to q[0]
	double ang_q0 = K_Cartesian2D::PseudoAng(x[q0], y[q0], cx, cy);
	//get pseudoangles of all nodes anticlockwise around c from q[0]
	std::vector<double> angle(n, 0.0);
	for (int i = 0; i < n; i++)
		angle[i] = K_Cartesian2D::PseudoAngFromP(x[i], y[i], cx, cy, ang_q0);
	angle[q0] = 0.0;
	//start boundary lists. Linked lists ie boundaryNext[4]=9 means 9 just after 4 if travelling anticlockwise around boundary. Note q0 will always be on boundary
	std::vector<int> boundaryNext(n, -1);
	std::vector<int> boundaryPrev(n, -1);
	if (K_Cartesian2D::Left(x[q[0]], y[q[0]], x[q[1]], y[q[1]], x[q[2]], y[q[2]]) == 1) //if q2 to the left of an observer at q0 facing q1, then anticlockwise order is 012
	{
		boundaryNext[q[0]] = q[1];
		boundaryNext[q[1]] = q[2];
		boundaryNext[q[2]] = q[0];
		boundaryPrev[q[0]] = q[2];
		boundaryPrev[q[2]] = q[1];
		boundaryPrev[q[1]] = q[0];
	}
	else //if q2 to the right of an observer at q0 facing q1, then anticlockwise order is 021
	{
		boundaryNext[q[0]] = q[2];
		boundaryNext[q[2]] = q[1];
		boundaryNext[q[1]] = q[0];
		boundaryPrev[q[0]] = q[1];
		boundaryPrev[q[1]] = q[2];
		boundaryPrev[q[2]] = q[0];
	}


	//assign first triangle to q0,q1,q2
	randTri[q[0]] = 0;
	randTri[q[1]] = 0;
	randTri[q[2]] = 0;


	//for nodes q[3], ..., q[n-1]
	for (k = 3; k < n; k++)
	{
		qk = q[k];
		//find adjacent pair qi and qj on boundary such that ang(qi) < ang(qk) <= ang(qj)
		int qj = q0;
		while (angle[qj] < angle[qk])
		{
			qj = boundaryNext[qj];
			if (qj == q0)
				break; //angle[q0]=0 not 4, so break
		}
		int qi = boundaryPrev[qj];


		//insert qk between qi and qj to make new triangle
		//update boundary points
		K_Thiessen::InsertBoundaryPoint(qi, qk, boundaryNext, boundaryPrev);
		//add triangle such that edge 2 is the non-boundary edge
		TriPts[3 * numTriangles + 0] = qi; //set TriPts[numTriangles, 0] to qi. Go anticlockwise:
		TriPts[3 * numTriangles + 1] = qk; //set TriPts[numTriangles, 1] to qk.
		TriPts[3 * numTriangles + 2] = qj; //set TriPts[numTriangles, 2] to qj. Edge 0 joins qi to qk, edge 1 joins qk to qj, edge 2 joins qj to qi
		//assign this triangle to q[k]
		randTri[qk] = numTriangles;
		//record adjacency between old triangle and new triangle
		//set Tij to already existing triangle joining qi and qj, and find its edge joining qi and qj
		int Tij = K_Thiessen::TriWithEdge(qi, qj, randTri[qi], TriAdj, TriPts);
		int edgeTij = K_Thiessen::GetEdge(Tij, qi, qj, TriPts);
		TriAdj[3 * numTriangles + 2] = Tij; //set TriAdj[numTriangles, 2] to Tij since edge 2 joins qi and qj
		TriAdj[3 * Tij + edgeTij] = numTriangles; //set TriAdj[Tij, edgeTij] to current triangle being added
		numTriangles++;


		//optimise new triangle with Tij
		K_Thiessen::Optimise(numTriangles - 1, Tij, x, y, randTri, TriAdj, TriPts);


		//join qk to boundaryNext[qj], then set qj=boundaryNext[qj] and repeat, until we can't
		int qjj = boundaryNext[qj];
		while (K_Cartesian2D::Left(x[qk], y[qk], x[qjj], y[qjj], x[qj], y[qj]) == 1) //while qj to the left of an observer at qk facing qjj
		{
			//add new triangle such that edge 1 lies on boundary
			TriPts[numTriangles * 3 + 0] = qj; //go anticlockwise
			TriPts[numTriangles * 3 + 1] = qk; //edge 0 joins qj to qk, edge 1 joins qk to qjj, edge 2 joins qjj to qj
			TriPts[numTriangles * 3 + 2] = qjj;
			//remove qj from boundary
			K_Thiessen::DeleteBoundaryPoint(qj, boundaryNext, boundaryPrev);
			//set Tkj to already existing triangle joining qk and qj, and find its edge joining qk and qj
			int Tkj = K_Thiessen::TriWithEdge(qk, qj, randTri[qk], TriAdj, TriPts);
			int edgeTkj = K_Thiessen::GetEdge(Tkj, qk, qj, TriPts);
			//set Tjj to already existing triangle joining qj and qjj, and find its edge joining qj and qjj
			int Tjj = K_Thiessen::TriWithEdge(qj, qjj, randTri[qj], TriAdj, TriPts);
			int edgeTjj = K_Thiessen::GetEdge(Tjj, qj, qjj, TriPts);
			//update adjacencies
			TriAdj[Tkj * 3 + edgeTkj] = numTriangles;
			TriAdj[Tjj * 3 + edgeTjj] = numTriangles;
			TriAdj[numTriangles * 3 + 0] = Tkj; //edge 0 of new triangle joins qj and qk
			TriAdj[numTriangles * 3 + 2] = Tjj; //edge 2 of new triangle joins qjj and qj
			numTriangles++;

			//optimise new triangle with Tjj
			K_Thiessen::Optimise(numTriangles - 1, Tjj, x, y, randTri, TriAdj, TriPts);

			//set qj to qjj and repeat
			qj = qjj;
			qjj = boundaryNext[qj];
		}//finished joining qk to points later in boundary


		//join qk to boundaryPrev[qi], then set qi=boundaryPrev[qi] and repeat, until we can't
		int qii = boundaryPrev[qi];
		while (K_Cartesian2D::Left(x[qk], y[qk], x[qii], y[qii], x[qi], y[qi]) == -1) //while qi to the right of an observer at qk facing qii
		{
			//add new triangle such that edge 0 lies on boundary
			TriPts[numTriangles * 3 + 0] = qii; //go anticlockwise
			TriPts[numTriangles * 3 + 1] = qk; //edge 0 joins qii to qk, edge 1 joins qk to qi, edge 2 joins qi to qii
			TriPts[numTriangles * 3 + 2] = qi;
			//remove qi from boundary
			K_Thiessen::DeleteBoundaryPoint(qi, boundaryNext, boundaryPrev);
			//set Tik to already existing triangle joining qi and qk, and find its edge joining qi and qk
			int Tik = K_Thiessen::TriWithEdge(qi, qk, randTri[qi], TriAdj, TriPts);
			int edgeTik = K_Thiessen::GetEdge(Tik, qi, qk, TriPts);
			//set Tii to already existing triangle joining qi and qii, and find its edge joining qi and qii
			int Tii = K_Thiessen::TriWithEdge(qi, qii, randTri[qi], TriAdj, TriPts);
			int edgeTii = K_Thiessen::GetEdge(Tii, qi, qii, TriPts);
			//update adjacencies
			TriAdj[Tik * 3 + edgeTik] = numTriangles;
			TriAdj[Tii * 3 + edgeTii] = numTriangles;
			TriAdj[numTriangles * 3 + 1] = Tik; //edge 1 of new triangle joins qk and qi
			TriAdj[numTriangles * 3 + 2] = Tii; //edge 2 of new triangle joins qi and qii
			numTriangles++;

			//optimise new triangle with Tii
			K_Thiessen::Optimise(numTriangles - 1, Tii, x, y, randTri, TriAdj, TriPts);

			//set qi to qii and repeat
			qi = qii;
			qii = boundaryPrev[qi];
		}//finished joining qk to points earlier in boundary
	}//move to k+1
}//end function


//insert boundary point qk just after qi
//boundaryNext and boundaryPrev are linked lists, for example
//boundaryNext[11] = 17
//boundaryPrev[17] = 11
//means 17 is immediately after node 11 if you travel around the boundary anticlockwise
void K_Thiessen::InsertBoundaryPoint(int qi, int qk, std::vector<int> & boundaryNext, std::vector<int> & boundaryPrev)
{
	//get point qj currently after qi
	int qj = boundaryNext[qi];
	//move qk in-between qi and qj by updating boudary node lists
	boundaryNext[qi] = qk;
	boundaryNext[qk] = qj;
	boundaryPrev[qj] = qk;
	boundaryPrev[qk] = qi;
}


//delete boundary point qi
//boundaryNext and boundaryPrev are linked lists, for example
//boundaryNext[11] = 17
//boundaryPrev[17] = 11
//means 17 is immediately after node 11 if you travel around the boundary anticlockwise
void K_Thiessen::DeleteBoundaryPoint(int qi, std::vector<int> & boundaryNext, std::vector<int> & boundaryPrev)
{
	//get previous and next points from qi
	int prev = boundaryPrev[qi];
	int next = boundaryNext[qi];
	boundaryNext[prev] = next;
	boundaryPrev[next] = prev;
	boundaryNext[qi] = -1;
	boundaryPrev[qi] = -1;
}


//returns index of edge of triangle joining qi and qj such that TriAdj[triangle,edge] = other triangle sharing both qi and qj
int K_Thiessen::GetEdge(int triangle, int qi, int qj, std::vector<int> & TriPts)
{
	int t3 = triangle * 3;
	for (int pt = 0; pt < 3; pt++) //pt=0,1,2
	{
		if (TriPts[t3 + pt] != qi && TriPts[t3 + pt] != qj) //if pt not qi or qj
		{
			int edge = pt + 1;
			if (edge == 3)
				edge = 0;
			return edge;
		}
	}
	//following line should never be reached
	return -1;
}


//sets t to triangle containing p = (x,y)
//sets edge to that of t facing p
//this function will only work if boundary is never concave, all internal angles must be <= 180
//returns true if t contains p, false if no triangle does (however, edge of t will still be on boundary facing p)
//quicker if t enters function as a triangle near to p
bool K_Thiessen::TRFIND(int& edge, int& t, double x, double y, std::vector<double> & X, std::vector<double> & Y, std::vector<int> & TriAdj, std::vector<int> & TriPts)
{
	//start with triangle t
	if (t < 0)
		t = 0;
	int t3 = t * 3; //so that TriPts[t,i] = TriPts[t3+i]
	//set next variable to true or false depending on whether p is in t, and set edge to the one facing p
	//edge will be 0 if p facing edge 01, 1 if p facing 12, 2 if p facing 20
	bool p_in_t = K_Cartesian2D::PtFromTriangle(x, y, edge, X[TriPts[t3 + 0]], Y[TriPts[t3 + 0]], X[TriPts[t3 + 1]], Y[TriPts[t3 + 1]], X[TriPts[t3 + 2]], Y[TriPts[t3 + 2]]);
	//while p outside t
	while (!p_in_t)
	{
		//if no triangle exists on other side of edge, return false
		if (TriAdj[t3 + edge] < 0)
			return false;
		//otherwise replace t with that triangle
		t = TriAdj[t3 + edge];
		t3 = t * 3;
		//set next variable to true or false depending on whether p is in t, and set edge to the one facing p
		p_in_t = K_Cartesian2D::PtFromTriangle(x, y, edge, X[TriPts[t3 + 0]], Y[TriPts[t3 + 0]], X[TriPts[t3 + 1]], Y[TriPts[t3 + 1]], X[TriPts[t3 + 2]], Y[TriPts[t3 + 2]]);
	}
	return true;
}


//function that returns a triangle t with edge connecting points a and b
//If a and b are on the boundary of the triangulation, there is only one possible output
//Ta is any triangle connected to node a
int K_Thiessen::TriWithEdge(int a, int b, int Ta, std::vector<int> & TriAdj, std::vector<int> & TriPts)
{
	int t = Ta;
	int i;
	int nextEdge = 0;

	//check triangles in clockwise order around a, starting from Ta, until we hit boundary or b
	while (t > -1) //while t is a triangle and not -1
	{
		for (i = 0; i < 3; i++) //check nodes 0,1,2 of triangle t
		{
			if (TriPts[t * 3 + i] == b) //if TriPts[t,i]=b then we found b!
				return t;
			if (TriPts[t * 3 + i] == a) //TriPts[t,i] = a
				nextEdge = i; //next triangle to check is other side of this edge
		}
		t = TriAdj[t * 3 + nextEdge]; //next triangle = TriAdj[t,nextEdge]
	}

	//check triangles in anticlockwise order around a, starting with Ta, until we hit boundary or b
	t = Ta;
	while (t > -1) //while t is a triangle, not -1
	{
		for (i = 0; i < 3; i++) //check all nodes of triangle t
		{
			if (TriPts[t * 3 + i] == b) //if TriPts[t,i] = b we found b!
				return t;
			if (TriPts[t * 3 + i] == a) //TriPts[t,i] = a
			{
				nextEdge = i - 1; //next triangle to check is other side of this edge
				if (nextEdge < 0)
					nextEdge = 2;
			}
		}
		t = TriAdj[t * 3 + nextEdge]; //TriAdj[t,nextEdge]
	}

	//b not connected to a
	return -1;
}

/*given two adjacent triangles T1=ABC and T2=ACD, optimise T1 with T2, i.e. see if it would be better to use T1=ABD and T2=BCD

		B____A                B____A
		|   /|                |\   |
		|T1/ |                | \T1|
 keep   | /T2|  or change to  |T2\ |
		|/___|                |___\|
		C    D                C    D

If the answer is yes, update randTri, TriAdj and TriPts.

We use these names for the surrounding triangles:

		 triAB
		  /\
		 /  \
		B____A
	   /|    |\
triBC / |    | \  triAD
	  \ |    | /
	   \|____|/
		C    D
		 \  /
		  \/
		 triCD

Order of T1 and T2 matters: At the end of the function, we apply the function again to the edges opposite B, ie optimise T1 with triAD, and T2 with triCD*/
void K_Thiessen::Optimise(int T1, int T2, std::vector<double> & X, std::vector<double> & Y, std::vector<int> & randTri, std::vector<int> & TriAdj, std::vector<int> & TriPts)
{
	//let T1=ABC and T2=ACD, travelling anticlockwise as per diagram above
	int A, B, C, D;
	K_Thiessen::GetABCD(T1, T2, A, B, C, D, TriAdj, TriPts);

	if (K_Thiessen::ShouldSwap(X[A], Y[A], X[B], Y[B], X[C], Y[C], X[D], Y[D]))
	{
		//get triangles surrounding ABCD
		int triBC = TriAdj[T1 * 3 + K_Thiessen::GetEdge(T1, B, C, TriPts)];
		int triCD = TriAdj[T2 * 3 + K_Thiessen::GetEdge(T2, C, D, TriPts)];
		int triAD = TriAdj[T2 * 3 + K_Thiessen::GetEdge(T2, A, D, TriPts)];

		//change points of T1
		K_Thiessen::ReplacePt(TriPts, T1, C, D); //replace C with D

		//changes points of T2
		K_Thiessen::ReplacePt(TriPts, T2, A, B); //replace A with B

		//change adjacencies of T1
		TriAdj[T1 * 3 + K_Thiessen::GetEdge(T1, B, D, TriPts)] = T2; //T2 now opposite T1's new edge BC
		TriAdj[T1 * 3 + K_Thiessen::GetEdge(T1, A, D, TriPts)] = triAD; //triAD now opposite T1's new edge AD

		//change adjacencies of T2
		TriAdj[T2 * 3 + K_Thiessen::GetEdge(T2, B, C, TriPts)] = triBC; //triBC now opposite T2's new edge BC
		TriAdj[T2 * 3 + K_Thiessen::GetEdge(T2, B, D, TriPts)] = T1; //T1 now opposite T2's new edge BD

		//changes random assigned triangle for A and C
		randTri[A] = T1; //may have been T2 before
		randTri[C] = T2; //may have been T1 before

		//change adjacencies of triBC
		if (triBC > -1) //if triBC exists
			TriAdj[triBC * 3 + K_Thiessen::GetEdge(triBC, B, C, TriPts)] = T2; //T2 now adjacent to triBC

		//change adjacencies of triAD
		if (triAD > -1) //if triAD exists
			TriAdj[triAD * 3 + K_Thiessen::GetEdge(triAD, A, D, TriPts)] = T1; //T1 now adjacent to triAD

		//having updated all tables, reapply function to T1 and T2 using B as "origin"
		if (triCD > -1) //if triCD exists
			K_Thiessen::Optimise(T2, triCD, X, Y, randTri, TriAdj, TriPts);
		//line above won't change T1 and triAD, so don't need to recalculate them
		if (triAD > -1) //if triAD exists
			K_Thiessen::Optimise(T1, triAD, X, Y, randTri, TriAdj, TriPts);
	}
}


/*finds vertices A,B,C,D such that ABCD is the quadrilateral formed by two adjacent triangles T1=ABC and T2=ACD (travelling anticlockwise):
B____A
|   /|
|T1/ |
| /T2|
|/___|
C    D
*/
void K_Thiessen::GetABCD(int T1, int T2, int& A, int& B, int& C, int& D, std::vector<int> & TriAdj, std::vector<int> & TriPts)
{
	//get common edge AC according to T1 and T2:
	int edge1 = 0;
	int edge2 = 0;
	for (int i = 0; i < 3; i++) //i=0,1,2
	{
		if (TriAdj[T1 * 3 + i] == T2) //TriAdj[T1,i]
			edge1 = i;
		if (TriAdj[T2 * 3 + i] == T1) //TriAdj[T2,i]
			edge2 = i;
	}
	C = TriPts[T1 * 3 + edge1]; //TriPts[T1, edge1]
	edge1++;
	if (edge1 > 2)
		edge1 = 0;
	A = TriPts[T1 * 3 + edge1]; //TriPts[T1, edge1]
	edge1++;
	if (edge1 > 2)
		edge1 = 0;
	B = TriPts[T1 * 3 + edge1]; //TriPts[T1, edge1]
	edge2--;
	if (edge2 < 0)
		edge2 = 2;
	D = TriPts[T2 * 3 + edge2]; //TriPts[T2, edge2]
}


//This function replaces the value oldVal with newVal in triangle t's row in TriPts
void K_Thiessen::ReplacePt(std::vector<int> & TriPts, int t, int oldVal, int newVal)
{
	int t3 = t * 3; //so that TriPts[t3+i] = TriPts[t,i]
	for (int i = 0; i < 3; i++)
	{
		if (TriPts[t3 + i] == oldVal) //TriPts[t, i]
		{
			TriPts[t3 + i] = newVal; //TriPts[t, i]
			return;
		}
	}
}


/*function that returns false if left triangulation preferable, true if right preferable
B___A      B___A
|  /|      |\  |
| / |  or  | \ |
|/__|      |__\|
C   D      C   D
false      true
Based on section III from https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770025881.pdf
*/
bool K_Thiessen::ShouldSwap(double xA, double yA, double xB, double yB, double xC, double yC, double xD, double yD)
{
	//if ABCD ever concave (internal angle >= 180), return false
	if (K_Cartesian2D::Left(xA, yA, xB, yB, xC, yC) < 1) //if C not left of an observer at A facing B
		return false;
	if (K_Cartesian2D::Left(xB, yB, xC, yC, xD, yD) < 1) //if D not left of an observer at B facing C
		return false;
	if (K_Cartesian2D::Left(xC, yC, xD, yD, xA, yA) < 1) //if A not left of an observer at C facing D
		return false;
	if (K_Cartesian2D::Left(xD, yD, xA, yA, xB, yB) < 1) //if B not left of an observer at D facing A
		return false;
	//only reaches here if ABCD strictly convex

	//calculate origin and radius of circle passing through A, B and C
	double x0, y0, rSquared; //circle has origin (x0,y0) and squared radius rSquared
	bool circleFinite = K_Cartesian2D::Circle3Pts(xA, yA, xB, yB, xC, yC, x0, y0, rSquared);
	//if A,B,C all in a line, return true unless all 4 are in a line
	if (!circleFinite)
	{
		//see if other circle defined by BCD finite
		circleFinite = K_Cartesian2D::Circle3Pts(xD, yD, xB, yB, xC, yC, x0, y0, rSquared);
		if (!circleFinite) //ABCD colinear
			return false;
		return true; //only ABC colinear
	}
	//squared distance of D from circle origin
	double dSquared = K_Cartesian2D::DistSquared(x0, y0, xD, yD);
	//if D inside circle
	if (dSquared < rSquared)
		return true;
	return false;
}



//function that returns index of node p in triangle t
int K_Thiessen::PtIndex(int p, int t, std::vector<int> & TriPts)
{
	int row = 3 * t; //row of TriPts corressponding to Tp
	//check each index, return if right one
	for (int i = 0; i < 3; i++)
	{
		if (TriPts[row + i] == p)
			return i;
	}
	return -1; //if this point reached, p not in t
}