#include <vector>
#include "KUtils.h"
#include "K_qshep2.h"

/*Functions that find Q(x,y) given F[k] at kth node (X[k],Y[k])
such that Q is twice-differentiable and coincides with all (X[k],Y[k],F[k])
X, Y, F and R have length n
A is 5 by n
R such that exactly Nw nodes are within distance sqrt(R[k]) of kth node
Returns Q = (w0q0 + w1*q1 + w2*q2 + ... + wn*qn) / (w0 + w1 + w2 + ... + wn)
where
qk = A[0,k] * (x-X[k])^2
   + A[1,k] * (x-X[k])*(y-Y[k])
   + A[2,k] * (y-Y[k])^2
   + A[3,k] * (x-X[k])
   + A[4,k] * (y-Y[k])
   + F[k]
wk = ( (r[k]-d) / (r[k]*d) )^2 where d = euclidean distance from (x,y) to kth node
*/

int maxInt2(const int a, const int b)
{
	if (a < b)
		return b;
	return a;
}

int minInt2(const int a, const int b)
{
	if (a > b)
		return b;
	return a;
}

double minDouble2(const double a, const double b)
{
	if (a > b)
		return b;
	return a;
}

double minDouble3(const double a, const double b, const double c)
{
	return minDouble2(minDouble2(a, b), c);
}

double minDouble4(const double a, const double b, const double c, const double d)
{
	return minDouble3(minDouble2(a, b), c, d);
}

double minDouble5(const double a, const double b, const double c, const double d, const double e)
{
	return minDouble4(minDouble2(a, b), c, d, e);
}

//This function inputs n points (represented by vectors x, y and f, all length n) and outputs values (including vectors fills lcell, lnext, rsq and a) which can be used to interpolate f for any (x,y)
void K_qshep2::qshep2(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, int nq, int nw, int nr,
	std::vector<int>& lcell, std::vector<int>& lnext, double* xmin, double* ymin, double* dx, double* dy, double* rmax, std::vector<double>& rsq, std::vector<double>& a, int* ier)
{

	/*
	THIS SUBROUTINE COMPUTES A SET OF PARAMETERS A AND RSQ DEFINING A SMOOTH
	(ONCE CONTINUOUSLY DIFFERENTIABLE) BI - VARIATE FUNCTION Q(X, Y) WHICH
	INTERPOLATES DATA VALUES F AT SCATTERED NODES(X, Y).THE INTERPOLANT Q MAY
	BE EVALUATED AT AN ARBITRARY POINT BY FUNCTION QS2VAL, AND ITS FIRST
	DERIVATIVES ARE COMPUTED BY SUBROUTINE QS2GRD.

	THE INTERPOLATION SCHEME IS A MODIFIED QUADRATIC SHEPARD METHOD--

	Q = (W(1)*Q(1) + W(2)*Q(2) + .. + W(N)*Q(N)) / (W(1) + W(2) + .. + W(N))

	FOR BIVARIATE FUNCTIONS W(K) AND Q(K).THE NODAL FUNCTIONS ARE GIVEN BY

	Q(K)(X, Y) = A(1, K)*(X - X(K))**2 + A(2, K)*(X - X(K))*(Y - Y(K))
	+ A(3, K)*(Y - Y(K))**2 + A(4, K)*(X - X(K))
	+ A(5, K)*(Y - Y(K)) + F(K) .

	THUS, Q(K) IS A QUADRATIC FUNCTION WHICH INTERPOLATES THE DATA VALUE AT
	NODE K.ITS COEFFICIENTS A(, K) ARE OBTAINED BY A WEIGHTED LEAST SQUARES FIT
	TO THE CLOSEST NQ DATA POINTS WITH WEIGHTS SIMILAR TO W(K).
	NOTE THAT THE RADIUS OF INFLUENCE FOR THE LEAST SQUARES FIT IS FIXED FOR
	EACH K, BUT VARIES WITH K.

	THE WEIGHTS ARE TAKEN TO BE W(K)(X, Y) = ((R(K) - D(K)) + / R(K)*D(K))**2

	WHERE(R(K) - D(K)) + = 0 IF R(K) <= D(K) AND D(K)(X, Y) IS THE EUCLIDEAN
	DISTANCE BETWEEN(X, Y) AND(X(K), Y(K)).THE RADIUS OF INFLUENCE R(K) VARIES
	WITH K AND IS CHOSEN SO THAT NW NODES ARE WITHIN THE RADIUS.
	NOTE THAT W(K) IS NOT DEFINED AT NODE(X(K), Y(K)), BUT Q(X, Y) HAS LIMIT F(K)
	AS(X, Y) APPROACHES(X(K), Y(K)).

	ON INPUT--

	N = NUMBER OF NODES AND ASSOCIATED DATA VALUES. N >= 6.

	X, Y = ARRAYS OF LENGTH N CONTAINING THE CARTESIAN COORDINATES OF THE NODES.

	F = ARRAY OF LENGTH N CONTAINING THE DATA VALUES IN ONE - TO - ONE CORRESPONDENCE WITH THE NODES.

	NQ = NUMBER OF DATA POINTS TO BE USED IN THE LEAST SQUARES FIT FOR
	COEFFICIENTS DEFINING THE NODAL FUNCTIONS Q(K). A HIGHLY RECOMMENDED VALUE IS NQ = 13. 5 <= NQ <= MIN(40, N - 1).

	NW = NUMBER OF NODES WITHIN(AND DEFINING) THE RADII OF INFLUENCE
	R(K) WHICH ENTER INTO THE WEIGHTS W(K).FOR N SUFFICIENTLY LARGE, A RECOMMENDED VALUE IS NW = 19. 1 <= NW <= MIN(40, N - 1).

	NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID DEFINED IN SUBROUTINE STORE2.A RECTANGLE CONTAINING THE NODES IS
	PARTITIONED INTO CELLS IN ORDER TO INCREASE SEARCH EFFICIENCY. NR = SQRT(N / 3) IS RECOMMENDED. NR >= 1.

	THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

	LCELL = ARRAY OF LENGTH >= NR * *2.

	LNEXT = ARRAY OF LENGTH >= N.

	RSQ = ARRAY OF LENGTH >= N.

	A = ARRAY OF LENGTH >= 5N.

	ON OUTPUT--

	LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.REFER TO STORE2.

	LNEXT = ARRAY OF LENGTH N CONTAINING NEXT - NODE INDICES.REFER TO STORE2.

	XMIN, YMIN, DX, DY = MINIMUM NODAL COORDINATES AND CELL DIMENSIONS.REFER TO STORE2.

	RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ-- MAXIMUM RADIUS R(K).

	RSQ = ARRAY CONTAINING THE SQUARES OF THE RADII R(K) WHICH ENTER INTO THE WEIGHTS W(K).

	A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR QUADRATIC NODAL FUNCTION Q(K) IN COLUMN K.

	NOTE THAT THE ABOVE OUTPUT PARAMETERS ARE NOT DEFINED UNLESS IER = 0.

	IER = ERROR INDICATOR--
	IER = 0 IF NO ERRORS WERE ENCOUNTERED.
	IER = 1 IF N, NQ, NW, OR NR IS OUT OF RANGE.
	IER = 2 IF DUPLICATE NODES WERE ENCOUNTERED.
	IER = 3 IF ALL NODES ARE COLLINEAR.

	MODULES REQUIRED BY QSHEP2 -- GETNP2, GIVENS, ROTATE, SETUP2, STORE2*/

	int ierr, irow, lnp, np;
	std::vector<int> npts(40);

	double av, c, ddx, ddy, dmin, fk, rq, rs, rsmx, rsold, rws, s, sum, t, xk, xmn, yk, ymn;

	std::vector<double> b(36); //b 6 by 6 array

	int neq = -1;
	double avsq = -1.0;
	double dtol = 0.01;
	double rtol = 0.00001;
	double sf = 1.0;

	//create boolean array <marked> of length n instead of using lnext[i]<0 to show that node i is marked
	std::vector<bool> marked(n, false);

	int nn = n;
	int nnq = nq;
	int nnw = nw;
	int nnr = nr;
	int nqwmax = maxInt2(nnq, nnw);
	int lmax = minInt2(40, nn - 1);
	if (5 <= nnq && 1 <= nnw && nqwmax <= lmax && nnr >= 1)
	{
		//CREATE THE CELL DATA STRUCTURE, AND INITIALIZE RSMX.

		for (int i = 0; i < n; i++) //lnext should be full of -1s before store2
			lnext[i] = -1;
		for (int i = 0; i < nr * nr; i++) //lcell should be full of zeros before store2
			lcell[i] = -1;
		store2(nn, x, y, nnr, lcell, lnext, &xmn, &ymn, &ddx, &ddy, &ierr); //x length n, y length n, lcell length nr^2, lnext length n.

		if (ierr != 0)
		{
			//NO UNIQUE SOLUTION DUE TO COLLINEAR NODES.
			*xmin = xmn;
			*ymin = ymn;
			*dx = ddx;
			*dy = ddy;
			*ier = 3;
			return;
		}
		rsmx = 0.0;

		//OUTER LOOP ON NODE K

		for (int k = 0; k < nn; k++)
		{
			xk = x[k];
			yk = y[k];
			fk = f[k];

			//MARK NODE K TO EXCLUDE IT FROM THE SEARCH FOR NEAREST NEIGHBORS.

			marked[k] = !marked[k];

			//INITIALIZE FOR LOOP ON NPTS.

			rs = 0.0;
			sum = 0.0;
			rws = 0.0;
			rq = 0.0;
			lnp = 0;

			//COMPUTE NPTS, LNP, RWS, NEQ, RQ, AND AVSQ.
			while (true)
			{
				sum += rs;
				if (lnp != lmax)
				{
					lnp = lnp + 1;
					rsold = rs;
					getnp2(xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, &np, &rs, marked);

					if (KUtils::AreDoublesEqual(sqrt(rs), 0.0))
					{
						//DUPLICATE NODES WERE ENCOUNTERED BY GETNP2.
						*ier = 2;
						return;
					}
					npts[lnp - 1] = np;
					if ((rs - rsold) / rs < rtol)
						continue; // GO TO 10
					if ( KUtils::AreDoublesEqual(rws , 0.0) && lnp > nnw)
						rws = rs;
					if ( KUtils::AreDoublesEqual(rq , 0.0) && lnp > nnq)
					{
						/*RQ = 0 (NOT YET COMPUTED) AND LNP > NQ.RQ = SQRT(RS) IS SUFFICIENTLY LARGE TO(STRICTLY) INCLUDE
						NQ NODES.THE LEAST SQUARES FIT WILL INCLUDE NEQ = LNP - 1 EQUATIONS FOR 5 <= NQ <= NEQ < LMAX <= N - 1.*/
						neq = lnp - 1;
						rq = sqrt(rs);
						avsq = sum / static_cast<double>(neq);
					}

					//BOTTOM OF LOOP -- TEST FOR TERMINATION.

					if (lnp > nqwmax)
						break;
					continue;
				}

				/*ALL LMAX NODES ARE INCLUDED IN NPTS.RWS AND / OR RQ**2 IS
				(ARBITRARILY)TAKEN TO BE 10 PERCENT LARGER THAN THE DISTANCE RS TO THE LAST NODE INCLUDED.*/

				if (KUtils::AreDoublesEqual(rws, 0.0))
					rws = 1.1 * rs;
				if (KUtils::AreDoublesEqual(rq, 0.0))
				{
					neq = lmax;
					rq = sqrt(1.1 * rs);
					avsq = sum / static_cast<double>(neq);
				}
				break;
			}

			//STORE RSQ(K), UPDATE RSMX IF NECESSARY, AND COMPUTE AV.
			rsq[k] = rws;
			if (rws > rsmx)
				rsmx = rws;
			av = sqrt(avsq);

			/*SET UP THE AUGMENTED REGRESSION MATRIX(TRANSPOSED) AS THE COLUMNS OF B,
			AND ZERO OUT THE LOWER TRIANGLE(UPPER TRIANGLE OF B) WITH GIVENS ROTATIONS -- QR DECOMPOSITION WITH ORTHOGONAL MATRIX Q NOT STORED.*/

			int i = 0;
			while (true)
			{
				i++;
				np = npts[i - 1];
				irow = minInt2(i, 6);
				setup2(xk, yk, fk, x[np], y[np], f[np], av, avsq, rq, b, irow - 1);

				if (i == 1)
					continue;
				int irm1 = irow - 1;
				for (int j = 0; j < irm1; j++)
				{
					givens(b, j, irow - 1, &c, &s);
					rotate(5 - j, c, s, b, j + 1, j, irow - 1);
				}
				if (i < neq)
					continue;

				//TEST THE SYSTEM FOR ILL - CONDITIONING.

				dmin = minDouble5(abs(b[0 * 6 + 0]), abs(b[1 * 6 + 1]), abs(b[2 * 6 + 2]), abs(b[3 * 6 + 3]), abs(b[4 * 6 + 4]));
				if (dmin * rq < dtol)
				{
					if (neq != lmax)
					{
						//INCREASE RQ AND ADD ANOTHER EQUATION TO THE SYSTEM TO IMPROVE THE CONDITIONING.THE NUMBER OF NPTS ELEMENTS IS ALSO INCREASED IF NECESSARY.

						bool goto30 = false;
						while (true)
						{
							rsold = rs;
							neq = neq + 1;
							if (neq != lmax)
							{
								if (neq != lnp)
								{
									//NEQ < LNP
									np = npts[neq];
									rs = (x[np] - xk) * (x[np] - xk) + (y[np] - yk) * (y[np] - yk);
									if ((rs - rsold) / rs < rtol)
										continue;
									rq = sqrt(rs);
									goto30 = true;
									break;
								}

								//ADD AN ELEMENT TO NPTS.
								lnp = lnp + 1;
								getnp2(xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, &np, &rs, marked);
								if (np == -1)
								{
									//DUPLICATE NODES WERE ENCOUNTERED BY GETNP2.
									*ier = 2;
									return;
								}
								npts[lnp - 1] = np;
								if ((rs - rsold) / rs < rtol)
									continue;
								rq = sqrt(rs);
								goto30 = true;
								break;
							}
							break;
						}
						if (goto30)
							continue;

						rq = sqrt(1.1 * rs);
						continue;
					}

					//STABILIZE THE SYSTEM BY DAMPING SECOND PARTIALS -- ADD MULTIPLES OF THE FIRST THREE UNIT VECTORS TO THE FIRST THREE EQUATIONS.

					for (i = 0; i < 3; i++)
					{
						b[6 * i + 5] = sf;
						for (int j = i + 1; j < 6; j++)
							b[6 * j + 5] = 0.0;
						for (int j = i; j < 5; j++)
						{
							givens(b, j, 5, &c, &s);
							rotate(6 - j - 1, c, s, b, j + 1, j, 5);
						}
					}

					//TEST THE STABILIZED SYSTEM FOR ILL - CONDITIONING.
					dmin = minDouble5(abs(b[0 * 6 + 0]), abs(b[1 * 6 + 1]), abs(b[2 * 6 + 2]), abs(b[3 * 6 + 3]), abs(b[4 * 6 + 4]));
					if (dmin * rq < dtol)
					{	//NO UNIQUE SOLUTION DUE TO COLLINEAR NODES.
						*xmin = xmn;//150
						*ymin = ymn;
						*dx = ddx;
						*dy = ddy;
						*ier = 3;
						return;
					}
				}
				break;
			}
			//SOLVE THE 5 BY 5 TRIANGULAR SYSTEM FOR THE COEFFICIENTS
			for (int ib = 0; ib < 5; ib++)//DO  ib = 1, 5
			{
				i = 5 - ib;
				t = 0.0;
				if (i != 5)
				{
					for (int j = i; j < 5; j++) //DO j=ip1,5
						t = t + b[6 * j + (i - 1)] * a[j * nn + k];	//t = t + b(j, i) * a(j, k)
				}
				a[(i - 1) * nn + k] = (b[5 * 6 + (i - 1)] - t) / b[(i - 1) * 6 + (i - 1)]; //a(i, k) = (b(6, i) - t) / b(i, i)
			}

			//SCALE THE COEFFICIENTS TO ADJUST FOR THE COLUMN SCALING.
			for (i = 1; i <= 3; i++)//DO  i = 1, 3
				a[(i - 1) * nn + k] = a[(i - 1) * nn + k] / avsq; //a(i, k) = a(i, k) / avsq
			a[3 * nn + k] /= av; //a(4, k) = a(4, k) / av
			a[4 * nn + k] /= av; //a(5, k) = a(5, k) / av

			//UNMARK K AND THE ELEMENTS OF NPTS.
			marked[k] = !marked[k];
			for (i = 0; i < lnp; i++)
			{
				np = npts[i];
				marked[np] = !marked[np];
			}
		}

		//NO ERRORS ENCOUNTERED.

		*xmin = xmn;
		*ymin = ymn;
		*dx = ddx;
		*dy = ddy;
		*rmax = sqrt(rsmx);
		*ier = 0;
		return;
	}

	//N, NQ, NW, OR NR IS OUT OF RANGE.
	*ier = 1;
}

	/*THIS FUNCTION RETURNS THE VALUE Q(PX, PY) WHERE Q IS THE WEIGHTED SUM
	OF QUADRATIC NODAL FUNCTIONS DEFINED IN SUBROUTINE QSHEP2.
	QS2GRD MAY BE CALLED TO COMPUTE A GRADIENT OF Q ALONG WITH THE VALUE, AND / OR TO TEST FOR ERRORS.

	ON INPUT--

	PX, PY = CARTESIAN COORDINATES OF THE POINT P AT WHICH Q IS TO BE EVALUATED.

	N = NUMBER OF NODES AND DATA VALUES DEFINING Q. N >= 6.

	X, Y, F = ARRAYS OF LENGTH N CONTAINING THE NODES AND DATA VALUES INTERPOLATED BY Q.

	NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID. REFER TO STORE2.NR >= 1.

	LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS. REFER TO STORE2.

	LNEXT = ARRAY OF LENGTH N CONTAINING NEXT - NODE INDICES. REFER TO STORE2.

	XMIN, YMIN, DX, DY = MINIMUM NODAL COORDINATES AND CELL DIMENSIONS. DX AND DY MUST BE POSITIVE.REFER TO STORE2.

	RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ -- MAXIMUM RADIUS.

	RSQ = ARRAY OF LENGTH N CONTAINING THE SQUARED RADII WHICH ENTER INTO THE WEIGHTS DEFINING Q.

	A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR THE NODAL FUNCTIONS DEFINING Q.

	INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.THE PARAMETERS OTHER
	THAN PX AND PY SHOULD BE INPUT UNALTERED FROM THEIR VALUES ON OUTPUT FROM
	QSHEP2.THIS FUNCTION SHOULD NOT BE CALLED IF A NONZERO ERROR FLAG WAS
	RETURNED BY QSHEP2.

	ON OUTPUT--

	QS2VAL = FUNCTION VALUE Q(PX, PY) UNLESS N, NR, DX, DY, OR RMAX IS INVALID, IN WHICH CASE NO VALUE IS RETURNED.

	MODULES REQUIRED BY QS2VAL -- NONE

	INTRINSIC FUNCTIONS CALLED BY QS2VAL -- IFIX, SQRT*/

double K_qshep2::qs2val(double px, double py, int n, std::vector<double> & x, std::vector<double> & y, std::vector<double> & f,
	int nr, std::vector<int> & lcell, std::vector<int> & lnext, double xmin, double ymin, double dx, double dy, double rmax, std::vector<double> & rsq, std::vector<double> & a)
{
	const double xp = px;
	const double yp = py;
	if (n < 6 || nr < 1 || dx <= 0.0 || dy <= 0.0 || rmax < 0.0)
		return 0.0;

	//SET IMIN, IMAX, JMIN, AND JMAX TO CELL INDICES DEFINING THE RANGE OF THE SEARCH FOR NODES WHOSE RADII INCLUDE P.
	//THE CELLS WHICH MUST BE SEARCHED ARE THOSE INTERSECTED BY(OR CONTAINED IN) A CIRCLE OF RADIUS RMAX CENTERED AT P.

	int imin = static_cast<int>((xp - xmin - rmax) / dx);
	int imax = static_cast<int>((xp - xmin + rmax) / dx);
	if (imin < 0)
		imin = 0;
	if (imax > nr - 1)
		imax = nr - 1;
	int jmin = static_cast<int>((yp - ymin - rmax) / dy);
	int jmax = static_cast<int>((yp - ymin + rmax) / dy);
	if (jmin < 0)
		jmin = 0;
	if (jmax > nr - 1)
		jmax = nr - 1;

	//THE FOLLOWING IS A TEST FOR NO CELLS WITHIN THE CIRCLE OF RADIUS RMAX.

	if (imin + 1 <= imax + 1 && jmin + 1 <= jmax + 1)
	{
		//ACCUMULATE WEIGHT VALUES IN SW AND WEIGHTED NODAL FUNCTION VALUES IN SWQ.
		//THE WEIGHTS ARE W(K) = ((R - D) + / (R*D))**2 FOR R**2 = RSQ(K) AND D = DISTANCE BETWEEN P AND NODE K.

		double sw = 0.0;
		double swq = 0.0;

		//OUTER LOOP ON CELLS(I, J).
		for (int j = jmin; j <= jmax; j++)
		{
			for (int i = imin; i <= imax; i++)
			{
				int k = lcell[nr * i + j]; //lcell(i, j)
				if (k + 1 != 0)
				{
					//INNER LOOP ON NODES K.
					while (true)
					{
						const double delx = xp - x[k];
						const double dely = yp - y[k];
						const double dxsq = delx * delx;
						const double dysq = dely * dely;
						const double ds = dxsq + dysq;
						const double rs = rsq[k];
						if (ds < rs)
						{
							if (KUtils::AreDoublesEqual(sqrt(ds) , 0.0))
							{
								//(PX, PY) = (X(K), Y(K))
								return f[k];
							}
							const double rds = rs * ds;
							const double rd = sqrt(rds);
							const double w = (rs + ds - rd - rd) / rds;
							sw = sw + w;
							swq = swq + w * (a[k] * dxsq + a[n + k] * delx * dely + a[2 * n + k] * dysq + a[3 * n + k] * delx + a[4 * n + k] * dely + f[k]);
						}

						//BOTTOM OF LOOP ON NODES IN CELL(I, J).

						const int kp = k;
						k = lnext[kp];
						if (k != kp)
							continue;
						break;
					}
				}
			}
		}

		//SW = 0 IFF P IS NOT WITHIN THE RADIUS R(K) FOR ANY NODE K.

		if (KUtils::AreDoublesEqual(sw , 0.0))
			return 0.0; //ALL WEIGHTS ARE 0 AT P

		return swq / sw;
	}
	return 0.0;
}

/*THIS SUBROUTINE COMPUTES THE VALUE AND GRADIENT AT(PX, PY) OF THE
INTERPOLATORY FUNCTION Q DEFINED IN SUBROUTINE QSHEP2.
Q(X, Y) IS A WEIGHTED SUM OF QUADRATIC NODAL FUNCTIONS.

ON INPUT--

PX, PY = CARTESIAN COORDINATES OF THE POINT AT WHICH
Q AND ITS PARTIALS ARE TO BE EVALUATED.

N = NUMBER OF NODES AND DATA VALUES DEFINING Q.
N >= 6.

X, Y, F = ARRAYS OF LENGTH N CONTAINING THE NODES AND
DATA VALUES INTERPOLATED BY Q.

NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID.
REFER TO STORE2.NR >= 1.

LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.
REFER TO STORE2.

LNEXT = ARRAY OF LENGTH N CONTAINING NEXT - NODE INDICES.
REFER TO STORE2.

XMIN, YMIN, DX, DY = MINIMUM NODAL COORDINATES AND CELL DIMENSIONS.
DX AND DY MUST BE POSITIVE.REFER TO STORE2.

RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ -- MAXIMUM RADIUS.

RSQ = ARRAY OF LENGTH N CONTAINING THE SQUARED RADII
WHICH ENTER INTO THE WEIGHTS DEFINING Q.

A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR THE
NODAL FUNCTIONS DEFINING Q.

INPUT PARAMETERS ARE NOT ALTERED BY THIS SUBROUTINE.
THE PARAMETERS OTHER THAN PX AND PY SHOULD BE INPUT UNALTERED FROM THEIR
VALUES ON OUTPUT FROM QSHEP2.THIS SUBROUTINE SHOULD NOT BE CALLED IF A
NONZERO ERROR FLAG WAS RETURNED BY QSHEP2.

ON OUTPUT--

Q = VALUE OF Q AT(PX, PY) UNLESS IER.EQ. 1, IN
WHICH CASE NO VALUES ARE RETURNED.

QX, QY = FIRST PARTIAL DERIVATIVES OF Q AT(PX, PY) UNLESS IER.EQ. 1.

IER = ERROR INDICATOR
IER = 0 IF NO ERRORS WERE ENCOUNTERED.
IER = 1 IF N, NR, DX, DY OR RMAX IS INVALID.
IER = 2 IF NO ERRORS WERE ENCOUNTERED BUT(PX, PY) IS NOT WITHIN
THE RADIUS R(K) FOR ANY NODE K(AND THUS Q = QX = QY = 0)*/
//SUBROUTINE qs2grd(px, py, n, x, y, f, nr, lcell, lnext, xmin, ymin, dx, dy, rmax, rsq, a, q, qx, qy, ier)
void K_qshep2::qs2grd(double px, double py, int n, std::vector<double> & x, std::vector<double> & y, std::vector<double> & f, int nr, std::vector<int> & lcell, std::vector<int> & lnext,
	double xmin, double ymin, double dx, double dy, double rmax, std::vector<double> & rsq, std::vector<double> & a, double* q, double* qx, double* qy, int* ier)
{
	double delx, dely, ds, dxsq, dysq, qk, qkx, qky, rd, rds, rs, sw, swq, swqx, swqy, sws, swx, swy, t, w, wx, wy, xp, yp;
	int i, imax, imin, j, jmax, jmin, k, kp;

	xp = px;
	yp = py;
	if (n >= 6 && nr >= 1 && dx > 0.0 && dy > 0.0 && rmax >= 0.0)
	{
		//SET IMIN, IMAX, JMIN, AND JMAX TO CELL INDICES DEFINING P.
		//THE RANGE OF THE SEARCH FOR NODES WHOSE RADII INCLUDE THE CELLS WHICH MUST
		//BE SEARCHED ARE THOSE INTERSECTED BY(OR CONTAINED IN) A CIRCLE OF RADIUS
		//RMAX CENTERED AT P.
		imin = static_cast<int> ((xp - xmin - rmax) / dx);
		imax = static_cast<int> ((xp - xmin + rmax) / dx);
		if (imin < 0)//	IF(imin < 1) imin = 1
			imin = 0;
		if (imax > nr - 1)//IF(imax > nr) imax = nr
			imax = nr - 1;
		jmin = static_cast<int> ((yp - ymin - rmax) / dy);
		jmax = static_cast<int> ((yp - ymin + rmax) / dy);
		if (jmin < 0) //IF(jmin < 1) jmin = 1
			jmin = 0;
		if (jmax > nr - 1) //IF(jmax > nr) jmax = nr
			jmax = nr - 1;

		//THE FOLLOWING IS A TEST FOR NO CELLS WITHIN THE CIRCLE OF RADIUS RMAX.
		if (imin > imax || jmin > jmax)
		{
			//NO CELLS CONTAIN A POINT WITHIN RMAX OF P, OR
			//SW = 0 AND THUS DS >= RSQ(K) FOR ALL K.
			*q = 0.0;// 50 q = 0.
			*qx = 0.0;
			*qy = 0.0;
			*ier = 2;
			return;//goto 50
		}
		//Q = SWQ / SW = SUM(W(K)*Q(K)) / SUM(W(K)) WHERE THE SUM IS FROM K = 1 TO N,
		//Q(K) IS THE QUADRATIC NODAL FUNCTION, AND W(K) = ((R - D) + / (R*D))**2 FOR
		//RADIUS R(K) AND DISTANCE D(K).
		//THUS
		//QX = (SWQX*SW - SWQ * SWX) / SW * *2  AND
		//QY = (SWQY*SW - SWQ * SWY) / SW * *2
		//WHERE SWQX AND SWX ARE PARTIAL DERIVATIVES WITH RESPECT TO X OF SWQ
		//AND SW, RESPECTIVELY.SWQY AND SWY ARE DEFINED SIMILARLY.

		sw = 0.0;
		swx = 0.0;
		swy = 0.0;
		swq = 0.0;
		swqx = 0.0;
		swqy = 0.0;

		//OUTER LOOP ON CELLS(I, J).
		for (j = jmin; j <= jmax; j++)
		{
			for (i = imin; i <= imax; i++)
			{
				k = lcell[i * nr + j]; //lcell(i, j)
				if (k != -1)
				{
					//INNER LOOP ON NODES K.
					while (true)
					{
						delx = xp - x[k];
						dely = yp - y[k];
						dxsq = delx * delx;
						dysq = dely * dely;
						ds = dxsq + dysq;
						rs = rsq[k];
						if (ds < rs)
						{
							if (KUtils::AreDoublesEqual(sqrt(ds) , 0.0))
							{
								//(PX, PY) = (X(K), Y(K))
								*q = f[k];
								*qx = a[3 * n + k];
								*qy = a[4 * n + k];
								*ier = 0;
								return;
							}
							rds = rs * ds;
							rd = sqrt(rds);
							w = (rs + ds - rd - rd) / rds;
							t = 2.0 * (rd - rs) / (ds * rds);
							wx = delx * t;
							wy = dely * t;
							qkx = 2.0 * a[k] * delx + a[n + k] * dely;
							qky = a[n + k] * delx + 2. * a[2 * n + k] * dely;
							qk = (qkx * delx + qky * dely) / 2.0;
							qkx = qkx + a[3 * n + k];
							qky = qky + a[4 * n + k];
							qk = qk + a[3 * n + k] * delx + a[4 * n + k] * dely + f[k];
							sw = sw + w;
							swx = swx + wx;
							swy = swy + wy;
							swq = swq + w * qk;
							swqx = swqx + wx * qk + w * qkx;
							swqy = swqy + wy * qk + w * qky;
						}

						//BOTTOM OF LOOP ON NODES IN CELL(I, J).
						kp = k;
						k = lnext[kp];
						if (k != kp)
							continue;
						break;
					}
				}
			}
		}

		//SW = 0 IFF P IS NOT WITHIN THE RADIUS R(K) FOR ANY NODE K.
		if (KUtils::AreDoublesEqual(sw , 0.0))
		{
			//NO CELLS CONTAIN A POINT WITHIN RMAX OF P, OR
			//SW = 0 AND THUS DS >= RSQ(K) FOR ALL K.
			*q = 0.0;
			*qx = 0.0;
			*qy = 0.0;
			*ier = 2;
			return;
		}
		*q = swq / sw;
		sws = sw * sw;
		*qx = (swqx * sw - swq * swx) / sws;
		*qy = (swqy * sw - swq * swy) / sws;
		*ier = 0;
		return;


	}

	//INVALID INPUT PARAMETER.
	*ier = 1;
}


	/*GIVEN A SET OF N NODES AND THE DATA STRUCTURE DEFINED IN SUBROUTINE STORE2,
	THIS SUBROUTINE USES THE CELL METHOD TO FIND THE CLOSEST UNMARKED NODE NP TO
	A SPECIFIED POINT P.
	NP IS THEN MARKED BY SETTING LNEXT(NP) TO - LNEXT(NP).  (A NODE IS MARKED IF AND ONLY IF THE CORRESPONDING LNEXT ELEMENT IS NEGATIVE.THE ABSOLUTE
	VALUES OF LNEXT ELEMENTS, HOWEVER, MUST BE PRESERVED.)   THUS, THE CLOSEST
	M NODES TO P MAY BE DETERMINED BY A SEQUENCE OF M CALLS TO THIS ROUTINE.
	NOTE THAT IF THE NEAREST NEIGHBOR TO NODE K IS TO BE DETERMINED(PX = X(K) AND PY = Y(K)), THEN K SHOULD BE MARKED BEFORE THE CALL TO THIS ROUTINE.

	THE SEARCH IS BEGUN IN THE CELL CONTAINING(OR CLOSEST TO) P AND PROCEEDS
	OUTWARD IN RECTANGULAR LAYERS UNTIL ALL CELLS WHICH CONTAIN POINTS WITHIN
	DISTANCE R OF P HAVE BEEN SEARCHED, WHERE R IS THE DISTANCE FROM P TO THE
	FIRST UNMARKED NODE ENCOUNTERED(INFINITE IF NO UNMARKED NODES ARE PRESENT).
	ON INPUT--
	PX, PY = CARTESIAN COORDINATES OF THE POINT P WHOSE
	NEAREST UNMARKED NEIGHBOR IS TO BE FOUND.
	X, Y = ARRAYS OF LENGTH N, FOR N >= 2, CONTAINING
	THE CARTESIAN COORDINATES OF THE NODES.
	NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID.
	NR >= 1.
	LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED WITH CELLS.
	LNEXT = ARRAY OF LENGTH N CONTAINING NEXT - NODE INDICES
	(OR THEIR NEGATIVES).
	XMIN, YMIN, DX, DY = MINIMUM NODAL COORDINATES AND CELL DIMENSIONS.
	DX AND DY MUST BE POSITIVE.
	INPUT PARAMETERS OTHER THAN LNEXT ARE NOT ALTERED BY THIS ROUTINE.WITH
	THE EXCEPTION OF(PX, PY) AND THE SIGNS OF LNEXT ELEMENTS, THESE PARAMETERS
	SHOULD BE UNALTERED FROM THEIR VALUES ON OUTPUT FROM SUBROUTINE STORE2.

	ON OUTPUT--

	NP = INDEX(FOR X AND Y) OF THE NEAREST UNMARKED NODE TO P,
	OR 0 IF ALL NODES ARE MARKED OR NR < 1 OR DX <= 0 OR
	DY <= 0.
	LNEXT(NP) < 0 IF NP.NE. 0.

	DSQ = SQUARED EUCLIDEAN DISTANCE BETWEEN P AND NODE
	NP, OR 0 IF NP = 0*/
void K_qshep2::getnp2(double px, double py, std::vector<double> & x, std::vector<double> & y, int nr, std::vector<int> & lcell, std::vector<int> & lnext,
	double xmin, double ymin, double dx, double dy, int* np, double* dsq, std::vector<bool> & marked)//x length n, y length n, lcell length nr^2, lnext length n. Extra input, marked, of length n
{
	int lmin = -1;
	double rsmin = -1.0;
	const double xp = px;
	const double yp = py;

	//TEST FOR INVALID INPUT PARAMETERS.
	if (nr >= 1 && dx > 0 && dy > 0)
	{
		/*!INITIALIZE PARAMETERS--
		FIRST = TRUE IFF THE FIRST UNMARKED NODE HAS YET TO BE ENCOUNTERED,
		IMIN, IMAX, JMIN, JMAX = CELL INDICES DEFINING THE RANGE OF THE SEARCH,
		DELX, DELY = PX - XMIN AND PY - YMIN,
		I0, J0 = CELL CONTAINING OR CLOSEST TO P,
		I1, I2, J1, J2 = CELL INDICES OF THE LAYER WHOSE INTERSECTION WITH THE RANGE
		DEFINED BY IMIN, ..., JMAX IS CURRENTLY BEING SEARCHED.*/

		bool first = true;
		int imin = 0;
		int imax = nr - 1;
		int jmin = 0;
		int jmax = nr - 1;
		const double delx = xp - xmin;
		const double dely = yp - ymin;
		int i0 = static_cast<int>(delx / dx);
		if (i0 + 1 < 1)
			i0 = 0;
		if (i0 + 1 > nr)
			i0 = nr - 1;
		int j0 = static_cast<int>(dely / dy);
		if (j0 < 0)
			j0 = 0;
		if (j0 > nr - 1)
			j0 = nr - 1;
		int i1 = i0;
		int i2 = i0;
		int j1 = j0;
		int j2 = j0;

		//OUTER LOOP ON LAYERS, INNER LOOP ON LAYER CELLS, EXCLUDING
		//THOSE OUTSIDE THE RANGE(IMIN, IMAX) X(JMIN, JMAX).

		//fortran stuff
		bool goto10 = true;
		while (goto10)
		{
			goto10 = false;
			int i = -1;
			for (int j = j1; j <= j2; j++)
			{
				if (j > jmax)
					break;
				if (j >= jmin)
				{
					for (i = i1; i <= i2; i++)
					{
						if (i > imax)
							break;
						if (i >= imin)
						{
							if (j == j1 || j == j2 || i == i1 || i == i2)
							{
								//SEARCH CELL(I, J) FOR UNMARKED NODES L.
								int l = lcell[i * nr + j]; //lcell(i, j)
								if (l != -1)
								{
									//LOOP ON NODES IN CELL(I, J).
									bool goto20 = true;
									while (goto20)
									{
										goto20 = false;
										const int ln = lnext[l];
										if (!marked[l])
										{
											//NODE L IS NOT MARKED.
											const double rsq = (x[l] - xp) * (x[l] - xp) + (y[l] - yp) * (y[l] - yp);
											if (first)
											{
												/*NODE L IS THE FIRST UNMARKED NEIGHBOR OF P ENCOUNTERED.
												INITIALIZE LMIN TO THE CURRENT CANDIDATE FOR NP, AND RSMIN TO THE
												SQUARED DISTANCE FROM P TO LMIN.IMIN, IMAX, JMIN, AND JMAX ARE UPDATED
												TO DEFINE THE SMALLEST RECTANGLE CONTAINING A CIRCLE OF RADIUS R =
												SQRT(RSMIN) CENTERED AT P, AND CONTAINED IN(1, NR) X(1, NR) (EXCEPT
												THAT, IF P IS OUTSIDE THE RECTANGLE DEFINED BY THE NODES, IT IS POSSIBLE
												THAT IMIN > NR, IMAX < 1, JMIN > NR, OR JMAX < 1).
												FIRST IS RESET TO FALSE.*/

												lmin = l;
												rsmin = rsq;
												const double r = sqrt(rsmin);
												imin = static_cast<int>((delx - r) / dx);
												if (imin < 0)
													imin = 0;
												imax = static_cast<int>((delx + r) / dx);
												if (imax + 1 > nr)
													imax = nr - 1;
												jmin = static_cast<int>((dely - r) / dy);
												if (jmin < 0)
													jmin = 0;
												jmax = static_cast<int>((dely + r) / dy);
												if (jmax > nr - 1)
													jmax = nr - 1;
												first = false;
											}
											else
											{
												//TEST FOR NODE L CLOSER THAN LMIN TO P.
												if (rsq < rsmin)
												{
													//UPDATE LMIN AND RSMIN.
													lmin = l;
													rsmin = rsq;
												}
											}
										}

										//TEST FOR TERMINATION OF LOOP ON NODES IN CELL(I, J).
										if (abs(ln) != l)
										{
											l = abs(ln);
											goto20 = true;
										}
									}//end while goto20
								}//END IF
							}//END IF
						}//END IF
					}//END DO
				}//END IF
			}//END DO loop40

			//TEST FOR TERMINATION OF LOOP ON CELL LAYERS.

			if (i1 > imin || i2  < imax || j1  > jmin || j2 < jmax)
			{
				i1--;
				i2++;
				j1--;
				j2++;
				goto10 = true;
			}
		}//end while goto10
		//UNLESS NO UNMARKED NODES WERE ENCOUNTERED, LMIN IS THE CLOSEST UNMARKED NODE TO P.

		if (!first)
		{
			*np = lmin;
			*dsq = rsmin;
			marked[lmin] = !marked[lmin];
			return;
		}
	}//END IF

	//ERROR -- NR, DX, OR DY IS INVALID OR ALL NODES ARE MARKED.
	*np = -1;
	*dsq = 0.0;
}


/*THIS ROUTINE CONSTRUCTS THE GIVENS PLANE ROTATION--
(C  S)
G = () WHERE C*C + S * S = 1 --WHICH ZEROS THE SECOND ENTRY OF THE
(-S  C)
2 - VECTOR(A B) - TRANSPOSE.A CALL TO GIVENS IS NORMALLY FOLLOWED BY A CALL
TO ROTATE WHICH APPLIES THE TRANSFORMATION TO A 2 BY N MATRIX.
THIS ROUTINE WAS TAKEN FROM LINPACK.

ON INPUT--
A, B = COMPONENTS OF THE 2 - VECTOR TO BE ROTATED.

ON OUTPUT--
A = VALUE OVERWRITTEN BY R = +/ -SQRT(A*A + B * B)
B = VALUE OVERWRITTEN BY A VALUE Z WHICH ALLOWS C
AND S TO BE RECOVERED AS FOLLOWS--
C = SQRT(1 - Z * Z), S = Z     IF ABS(Z) <= 1.
C = 1 / Z, S = SQRT(1 - C * C) IF ABS(Z) > 1.
C = +/ -(A / R)
S = +/ -(B / R)*/
void K_qshep2::givens(std::vector<double> & B, int j, int irow, double* c, double* s)
{
	double u, v;
	double a = B[j * 6 + j];
	double b = B[j * 6 + irow];

	const double aa = a;
	const double bb = b;
	if (abs(aa) > abs(bb))
	{
		//ABS(A) > ABS(B)
		u = aa + aa;
		v = bb / u;
		const double r = sqrt(0.25 + v * v) * u;
		*c = aa / r;
		*s = v * (*c + *c);

		//NOTE THAT R HAS THE SIGN OF A, C > 0, AND S HAS SIGN(A)*SIGN(B).
		b = *s;
		a = r;

		B[j * 6 + j] = a;
		B[j * 6 + irow] = b;
		return;
	}

	//ABS(A) <= ABS(B)

	if (!KUtils::AreDoublesEqual(bb, 0.0))
	{
		u = bb + bb;
		v = aa / u;

		//STORE R IN A.
		a = sqrt(0.25 + v * v) * u;
		*s = bb / a;
		*c = v * (*s + *s);

		//NOTE THAT R HAS THE SIGN OF B, S > 0, AND C HAS SIGN(A)*SIGN(B).
		b = 1.0;
		if (!KUtils::AreDoublesEqual(*c, 0.0))
			b = 1.0 / *c;

		B[j * 6 + j] = a;
		B[j * 6 + irow] = b;
		return;
	}

	//A = B = 0.
	*c = 1.0;
	*s = 0.0;
}

/*THIS ROUTINE APPLIES THE GIVENS ROTATION
   C  S
  -S  C
TO THE 2 BY N MATRIX
(X(1) ... X(N))
(Y(1) ... Y(N))

ON INPUT--
N = NUMBER OF COLUMNS TO BE ROTATED.
C, S = ELEMENTS OF THE GIVENS ROTATION.THESE MAY BE DETERMINED BY SUBROUTINE GIVENS.
X, Y = ARRAYS OF LENGTH >= N CONTAINING THE VECTORS TO BE ROTATED.
PARAMETERS N, C, AND S ARE NOT ALTERED BY THIS ROUTINE.

ON OUTPUT--
X, Y = ROTATED VECTORS*/

void K_qshep2::rotate(int n, double c, double s, std::vector<double> & b, int jp1, int j, int irow)
{
	if (n <= 0 || (KUtils::AreDoublesEqual(c, 1.0) && KUtils::AreDoublesEqual(s, 0.0)))
		return;

	for (int i = 0; i < n; i++)
	{
		const double xi = b[(jp1 + i) * 6 + j];
		const double yi = b[(jp1 + i) * 6 + irow];
		b[(jp1 + i) * 6 + j] = c * xi + s * yi;
		b[(jp1 + i) * 6 + irow] = -s * xi + c * yi;
	}
}

/*THIS ROUTINE SETS UP THE I - TH ROW OF AN AUGMENTED REGRESSION MATRIX FOR A
WEIGHTED LEAST - SQUARES FIT OF A QUADRATIC FUNCTION Q(X, Y) TO A SET OF DATA
VALUES F, WHERE Q(XK, YK) = FK.THE FIRST 3 COLUMNS(QUADRATIC TERMS) ARE
SCALED BY 1 / S2 AND THE FOURTH AND FIFTH COLUMNS(LINEAR TERMS) ARE SCALED
BY 1 / S1.THE WEIGHT IS(R - D) / (R*D) IF R > D AND 0 IF R <= D, WHERE D
IS THE DISTANCE BETWEEN NODES I AND K.

ON INPUT--
XK, YK, FK = COORDINATES AND DATA VALUE AT NODE K -- INTERPOLATED BY Q.
XI, YI, FI = COORDINATES AND DATA VALUE AT NODE I.
S1, S2 = RECIPROCALS OF THE SCALE FACTORS.
R = RADIUS OF INFLUENCE ABOUT NODE K DEFINING THE WEIGHT.
ROW = ARRAY OF LENGTH 6.
INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

ON OUTPUT--
ROW = VECTOR CONTAINING A ROW OF THE AUGMENTED REGRESSION MATRIX.*/
void K_qshep2::setup2(double xk, double yk, double fk, double xi, double yi, double fi, double s1, double s2, double r, std::vector<double> & M, int irow)
{
	const double dx = xi - xk;
	const double dy = yi - yk;
	const double dxsq = dx * dx;
	const double dysq = dy * dy;
	const double d = sqrt(dxsq + dysq);
	if (d > 0 && d < r)
	{
		const double w = (r - d) / (r * d);
		const double w1 = w / s1;
		const double w2 = w / s2;
		M[irow] = dxsq * w2;
		M[irow + 6] = dx * dy * w2;
		M[irow + 12] = dysq * w2;
		M[irow + 18] = dx * w1;
		M[irow + 24] = dy * w1;
		M[irow + 30] = (fi - fk) * w;
		return;
	}

	//NODES K AND I COINCIDE OR NODE I IS OUTSIDE OF THE RADIUS OF INFLUENCE. SET ROW TO THE ZERO VECTOR.
	for (int i = 0; i < 6; i++)
		M[irow + i * 6] = 0.0;
}


/*GIVEN A SET OF N ARBITRARILY DISTRIBUTED NODES IN THE PLANE, THIS
SUBROUTINE CREATES A DATA STRUCTURE FOR A CELL - BASED METHOD OF SOLVING
CLOSEST - POINT PROBLEMS.THE SMALLEST RECTANGLE CONTAINING THE NODES IS
PARTITIONED INTO AN NR BY NR UNIFORM GRID OF CELLS, AND NODES ARE
ASSOCIATED WITH CELLS.IN PARTICULAR, THE DATA STRUCTURE STORES THE
INDICES OF THE NODES CONTAINED IN EACH CELL.
FOR A UNIFORM RANDOM DISTRIBUTION OF NODES, THE NEAREST NODE TO AN
ARBITRARY POINT CAN BE DETERMINED IN CONSTANT EXPECTED TIME.

ON INPUT--
N = NUMBER OF NODES.N >= 2.
X, Y = ARRAYS OF LENGTH N CONTAINING THE CARTESIAN COORDINATES OF THE NODES.
NR = NUMBER OF ROWS AND COLUMNS IN THE GRID.THE CELL DENSITY
(AVERAGE NUMBER OF NODES PER CELL) IS D = N / (NR**2).
A RECOMMENDED VALUE, BASED ON EMPIRICAL EVIDENCE,
IS D = 3 --NR = SQRT(N / 3).NR >= 1.

THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.

LCELL = ARRAY OF LENGTH >= NR * *2.
LNEXT = ARRAY OF LENGTH >= N.
ON OUTPUT--
LCELL = NR BY NR CELL ARRAY SUCH THAT LCELL(I, J) CONTAINS THE INDEX
(FOR X AND Y) OF THE FIRST NODE(NODE WITH SMALLEST INDEX) IN
CELL(I, J), OR LCELL(I, J) = 0 IF NO NODES ARE CONTAINED IN
THE CELL.THE UPPER RIGHT CORNER OF CELL(I, J) HAS
COORDINATES(XMIN + I * DX, YMIN + J * DY).
LCELL IS NOT DEFINED IF IER.NE. 0.
LNEXT = ARRAY OF NEXT - NODE INDICES SUCH THAT LNEXT(K) CONTAINS THE
INDEX OF THE NEXT NODE IN THE CELL WHICH CONTAINS NODE K, OR
LNEXT(K) = K IF K IS THE LAST NODE IN THE CELL FOR
K = 1, ..., N.  (THE NODES CONTAINED IN A CELL ARE ORDERED BY THEIR INDICES.)
IF, FOR EXAMPLE, CELL(I, J) CONTAINS NODES 2, 3, AND 5 (AND NO OTHERS), THEN LCELL(I, J) = 2, LNEXT(2) = 3, LNEXT(3) = 5,
AND LNEXT(5) = 5.  LNEXT IS NOT DEFINED IF IER.NE. 0.
XMIN, YMIN = CARTESIAN COORDINATES OF THE LOWER LEFT CORNER OF THE
RECTANGLE DEFINED BY THE NODES(SMALLEST NODAL COORDINATES) UNLESS IER = 1.  THE UPPER RIGHT CORNER IS
(XMAX, YMAX) FOR XMAX = XMIN + NR * DX AND
YMAX = YMIN + NR * DY.

DX, DY = DIMENSIONS OF THE CELLS UNLESS IER = 1.  DX = (XMAX - XMIN) / NR
AND DY = (YMAX - YMIN) / NR WHERE XMIN, XMAX, YMIN, AND YMAX ARE
THE EXTREMA OF X AND Y.

IER = ERROR INDICATOR--
IER = 0 IF NO ERRORS WERE ENCOUNTERED.
IER = 1 IF N < 2 OR NR < 1.
IER = 2 IF DX = 0 OR DY = 0.*/

void K_qshep2::store2(int n, std::vector<double> & x, std::vector<double> & y, int nr, std::vector<int> & lcell, std::vector<int> & lnext,
	double* xmin, double* ymin, double* dx, double* dy, int* ier) //x length n, y length n, lcell length nr^2, lnext length n
{
	int i, k;

	const int nn = n;
	const int nnr = nr;
	if (nn >= 2 && nnr >= 1)
	{
		//COMPUTE THE DIMENSIONS OF THE RECTANGLE CONTAINING THE NODES.
		double xmn = x[0];
		double xmx = xmn; //will be maximum value in x
		double ymn = y[0];
		double ymx = ymn; //will be maximum value in y
		for (k = 1; k < nn; k++)
		{
			if (x[k] < xmn) xmn = x[k];
			if (x[k] > xmx) xmx = x[k];
			if (y[k] < ymn) ymn = y[k];
			if (y[k] > ymx) ymx = y[k];
		}
		*xmin = xmn; //minimum value in x
		*ymin = ymn; //minimum value in y

		//COMPUTE CELL DIMENSIONS AND TEST FOR ZERO AREA.
		const double delx = (xmx - xmn) / nnr; //span of x divided by nr = vertical range of a cell of lcell
		const double dely = (ymx - ymn) / nnr; //span of y divided by nr = horizontal range of a cell of lcell
		*dx = delx;
		*dy = dely;
		if (KUtils::AreDoublesEqual(delx, 0.0) || KUtils::AreDoublesEqual(dely, 0.0)) //IF(delx == 0..OR.dely == 0.) GO TO 50
		{
			*ier = 2;
			return;
		}

		//INITIALIZE LCELL.
		for (i = 0; i < nnr * nnr; i++)
			lcell[i] = -1;

		//LOOP ON NODES, STORING INDICES IN LCELL AND LNEXT.
		for (k = 0; k < nn; k++)
		{
			const int kb = nn - k - 1;
			//want kb to be nn-1,nn-2,... as k goes 0,1,...
			//find entry of lcell that (x[kb],y[kb]) should land in
			i = static_cast<int>((x[kb] - xmn) / delx);
			if (i >= nnr)
				i = nnr - 1;
			int j = static_cast<int>((y[kb] - ymn) / dely);
			if (j >= nnr)
				j = nnr - 1;
			//find out if a (x[l],y[l]) (where l>kb) already landed in that entry. If it did, save it in lnext
			const int l = lcell[i * nnr + j]; //l = lcell(i, j)
			lnext[kb] = l;
			if (l == -1) { lnext[kb] = kb; }
			lcell[i * nnr + j] = kb;
		}

		//NO ERRORS ENCOUNTERED
		*ier = 0;
		return;
	}

	//INVALID INPUT PARAMETER
	*ier = 1;
}


//Function that changes values (including a and rsq) so that e01shf() can use them
//Before this function is called, removeDuplicates() or some method must be used to ensure there are no duplicate coordinates (x[i],y[i]) = (x[j],y[j])
void K_qshep2::e01sgf(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, int nq, int nw, int nr, std::vector<int> & lcell, std::vector<int> & lnext,
	double& xmin, double& ymin, double& dx, double& dy, double& rmax, std::vector<double> & rsq, std::vector<double> & a, int& ifail)
{
	//nw = number of data points inside each radius of influence RW
	int nw_ = nw;
	if (nw_ <= 0)
		nw_ = minInt2(19, n - 1);
	//nq = number of data points to be used in the least squares fit
	int nq_ = nq;
	if (nq_ <= 0)
		nq_ = minInt2(13, n - 1);
	//nr = number of data
	const int nr_ = nr; //changing nr makes no difference to results
	if (nr < 2)
		nr = 2;
	//apply QSHEP2 algorithm
	qshep2(n, x, y, f, nq_, nw_, nr_, lcell, lnext, &xmin, &ymin, &dx, &dy, &rmax, rsq, a, &ifail);
}



//X, Y, F vectors of length n representing n points (x,y) each assigned a value f
//U, V, Q, QX, QY double arrays. Function interpolates value of f, df/dx and df/dy at the points from U and V, stores answers in Q, QX, QY
void K_qshep2::e01shf(int n, std::vector<double>& x, std::vector<double>& y, std::vector<double>& f, int nr, std::vector<int> & lcell, std::vector<int> & lnext, double xmin, double ymin,
	double dx, double dy, double rmax, std::vector<double> & rsq, std::vector<double> & a, int lenUV, const double* U, const double* V, double* Q, double* QX, double* QY, int& IFAIL)
{
	//for each (u,v)
	double q, qx, qy;
	IFAIL = 0;
	int ifailTemp = 0;
	for (int i = 0; i < lenUV; i++)
	{
		qs2grd(U[i], V[i], n, x, y, f, nr, lcell, lnext, xmin, ymin, dx, dy, rmax, rsq, a, &q, &qx, &qy, &ifailTemp);
		Q[i] = q;
		QX[i] = qx;
		QY[i] = qy;
		//if we only want Q, not QX or QY, then can do Q[i] = K_qshep2::qs2val(U[i], V[i], n, xVec, yVec, fVec, nr, lcell, lnext, xmin, ymin, dx, dy, rmax, rsq, a)
		if (ifailTemp == 1)
			IFAIL = 2; //N, NR, DX, DY or RMAX is invalid
		if (ifailTemp == 2)
			IFAIL = 3; //(u,v) not within the radius R(k) of any node
	}
}

//removes duplicate coordinates (x,y) where "duplicate" means within dist of each other
//In: X, Y, F (length n)
//Out: newX, newY, newF (at start, length n, full of -1s. But number of sensible values returned by function)
int K_qshep2::removeDuplicates(std::vector<double>& X, std::vector<double>& Y, std::vector<double>& F, std::vector<double>& newX, std::vector<double>& newY,
	std::vector<double>& newF, int n, double dist)
{
	//min and max values
	double xMin = X[0];
	double xMax = X[0];
	double yMin = Y[0];
	double yMax = Y[0];
	for (int i = 0; i < n; i++)
	{
		if (X[i] < xMin)
			xMin = X[i];
		if (X[i] > xMax)
			xMax = X[i];
		if (Y[i] < yMin)
			yMin = Y[i];
		if (Y[i] > yMax)
			yMax = Y[i];
	}

	//make a bit wider, just in case
	xMin -= dist;
	xMax += dist;
	yMin -= dist;
	yMax += dist;
	const double xSpan = xMax - xMin;
	const double ySpan = yMax - yMin;

	//dimensions of grid, such that each cell has height dx>dist and width dy>dist
	//xMin will be the top, xMax the bottom, yMin the left, yMax the right
	int rows = static_cast<int>((sqrt(n)) + 1); //default
	int cols = rows; //default
	//ensure height per row > dist
	if (xSpan / rows <= dist) //if height per row <= dist, increase height by decreasing number of rows:
		rows = static_cast<int> (xSpan / dist) - 1; //now height per row > dist
	if (rows < 1)
		rows = 1;
	//ensure width per column > dist
	if (ySpan / cols <= dist) //if width per column <= dist, increase width by decreasing number of columns:
		cols = static_cast<int> (ySpan / dist) - 1; //now width per column > dist
	if (cols < 1)
		cols = 1;
	const double dx = xSpan / rows; //height per row
	const double dy = ySpan / cols; //width per column

	//initiate grid such that grid[a,b]=i if the point (X[i],Y[i]) satisfies the following
	//   a <= (X[i]-xMin)/dx < a+1 and
	//   b <= (Y[i]-yMin)/dy < b+1 and
	//   next[i] = j means (X[j],Y[j]) lies in same cell of grid as (X[i],Y[i]) (unless j=-1)
	std::vector<int> grid(rows * cols, -1); //m by m, all -1s
	std::vector<int> next(n, -1); //length n, all -1s
	int lengthNew = 0; //will be length of output vectors
	const double d2 = dist * dist;
	for (int i = 0; i < n; i++)
	{
		bool iDuplicate = false;
		const int a = static_cast<int>((X[i] - xMin) / dx);
		const int b = static_cast<int>((Y[i] - yMin) / dy);
		//loop through 9 cells around grid[a,b]
		for (int aa = maxInt2(a - 1, 0); aa <= minInt2(a + 1, rows - 1); aa++) //want aa = a-1, a, a+1
		{
			for (int bb = maxInt2(b - 1, 0); bb <= minInt2(b + 1, cols - 1); bb++) //want bb = b-1, b, b+1
			{
				//loop through points (j) in grid[aa,bb]
				int j = grid[aa * cols + bb];
				while (j != -1)
				{
					//if point j within dist of point i, skip i and move to i+1
					const double diffx = X[i] - X[j];
					const double diffy = Y[i] - Y[j];
					if (diffx * diffx + diffy * diffy <= d2)
					{
						iDuplicate = true;
						break;
					}
					//if point j further than dist from point i, check next j
					j = next[j];
				}
				if (iDuplicate)
					break; //break out of bb loop
			}
			if (iDuplicate)
				break; //break out of aa loop
		}
		if (iDuplicate)
			continue; //move on to i+1
		//if this point is reached, i is not a duplicate, update:
		next[i] = grid[a * cols + b]; //grid[a,b]
		grid[a * cols + b] = i;
		newX[lengthNew] = X[i];
		newY[lengthNew] = Y[i];
		newF[lengthNew] = F[i];
		lengthNew++;
	}
	return lengthNew;
}