#include "matvec_functions.h"

#include <string>

/**
 * pretty printing
 */
void prettyPrintArray(const std::string &var, Uint1dArray& array)
{
	printf("%s = ", var.c_str());
	array.print();
}

/**
 * pretty printing
 */
void prettyPrintArray(const std::string &var, Float1dArray& array)
{
	printf("%s = ", var.c_str());
	array.print();
}

/**
 * Test set if empty
 */
bool isempty(Uint1dArray &a)
{
	return (a.x == 0);
}

/**
 * copy (replace
 */
int copyArray(Uint1dArray &from, Uint1dArray &to)
{
	if (from.x > to.max)
		return -1;

	memcpy(to.data, from.data, from.size());
	to.shrink(from.x);
	return 0;
}

/**
 * copy (replace
 */
int copyArray(Float1dArray &from, Float1dArray &to)
{
	if (from.x > to.max)
		return -1;

	memcpy(to.data, from.data, from.size());
	to.shrink(from.x);
	return 0;
}

/**
 * copy (replace
 */
int copyArray(Float2dArray &from, Float2dArray &to)
{
	if (from.x > to.x_max)
		return -1;

	if (from.y > to.y_max)
		return -1;

	memcpy(to.data, from.data, from.size());
	to.shrink(from.x, from.y);
	return 0;
}

/**
 * copy (replace
 */
int copyArray(Double1dArray &from, Double1dArray &to)
{
	if (from.x > to.max)
		return -1;

	memcpy(to.data, from.data, from.size());
	to.shrink(from.x);
	return 0;
}

/**
 * Extract row
 */
void extractRow(Uint2dArray &from, const unsigned int &row, const unsigned int &length, Uint1dArray &to)
{
	memcpy(to.data, &(from.data[row*from.y]), length*sizeof(unsigned int));
	to.shrink(length);
}


/**
 * 2d to 1d
 */
void convert2dTo1d(Uint2dArray &a, Uint1dArray &b)
{
	if (b.length() == a.length())
		memcpy(b.data, a.data, a.size());
}

/**
 * 2d to 1d
 */
void convert2dTo1d(Float2dArray &a, Float1dArray &b)
{
	if (b.length() == a.length())
		memcpy(b.data, a.data, a.size());
}

void convert2dTo1dFloatToUint(Float2dArray &a, Uint1dArray &b)
{
	for (int i = 0; i < a.x; ++i)
	{
		for (int j = 0; j < a.y; ++j)
		{
			b[i*a.y+j] = (unsigned int)(a(i, j));
		}
	}
}

/**
 * Transpose matrix: m -> m'
 */
void transpose(Float2dArray &m, Float2dArray& m_trans)
{
	for (int i = 0; i < m.x; ++i)
	{
		for (int j = 0; j < m.y; ++j)
		{
			m_trans(j, i) = m(i, j);
		}
	}
}


/**
 * Transpose matrix: m -> m'
 */
void transpose(Double2dArray &m, Double2dArray& m_trans)
{
	for (int i = 0; i < m.x; ++i)
	{
		for (int j = 0; j < m.y; ++j)
		{
			m_trans(j, i) = m(i, j);
		}
	}
}

/**
 * Transpose matrix: m -> m'
 */
void transposeFloatToUint(Float2dArray &m, Uint2dArray& m_trans)
{
	for (int i = 0; i < m.x; ++i)
	{
		for (int j = 0; j < m.y; ++j)
		{
			m_trans(j, i) = (unsigned int)(m(i, j));
		}
	}
}

/**
 * Standard matrix multiply
 */
void multiplyMatrix(Float2dArray &a, Float2dArray &b, Float2dArray &c)
{
	memset(c.data, 0, c.length() * sizeof(float));
	for (int row = 0; row < a.x; ++row)
	{
		for (int col = 0; col < b.y; ++col)
		{
			for (int inner = 0; inner < a.y; ++inner)
			{
				c(row, col) += a(row, inner) * b(inner, col);
			}
		}
	}
}

/**
 * Standard matrix multiply
 */
void multiplyMatrix(Double2dArray &a, Double2dArray &b, Double2dArray &c)
{
	memset(c.data, 0, c.length() * sizeof(double));
	for (int row = 0; row < a.x; ++row)
	{
		for (int col = 0; col < b.y; ++col)
		{
			for (int inner = 0; inner < a.y; ++inner)
			{
				c(row, col) += a(row, inner) * b(inner, col);
			}
		}
	}
}


/**
 * Standard scalar divide (c = a./b)
 */
int scalarDivideArray(Float1dArray &a, Float1dArray &b, Float1dArray &c)
{
	for (int col = 0; col < a.x; ++col)
	{
		c[col] = a[col] / b[col];
	}
}


/**
 * 1D mean
 */
unsigned int max(Uint1dArray& p)
{
	unsigned int max = 0;
	for (int i = 0; i < p.x; ++i)
	{
		if (p[i] >= max)
		{
			max = p[i];
		}
	}
	return max;
}

/**
 * 1D mean
 */
float mean(Double1dArray& p)
{
	if (p.x == 0)
		return 0.0f;

	double sum = 0.0f;
	for (int i = 0; i < p.x; ++i)
		sum += p[i];
	return sum / p.x;
}

/**
 * 1D mean
 */
float mean(Float1dArray& p)
{
	if (p.x == 0)
		return 0.0f;

	float sum = 0.0f;
	for (int i = 0; i < p.x; ++i)
		sum += p[i];
	return sum / p.x;
}

/**
 * 2d mean (row-wise)
 */
void mean(Float2dArray& in, Float1dArray &out)
{
	for (int xi = 0; xi < in.x; ++xi)
	{
		float sum = 0.0f;
		for (int yi = 0; yi < in.y; ++yi)
		{
			sum += in(xi, yi);
		}
		out[xi] = sum / (in.y);
	}
}

/**
 * Standard Deviation 1d
 */
float stddev(Float1dArray &p, const float arraymean)
{
	float sum = 0.0f;
	for (int i = 0; i < p.x; ++i)
		sum += pow(p[i] - arraymean, 2);
	return sqrt(sum / (p.x - 1));
}


/**
 * Standard Deviation 1d
 */
float stddev(Double1dArray &p, const double arraymean)
{
	double sum = 0.0f;
	for (int i = 0; i < p.x; ++i)
		sum += pow(p[i] - arraymean, 2);
	return sqrt(sum / (p.x - 1));
}

/**
 * Standard Deviation 2d  (row-wise)
 */
void stddev(Float2dArray &in, Float1dArray &arraymean, Float1dArray &out)
{
	for (int xi = 0; xi < in.x; ++xi)
	{
		float sum = 0.0f;
		for (int yi = 0; yi < in.y; ++yi)
		{
			sum += pow(in(xi, yi) - arraymean[xi], 2);
		}
		out[xi] = sqrt(sum / (in.y - 1));
	}
}

/**
 * Variance 1d
 */
float var(Float1dArray &p, const float arraymean)
{
	float sum = 0.0f;
	for (int i = 0; i < p.x; ++i)
		sum += pow(p[i] - arraymean, 2);
	return sum / (p.x - 1);
}

/**
 * Variance 2d  (row-wise)
 */
void var(Float2dArray& in, Float1dArray &arraymean, Float1dArray &out)
{
	for (int xi = 0; xi < in.x; ++xi)
	{
		float sum = 0.0f;
		for (int yi = 0; yi < in.y; ++yi)
		{
			sum += pow(in(xi, yi) - arraymean[xi], 2);
		}
		out[xi] = sum / (in.y - 1);
	}
}

/**
 * Array subtraction
 */
void subtract(Float1dArray &x, Float1dArray &y, Float1dArray &difference)
{
	for (int i = 0; i < x.x; ++i)
	{
		difference[i] = x[i] - y[i];
	}
}


/**
 * Array subtraction
 */
void subtract(Double1dArray &x, Double1dArray &y, Double1dArray &difference)
{
	for (int i = 0; i < x.x; ++i)
	{
		difference[i] = x[i] - y[i];
	}
}


/*
 * svdcomp - SVD decomposition routine.
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
 */

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))

static double PYTHAG(double a, double b) {
	double at = fabs(a), bt = fabs(b), ct, result;

	if (at > bt) {
		ct = bt / at;
		result = at * sqrt(1.0 + ct * ct);
	} else if (bt > 0.0) {
		ct = at / bt;
		result = bt * sqrt(1.0 + ct * ct);
	} else
		result = 0.0;
	return (result);
}

int dsvd(Double2dArray &a, Double1dArray &w, Double2dArray &v) {
	int flag, i, its, j, jj, k, l, nm;
	double c, f, h, s, x, y, z;
	double anorm = 0.0, g = 0.0, scale = 0.0;
	double *rv1;
	int m = a.x;
	int n = a.y;

	if (m < n) {
		fprintf(stderr, "#rows must be > #cols \n");
		return (0);
	}

	rv1 = (double *) malloc((unsigned int) n * sizeof(double));

	/* Householder reduction to bidiagonal form */
	for (i = 0; i < n; i++) {
		/* left-hand reduction */
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i < m) {
			for (k = i; k < m; k++)
				scale += fabs((double) a(k, i));
			if (scale) {
				for (k = i; k < m; k++) {
					a(k, i) = (double) ((double) a(k, i) / scale);
					s += ((double) a(k, i) * (double) a(k, i));
				}
				f = (double) a(i, i);
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a(i, i) = (double) (f - g);
				if (i != n - 1) {
					for (j = l; j < n; j++) {
						for (s = 0.0, k = i; k < m; k++)
							s += ((double) a(k, i) * (double) a(k, j));
						f = s / h;
						for (k = i; k < m; k++)
							a(k, j) += (double) (f * (double) a(k, i));
					}
				}
				for (k = i; k < m; k++)
					a(k, i) = (double) ((double) a(k, i) * scale);
			}
		}
		w[i] = (double) (scale * g);

		/* right-hand reduction */
		g = s = scale = 0.0;
		if (i < m && i != n - 1) {
			for (k = l; k < n; k++)
				scale += fabs((double) a(i, k));
			if (scale) {
				for (k = l; k < n; k++) {
					a(i, k) = (double) ((double) a(i, k) / scale);
					s += ((double) a(i, k) * (double) a(i, k));
				}
				f = (double) a(i, l);
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a(i, l) = (double) (f - g);
				for (k = l; k < n; k++)
					rv1[k] = (double) a(i, k) / h;
				if (i != m - 1) {
					for (j = l; j < m; j++) {
						for (s = 0.0, k = l; k < n; k++)
							s += ((double) a(j, k) * (double) a(i, k));
						for (k = l; k < n; k++)
							a(j, k) += (double) (s * rv1[k]);
					}
				}
				for (k = l; k < n; k++)
					a(i, k) = (double) ((double) a(i, k) * scale);
			}
		}
		anorm = fmax(anorm, (fabs((double) w[i]) + fabs(rv1[i])));
	}

	/* accumulate the right-hand transformation */
	for (i = n - 1; i >= 0; i--) {
		if (i < n - 1) {
			if (g) {
				for (j = l; j < n; j++)
					v(j, i) = (double) (((double) a(i, j) / (double) a(i, l))
							/ g);
				/* double division to avoid underflow */
				for (j = l; j < n; j++) {
					for (s = 0.0, k = l; k < n; k++)
						s += ((double) a(i, k) * (double) v(k, j));
					for (k = l; k < n; k++)
						v(k, j) += (double) (s * (double) v(k, i));
				}
			}
			for (j = l; j < n; j++)
				v(i, j) = v(j, i) = 0.0;
		}
		v(i, i) = 1.0;
		g = rv1[i];
		l = i;
	}

	/* accumulate the left-hand transformation */
	for (i = n - 1; i >= 0; i--) {
		l = i + 1;
		g = (double) w[i];
		if (i < n - 1)
			for (j = l; j < n; j++)
				a(i, j) = 0.0;
		if (g) {
			g = 1.0 / g;
			if (i != n - 1) {
				for (j = l; j < n; j++) {
					for (s = 0.0, k = l; k < m; k++)
						s += ((double) a(k, i) * (double) a(k, j));
					f = (s / (double) a(i, i)) * g;
					for (k = i; k < m; k++)
						a(k, j) += (double) (f * (double) a(k, i));
				}
			}
			for (j = i; j < m; j++)
				a(j, i) = (double) ((double) a(j, i) * g);
		} else {
			for (j = i; j < m; j++)
				a(j, i) = 0.0;
		}
		++a(i, i);
	}

	/* diagonalize the bidiagonal form */
	for (k = n - 1; k >= 0; k--) { /* loop over singular values */
		for (its = 0; its < 30; its++) { /* loop over allowed iterations */
			flag = 1;
			for (l = k; l >= 0; l--) { /* test for splitting */
				nm = l - 1;
				if (fabs(rv1[l]) + anorm == anorm) {
					flag = 0;
					break;
				}
				if (fabs((double) w[nm]) + anorm == anorm)
					break;
			}
			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * rv1[i];
					if (fabs(f) + anorm != anorm) {
						g = (double) w[i];
						h = PYTHAG(f, g);
						w[i] = (double) h;
						h = 1.0 / h;
						c = g * h;
						s = (-f * h);
						for (j = 0; j < m; j++) {
							y = (double) a(j, nm);
							z = (double) a(j, i);
							a(j, nm) = (double) (y * c + z * s);
							a(j, i) = (double) (z * c - y * s);
						}
					}
				}
			}
			z = (double) w[k];
			if (l == k) { /* convergence */
				if (z < 0.0) { /* make singular value nonnegative */
					w[k] = (double) (-z);
					for (j = 0; j < n; j++)
						v(j, k) = (-v(j, k));
				}
				break;
			}
			if (its >= 30) {
				free((void*) rv1);
				fprintf(stderr, "No convergence after 30,000! iterations \n");
				return (0);
			}

			/* shift from bottom 2 x 2 minor */
			x = (double) w[l];
			nm = k - 1;
			y = (double) w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = PYTHAG(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f)))- h)) / x;

			/* next QR transformation */
			c = s = 1.0;
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = (double) w[i];
				h = s * g;
				g = c * g;
				z = PYTHAG(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for (jj = 0; jj < n; jj++) {
					x = (double) v(jj, j);
					z = (double) v(jj, i);
					v(jj, j) = (double) (x * c + z * s);
					v(jj, i) = (double) (z * c - x * s);
				}
				z = PYTHAG(f, h);
				w[j] = (double) z;
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = (c * g) + (s * y);
				x = (c * y) - (s * g);
				for (jj = 0; jj < m; jj++) {
					y = (double) a(jj, j);
					z = (double) a(jj, i);
					a(jj, j) = (double) (y * c + z * s);
					a(jj, i) = (double) (z * c - y * s);
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = (double) x;
		}
	}
	free((void*) rv1);
	return (1);
}

/** END svdcomp */

/**
 * Use SVD for pseudo inverse of matrix a
 */
void pinv(Double2dArray &a, Double2dArray &a_pinv) {

	// Matlab code
	//	[m,n]=size(A);
	//	[U,S,V]=svd(A);
	Double1dArray s(a.x);
	Double2dArray v(a.x, a.y);

	// SVD
	if (dsvd(a, s, v) != 1)
		printf("SVD failed!!!!\n");

	// Take the reciprocal of s and convert to 2D
	Double2dArray src(a.x, a.y);
	memset(src.data, 0, src.size());
	for (int i = 0; i < s.x; ++i) {
		if (s[i] == 0)
			src(i, i) = 0;
		else
			src(i, i) = 1/s[i];
	}

	// Transpose u
	Double2dArray u_t(a.x, a.y);
	transpose(a, u_t);

	// A_pseu=V*SRc*U'
	multiplyMatrix(v, src, a);
	multiplyMatrix(a, u_t, a_pinv);

}
