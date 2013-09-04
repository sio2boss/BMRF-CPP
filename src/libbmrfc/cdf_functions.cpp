// tcdf and dependent functions taken from:
// 	http://www-stat.stanford.edu/~serban/gxna/src/cdf.cpp

#include "cdf_functions.h"
#include <math.h>

#define  A1  (-3.969683028665376e+01)
#define  A2   2.209460984245205e+02
#define  A3  (-2.759285104469687e+02)
#define  A4   1.383577518672690e+02
#define  A5  (-3.066479806614716e+01)
#define  A6   2.506628277459239e+00

#define  B1  (-5.447609879822406e+01)
#define  B2   1.615858368580409e+02
#define  B3  (-1.556989798598866e+02)
#define  B4   6.680131188771972e+01
#define  B5  (-1.328068155288572e+01)

#define  C1  (-7.784894002430293e-03)
#define  C2  (-3.223964580411365e-01)
#define  C3  (-2.400758277161838e+00)
#define  C4  (-2.549732539343734e+00)
#define  C5   4.374664141464968e+00
#define  C6   2.938163982698783e+00

#define  D1   7.784695709041462e-03
#define  D2   3.224671290700398e-01
#define  D3   2.445134137142996e+00
#define  D4   3.754408661907416e+00

#define P_LOW   0.02425
/* P_high = 1 - p_low*/
#define P_HIGH  0.97575

long double normsinv(long double p)
{
	long double x;
	long double q, r, u, e;
	if ((0.0f < p) && (p < P_LOW))
	{
		q = sqrt(-2.0f * log(p));
		x = (((((C1 * q + C2)*q+C3)*q+C4)*q+C5)*q+C6)
				/ ((((D1 * q + D2) * q + D3) * q + D4) * q + 1);
	}
	else
	{
		if ((P_LOW <= p) && (p <= P_HIGH))
		{
			q = p - 0.5f;
			r = q * q;
			x = (((((A1 * r + A2) * r + A3)*r+A4)*r+A5)*r+A6) * q
					/ (((((B1 * r + B2) * r + B3)*r+B4)*r+B5)*r+1);
		}
		else
		{
			if ((P_HIGH < p) && (p < 1))
			{
				q = sqrt(-2.0f * log(1.0f - p));
				x = -(((((C1 * q + C2)*q+C3)*q+C4)*q+C5)*q+C6)
						/ ((((D1 * q + D2) * q + D3) * q + D4) * q + 1.0f);
			}
		}
	}

	// If you are compiling this under UNIX OR LINUX, you may uncomment this block for better accuracy.
	if ((0.0f < p) && (p < 1.0f)) {
		e = 0.5f * erfc(double(-x / sqrt(2.0f))) - p;
		u = e * sqrt(2 * M_PI) * exp(x * x / 2.0f);
		x = x - u / (1.0f + x * u / 2);
	}

	return x;
}

const double SQRT2PI = 2.5066282746;
float logGamma(double x)
{
	const double c[8] =
	{ 676.5203681218851, -1259.1392167224028, 771.32342877765313,
			-176.61502916214059, 12.507343278686905, -0.13857109526572012,
			9.9843695780195716e-6, 1.5056327351493116e-7 };
	double sum = 0.99999999999980993;
	double y = x;
	for (int j = 0; j < 8; j++)
		sum += c[j] / ++y;
	return log(SQRT2PI * sum / x) - (x + 7.5) + (x + 0.5) * log(x + 7.5);
}

/**
 * helper function for incomplete beta
 * computes continued fraction
 * source: Numerical Recipes in C
 */
double betaContFrac(double a, double b, double x)
{
	const int MAXIT = 1000;
	const double EPS = 3e-7;
	const double FPMIN = 1e-30;
	double qab = a + b;
	double qap = a + 1;
	double qam = a - 1;
	double c = 1;
	double d = 1 - qab * x / qap;
	if (fabs(d) < FPMIN)
		d = FPMIN;
	d = 1 / d;
	double h = d;
	int m;
	for (m = 1; m <= MAXIT; m++)
	{
		int m2 = 2 * m;
		double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
		d = 1 + aa * d;
		if (fabs(d) < FPMIN)
			d = FPMIN;
		c = 1 + aa / c;
		if (fabs(c) < FPMIN)
			c = FPMIN;
		d = 1 / d;
		h *= (d * c);
		aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
		d = 1 + aa * d;
		if (fabs(d) < FPMIN)
			d = FPMIN;
		c = 1 + aa / c;
		if (fabs(c) < FPMIN)
			c = FPMIN;
		d = 1 / d;
		double del = d * c;
		h *= del;
		if (fabs(del - 1) < EPS)
			break;
	}
	if (m > MAXIT)
	{
		printf("betaContFrac: too many iterations\n");
	}
	return h;
}

/**
 *
 */
double betaInc(double a, double b, double x)
{
	if (x == 0)
		return 0;
	else if (x == 1)
		return 1;
	else
	{
		double logBeta = logGamma(a + b) - logGamma(a) - logGamma(b) + a * log(x)
				+ b * log(1 - x);
		if (x < (a + 1) / (a + b + 2))
			return exp(logBeta) * betaContFrac(a, b, x) / a;
		else
			return 1 - exp(logBeta) * betaContFrac(b, a, 1 - x) / b;
	}
}

/**
 * T-test CDF
 */
void tcdf(Float1dArray& in, double dfe, Double1dArray& out)
{
	for (int i = 0; i < in.x; ++i)
	{
		out[i] = 0.5 * betaInc(dfe / 2, 0.5, dfe / (dfe + pow(double(in[i]), 2)));
	}
}


/**
 * Hardcoded 100 value CDF
 */
void make_cdf(Uint1dArray &freq, Float1dArray &cdf_array, const unsigned int &bootstraps) {

	Float1dArray confidence_obs(freq.x);
	for (int i=0; i<freq.x; ++i)
		confidence_obs[i] = float(freq[i])/float(bootstraps);

	for (int i=0; i < cdf_array.x; ++i)
	{
		unsigned int prob_sum = 0;
		float prob = float(i+1) / float(cdf_array.x);
		for (int j=0; j < confidence_obs.x; ++j)
		{
			if (confidence_obs[j] >= prob)
				prob_sum++;
		}
		cdf_array[i] = float(prob_sum) / float(confidence_obs.x);
	}
}

struct Erf {
	static const int ncof=28;

	inline double erf(double x) {
		if (x >=0.) return 1.0 - erfccheb(x);
		else return erfccheb(-x) - 1.0;
	}

	inline double erfc(long double x) {
		if (x >= 0.) return erfccheb(x);
		else return 2.0 - erfccheb(-x);
	}

	double erfccheb(double z){
		int j;
		double t,ty,tmp,d=0.,dd=0.;
		if (z < 0.) throw("erfccheb requires nonnegative argument");
		t = 2./(2.+z);
		ty = 4.*t - 2.;
		for (j=ncof-1;j>0;j--) {
			tmp = d;
			d = ty*d - dd + cof[j];
			dd = tmp;
		}
		return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
	}

	double inverfc(double p) {
		double x,err,t,pp;
		if (p >= 2.0) return -100.;
		if (p <= 0.0) return 100.;
		pp = (p < 1.0)? p : 2. - p;
		t = sqrt(-2.*log(pp/2.));
		x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
		for (int j=0;j<2;j++) {
			err = erfc(x) - pp;
			x += err/(1.12837916709551257*exp(-1.0f*x*x)-x*err);
		}
		return (p < 1.0? x : -x);
	}

	inline double inverf(double p) {return inverfc(1.-p);}

};

extern "C" const double cof[] = {-1.3026537197817094, 6.4196979235649026e-1,
	1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
	3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
	-1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
	6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
	9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
	-1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};

struct Normaldist : Erf {
	double mu, sig;
	Normaldist(double mmu = 0., double ssig = 1.) : mu(mmu), sig(ssig) {
		if (sig <= 0.) throw("bad sig in Normaldist");
	}
	double p(double x) {
		double z = (x-mu)/sig;
		return (0.398942280401432678/sig)*exp(-0.5*z*z);
	}
	double cdf(double x) {
		return 0.5*erfc(-0.707106781186547524*(x-mu)/sig);
	}
	double invcdf(double p) {
		if (p <= 0. || p >= 1.) throw("bad p in Normaldist");
		return -1.41421356237309505*sig*inverfc(2.*p)+mu;
	}
};

/**
 * Inverse CDF http://home.online.no/~pjacklam/notes/invnorm/#C
 */
void icdf_normal(Double1dArray& in, Double1dArray& out, float mu, float sigma)
{
	Normaldist nd(mu, sigma);
	for (int i = 0; i < in.x; ++i)
	{
		out[i] = nd.invcdf(in[i]);
	}
}
