/*
 * Example - Numerical Solution Of Non-Singular Fredholm Integral Equations Of
 *              The Second Kind, Sastry, 1973.
 *
 */

#include <math.h>

static double d = +1.0;

double
getLower() {
    return -1.0;
}

double
getUpper() {
	return +1.0;
}

double
getLambda() {
	return +1.0 / M_PI;
}

double
g( double x ) {
	return +1.0;
}

double
K( double x, double s ) {
	return d / ( d * d + ( x - s ) * ( x - s ) );
}

extern	double	ChebyshevEval( double, double [], int, double, double );

double
f( double x ) {
    double a[] = {
        +1.4151850 * 0.5,
        +0.0,
        +0.0493851,
        +0.0,
        -0.0010481,
        +0,0,
        -0.0002310,
        +0.0,
        +0.0000391
    };
#define N   (sizeof(a)/sizeof(a[0]))

    return ChebyshevEval( x, a, N, getLower(), getUpper() );
}

int
fredholm() {
    return 2;
}