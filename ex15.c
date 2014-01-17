/*
 * int [0,1] exp( x * y ) * f(y) dy = ( exp( x + 1 ) - 1 ) / ( x + 1 )
 *
 * f(x) = exp( x )
 */

#include <math.h>

static double D = 0.5;

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
	return -D / M_PI;
}

double
g( double x ) {
	return 1.0;
}

double
K( double x, double y ) {
    return 1.0 / ( D * D + ( x - y ) * ( x - y ) );
}

extern	double	ChebyshevEval( double, double [], int, double, double );

double
f( double x ) {
    double a[] = {
        +1.774444805,
        +0.0,
        -0.1400430745,
        +0.0,
        +0.004962498177,
        +0,0,
        +0.0003746723341,
        +0.0,
        -0.00004368650886
    };
#define N   (sizeof(a)/sizeof(a[0]))

    a[0]    *= 2.0;

    return ChebyshevEval( x, a, N, -1.0, +1.0 );
}

int
fredholm() {
    return 2;
}