/*
 * Example - Love's equation
 *
 *  f(x) - 1/pi int [-1,+1] f(y) * d / ( d ^ 2 + ( x - y ) ^2 ) dy = 1
 *
 * http://www.mathematica-journal.com/issue/v8i4/tricks/contents/Tricks8-4_2.html
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
	return -1.0 / M_PI;
}

double
g( double x ) {
	return +1.0;
}

double
K( double x, double y ) {
    return d / ( d * d + ( x - y ) * ( x - y ) );
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

    return ChebyshevEval( x, a, N, getLower(), getUpper() );
}

int
fredholm() {
    return 2;
}
