/*
 * f(x) - int[0,1] |x-y| * f(y) dy = -( 2 * x^3 - 9 * x + 2 ) / 6
 *
 * f(x) = x
 */

#include <math.h>

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
	return 1.0 / M_PI;
}

double
g( double x ) {
	return +1.0;
}

double
K( double x, double y ) {
	double	diff	= x - y;

	return 1.0 + diff * diff;;
}

extern double ChebyshevEval( double, double [], int, double, double );

double
f( double x ) {
	static double a[15];;

	a[0]	= 0.708;
	a[1]	= 0.0;
	a[2]	= 0.049;

	return ChebyshevEval( x, a, 3, getLower(), getUpper() );
}

int
fredholm() {
    return 2;
}
