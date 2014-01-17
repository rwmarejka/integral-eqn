/*
 * f(x) - int[0,1] |x-y| * f(y) dy = -( 2 * x^3 - 9 * x + 2 ) / 6
 *
 * f(x) = x
 */

#include <math.h>

double
getLower() {
	return 0.0;
}

double
getUpper() {
	return 1.0;
}

double
getLambda() {
	return -1.0;
}

double
g( double x ) {
	return -( 2.0 * x * x * x - 9.0 * x + 2.0 ) / 6.0;
}

double
K( double x, double y ) {
	return fabs( x - y );
}

double
f( double x ) {
	return x;
}

int
fredholm() {
    return 2;
}