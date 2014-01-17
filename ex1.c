/*
 * f(x) + int [0,pi] cos(x) * cos(y) * f(y) dy = (1 + pi/2) * cos(x);
 *
 * f(x) = cos(x)
 *
 */

#include <math.h>

double
getLambda() {
	return 1.0;
}

double
getUpper() {
	return M_PI;
}

double
getLower() {
	return 0.0;
}

double
g( double x ) {
	return ( 1 + M_PI / 2.0 ) * cos( x );
}

double
K( double x, double y ) {
	return cos( x ) * cos( y );
}

double
f( double x ) {
	return cos( x );
}

int
fredholm() {
    return 2;
}