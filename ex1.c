/*
 * QA Test 1 - lambda = 1, range [0,pi], use trig functions
 *
 */

#include <math.h>

double
getLower() {
	return 0.0;
}

double
getUpper() {
	return M_PI;
}

double
getLambda() {
	return +1.0;
}

double
g( double x ) {
	return ( +1.0 + M_PI / 2.0 ) * cos( x );
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