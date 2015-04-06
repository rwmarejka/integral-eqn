/*
 * Example - Solving Fredholm Equations of the Second Kind in MATLAB, pg 10
 *
 */

#include <math.h>

double
getLower() {
	return 0.0;
}

double
getUpper() {
	return +1.0;
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
	return fabs( y - x );
}

double
f( double x ) {
	return x;
}

int
fredholm() {
    return 2;
}