/*
 * Love's equation
 *
 *	f(x) -1 int [-1,+1] f(y) / ( 1 + ( x - y )^2 ) dy = -atan( x - 1 ) + atan( x + 1 ) + 1
 *
 *	f(x) = 1
 */

#include <math.h>

double
getLambda() {
	return +1.0;
}

double
getUpper() {
	return +1.0;
}

double
getLower() {
	return -1.0;
}

double
g( double x ) {
	return -atan( x - 1.0 ) + atan( x + 1.0 ) + 1.0;
}

double
K( double x, double y ) {
	double	d	= x - y;

	return 1.0 / ( 1.0 + d * d );
}

double
f( double x ) {
	return 1.0;
}
