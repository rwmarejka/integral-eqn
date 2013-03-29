/*
 *	f(x) - (1/2) * int [-1,+1] sin( pi * x * y ) * f(y) dy = x ^ 3
 *
 *	f(x) = a5 * x^5 + a3 * x^3 + a1 * x
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
	return -0.5;
}

double
g( double x ) {
	return x * x * x;
}

double
K( double x, double y ) {
	return sin( M_PI * x * y );
}

double
f( double x ) {
	double	x2	= x * x;

	return ( ( ( 0.014041 ) * x2 + 0.847750 ) * x2 + 0.565420 ) * x;
}
