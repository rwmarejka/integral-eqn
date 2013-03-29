/*
 * f(x) - int [0,2*pi] sin( x + 2 * y ) * f(y) dy = x
 *
 * f(x) = x - pi * cos(x) 
 */

#include <math.h>

double
getLambda() {
	return -1.0;
}

double
getUpper() {
	return 2.0 * M_PI;
}

double
getLower() {
	return 0.0;
}

double
g( double x ) {
	return x;
}

double
K( double x, double y ) {
	return sin( x + 2.0 * y );;
}

double
f( double x ) {
	return x - M_PI * cos( x );
}
