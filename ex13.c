/*
 * Example - Numerical Soution of Intergal Equations, pg 104
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
	return +1.0;
}

double
g( double x ) {
	return ( exp( x + 1.0 ) - 1.0 ) / ( x + 1.0 );
}

double
K( double x, double y ) {
	return exp( x * y );
}

double
f( double x ) {
    return exp( x );
}

int
fredholm() {
    return 1;
}