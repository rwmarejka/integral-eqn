/*
 * Example 21 - The Mathematica Journal 9:2 2004, Integral Equations by Stan Richardson
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
	return -0.5;
}

double
g( double x ) {
	return exp( -x ) - 0.5 + 0.5 * exp( -( x + 1.0 ) );
}

double
K( double x, double y ) {
	return ( x + 1.0 ) * exp( - x * y );
}

double
f( double x ) {
	return exp( -x );
}

int
fredholm() {
    return 2;
}