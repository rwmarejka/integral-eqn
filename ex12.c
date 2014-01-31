/*
 * Example - Numerical Soution of Intergal Equations, pg 98
 *
 * f(x) - (1/4) int [0,pi/2] x y f(y) dy = sin( x ) - x / 4
 *
 * f(x) = sin( x )
 */

#include <math.h>

double
getLower() {
	return 0.0;
}

double
getUpper() {
	return M_PI / 2.0;
}

double
getLambda() {
	return -0.25;
}

double
g( double x ) {
	return sin( x ) - 0.25 * x;
}

double
K( double x, double y ) {
    return x * y;
}

double
f( double x ) {
    return sin( x );
}

int
fredholm() {
    return 2;
}
