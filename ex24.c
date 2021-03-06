/*
 * Example - The Classical Theory of Integral Equations: A Concise Treatment
 *             Chapter 2 -  Feedholm Intergal Equations of the Second Kind (General Kernel), pg 42
 *
 */

#include <math.h>

double
getLower() {
	return 0.0;
}

double
getUpper() {
	return +2.0 * M_PI;
}

double
getLambda() {
	return -1.0;
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

int
fredholm() {
    return 2;
}