/*
 * Example - The Classical Theory of Integral Equations: A Concise Treatment
 *             Chapter 2 -  Feedholm Intergal Equations of the Second Kind (General Kernel), pg 71
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
	return cos( x );
}

double
K( double x, double y ) {
	return cos( x * y );
}

double
f( double x ) {
    double x2 = x * x;
    double x4 = x2 * x2;
    double y;

    y =   (   11532090.0 /        6397711.0 )
        - (    7944195.0 /       12795422.0 ) * x2
        + (     607005.0 /       12795422.0 ) * x4
        + ( 1620362251.0 /  1064067293520.0 ) * x4 * x2
        + ( 1473509027.0 / 55331499263040.0 ) * x4 * x4;

	return y;
}

int
fredholm() {
    return 2;
}