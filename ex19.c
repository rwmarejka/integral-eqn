/*
 * Example - The Classical Theory of Integral Equations: A Concise Treatment
 *             Chapter 2 -  Feedholm Intergal Equations of the Second Kind (General Kernel), pg 81
 *
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
	return sin( 0.5 * M_PI * x * y );
}

double
f( double x ) {
    double x2 = x * x;
    double x4 = x2 * x2;
    double a[] = {
        +0.0,
        +0.565420,
        +0.0,
        +0.847750,
        +0.0,
        +0.014041,
    };
#define NTERMS (sizeof(a)/sizeof(a[0]))

    double y = a[1] * x + a[3] * x2 * x + a[5] * x4 * x;

	return y;
}

int
fredholm() {
    return 2;
}