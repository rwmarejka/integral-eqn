/*
 * Example - The Classical Theory of Integral Equations: A Concise Treatment
 *             Chapter 2 -  Feedholm Intergal Equations of the Second Kind (General Kernel), pg 77
 *
 */

#include <math.h>

extern double PowerEval( double, double [], unsigned );

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
	return sin( M_PI * x );
}

double
K( double x, double y ) {
	return +1.0 / ( 1.0 + x + y );
}

double
f( double x ) {
    double a[] = {
        +0.297340,
        +2.932409,
        +0.124667,
        -5.092254,
        -0.557749,
        +3.710938,
        -1.240580
    };
    
    return PowerEval( x, a, sizeof(a[0]) );
}

int
fredholm() {
    return 2;
}
