/*
 * Example - The Classical Theory of Integral Equations: A Concise Treatment
 *             Chapter 2 -  Feedholm Intergal Equations of the Second Kind (General Kernel), pg 71
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
	return cos( x );
}

double
K( double x, double y ) {
	return cos( x * y );
}

double
f( double x ) {  
    double a[] = {
        +(   11532090.0 /        6397711.0 ),
        +0.0,
        -(    7944195.0 /       12795422.0 ),
        +0.0,
        +(     607005.0 /       12795422.0 ),
        +0.0,
        +( 1620362251.0 /  1064067293520.0 ),
        +0.0,
        +( 1473509027.0 / 55331499263040.0 )
    };
    
    return PowerEval( x, a, sizeof(a)/sizeof(a[0]) );
}

int
fredholm() {
    return 2;
}
