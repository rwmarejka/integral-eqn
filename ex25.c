/*
 * Example - The Classical Theory of Integral Equations: A Concise Treatment
 *             Chapter 2 -  Feedholm Intergal Equations of the Second Kind (General Kernel), pg 43
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
	return -0.1;
}

double
g( double x ) {
	return x * x * x * x;
}

double
K( double x, double y ) {
    return exp( x * y );
}

double
f( double x ) { 
    double a[] = {
        +0.236282,
        +0.187227,
        +0.078719,
        +0.022725,
        +0.005017,
    };
    
    return x * x * x * x + 0.1 * PowerEval( x, a, sizeof(a)/sizeof(a[0]) );
}

int
fredholm() {
    return 2;
}
