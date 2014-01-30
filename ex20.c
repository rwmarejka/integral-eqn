/*
 * Example 6 - The Classical Theory of Integral Equations: A Concise Treatment
 *             Chapter 2 -  Feedholm Intergal Equations of the Second Kind (General Kernel), pg 43
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
	return -0.1;
}

double
g( double x ) {
	return x * x * x * x;
}

double
K( double x, double y ) {
	return
    1.0
    + x             * y
    + x * x         * y * y         / ( 1.0 * 2.0             )
    + x * x * x     * y * y * y     / ( 1.0 * 2.0 * 3.0       )
    + x * x * x * x * y * y * y * y / ( 1.0 * 2.0 * 3.0 * 4.0 );
}

double
f( double x ) {
    double x2 = x * x;
    double x4 = x2 * x2;
    double a[] = {
        +0.236282,
        +0.187227,
        +0.078719,
        +0.022725,
        +0.005017,
    };
#define NTERMS (sizeof(a)/sizeof(a[0]))

    double y = x4 + 0.1 * ( a[0] + a[1] * x + a[2] * x2 + a[3] * x2 * x + a[4] * x4 );

	return y;
}

int
fredholm() {
    return 2;
}