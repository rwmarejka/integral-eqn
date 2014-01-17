/*
 *	f(x) - int [0,1] exp( x * y ) * f(y) dy = x ^ 4
 *
 *	f(x) = x ^ 4 + 1/10 * tau( x )
 */

#include <math.h>

double
getLower() {
	return 0.0;
}

double
getUpper() {
	return 1.0;
}

double
getLambda() {
	return -0.1;
}

double
g( double x ) {
	double	x2	= x * x;
	return x2 * x2;
}

double
K( double x, double y ) {
	return exp( x * y );
}

static double
tau( double x ) {
	return
	( ( ( ( 0.005017 ) * x + 0.022725 ) * x + 0.078719 ) * x + 0.187227 ) * x + 0.236282;
}

double
f( double x ) {
	double	x2	= x * x;

	return x2 * x2 + 0.1 * tau( x );
}

int
fredholm() {
    return 2;
}
