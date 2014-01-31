/*
 * Example 23 - Solving Fredholm Equations of the Second Kind in MATLAB, pg 10
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
	return -1.0;
}

double
g( double x ) {
	return ( 1.0 - 1.0 / ( M_PI * M_PI ) ) * sin( M_PI * x );
}

double
K( double s, double t ) {
    double z;

    if ( s <= t )
        z   = s * ( 1.0 - t );
    else
        z   = t * ( 1.0 - s );

	return z;
}

double
f( double x ) {
	return sin( M_PI * x );
}

int
fredholm() {
    return 2;
}