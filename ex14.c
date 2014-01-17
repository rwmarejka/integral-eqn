/*
 * int [0,1] exp( x * y ) * f(y) dy = ( exp( x + 1 ) - 1 ) / ( x + 1 )
 *
 * f(x) = exp( x )
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
	return 1.0;
}

double
g( double x ) {
	return x * x * x;
}

double
K( double x, double y ) {
    if ( x <= y ) {
        return x * ( 1.0 - y );
    } else {
        return y * ( 1.0 - x );
    }
}

double
f( double x ) {
    return 7.0 * exp( 1.0 ) * ( exp( x ) - exp( -x ) ) / ( exp( 2.0 ) - 1.0 ) - 6.0 * x;
}

int
fredholm() {
    return 2;
}