/*
 * f(x) - int[0,1] K(x,y) f(y) dy = ( 1 - 1  / pi^2 ) * sin( pi x )
 *
 * K(x,y) = ( x <= y )? 1 - y : ( 1 - x )
 *
 * f(x) = sin( pi * x )
 */

#include <math.h>

double
getLambda() {
	return -1.0;
}

double
getUpper() {
	return +1.0;
}

double
getLower() {
	return 0.0;
}

double
g( double x ) {
	return ( 1.0 - 1.0 / ( M_PI * M_PI ) ) * sin ( M_PI * x );
}

double
K( double x, double y ) {
	if ( x <= y )
		return x * ( 1.0 - y );
	else /* y <= x */
		return y * ( 1.0 - x );
}

double
f( double x ) {
	return sin( x * M_PI );
}
