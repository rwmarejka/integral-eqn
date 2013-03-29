/*
 *
 * f(x) - (1/2) int[0,1] ( 1/( 1 + x + y ) ) * f(y) dy = sin( pi * x )
 *
 * f(x) = sum(j=0,6) a[j] * x^j
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
	return -0.5;
}

double
g( double x ) {
	return sin( M_PI * x );
}

double
K( double x, double y ) {
	return 1.0 / ( 1.0 + x + y );
}

double
f( double x ) {
	return
	( ( ( ( ( ( -1.240580 ) * x + 3.710938 ) * x + -0.557749 ) * x + -5.092254 ) * x + 0.124667 ) * x + 2.932409 ) * x + 0.297340;
}
