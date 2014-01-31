/*
 * Example - Numerical Analysis Of Love's Integral Equation, Pilkington
 *
 * http://www.math.oregonstate.edu/~math_reu/proceedings/REU_Proceedings/Proceedings1993/1993Pilkington.pdf
 */

#include <math.h>

static double d = -1.0;

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
	return -d / M_PI;
}

double
g( double x ) {
	return +1.0 + ( atan( 1.0 + x ) + atan( 1.0 - x ) ) / M_PI;
}

double
K( double x, double y ) {
	return +1.0 / ( d * d + ( x - y ) * ( x - y ) );
}

double
f( double x ) {
	return +1.0;
}

int
fredholm() {
    return 2;
}
