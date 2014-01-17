/*
 * Love's equation
 *
 *	f(x) - d/pi * int [-1,+1] f(y) / ( d^2 + ( x - y )^2 ) dy = atan( 1 + x ) + atan( 1 - x )
 *
 *	f(x) = 1
 *
 * http://www.math.oregonstate.edu/~math_reu/proceedings/REU_Proceedings/Proceedings1993/1993Pilkington.pdf
 */

#include <math.h>

static double D = -1.0;

double
getLambda() {
	return -D / M_PI;
}

double
getUpper() {
	return +1.0;
}

double
getLower() {
	return -1.0;
}

double
g( double x ) {
	return 1.0 + ( atan( 1.0 + x ) + atan( 1.0 - x ) ) / M_PI;
}

double
K( double x, double y ) {
	return 1.0 / ( D * D + ( x - y ) * ( x - y ) );
}

double
f( double x ) {
	return 1.0;
}

int
fredholm() {
    return 2;
}
