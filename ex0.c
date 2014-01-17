/*
 * f(x) + int [-1,1] f(y) dy = 3
 *
 * f(x) = 1
 *
 */

#include <math.h>

double
getLambda() {
	return 1.0;
}

double
getUpper() {
	return 1.0;
}

double
getLower() {
	return -1.0;
}

double
g( double x ) {
	return 3.0;
}

double
K( double x, double y ) {
	return 1.0;
}

double
f( double x ) {
	return 1.0;
}

int
fredholm() {
    return 2;
}