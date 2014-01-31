/*
 * QA Test 0 - everything is 1, except g(x)=3
 *
 */

#include <math.h>

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
	return +1.0;
}

double
g( double x ) {
	return +3.0;
}

double
K( double x, double y ) {
	return +1.0;
}

double
f( double x ) {
	return +1.0;
}

int
fredholm() {
    return 2;
}