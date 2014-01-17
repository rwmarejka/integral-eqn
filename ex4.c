/*
 *
 * Love's equation
 *
 *	f(x) + 1/pi int [-1,+1] f(y) * d / ( d ^ 2 + ( x - y ) ^2 ) dy = 1
 *
 * http://www.mathematica-journal.com/issue/v8i4/tricks/contents/Tricks8-4_2.html
 */

#include <math.h>

#define	DISTANCE	(+1.0)

double
getLambda() {
	return -1.0 / M_PI;
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
	return +1.0;
}

double
K( double x, double y ) {
	double	delta	= x - y;

	return( DISTANCE / ( DISTANCE * DISTANCE + delta * delta ) );
}

static double
T0( double x ) {
	return 1.0;
}

static double
T2( double x ) {
	return 2.0 * x * x - 1.0;
}

static double
T4( double x ) {
	double	x2	= x * x;

	return 8.0 * x2 * x2 - 8.0 * x2 + 1.0;
}

static double
T6( double x ) {
	double	x2	= x * x;

	return 32.0 * x2 * x2 * x2 - 48.0 * x2 * x2 + 18.0 * x2  - 1.0;
}

static double
T8( double x ) {
	double	x2	= x * x;

	return 128.0 * x2 * x2 * x2 * x2 - 256.0 * x2 * x2 * x2 + 160.0 * x2 * x2 - 32.0 * x2 + 1.0;
}


/*
 * Note: f(x) only valid at DISTANCE = 1
 */

double
f( double x ) {
	return    1.774444805      * T0( x )
		- 0.1400430745     * T2( x )
		+ 0.004962498177   * T4( x )
		+ 0.0003746723341  * T6( x )
		+ 0.00004368650886 * T8( x );
}

int
fredholm() {
    return 2;
}
