/*
% cc -g -Xc -o int int.c -lm
 *
 * INTEGRAL - generate an approximate solution to a Fredholm integral
 *		equation of the second kind using Chebyshev polynominals.
 *
 * written: rwmarejka@mac.com
 */

/* Feature Test Macros		*/

#define	_POSIX_SOURCE

/* Include Files            */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/* Constants and Macros		*/

#define NTERMS  50

#if !defined(M_PI)
#	define	M_PI	(3.1415926535897932384626433)
#endif

/*
 * XK - translate T sub n ( x sub k ) to the range (a,b) where x sub k = cos ( k * pi / n )
 */

#define XK(k,n,a,b)     (0.5*((b)+(a)+((b)-(a))*cos((k)*M_PI/(n))))

/* External References		*/

extern	void	Chebyshev2Coeff( double [NTERMS][NTERMS], unsigned, double, double, double (*)( double, double ) );
extern  double  Chebyshev2Eval( double, double, double [NTERMS][NTERMS], unsigned, double, double );

/* External Declarations	*/

double  k( double x, double y ) {
#if 0
    return cos( x * y );
#else
    return exp( x );

    return 1                                                                        // T0(x) * T0(y)
            + y                                                                     // T0(x) * T1(y)
            + x                                                                     // T1(x) * T0(y)
            + x * y                                                                 // T1(x) * T1(y)
            + x                 * ( 2 * y * y - 1 )                                 // T1(x) * T2(y)
            + ( 2 * x * x - 1 ) * ( 2 * y * y - 1 )                                 // T2(x) * T2(y)
            + ( 2 * x * x - 1 ) * ( 4 * y * y * y - 3 * y )                         // T2(x) * T3(y)
            + ( 4 * x * x * x - 3 * x )                                             // T3(x) * T0(y)
            + ( 8 * x * x * x * x - 8 * x * x + 1 )                                 // T4(x) * T0(y)
            + ( 8 * x * x * x * x - 8 * x * x + 1 ) * y                             // T4(x) * T1(y)
            + ( 8 * x * x * x * x - 8 * x * x + 1 ) * ( 2 * y * y - 1 )             // T4(x) * T2(y)
            + ( 8 * x * x * x * x - 8 * x * x + 1 ) * ( 4 * y * y * y - 3 * y )     // T4(x) * T3(y)
        ;
//    return sin( x ) * cos( y );
//    return cos( 10.0 * x * y * y ) + exp( -( x * x ) );
#endif
}

/*	main - kiss of life		*/

int
main( int argc, char *argv[] ) {
    double  a[NTERMS][NTERMS];
    double  lower   = -1.0;
    double  upper   = +1.0;
    unsigned    n   = 6;

    Chebyshev2Coeff( a, n, lower, upper, k );

    {
        unsigned    i   = 0;

        for ( ; i < n; ++i ) {
            unsigned    j = 0;
            
            for ( ; j < n; ++j )
                printf( " %+11.9lf", a[i][j] );
            
            putchar( '\n' );
        }
        putchar( '\n' );
    }
    
    {
        double  x   = 0.0;
        double  dx  = 0.1;
        
        for ( x = 0.0; x <= 1.0; x += dx ) {
            double  y   = k( x, x );
            double  y1  = Chebyshev2Eval( x, x, a, n, lower, upper );
            double  dy  = y - y1;
            
            printf( "\t%.2lf  %+.9lf  %+.9lf  %+.9lf\n", x, y, y1, dy );
        }
    }
}

/*    ChebyshevEval - evaluate a Chebyshev series at a point        */

double ChebyshevEval( double x, double a[NTERMS], unsigned n, double lower, double upper ) {
    double    c0    = 0.0;
    double    c1    = 0.0;
    double    c2    = 0.0;
    double    mul   = 2.0 * ( 2.0 * x - ( upper + lower ) ) / ( upper - lower );
    int       degree    = n - 1;

    assert( n > 0 );

//    for ( ; degree >= 0; --degree ) {
    while ( n-- ) {
        c2    = c1;
        c1    = c0;
        c0    = mul * c1 - c2 + a[n];
    }

//    return 0.5 * ( c0 - c2 );
    return( 0.5 * ( c0 - c2 + a[0] ) );
}

/*	Chebyshev2Coeff - generate Chebyshev coefficients to K(x,y)	*/

void Chebyshev2Coeff( double a[NTERMS][NTERMS], unsigned n, double lower, double upper, double (*f)( double, double ) ) {
	unsigned    i      = 0;
	unsigned    np1    = n--;
    double      factor = 4.0 / ( n * n );
    double      X[NTERMS];
    double      K[NTERMS][NTERMS];

	for ( ; i <= n; ++i )               /* compute the abscissa     */
		X[i] = XK( i, n, lower, upper );

	for ( i=0; i <= n; ++i ) {			/* compute K(x,y)           */
		unsigned j  = 0;
        double   Xi = X[i];

		for ( ; j <= n; ++j )
			K[i][j]	= (*f)( Xi, X[j] );

        K[i][0] *= 0.5;
		K[i][n]	*= 0.5;
	}

	for ( i=0; i <= n; ++i ) {          /* compute K matrix         */
		unsigned j  = 0;
		double	 Xi	= X[i];

		for ( ; j <= n; ++j ) {
			unsigned k  = 0;
			double	 b[NTERMS];
            double   Xj = X[j];

			for ( ; k <= n; ++k )
				b[k]	 = ChebyshevEval( Xj, K[k], np1, lower, upper );

            b[0]    *= 0.5;             /* half first and last term         */
			b[n]	*= 0.5;
			a[i][j]	 = factor * ChebyshevEval( Xi, b, np1, lower, upper );
		}
	}
#if 1
	for ( i=0; i <= n; ++i ) {
		a[0][i]	*= 0.5;
		a[i][0]	*= 0.5;
		a[n][i]	*= 0.5;
		a[i][n]	*= 0.5;
	}
#else
    a[0][0] *= 0.25;
    
    for ( i=1; i <= n; ++i )
        a[1][i] *= 0.5;

    for ( i=2; i <= n; ++i )
        a[i][0] *= 0.5;
#endif
	return;
}

/*	Chebyshev2Eval - evaluate a two-dimensional Chebyshev series at a point		*/

double Chebyshev2Eval( double x, double y, double K[NTERMS][NTERMS], unsigned n, double lower, double upper ) {
    unsigned    i   = 0;
    double      b[NTERMS];

    for ( ; i < n; ++i )
        b[i]     = ChebyshevEval( y, K[i], n, lower, upper );

    return ChebyshevEval( x, b, n, lower, upper );
}
