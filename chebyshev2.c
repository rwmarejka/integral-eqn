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

#include "polynomial.h"
#include "matrix.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/* Constants and Macros		*/

#define NTERMS  50

#if !defined(M_PI)
#	define	M_PI	(3.1415926535897932384626433)
#endif

/* External References		*/

extern	matrix_t	*Chebyshev2Coeff( unsigned, double, double, double (*)( double, double ) );
extern  double      Chebyshev2Eval( double, double, matrix_t *, double, double );

/* External Declarations	*/

/* toInterval - convert [ -1, +1 ] -> [ lower, upper ] */

static double toInterval( double t, double lower, double upper ) {
    return 0.5 * ( ( upper - lower ) * t + ( upper + lower ) );
}

/*
 * fromInterval - convert [ lower, upper ] to [ -1, +1 ]
 */

static double fromInterval( double x, double lower, double upper ) {
    return ( 2.0 * x - ( upper + lower ) ) / ( upper - lower );
}

/*    ChebyshevEval - evaluate a Chebyshev series at a point        */

static double ChebyshevEval( double x, double *a, unsigned n, double lower, double upper ) {
    double    c0    = 0.0;
    double    c1    = 0.0;
    double    c2    = 0.0;
    double    mul   = 2.0 * ( 2.0 * x - ( upper + lower ) ) / ( upper - lower );

    assert( n > 0 );

    while ( n-- ) {
        c2    = c1;
        c1    = c0;
        c0    = mul * c1 - c2 + a[n];
    }

    return 0.5 * ( c0 - c2 );
}


/*    Chebyshev2Eval - evaluate a two-dimensional Chebyshev series at a point        */

double Chebyshev2Eval( double x, double y, matrix_t *K, double lower, double upper ) {
    unsigned    i;
    unsigned    n   = K->nrow;
    double      b[NTERMS];

    assert( K->nrow == K->ncol );

    for ( i = 0; i < n; ++i )
        b[i]     = ChebyshevEval( y, K->m[i], n, lower, upper );

    return ChebyshevEval( x, b, n, lower, upper );
}


/*	Chebyshev2Coeff - generate Chebyshev coefficients to K(x,y)	*/

matrix_t *Chebyshev2Coeff( unsigned n, double lower, double upper, double (*f)( double, double ) ) {
	unsigned    i;
	unsigned    np1;
    double      factor;
    double      um;
    vector_t    *X      = vector_alloc( n );
    matrix_t    *K      = matrix_alloc( n, n );
    matrix_t    *a      = matrix_alloc( n, n );
    double     **ap     = matrix_getData( a );
    double     **Kp     = matrix_getData( K );

    np1     = n--;
    factor  = 4.0 / ( n * n );
    um      = M_PI / n;

	for ( i = 0; i <= n; ++i )               /* compute the abscissa     */
        vector_set( X, i, toInterval( cos( i * um ), lower, upper ) );

	for ( i = 0; i <= n; ++i ) {			/* compute K(x,y)           */
		unsigned j;
        double   Xi = vector_get( X, i );

		for ( j = 0; j <= n; ++j )
			Kp[i][j]	= (*f)( Xi, vector_get( X, j ) );
        
        Kp[i][n]    *= 0.5;                 // half the last column in the matrix
	}

	for ( i = 0; i <= n; ++i ) {          /* compute K matrix         */
		unsigned j;
        double   Xi = vector_get( X, i );

		for ( j = 0; j <= n; ++j ) {
			unsigned k;
			double	 b[NTERMS];
            double   Xj = vector_get( X, j );

			for ( k = 0; k <= n; ++k )
				b[k]	 = ChebyshevEval( Xj, Kp[k], np1, lower, upper );

			b[n]        *= 0.5;         // half the last term
			ap[i][j]    = factor * ChebyshevEval( Xi, b, np1, lower, upper );
		}
	}
#if 0
	for ( i = 0; i <= n; ++i ) {
		ap[0][i]	*= 0.5;
		ap[i][0]	*= 0.5;
		ap[n][i]	*= 0.5;
		ap[i][n]	*= 0.5;
	}
#elif 0
    ap[0][0] *= 0.25;
    
    for ( i=1; i <= n; ++i )
        ap[1][i] *= 0.5;

    for ( i=2; i <= n; ++i )
        ap[i][0] *= 0.5;
#endif
    
    vector_free( X );
    matrix_free( K );

	return a;
}

#if defined(TEST)

double  k( double x, double y ) {
//  return cos( x );
//  return exp( y );
  return cos( x ) * exp( y );
//  return cos( 10.0 * x * y * y ) + exp( -( x * x ) );
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
}

/*    main - kiss of life        */

int
main( int argc, char *argv[] ) {
    double  lower   = -1.0;
    double  upper   = +1.0;
    unsigned    n   = 8;
    matrix_t    *a;
    double      **ap;

    a   = Chebyshev2Coeff( n, lower, upper, k );
    ap  = matrix_getData( a );

    matrix_print( stdout, "A[][] =", a, "%+.6lf" );
    putchar( '\n' );
    
    {
        double  x;
        double  dx  = 0.125;
        
        for ( x = -1.0; x <= 1.0; x += dx ) {
            double  xpt = toInterval( x, lower, upper );
            double  y   = k( xpt, xpt );
            double  y1  = Chebyshev2Eval( xpt, xpt, a, lower, upper );
            double  dy  = y - y1;
            
            printf( "\t%+5.3lf  %+.6lf  %+.6lf  %+.6lf\n", x, y, y1, dy );
        }
    }
    
    matrix_free( a );
}

#endif
