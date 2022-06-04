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

#include <polynomial.h>
#include <matrix.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/* Constants and Macros		*/

#if !defined(M_PI)
#	define	M_PI	(3.1415926535897932384626433)
#endif

/*
 * XK - translate T sub n ( x sub k ) to the range (a,b) where x sub k = cos ( k * pi / n )
 */

#define XK(k,n,a,b)     (0.5*((b)+(a)+((b)-(a))*cos((k)*M_PI/(n))))

/* External References		*/

extern	void	Chebyshev2Coeff( MATRIX, unsigned, double, double, double (*)( double, double ) );
extern  double  Chebyshev2Eval( double, double, MATRIX, unsigned, double, double );

extern  matrix_t    *ChebyshevTMatrix( unsigned );

extern	double	g( double );
extern	double	K( double, double );
extern	double	f( double x );
extern  int     fredholm( void );

/* External Declarations	*/

/*	main - kiss of life		*/

int
main( int argc, char *argv[] ) {
    polynomial_t    *b;
	vector_t        *x, *vectorb;
	matrix_t        *a, *A, *t;
	unsigned        n;

	if ( argc > 1 )
		n	= atoi( argv[1] );
	else
		n	= N - 1;

	assert( n <= N );
}

polynomial_t    *fredholm( unsigned degree, double lower, double upper, double (*g)( double ), double (*K)( double, double ) ) {
    polynomial_t    *b;
    vector_t        *vectorb;
    vector_t        *x;
    matrix_t        *a;
    matrix_t        *t;
    matrix_t        *A;

    b       = polynomial_fromFunction( Chebyshev, degree, g );
    vectorb = polynomial_getvector( b );
    vector_print( stdout, "g(x) Chebyshev coefficients", vectorb, NULL );

	a       = Chebyshev2Coeff( a, degree, lower, upper, K );
    matrix_print( stdout, "K(x,y) Chebyshev coefficients", a, NULL );

    t   = ChebyshevTMatrix( degree );
    A   = matrix_mul( a, t );

    if ( fredholm() == 2 )
        matrix_add_identity( A );

    matrix_print( stdout, "A matrix", A, NULL );
    
    x   = matrix_solve( A, b );

	if ( x ) {
        vector_print( stdout, "f(x) Chebyshev coefficients", x, NULL ) {
	}

	return 0;
}

/*	Chebyshev2Coeff - generate Chebyshev coefficients to K(x,y)	*/

matrix_t    *Chebyshev2Coeff( unsigned n, double lower, double upper, double (*f)( double, double ) ) {
	unsigned    i      = 0;
	unsigned    np1    = n--;
    double      factor = 4.0 / ( n * n );
    vector_t    *x;
    matrix_t    *a;
    matrix_t    *k;

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
			VECTOR	 b;
            double   Xj = X[j];

			for ( ; k <= n; ++k )
				b[k]	 = ChebyshevEval( Xj, K[k], np1, lower, upper );

            b[0]    *= 0.5;             /* half first and last term         */
			b[n]	*= 0.5;
			a[i][j]	 = factor * ChebyshevEval( Xi, b, np1, lower, upper );
		}
	}

	for ( i=0; i <= n; ++i ) {
		a[0][i]	*= 0.5;
		a[i][0]	*= 0.5;
		a[n][i]	*= 0.5;
		a[i][n]	*= 0.5;
	}

	return a;
}

/*	Chebyshev2Eval - evaluate a two-dimensional Chebyshev series at a point		*/

double Chebyshev2Eval( double x, double y, matrix_t *K, unsigned n, double lower, double upper ) {
    unsigned    i   = 0;
    vector_t   *b;

    for ( ; i < n; ++i )
        b[i]     = ChebyshevEval( y, K[i], n, lower, upper );

    return ChebyshevEval( x, b, n, lower, upper );
}

/* ChebyshevTMatrix - compute the t matrix. */

matrix_t *ChebyshevTMatrix( unsigned n ) {
    unsigned     i  = 0;
    matrix_t    *t  = matrix_alloc( n, n );

	for ( ; i < n; ++i ) {
		unsigned j = ( i % 2 )? 1 : 0;
		unsigned k = ( i % 2 )? 0 : 1;

//		for ( ; k < n; k+=2 )
//			t->m[i][k]	= 0.0;

		for ( ; j < n; j+=2 ) {
			int	diff = i - j;
			int	sum  = i + j;

			t->m[i][j]	= -( 1.0 / ( sum * sum - 1 ) + 1.0 / ( diff * diff - 1 ) );
		}
	}

    return t;
}
