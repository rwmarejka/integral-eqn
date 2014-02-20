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

#define	TOLERANCE    1.0e-15			/* machine tolerance	*/

#if !defined(M_PI)
#	define	M_PI	(3.1415926535897932384626433)
#endif

#define	N	25                          /* size of matrix and vectors   */
#define	M	 5                          /* number of sample points      */

/*
 * XK - translate T sub n ( x sub k ) to the range (a,b) where x sub k = cos ( k * pi / n )
 */

#define XK(k,n,a,b)     (0.5*((b)+(a)+((b)-(a))*cos((k)*M_PI/(n))))
#define	_exchange(a,b)	{double _T = (a); (a)=(b); (b)=_T;}

/* Data Declarations		*/

typedef double	MATRIX[N][N];
typedef double	VECTOR[N];

typedef struct {
	double	x;
	double	y;
	double	z;
} POINT;

/* External References		*/

extern	void	ChebyshevCoeff( VECTOR, unsigned, double, double, double (*)( double ) );
extern	double	ChebyshevEval( double, VECTOR, unsigned, double, double );
extern	void	Points( POINT [M], unsigned, VECTOR, unsigned, double, double );
extern	int     PointWrite( FILE *, char *, POINT [M], unsigned, double (*)( double ) );
extern	int     VectorWrite( FILE *, char *, VECTOR, unsigned );

extern	void	Chebyshev2Coeff( MATRIX, unsigned, double, double, double (*)( double, double ) );
extern  double  Chebyshev2Eval( double, double, MATRIX, unsigned, double, double );
extern  void    Points2( POINT [M][M], unsigned, MATRIX, unsigned, double, double );
extern  int     Point2Write( FILE *, char *, POINT [M][M], unsigned, double (*)( double, double ) );

extern	void	Integrate( MATRIX, MATRIX, unsigned, double, double, double );

extern	double	MatrixSolve( unsigned, MATRIX, VECTOR, VECTOR );
extern	double	MatrixLUdecomp( unsigned, MATRIX, int [N] );
extern	double	MatrixLUdeterminant( unsigned, MATRIX, int [N] );
extern	void	MatrixLUsolve( unsigned, MATRIX, VECTOR, int [N] );
extern	int     MatrixWrite( FILE *, char *, MATRIX, unsigned, unsigned );

extern	double	getLower( void );
extern	double	getUpper( void );
extern	double	getLambda( void );
extern	double	g( double );
extern	double	K( double, double );
extern	double	f( double x );
extern  int     fredholm( void );

/* External Declarations	*/

/*	main - kiss of life		*/

int
main( int argc, char *argv[] ) {
	VECTOR	 x, b;
	MATRIX	 a, A;
	POINT	 pt[M];
	double	 upper   = getUpper();
	double	 lower   = getLower();
	double	 lambda  = getLambda();
    double   det;
	unsigned n;

	if ( argc > 1 )
		n	= atoi( argv[1] );
	else
		n	= N - 1;

	assert( n <= N );

	fprintf( stdout, "lambda: %9.6f, %9.6f <= x <= %9.6f\n\n", lambda, lower, upper );

	ChebyshevCoeff( b, n, lower, upper, g );
	VectorWrite( stdout, "g(x) Chebyshev coefficients", b, n );

	Chebyshev2Coeff( a, n, lower, upper, K );
	MatrixWrite( stdout, "K(x,y) Chebyshev coefficients", a, n, n );

    if ( 1 ) {
        POINT grid[M][M];

        Points2( grid, M, a, n, lower, upper );
        Point2Write( stdout, "K(x,y)", grid, M, K );
    }

	Integrate( a, A, n, lambda, lower, upper );
	MatrixWrite( stdout, "A matrix", A, n, n );

	det = MatrixSolve( n, A, x, b );

    fprintf( stdout, "\tdet|A| = %9.6f\n\n", det );

	if ( det != 0.0 ) {
		VectorWrite( stdout, "f(x) Chebyshev coefficients", x, n );
		Points( pt, M, x, n, lower, upper );
		PointWrite( stdout, "f(x)", pt, M, f );
	}

	return( 0 );
}

/*	ChebyshevCoeff - generate Chebyshev coefficients to f(x)	*/

	void
ChebyshevCoeff( VECTOR a, unsigned n, double lower, double upper, double (*f)( double ) ) {
    unsigned i      = 0;
    unsigned np1    = n--;
    double   factor = 2.0 / n;
    VECTOR   X;
    VECTOR   Y;

    for ( ; i <= n; ++i ) {             /* compute the abscissa and f(x)    */
        X[i]    = XK( i, n, lower, upper );
        Y[i]    = (*f)( X[i] );
    }

    Y[0]    *= 0.5;                     /* half first and last term         */
    Y[n]    *= 0.5;

    for ( i=0; i <= n; ++i )
        a[i]    = factor * ChebyshevEval( X[i], Y, np1, lower, upper );

    a[0]    *= 0.5;                     /* half first and last term         */
    a[n]    *= 0.5;

	return;
}

/*	Chebyshev2Coeff - generate Chebyshev coefficients to K(x,y)	*/

	void
Chebyshev2Coeff( MATRIX a, unsigned n, double lower, double upper, double (*f)( double, double ) ) {
	unsigned i      = 0;
	unsigned np1    = n--;
    double   factor = 4.0 / ( n * n );
    VECTOR   X;
	MATRIX	 K;

	for ( ; i <= n; ++i )               /* compute the abscissa     */
		X[i] = XK( i, n, lower, upper );

	for ( i=0; i <= n; ++i ) {			/* compute K(x,y)           */
		unsigned j  = 0;
        double   Xi = X[i];

		for ( ; j <= n; ++j )
			K[i][j]	= (*f)( Xi, X[j] );

        K[i][0] *= 0.5;                 /* half first and last term         */
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

	for ( i=0; i <= n; ++i ) {          /* half the first and last terms */
		a[0][i]	*= 0.5;
		a[i][0]	*= 0.5;
		a[n][i]	*= 0.5;
		a[i][n]	*= 0.5;
	}

	return;
}

/*	ChebyshevEval - evaluate a Chebyshev series at a point		*/

	double
ChebyshevEval( double x, VECTOR a, unsigned n, double lower, double upper ) {
	double	c0	= 0.0;
	double	c1	= 0.0;
	double	c2	= 0.0;
	double	mul	= 2.0 * ( 2.0 * x - ( upper + lower ) ) / ( upper - lower );

	while ( n-- ) {
		c2	= c1;
		c1	= c0;
		c0	= mul * c1 - c2 + a[n];
	}

	return( 0.5 * ( c0 - c2 + a[0] ) );
}

/*	Chebyshev2Eval - evaluate a two-dimensional Chebyshev series at a point		*/

double
Chebyshev2Eval( double x, double y, MATRIX K, unsigned n, double lower, double upper ) {
    unsigned i = 0;
    VECTOR   b;

    for ( ; i < n; ++i )
        b[i]     = ChebyshevEval( y, K[i], n, lower, upper );

    return ChebyshevEval( x, b, n, lower, upper );
}

/*	Integrate - computes the matrix [ I + lambda * k * t ]						*/

	void
Integrate( MATRIX b, MATRIX d, unsigned n, double lambda, double lower, double upper ) {
	unsigned i      = 0;
	double	 factor	= lambda * ( upper - lower ) / 2.0;
	MATRIX	 t;

	for ( ; i < n; ++i ) {                      /* compute t                    */
		unsigned j = ( i % 2 )? 1 : 0;
		unsigned k = ( i % 2 )? 0 : 1;

		for ( ; k < n; k+=2 )
			t[i][k]	= 0.0;

		for ( ; j < n; j+=2 ) {
			int	diff = i - j;
			int	sum  = i + j;

			t[i][j]	= -( 1.0 / ( sum * sum - 1 ) + 1.0 / ( diff * diff - 1 ) );
		}
	}

	for ( i=0; i < n; ++i ) {                   /* d = b * t                    */
		unsigned j   = 0;

		for ( ; j < n; ++j ) {
			unsigned k   = 1;
			double	 sum = b[i][0] * t[0][j];

			for ( ; k < n; ++k )
				sum	+= b[i][k] * t[k][j];

			d[i][j]	= factor * sum;
		}

        if ( fredholm() == 2 )
            d[i][i]	+= 1.0;                     /* d += I                       */
	}

	return;
}

/*	Points - evaluate the Chebyshev approximation to f(x) at M points in the interval. */

	void
Points( POINT pt[M], unsigned m, VECTOR a, unsigned n, double lower, double upper ) {
	unsigned i  = 0;
    double   x  = lower;
	double   dx = ( upper - lower ) / ( m - 1 );

	for ( ; i < m; ++i, ++pt ) {
        pt->x  = x;
		pt->y  = ChebyshevEval( x, a, n, lower, upper );
        x 	  += dx;
	}
	return;
}

/*	Points - evaluate the Chebyshev approximation to K(x,y) at M points in the interval. */

    void
Points2( POINT grid[M][M], unsigned m, MATRIX K, unsigned n, double lower, double upper ) {
    unsigned i  = 0;
    double   x  = lower;
    double   dx = ( upper - lower ) / ( m - 1 );

    for ( ; i < m; ++i ) {
        unsigned j = 0;
        double   y = lower;

        for ( ; j < m; ++j ) {
            grid[i][j].x = x;
            grid[i][j].y = y;
            grid[i][j].z = Chebyshev2Eval( x, y, K, n, lower, upper );
            y += dx;
        }
        x += dx;
    }
    return;
}

/*	MatrixSolve - solve a system of linear equations	*/

	double
MatrixSolve( unsigned n, MATRIX a, VECTOR x, VECTOR b ) {
	int     ipvt[N];	/* the pivot vector		*/
	double	status;		/* indicates singularity	*/
/*
 * Decompose a to upper triangular...
 */
	status	= MatrixLUdecomp( n, a, ipvt );

	if ( status != 0.0 ) {		/* check status of system	*/
		int i   = 0;

		for ( ; i < n; ++i )
			x[i]	= b[i];

		MatrixLUsolve( n, a, x, ipvt );

		status	= MatrixLUdeterminant( n, a, ipvt );
	}

	return( status );
}

/*	MatrixLUsolve - solve by back substitution.	*/

	void
MatrixLUsolve( unsigned n, MATRIX a, VECTOR b, int ipvt[N] ) {
	if ( n != 1 ) {			/* if not the trivial case	*/
        int    k  = 0;
        int    i;
        double t;
/*
 * Forward elimination...
 */
		for ( ; k < (n-1); ++k ) {
			int m = ipvt[k];

			_exchange( b[k], b[m] );
			t	= b[k];

			for ( i=k+1; i < n; ++i )
				b[i]	+= a[i][k] * t;
		}
/*
 * Back substitution...
 */
		for ( k=(n-1); k > 0; --k ) {
			b[k] /= a[k][k];
			t     = -b[k];

			for ( i=0; i < k; ++i )
				b[i]	+= a[i][k] * t;
		}
	}

	b[0]	/= a[0][0];			/* the trivial case	*/

	return;
}

/*	MatrixLUdecompDECOMP - reduce matrix to upper triangular.	*/

	double
MatrixLUdecomp( unsigned n, MATRIX a, int ipvt[N] ) {
    double  status;

	ipvt[n-1] = 1;

	if ( n != 1 ) {
        int    i;
        int    k     = 0;
        double t;
        double anorm = 0.0;
        double ynorm = 0.0;
        double znorm = 0.0;
        VECTOR work;
/*
 * Compute 1-norm of a
 */
		for ( ; k < n; ++k ) {
			t	= 0.0;
			for ( i=0; i < n; ++i )
				t	+= fabs( a[i][k] );

			if ( t > anorm )
				anorm	= t;
		}
/*
 * Gaussian elimination with partial pivoting
 */
		for ( k=0; k < (n-1); ++k ) {
            unsigned j;
            unsigned m = k;
/*
 * Find the pivot...
 */
			for ( i=(k+1); i < n; ++i )
				if ( fabs( a[i][k] ) > fabs( a[m][k] ) )
					m = i;

			ipvt[k]	= m;

			if ( m != k )
				ipvt[n-1]	= -ipvt[n-1];

			t	= a[m][k];
			_exchange( a[m][k], a[k][k] );

			if ( t != 0.0 )	/* skip step on zero pivot	*/
/*
 * Compute multipliers
 */
				for ( i=(k+1); i < n; ++i )
					a[i][k]	/= (-t);
/*
 * Interchange and eliminate by columns 
 */
			for ( j=(k+1); j < n; ++j ) {
				t	= a[m][j];
				_exchange( a[m][j], a[k][j] );

				if ( t != 0.0 )
					for ( i=(k+1); i < n; ++i )
						a[i][j]	+= a[i][k] * t;
			}
		}
/*
 *	Status	= ( 1 - norm of a ) * ( an estimate of 1 - norm of a-inverse )
 *
 * Estimate obtained by one step of inverse iteration for the small
 * singular vector. This involves solving two systems of equations,
 *
 *	( a-transpose ) * y = e and  a * z = y
 *
 * where e is a vector of +1 or -1 chosen to cause growth in y.
 *
 *	estimate = ( 1 - norm of z ) / ( 1 - norm of y )
 *
 * solve
 *	( a-transpose ) * y = e
 */
		for ( k=0; k < n; ++k ) {
            double ek;

			t	= 0.0;
			if ( k != 0 )
				for ( i=0; i < k; ++i )
					t	+= a[i][k] * work[i];

			ek	= ( t < 0.0 )? -1.0 : 1.0;
/*
 * Test for singularity...
 */
			if ( fabs( a[k][k] ) < TOLERANCE )
				return( 0 );

			work[k]	= -( ek + t ) / a[k][k];
		}

		for ( k=(n-2); k >= 0; --k ) {
            unsigned m;

			t	= 0.0;

			for ( i=(k+1); i < n; ++i )
				t	+= a[i][k] * work[k];

			work[k]	= t;
			m	= ipvt[k];

			if ( m != k )
				_exchange( work[k], work[m] );
		}

		for ( i=0; i < n; ++i )
			ynorm	+= fabs( work[i] );
/*
 * Solve: a * z = y
 */
		MatrixLUsolve( n, a, work, ipvt );

		for ( i=0; i < n; ++i )
			znorm	+= fabs( work[i] );
/*
 * Estimate condition
 */
		status	= anorm * znorm / ynorm;

		if ( status < 1.0 )
			status	= 1.0;

		return( status );
	} else {                    /* the trivial case		*/
		status	= 1.0;

		if ( fabs( a[0][0] ) > TOLERANCE )
			return( status );
		else                    /* singularity???		*/
			return( 0.0 );
	}
/* NOTREACHED */
}

/*
 * DETERMINANT - compute determinant of upper triangular matrix.
 * 
 */

	double
MatrixLUdeterminant( unsigned n, MATRIX a, int ipvt[N] ) {
    unsigned i = 0;
    double	 r = ipvt[n-1];

	for ( ; i < n; ++i )
		r	*= a[i][i];

	return( r );
}

/*	VectorWrite - write a vector to a file.				*/

	int
VectorWrite( FILE *fp, char *header, VECTOR a, unsigned n ) {
	unsigned i = 0;

	fprintf( fp, "%s:\n\n", header );

	for ( ; i < n; ++i )
		fprintf( fp, "\t[%3d]:  %9.6f\n", i, *a++ );

	fputc( '\n', fp );

	return( ferror( fp ) );
}

/*	MatrixWrite - write a matrix to a file.				*/

	int
MatrixWrite( FILE *fp, char *header, MATRIX a, unsigned n, unsigned m ) {
	unsigned i = 0;

	if ( ( n > 10 ) || ( m > 10 ) )
		return( ferror( fp ) );

	fprintf( fp, "%s:\n\n", header );

	for ( ; i < n; ++i ) {
		unsigned j = 0;

		fputc( '\t', fp );

		while ( j < m ) {
			fprintf( fp, "%9.6f", a[i][j] );

			if ( ++j != m )
				fputc( ' ', fp );
		}

		fputc( '\n', fp );
	}

	fputc( '\n', fp );

	return( ferror( fp ) );
}

/*	PointWrite - write point vector f(x) to a file.			*/

	int
PointWrite( FILE *fp, char *header, POINT pt[M], unsigned n, double (*f)( double ) ) {
	unsigned i = 0;

	fprintf( fp, "%s:\n\n", header );
    fprintf( fp, "\t x          f(x)       f*(x)      f-f*\n\n" );
//  fprintf( fp, "\t-1.000000   1.000000   0.997352   0.002648" );

	while ( i < n ) {
		double	y	= (*f)( pt->x );
		double	delta	= y - pt->y;

		fprintf( fp, "\t%9.6f", pt->x );
		fprintf( fp, "  %9.6f", y );
		fprintf( fp, "  %9.6f", pt->y );
		fprintf( fp, "  %9.6f", delta );

		fputc( '\n', fp );
		( ++i, ++pt );
	}

	fputc( '\n', fp );

	return( ferror( fp ) );
}

/*	PointWrite - write the point grid for K(x,y) to a file.			*/

int
Point2Write( FILE *fp, char *header, POINT pt[M][M], unsigned n, double (*f)( double, double ) ) {
	unsigned i = 0;

	fprintf( fp, "%s:\n\n", header );
    fprintf( fp, "\t x          y          K(x,y)     K*(x,y)    K-K*\n\n" );
//  fprintf( fp, "\t-1.000000   1.000000   0.997352   0.997352   0.002648" );

	for ( ; i < n; ++i ) {
        unsigned j = 0;

        for ( ; j < n; ++j ) {
            double	z	    = (*f)( pt[i][j].x, pt[i][j].y );
            double	delta	= z - pt[i][j].z;

            fprintf( fp, "\t%9.6f", pt[i][j].x );
            fprintf( fp, "  %9.6f", pt[i][j].y );
            fprintf( fp, "  %9.6f", z );
            fprintf( fp, "  %9.6f", pt[i][j].z );
            fprintf( fp, "  %9.6f", delta );

            fputc( '\n', fp );
        }
	}

	fputc( '\n', fp );
    
	return( ferror( fp ) );
}
