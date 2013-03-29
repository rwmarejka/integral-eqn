/*
% cc -g -Xc -o int int.c -lm
 *
 * INTEGRAL - generate an approximate solution to a Fredholm integral
 *		equation of the second kind using Chenyshev polynominals.
 *
 * written: rwmarejka@mac.com
 */

/* Feature Test Macros		*/

#define	_POSIX_SOURCE

/* Include Files		*/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/* Constants and Macros		*/

#define	TOLERANCE    1.0e-15			/* machine tolerance	*/

#if !defined(M_PI)
#	define	M_PI	(3.1415926535897932384626433)
#endif

#define	N	25		/* size of matrix and vectors	*/
#define	M	11		/* number of sample points	*/

#define	RE(x,l,u)       (0.5*((u)+(l)+((u)-(l))*(x)))
#define	VALUE(i,n)      (cos(((i)*M_PI)/(n)))
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

extern	void	ChebyshevCoeff( VECTOR, int, double, double, double (*)( double ) );
extern	void	Chebyshev2Coeff( MATRIX, int, double, double, double (*)( double, double ) );
extern	double	ChebyshevEval( double, VECTOR, int, double, double );

extern	void	Integrate( MATRIX, MATRIX, int, double, double, double );
extern	void	Points( POINT [M], int, VECTOR, int, double, double );

extern	double	MatrixSolve( int, MATRIX, VECTOR, VECTOR );
extern	double	MatrixLUdecomp( int, MATRIX, int [N] );
extern	double	MatrixLUdeterminant( int, MATRIX, int [N] );
extern	void	MatrixLUsolve( int, MATRIX, VECTOR, int [N] );

extern	int     VectorWrite( FILE *, char *, VECTOR, int );
extern	int     MatrixWrite( FILE *, char *, MATRIX, int, int );
extern	int     PointWrite( FILE *, char *, POINT [M], int, double (*)( double ) );

extern	double	getLower( void );
extern	double	getUpper( void );
extern	double	getLambda( void );
extern	double	g( double );
extern	double	K( double, double );
extern	double	f( double x );

/* External Declarations	*/

/*	main - kiss of life		*/

int
main( int argc, char *argv[] ) {
	VECTOR	x, b;
	MATRIX	a, A;
	POINT	pt[M];
	double	upper;
	double	lower;
	double	lambda;
	int	n;

	lower	= getLower();
	upper	= getUpper();
	lambda	= getLambda();

	if ( argc > 1 )
		n	= atoi( argv[1] );
	else
		n	= ( N > 25 )? 25 : N;

	assert( n <= N );

	fprintf( stdout, "[a,b]:\t[%9.6f,%9.6f]\n", lower, upper );
	fprintf( stdout, "lambda:\t%9.6f\n\n", lambda );

	ChebyshevCoeff( b, n, lower, upper, g );
	b[0]	*= 0.5;
	VectorWrite( stdout, "Chebyshev Coefficients for g(x)", b, n );

	Chebyshev2Coeff( a, n, lower, upper, K );
	MatrixWrite( stdout, "Chebyshev Coefficients for K(x,y)", a, n, n );
	Integrate( a, A, n, lambda, lower, upper );
	MatrixWrite( stdout, "Matrix to be solved", A, n, n );

	if ( MatrixSolve( n, A, x, b ) ) {
		VectorWrite( stdout, "Chebyshev Coefficients for f(x)", x, n );
		x[0]	*= 2.0;
		Points( pt, M, x, n, lower, upper );
		PointWrite( stdout, "Approximate f(x)", pt, M, f );
	}
	return( 0 );
}

/*	ChebyshevCoeff - generate Chebyshev coefficients to f(x)	*/

	void
ChebyshevCoeff( VECTOR a, int n, double lower, double upper, double (*f)( double ) ) {
	VECTOR	b;
	double	um;
	double	fo;
	int	i;

	--n;
	um	= M_PI / n;
	fo	= 2.0  / n;

	for ( i=0; i <= n; ++i ) {
		double	p	= RE( cos( i * um ), lower, upper );

		b[i]	= (*f)( p );
	}

	for ( i=0; i <= n; ++i ) {
		int	j;
		double	sum	= 0.5 * ( b[0] + ( ( i % 2 )? -b[n] : b[n] ) );
		double	ui	= i * um;

		for ( j=1; j < n; ++j )
			sum	+= b[j] * cos( j * ui );

		a[i]	= fo * sum;		/* one coefficient	*/
	}

	a[n]	*= 0.5;			/* half last term		*/

	return;
}

/*	Chebyshev2Coeff - generate Chebyshev coefficients to K(x,y)	*/

	void
Chebyshev2Coeff( MATRIX a, int n, double lb, double ub, double (*f)( double, double ) ) {
	int	i;
	int	n2;
	int	np1;
	MATRIX	Kxy;

	np1	= n--;
	n2	= n * n;

	for ( i=0; i <= n; ++i ) {
		int	j;

		for ( j=0; j <= n; ++j )
			Kxy[i][j]	= (*f)( RE( VALUE( i, n ), lb, ub ), RE( VALUE( j, n ), lb, ub ) );

	}

	for ( i=0; i <= n; ++i ) {
		int	j;
		double	x	= RE( VALUE( i, n ), lb, ub );

		for ( j=0; j <= n; ++j ) {
			int	k;
			VECTOR	b;
			double	y	= RE( VALUE( j, n ), lb, ub );

			for ( k=0; k <= n; ++k ) {
				int	l;
				VECTOR	v;

				for ( l=0; l <= n; ++l )
					v[l]	= Kxy[k][l];

				v[n]	*= 0.5;		/* half last term	*/
				b[k]	 = ChebyshevEval( y, v, np1, lb, ub );
			}

			b[n]	*= 0.5;			/* half last term	*/
			a[i][j]	 = 4.0 * ChebyshevEval( x, b, np1, lb, ub ) / n2;
		}
	}

	for ( i=0; i <= n; ++i ) {		/* half the perimeter terms	*/
		a[0][i]	*= 0.5;
		a[i][0]	*= 0.5;
		a[n][i]	*= 0.5;
		a[i][n]	*= 0.5;
	}
	return;
}

/*	ChebyshevEval - evaluate a Chebyshev series at a point		*/

	double
ChebyshevEval( double x, VECTOR a, int n, double lower, double upper ) {
	double	c0	= 0.0;
	double	c1	= 0.0;
	double	c2	= 0.0;
	double	mul	= 2.0 * ( 2.0 * x - ( upper + lower ) ) / ( upper - lower );

	while ( n-- ) {
		c2	= c1;
		c1	= c0;
		c0	= mul * c1 - c2 + a[n];
	}

	return( 0.5 * ( c0 - c2 ) );
}

/*	Integrate							*/

	void
Integrate( MATRIX b, MATRIX d, int n, double lambda, double lower, double upper ) {
	int     i;
	double	factor	= lambda * ( upper - lower ) / 2.0;
	MATRIX	t;

	for ( i=0; i < n; ++i ) {
		int	j	= ( i % 2 )? 1 : 0;
		int	k	= ( i % 2 )? 0 : 1;

		for ( ; k < n; k+=2 )
			t[i][k]	= 0.0;

		for ( ; j < n; j+=2 ) {
			int	diff	= i - j;
			int	sum     = i + j;

			t[i][j]	= -( 1.0 / ( sum * sum - 1.0 ) + 1.0 / ( diff * diff - 1.0 ) );
		}
	}

	for ( i=0; i < n; ++i ) {
		int	j;

		for ( j=0; j < n; ++j ) {
			int	k;
			double	sum	= b[i][0] * t[0][j];

			for ( k=1; k < n; ++k )
				sum	+= b[i][k] * t[k][j];

			d[i][j]	= factor * sum;
		}

		d[i][i]	+= 1.0;
	}

	return;
}

/*	Points - generate N points from the Chebyshev coefficients to f(x).	*/

	void
Points( POINT pt[M], int m, VECTOR a, int n, double lower, double upper ) {
	int	i	= 0;
	double	dx	= ( upper - lower ) / ( m - 1 );

	for ( ; i < m; ++i, ++pt ) {
		pt->x	= lower + i * dx;
		pt->y	= ChebyshevEval( pt->x, a, n, lower, upper );
	}
	return;
}

/*	MatrixSolve - solve a system of linear equations	*/

	double
MatrixSolve( int n, MATRIX a, VECTOR x, VECTOR b ) {
	int	ipvt[N];	/* the pivot vector		*/
	double	status;		/* indicates singularity	*/
/*
 * Decompose a to upper triangular...
 */
	status	= MatrixLUdecomp( n, a, ipvt );

	if ( status != 0.0 ) {		/* check status of system	*/
		int	i;

		for ( i=0; i < n; ++i )
			x[i]	= b[i];

		MatrixLUsolve( n, a, x, ipvt );

		status	= MatrixLUdeterminant( n, a, ipvt );
	}

	return( status );
}

/*	MatrixLUsolve - solve by back substitution.	*/

	void
MatrixLUsolve( int n, MATRIX a, VECTOR b, int ipvt[N] ) {
	int	i, k, m;
	double	t;			/* avoid indexing		*/

	if ( n != 1 ) {			/* if not the trivial case	*/
/*
 * Forward elimination...
 */
		for ( k=0; k < (n-1); ++k ) {
			m	= ipvt[k];
			_exchange( b[k], b[m] );
			t	= b[k];

			for ( i=(k+1); i < n; ++i )
				b[i]	+= a[i][k] * t;
		}
/*
 * Back substitution...
 */
		for ( k=(n-1); k > 0; --k ) {
			b[k]	/= a[k][k];
			t	 = -b[k];

			for ( i=0; i < k; ++i )
				b[i]	+= a[i][k] * t;
		}
	}
	b[0]	/= a[0][0];			/* the trivial case	*/

	return;
}

/*	MatrixLUdecompDECOMP - reduce matrix to upper triangular.	*/

	double
MatrixLUdecomp( int n, MATRIX a, int ipvt[N] ) {
	int	i, j, k, m;
	double	t, ek, status;
	double	anorm	= 0.0;
	double	ynorm	= 0.0;
	double	znorm	= 0.0;
	VECTOR	work;

	ipvt[n-1] = 1;

	if ( n != 1 ) {
/*
 * Compute 1-norm of a
 */
		for ( j=0; j < n; ++j ) {
			t	= 0;
			for ( i=0; i < n; ++i )
				t	+= fabs( a[i][j] );

			if ( t > anorm )
				anorm	= t;
		}
/*
 * Gaussian elimination with partial pivoting
 */
		for ( k=0; k < (n-1); ++k ) {
			m	= k;
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

			if ( t != 0 )	/* skip step on zero pivot	*/
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

				if ( t != 0 )
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
			t	= 0;
			if ( k != 0 )
				for ( i=0; i < k; ++i )
					t	+= a[i][k] * work[i];

			ek	= ( t < 0 )? -1 : 1;
/*
 * Test for singularity...
 */
			if ( fabs( a[k][k] ) < TOLERANCE )
				return( 0 );

			work[k]	= -( ek + t ) / a[k][k];
		}
		for ( k=(n-2); k >= 0; --k ) {
			t	= 0;

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

		if ( status < 1 )
			status	= 1;

		return( status );
	}
	else {				/* the trivial case		*/
		status	= 1;

		if ( fabs( a[0][0] ) > TOLERANCE )
			return( status );
		else			/* singularity???		*/
			return( 0 );
	}
/* NOTREACHED */
}

/*
 * DETERMINANT - compute determinant of upper triangular matrix.
 * 
 */

	double
MatrixLUdeterminant( int n, MATRIX a, int ipvt[N] ) {
	register int	i	= 0;
	register double	r	= ipvt[n-1];

	for ( ; i < n; ++i )
		r	*= a[i][i];

	return( r );
}

/*	VectorWrite - write a vector to a file.				*/

	int
VectorWrite( FILE *fp, char *header, VECTOR a, int n ) {
	int	i;

	fprintf( fp, "%s:\n\n", header );

	for ( i=0; i < n; ++i )
		fprintf( fp, "\t[%3d]:  %9.6f\n", i, *a++ );

	fputc( '\n', fp );

	return( ferror( fp ) );
}

/*	MatrixWrite - write a matrix to a file.				*/

	int
MatrixWrite( FILE *fp, char *header, MATRIX a, int n, int m ) {
	int	i;

	if ( ( n > 10 ) || ( m > 10 ) )
		return( ferror( fp ) );

	fprintf( fp, "%s:\n\n", header );

	for ( i=0; i < n; ++i ) {
		int	j	= 0;

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

/*	PointWrite - write an array of (x,y) to a file.			*/

	int
PointWrite( FILE *fp, char *header, POINT pt[M], int n, double (*f)( double ) ) {
	int	i	= 0;

	fprintf( fp, "%s:\n\n", header );

	while ( i < n ) {
		double	y	= (*f)( pt->x );
		double	delta	= y - pt->y;

		fprintf( fp, "\t[%3d]:  %9.6f", i, pt->x );
		fprintf( fp, "  %9.6f", y );
		fprintf( fp, "  %9.6f", pt->y );
		fprintf( fp, "  %9.6f", delta );

		fputc( '\n', fp );
		( ++i, ++pt );
	}

	fputc( '\n', fp );

	return( ferror( fp ) );
}
