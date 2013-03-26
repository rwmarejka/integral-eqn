/*
%  cc -O -Xc -o int-st integral.c -lm
%t cc -O -Xc -o int-mt -D_REENTRANT integral.c -lm -lthread
 *
 * INTEGRAL - generate an approximate solution to a Fredholm integral
 *		equation of the second kind using Chenyshev polynominals.
 *
 * written: Richard.Marejka@Canada.Sun.COM
 */

static char	*SCCSid	= "@(#) integral.c 2.4 02/06/10 Richard Marejka";

/* Feature Test Macros		*/

#define	_POSIX_SOURCE
#define	DEBUG	0

/* Include Files		*/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#if defined(_REENTRANT)
#	include <thread.h>
#	include <synch.h>
#	include <signal.h>
#	include <string.h>
#endif

/* Constants and Macros		*/

#if defined(_REENTRANT)
#	define	PARALLEL	8
#	define	IO_lock()	(mutex_lock(&io_lock))
#	define	IO_unlock()	(mutex_unlock(&io_lock))
#endif

#define	TOLERANCE    1.0e-15			/* machine tolerance	*/

#if !defined(M_PI)
#	define	M_PI	(3.1415926535897932384626433)
#endif

#define	N	100		/* size of matrix and vectors	*/
#define	M	11		/* number of sample points	*/

#define	RE(x,l,u)	(0.5*((u)+(l)+((u)-(l))*(x)))
#define	VALUE(i,n)	(cos(((i)*M_PI)/(n)))
#define	_exchange(a,b)	{double _T = (a); (a)=(b); (b)=_T;}

/* Data Declarations		*/

typedef double	MATRIX[N][N];
typedef double	VECTOR[N];

typedef struct {
	double	x;
	double	y;
	double	z;
} POINT;

typedef struct {
	int	nrow;		/* number of rows		*/
	int	ncol;		/* number of columns		*/
	MATRIX	m;		/* the matrix			*/
} matrix_t;

#if defined(_REENTRANT)
typedef struct {
	matrix_t	*a;	/* the matrix			*/
	int		row;	/* starting row number		*/
	int		nrow;	/* number of rows		*/
	void		*opd;	/* operation private data	*/
} MatrixOpArg;
#endif

/* External References		*/

extern	double	g( double );
extern	double	K( double, double );

extern	void	ChebyshevCoeff( VECTOR, int, double, double, double (*)( double ) );
extern	void	Chebyshev2Coeff( MATRIX, int, double, double, double (*)( double, double ) );
extern	double	ChebyshevEval( double, VECTOR, int, double, double );
extern	double	Chebyshev2Eval( double, double, MATRIX, int, double, double );

extern	void	Integrate( MATRIX, MATRIX, int, double, double, double );
extern	int	SolveEquations( MATRIX, VECTOR, VECTOR, int );
extern	void	Points( POINT [M], int, VECTOR, int, double, double );
extern	void	Points2( POINT [M][M], int, MATRIX, int, double, double );

extern	double	MatrixSolve( int, MATRIX, VECTOR, VECTOR );
extern	double	MatrixLUdecomp( int, MATRIX, int [N] );
extern	double	MatrixLUdeterminant( int, MATRIX, int [N] );
extern	void	MatrixLUsolve( int, MATRIX, VECTOR, int [N] );

extern	int	VectorWrite( FILE *, char *, VECTOR, int );
extern	int	MatrixWrite( FILE *, char *, MATRIX, int, int );
extern	int	PointWrite( FILE *, char *, POINT [M], int );

#if defined(_REENTRANT)
	extern	int	MatrixOp( matrix_t *, void *(*)( void * ), void * );
	extern	void	*t_mop( MatrixOpArg * );
	extern	void	*t_matrix( void * );
	extern	void	*b_vector( void * );
	extern	void	*K_mop( MatrixOpArg * );
	extern	void	*K_matrix( void * );
	extern	void	*A_mop( MatrixOpArg * );
	extern	void	*A_matrix( void * );
	extern	void	*b_mop( MatrixOpArg * );
	extern	void	*b_matrix( void * );
#endif

/* External Declarations	*/

#if defined(_REENTRANT)
	mutex_t		io_lock;
	matrix_t	t;
	matrix_t	k;
	matrix_t	Kc;
	matrix_t	A;
	matrix_t	a;

	int	nCPU;
#else
	MATRIX	a, A;
#endif

double	lower	= 0.0;
double	upper	= M_PI;
double	lambda	= 1.0;
int	n	= 30;

VECTOR	x, b;
POINT	pt[M];

/*	main - kiss of life		*/

main( int argc, char *argv[] ) {
#if defined(_REENTRANT)
	int		s;
	thread_t	b_tid;
	thread_t	t_tid;
	thread_t	K_tid;
	thread_t	A_tid;
#endif
	if ( argc > 1 )
		n	= atoi( argv[1] );

	assert( n <= N );

#if defined(_REENTRANT)
	if ( argc > 2 )
		nCPU	= atoi( argv[2] );

	t.nrow	= n;
	t.ncol	= n;
	k.nrow	= n;
	k.ncol	= n;
	Kc.nrow	= n;
	Kc.ncol	= n;
	A.nrow	= n;
	A.ncol	= n;
	a.nrow	= n;
	a.ncol	= n;
#endif
#if defined(_REENTRANT)
	if ( s = thr_create( NULL, 0, b_vector, (void *) g, THR_BOUND, &b_tid ) ) {
		fprintf( stderr, "integral: thr_create: %s\n", strerror( s ) );
		exit( 1 );
	}

	if ( s = thr_create( NULL, 0, t_matrix, NULL, THR_BOUND, &t_tid ) ) {
		fprintf( stderr, "integral: thr_create: %s\n", strerror( s ) );
		exit( 1 );
	}

	if ( s = thr_create( NULL, 0, K_matrix, (void *) K, THR_BOUND, &K_tid ) ) {
		fprintf( stderr, "integral: thr_create: %s\n", strerror( s ) );
		exit( 1 );
	}

	thr_join( b_tid, NULL, NULL );
	thr_join( t_tid, NULL, NULL );
	thr_join( K_tid, NULL, NULL );

	if ( s = thr_create( NULL, 0, b_matrix, NULL, THR_BOUND, &b_tid ) ) {
		fprintf( stderr, "integral: thr_create: %s\n", strerror( s ) );
		exit( 1 );
	}

	thr_join( b_tid, NULL, NULL );

	if ( s = thr_create( NULL, 0, A_matrix, NULL, THR_BOUND, &A_tid ) ) {
		fprintf( stderr, "integral: thr_create: %s\n", strerror( s ) );
		exit( 1 );
	}

	thr_join( A_tid, NULL, NULL );

	if ( MatrixSolve( n, A.m, x, b ) ) {
		VectorWrite( stdout, "Chebyshev Coefficients for f(x)", x, n );
		x[0]	*= 2.0;
		Points( pt, M, x, n, lower, upper );
		PointWrite( stdout, "Approximate f(x)", pt, M );
	}

#else
	ChebyshevCoeff( b, n, lower, upper, g );
	b[0]	*= 0.5;
	VectorWrite( stdout, "Chebyshev Coefficients for g(x)", b, n );

	Chebyshev2Coeff( a, n, lower, upper, K );
#if DEBUG
	MatrixWrite( stdout, "Chebyshev Coefficients for K(x,y)", a, n, n );
#endif
	Integrate( a, A, n, lambda, lower, upper );
#if DEBUG
	MatrixWrite( stdout, "Matrix to be solved", A, n, n );
#endif
	if ( MatrixSolve( n, A, x, b ) ) {
		VectorWrite( stdout, "Chebyshev Coefficients for f(x)", x, n );
		x[0]	*= 2.0;
		Points( pt, M, x, n, lower, upper );
		PointWrite( stdout, "Approximate f(x)", pt, M );
	}
#endif
	return( 0 );
}

/*	g(x)								*/

	double
g( double x ) {
	return( ( 1 + M_PI / 2 ) * cos( x ) );
}

/*	K(x,y)								*/

	double
K( double x, double y ) {
	return( cos( x ) * cos( y ) );
}

/*	ChebyshevCoeff - generate Chebyshev coefficients to f(x)	*/

	void
ChebyshevCoeff( VECTOR a, int n, double lb, double ub, double (*f)( double ) ) {
	VECTOR	b;
	double	um;
	double	fo;
	int	i;

	--n;
	um	= M_PI / n;
	fo	= 2.0  / n;

	for ( i=0; i <= n; ++i ) {
		double	p	= RE( cos( i * um ), lb, ub );

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

#if !defined(_REENTRANT)

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
#endif

/*	ChebyshevEval - evaluate a Chebyshev series at a point		*/

	double
ChebyshevEval( double x, VECTOR a, int n, double lb, double ub ) {
	double	c0	= 0.0;
	double	c1	= 0.0;
	double	c2	= 0.0;
	double	mul	= 2.0 * ( 2.0 * x - ( ub + lb ) ) / ( ub - lb );

	while ( n-- ) {
		c2	= c1;
		c1	= c0;
		c0	= mul * c1 - c2 + a[n];
	}

	return( 0.5 * ( c0 - c2 ) );
}

/*	Integrate							*/

#if !defined(_REENTRANT)

	void
Integrate( MATRIX b, MATRIX d, int n, double lam, double lb, double ub ) {
	int	i;
	double	factor	= ( ub - lb ) / ( 2 * lam );
	MATRIX	t;

	for ( i=0; i < n; ++i ) {
		int	j	= ( i % 2 )? 1 : 0;
		int	k	= ( i % 2 )? 0 : 1;

		for ( ; k < n; k+=2 )
			t[i][k]	= 0.0;

		for ( ; j < n; j+=2 )
			t[i][j]	= 0.5 * ( 1.0/(i+j+1.0) - 1.0/(i+j-1.0) + 1.0/(i-j+1.0) - 1.0/(i-j-1.0) );
	}
#if DEBUG && 0
	MatrixWrite( stdout, "Integration Matrix", t, n, n );
#endif
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
#endif

/*	Points - generate N points from the Chebyshev coefficients to f(x).	*/

	void
Points( POINT pt[M], int m, VECTOR a, int n, double lb, double ub ) {
	int	i	= 0;
	double	dx	= ( ub - lb ) / ( m - 1 );

	for ( ; i < m; ++i, ++pt ) {
		pt->x	= lower + i * dx;
		pt->y	= ChebyshevEval( pt->x, a, n, lb, ub );
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
		fprintf( fp, "\t[%3d]:\t%18.12f\n", i, *a++ );

	fputc( '\n', fp );

	return( ferror( fp ) );
}

/*	MatrixWrite - write a matrix to a file.				*/

	int
MatrixWrite( FILE *fp, char *header, MATRIX a, int n, int m ) {
	int	i;

	fprintf( fp, "%s:\n\n", header );

	if ( ( n > 10 ) || ( m > 10 ) )
		return( ferror( fp ) );

	for ( i=0; i < n; ++i ) {
		int	j	= 0;

		fputc( '\t', fp );

		while ( j < m ) {
			fprintf( fp, "%10.6f", a[i][j] );

			if ( ++j != m )
				fputc( '\t', fp );
		}

		fputc( '\n', fp );
	}

	fputc( '\n', fp );

	return( ferror( fp ) );
}

/*	PointWrite - write an array of (x,y) to a file.			*/

	int
PointWrite( FILE *fp, char *header, POINT pt[M], int n ) {
	int	i	= 0;

	fprintf( fp, "%s:\n\n", header );

	while ( i < n ) {
		fprintf( fp, "\t[%3d]:\t%10.6f\t%18.12f\n", i, pt->x, pt->y );
		( ++i, ++pt );
	}

	fputc( '\n', fp );

	return( ferror( fp ) );
}

/*
 * T H R E A D S   O P E R A T I O N S
 */

#if defined(_REENTRANT)

/*	MatrixOp - execute a matrix operation in parallel	*/

	int
MatrixOp( matrix_t *a, void *(*op)( void * ), void *opd ) {
	int		i;
	int		j;
	int		rpt;
	int		mod;
	thread_t	tid[PARALLEL];
	MatrixOpArg	arg[PARALLEL];

	if ( nCPU == 0 ) {
		nCPU	= sysconf( _SC_NPROCESSORS_ONLN );

		assert( nCPU != -1 );
		assert( nCPU != 0 );

	}

	if ( nCPU > PARALLEL )
		nCPU	= PARALLEL;

	rpt	= a->nrow / nCPU;
	mod	= a->nrow % nCPU;

				/* fill the argument structures	*/
	for ( i=0, j=0; i < nCPU; ++i ) {
		arg[i].a	= a;
		arg[i].opd	= opd;
		arg[i].nrow	= rpt;
		arg[i].row	= j;

		if ( mod )
			( arg[i].nrow++, mod-- );

		j	+= arg[i].nrow;
	}
				/* create the threads		*/
	for ( i=0; i < nCPU; ++i ) {
		int	n;

		if ( n = thr_create( NULL, 0, op, &arg[i], THR_BOUND, &tid[i] ) ) {
			for ( ; i; i-- )
				thr_kill( tid[i], SIGKILL );

			return( n );
		}
	}

				/* wait for them to finish	*/
	for ( i=0; i < nCPU; ++i )
		thr_join( tid[i], NULL, (void **) NULL );

	return( 0 );
}

/*	t_mop - sub-operation for filling the t matrix			*/

	void *
t_mop( MatrixOpArg *arg ) {
	int	i;
	int	lim;
	MATRIX	*m	= &arg->a->m;

	lim	= ( i = arg->row ) + arg->nrow;

	for ( ; i < lim; ++i ) {
		int	j	= ( i % 2 )? 1 : 0;
		int	k	= ( i % 2 )? 0 : 1;

		for ( ; k < n; k+=2 )
			(*m)[i][k]	= 0.0;

		for ( ; j < n; j+=2 )
			(*m)[i][j]	= 0.5 * ( 1.0/(i+j+1.0) - 1.0/(i+j-1.0) + 1.0/(i-j+1.0) - 1.0/(i-j-1.0) );
	}

	return( NULL );
}

/*	t_matrix - fill in the t matrix					*/

	void *
t_matrix( void *arg ) {
	MatrixOp( &t, (void *(*)(void *)) t_mop, NULL );
#if DEBUG
	IO_lock();
	MatrixWrite( stdout, "Integration Matrix", t.m, t.nrow, t.ncol );
	IO_unlock();
#endif
	return( NULL );
}

/*	b_vector - generate Chebyshev coefficients for g(x)		*/

	void *
b_vector( void *arg ) {
	ChebyshevCoeff( b, n, lower, upper, (double (*)(double)) arg );
	b[0]	*= 0.5;
	IO_lock();
	VectorWrite( stdout, "Chebyshev Coefficients for g(x)", b, n );
	IO_unlock();

	return( NULL );
}

/*	K_mop - sub-operation for filling the K matrix			*/

	void *
K_mop( MatrixOpArg *arg ) {
	int	k;
	int	lim;
	int	n			= arg->a->nrow - 1;
	double	(*f)( double, double )	= (double (*)(double, double)) arg->opd;
	MATRIX	*m			= &arg->a->m;

	lim	= ( k = arg->row ) + arg->nrow;

	for ( ; k < lim; ++k ) {
		int	l;

		for ( l=0; l <= n; ++l )
			(*m)[k][l]	= (*f)( RE( VALUE( k, n ), lower, upper ), RE( VALUE( l, n ), lower, upper ) );

		(*m)[k][n]	*= 0.5;		/* half last term	*/
	}

	return( NULL );
}

/*	K_matrix - fill in the K matrix					*/

	void *
K_matrix( void *arg ) {
	MatrixOp( &k, (void *(*)(void *)) K_mop, arg );
#if DEBUG
	IO_lock();
	MatrixWrite( stdout, "K(x,y) Evaluation Matrix", k.m, k.nrow, k.ncol );
	IO_unlock();
#endif
	return( NULL );
}

/*	A_mop - sub-operation for computing A matrix			*/

	void *
A_mop( MatrixOpArg *arg ) {
	int	i;
	int	lim;
	int	n	= arg->a->nrow;
	double	factor	= ( upper - lower ) / ( 2 * lambda );
	MATRIX	*m	= &arg->a->m;

	lim	= ( i = arg->row ) + arg->nrow;

	for ( ; i < lim; ++i ) {
		int	j;

		for ( j=0; j < n; ++j ) {
			int	k;
			double	sum	= 0.0;

			for ( k=0; k < n; ++k )
				sum	+= Kc.m[i][k] * t.m[k][j];

			(*m)[i][j]	= factor * sum;
		}

		(*m)[i][i]	+= 1.0;
	}

	return( NULL );
}

/*	A_matrix - compute the A matrix to be solved			*/

	void *
A_matrix( void *arg ) {
	MatrixOp( &A, (void *(*)(void *)) A_mop, NULL );
#if DEBUG
	IO_lock();
	MatrixWrite( stdout, "Matrix to be solved", A.m, A.nrow, A.ncol );
	IO_unlock();
#endif
	return( NULL );
}

/*	b_mop - sub-operation						*/

	void *
b_mop( MatrixOpArg *arg ) {
	int	i;
	int	lim;
	int	n	= arg->a->nrow;
	MATRIX	*m	= &arg->a->m;
	MATRIX	*K	= &k.m;
	int	n2;
	int	np1;

	lim	= ( i = arg->row ) + arg->nrow;
	np1	= n--;
	n2	= n * n;

	for ( ; i < lim; ++i ) {
		int	j;
		double	x	= RE( VALUE( i, n ), lower, upper );

		for ( j=0; j <= n; ++j ) {
			int	k;
			VECTOR	b;
			double	y	= RE( VALUE( j, n ), lower, upper );

			for ( k=0; k <= n; ++k )
				b[k]	 = ChebyshevEval( y, (*K)[k], np1, lower, upper );

			b[n]	*= 0.5;			/* half last term	*/
			(*m)[i][j]	 = 4.0 * ChebyshevEval( x, b, np1, lower, upper ) / n2;
		}
	}

	return( NULL );
}

/*	b_matrix - compute the b matrix, Chebyshev coefficients for K(x,y)	*/

	void *
b_matrix( void *arg ) {
	int	i;
	int	n	= Kc.nrow - 1;

	MatrixOp( &Kc, (void *(*)(void *)) b_mop, NULL );

	for ( i=0; i <= n; ++i ) {		/* half the perimeter terms	*/
		Kc.m[0][i]	*= 0.5;
		Kc.m[i][0]	*= 0.5;
		Kc.m[n][i]	*= 0.5;
		Kc.m[i][n]	*= 0.5;
	}
#if DEBUG
	IO_lock();
	MatrixWrite( stdout, "Chebyshev Coefficients for K(x,y)", Kc.m, Kc.nrow, Kc.ncol );
	IO_unlock();
#endif
	return( NULL );
}
#endif	/* _REENTRANT	*/
