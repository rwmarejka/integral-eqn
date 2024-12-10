/*
% cc -DTEST=1 -I ../libm -o inteq int.c ../libm/matrix.o ../libm/chebyshev.o ../libm/combination.o -lm
 *
 * INTEGRAL - generate an approximate solution to a Fredholm integral
 *		equation of the first or second kind using Chebyshev polynominals.
 *
 * written: rwmarejka@mac.com
 */

/* Feature Test Macros		*/

#define	_POSIX_SOURCE

/* Include Files            */

#include "matrix.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

/* Constants and Macros		*/

/* External References		*/

extern  vector_t    *ChebyshevCoeff ( unsigned, double, double, double (*)( double         ) );
extern  matrix_t    *Chebyshev2Coeff( unsigned, double, double, double (*)( double, double ) );
extern  double      ChebyshevEval   ( double, vector_t *, unsigned, double, double );

extern  matrix_t    *ChebyshevTMatrix( unsigned );
extern  vector_t    *fredholm( unsigned, double, double, double (*)( double ), double (*)( double, double ), double, unsigned );

/* External Declarations	*/


/*
 * Fredholm - solve a Fredholm integral equation of the first or second kind.
 */

vector_t    *fredholm( unsigned degree, double lower, double upper, double (*g)( double ), double (*K)( double, double ), double lambda, unsigned fredholmType ) {
    vector_t    *b;         // g(x)
    vector_t    *x;         // f(x)
    matrix_t    *a;         // K(x,y)
    matrix_t    *t;         // t matrix
    matrix_t    *A;         // A = a x t + I

    b       = ChebyshevCoeff( degree + 1, lower, upper, g );                // g(x)
    vector_set( b, 0, vector_get( b, 0 ) * 0.5 );

    vector_print( stdout, "g(x)\n", b, NULL );

	a   = Chebyshev2Coeff( degree + 1, lower, upper, K );                   // K(x,y)
    {
        double      **ap = matrix_getData( a );
        unsigned    i;
        unsigned    n   = degree;

        for ( i = 0; i <= n; ++i ) {
            ap[0][i]    *= 0.5;
            ap[i][0]    *= 0.5;
            ap[n][i]    *= 0.5;
            ap[i][n]    *= 0.5;
        }
    }

    matrix_print( stdout, "K(x,y)\n", a, NULL );

    t   = ChebyshevTMatrix( degree + 1 );                                   // t matrix
#if defined(DEBUG)
    matrix_print( stdout, "t matrix\n", t, NULL );
#endif

    A   = matrix_mul( a, t );                                               // A = a x t
    
    if ( ( upper - lower ) != 2.0 )
        matrix_mul_constant( A, ( upper - lower ) * 0.5 );

    if ( lambda != 1.0 )
        matrix_mul_constant( A, lambda );

    if ( fredholmType == 2 ) {                                               // for Fredholm equaitons of the 2nd kind add f(x)
        matrix_add_identity( A );                                           // A = A + I
        matrix_print( stdout, "A = K x t + I\n", A, NULL );
    } else if ( fredholmType == 1 )
        matrix_print( stdout, "A = K x t\n",     A, NULL );
    
    x   = matrix_solve( A, b );                                             // solve  ( a x t + I ) X = b  for X -> f(x)

    vector_free( b );                                                       // release the working space
    matrix_free( a );
    matrix_free( t );
    matrix_free( A );

    if ( x )
        vector_set( x, 0, vector_get( x, 0 ) * 2.0 );
    
	return x;
}

/*
 * ChebyshevTMatrix - compute the t matrix.
 */

matrix_t    *ChebyshevTMatrix( unsigned n ) {
    int         i;
    matrix_t    *t  = matrix_alloc( n, n );

	for ( i = 0; i < n; ++i ) {
		int j = ( i % 2 )? 1 : 0;

		for ( ; j < n; j += 2 ) {
			int	diff = i - j;
			int	sum  = i + j;

			t->m[i][j]	= -( 1.0 / ( sum * sum - 1 ) + 1.0 / ( diff * diff - 1 ) );
		}
	}

    return t;
}


#if defined(TEST)

typedef enum  { REGULAR, ODD, EVEN } Function_Symmetry_t;

long double cheb_eval( double, double *, int, Function_Symmetry_t, double, double );

static  const double _pi       = 3.14159265358979323846L;  /* pi           */
static  const double _pi_div_2 = 1.57079632679489661923L;  /* pi / 2       */
static  const double _e        = 2.71828182845904523536L;  /* e            */
static  const double _4_div_pi = 1.27323954473516268615L;  /* 4 / pi       */
static  const double _2_div_pi = 0.63661977236758134308L;  /* 2 / pi       */


/*  QA Test 0 - everything is 1, except g(x)=3              */

double  ex00_f( double x           ) {   return 1.0;   }
double  ex00_g( double x           ) {   return 3.0;   }
double  ex00_K( double x, double y ) {   return 1.0;   }

/*  ex03 - Fox & Parker                                     */

double  ex03_f( double x           ) {  return 4.0 / 3.0 + 2.0 * x - x * x; }
double  ex03_g( double x           ) {  return -( x * x );                  }
double  ex03_K( double x, double y ) {  return  ( x + y );                  }

/*  ex04 -                                                  */

double  ex04_f( double x           ) {  return ( _e * exp( x ) ) / ( _e * _e +_e - 1.0 ); }
double  ex04_g( double x           ) {  return exp( x );                                  }
double  ex04_K( double x, double y ) {  return exp( x );                                  }

/* ex12 - Numerical Soution of Intergal Equations, pg 98    */

double  ex12_f( double x           ) {  return sin( x );            }
double  ex12_g( double x           ) {  return sin( x ) - 0.25 * x; }
double  ex12_K( double x, double y ) {  return x * y;               }

/* ex13 - Numerical Soution of Intergal Equations, pg 104   */

double  ex13_f( double x           ) {  return exp( x );                                }
double  ex13_g( double x           ) {  return ( exp( x + 1.0 ) - 1.0 ) / ( x + 1.0 );  }
double  ex13_K( double x, double y ) {  return exp( x * y );                            }

/* ex61 - Fox and Goodwin, Example 1 (13)                   */

double  ex61_f( double x           ) {  return sin( x );                }
double  ex61_g( double x           ) {  return -_2_div_pi * cos( x );   }
double  ex61_K( double x, double y ) {  return cos( x - y );            }

/* ex62 - Fox and Goodwin, Example 2 (17)                   */

static const double d   = 1.0;

double  ex62_f1( double x ) {       // d = +1
    static double a[]   = {
        +1.4151850,
        +0.0493851,
        -0.0010475,
        -0.0002327,
        +0.0000200,
        +0.0000010,
        -0.0000002,
    };
    static const unsigned nn    = sizeof( a ) / sizeof( a[0] );

    return cheb_eval( x, a, nn, EVEN, -1.0, +1.0 );
}

double  ex62_f2( double x ) {       // d = -1
    static double a[]   = {
        +3.5488889,
        -0.1400431,
        +0.0049620,
        +0.0003763,
        -0.0000437,
        -0.0000016,
        +0.0000005,
    };
    static const unsigned nn    = sizeof( a ) / sizeof( a[0] );

    return cheb_eval( x, a, nn, EVEN, -1.0, +1.0 );
}

double  ex62_g( double x           ) {  return 1.0;                                     }
double  ex62_K( double x, double y ) {  return d / ( d * d + ( x - y ) * ( x - y ) );   }

/* ex63 - Fox and Goodwin, Example 1 (49)                   */

double  ex63_f( double x           ) {  return 1.0;      }
double  ex63_g( double x           ) {  return 0.5 + x;  }
double  ex63_K( double x, double y ) {  return x + y;    }

/* ex71 - R. E. Scraton, Example 1. (10)                    */

double  ex71_f1( double x           ) { return cosh( x ) / ( 0.5 * 2.0 * sinh( 2.0 ) + ( 2.0 - 1.0 ) );     }
double  ex71_f2( double x           ) { return cosh( x ) / ( 0.5 * 0.5 * sinh( 2.0 ) + ( 0.5 - 1.0 ) );     }
double  ex71_g ( double x           ) { return -cosh( x );                                                  }
double  ex71_K ( double x, double y ) { return cosh( x + y );                                               }

/* ex72 - R. E. Scraton, Example 3. (13)                    */

double  ex72_f( double x           ) {  return cosh( x ) / cosh( 1.0 );                                 }
double  ex72_g( double x           ) {  return 1.0;                                                     }
double  ex72_K( double x, double y ) {  return ( y <= x )? 1.0 + y - x - x * y : 1.0 - y + x - x * y;   }

/* ex81 - Baker et al, Example 2 (39)                       */

double  ex81_f( double x           ) { return x;                                                   }
double  ex81_g( double x           ) { return ( pow( ( 1.0 + x * x ), 1.5 ) - x * x * x ) / 3.0;   }
double  ex81_K( double x, double y ) { return sqrt( x * x + y * y );                               }

/* ex91 - Elliott, Two Examples (54)                        */

static const double k       = +1.2;

double  ex91_f( double x           ) {
    static double a[] = {
        +3.2157560,
        +0.0089949,
        -0.1336355,
        +0.0605130,
        -0.0082906,
        -0.0037893,
        +0.0026980,
        -0.0006514,
        -0.0001003,
        +0.0001418,
        -0.0000491,
        +0.0000004,
        +0.0000076,
        -0.0000036,
        +0.0000005,
        +0.0000001,
    };
    static const unsigned nn    = sizeof( a ) / sizeof( a[0] );

    return cheb_eval( x, a, nn, ODD, -1.0, +1.0 );
}

double  ex91_g( double x           ) {
    double  c   = cos( _pi * x );
    double  s   = sin( _pi * x );

    return 2.0 * atan( k * s / ( k * k * ( c + c * c ) + s * s ) );
}

double  ex91_K( double x, double y ) { return k / ( ( k * k + 1.0 ) - ( k * k - 1.0 ) * cos( _pi * ( x + y ) ) );  }

/* ex101 - Rahbar and Hashemizadeh, Examples (13)"          */

double  ex101_f( double x           ) { return 3.0 * x ;    }
double  ex101_g( double x           ) { return x;           }
double  ex101_K( double x, double y ) { return x * y;       }

/* ex102 - Rahbar and Hashemizadeh, Examples (14)           */

double  ex102_f( double x           ) { return sin( x ) / cos( 1.0 );   }
double  ex102_g( double x           ) { return x;                       }
double  ex102_K( double x, double y ) { return ( x < y )? x : y;        }

/* ex103 - Rahbar and Hashemizadeh, Examples (16)           */

double  ex103_f( double x           ) { return exp( x );            }
double  ex103_g( double x           ) { return _e * x + 1.0;        }
double  ex103_K( double x, double y ) { return ( x < y )? x : y;    }

/* ex104 - Rahbar and Hashemizadeh, Examples (17)           */

double  ex104_f( double x           ) {
    static const double  e2 = _e * _e;
    static const double  c2 = 0.125 * ( e2 * e2 + 6.0 * e2 + 1 ) / ( e2 + 1.0 );
    static const double  c1 = c2 + 1.0 / ( e2 + 1.0 );
    double               ex = exp( x );

    return 0.5 * x * ex + c1 * ex + c2 * exp( -x );
}
double  ex104_g( double x           ) { return exp( x );        }
double  ex104_K( double x, double y ) { return fabs( x - y );   }

/* ex111 - Pilkington, Nystrom's method for solving... (II) */

double  ex111_f( double x           ) { return 1.0;                                                 }
double  ex111_g( double x           ) { return 1.0 + ( atan( 1.0 + x ) + atan( 1.0 - x ) ) / _pi;   }
double  ex111_K( double x, double y ) { return 1.0 / ( 1 + ( x - y ) * ( x - y ) );                 }

/* ex121 */

double  ex121_f( double x           ) { return exp( -x );                                   }
double  ex121_g( double x           ) { return exp( -x ) - 0.5 + 0.5 * exp( -( x + 1 ) );   }
double  ex121_K( double x, double y ) { return ( x + 1.0 ) * exp( -x * y );                 }



typedef struct {
    char        *name;
    int         FredholmType;
    double      lowerBound;
    double      upperBound;
    double      lambda;
    double      (*f)( double );
    double      (*K)( double, double );
    double      (*g)( double );
    unsigned    degree;
    double      dx;
    unsigned    execute;
    unsigned    output;
} IntegralEquation_t;

IntegralEquation_t  examples[] = {
    { "Test 0 - f(x) = 1, K(x,y) = 1, g(x) = 3, [-1,+1]",         2,  -1.0,   1.0,         1.0,        ex00_f,    ex00_K,   ex00_g,    2,  0.500,  0,  1  },
    { "ex03 - Fox & Parker, Chapter 6 (6)",                       2,  -1.0,   1.0,        -1.0,        ex03_f,    ex03_K,   ex03_g,    3,  0.500,  0,  1  },
    { "ex04 -             ",                                      2,  -1.0,   1.0,         1.0,        ex04_f,    ex04_K,   ex04_g,    4,  0.250,  0,  1  },
    { "ex12 - Numerical Solution of Integral Equations, pg 98",   2,   0.0,   _pi_div_2,  -0.25,       ex12_f,    ex12_K,   ex12_g,    5,  0.125,  0,  1  },
    { "ex13 - Numerical Solution of Integral Equations, pg 104",  1,   0.0,   1.0,         1.0,        ex13_f,    ex13_K,   ex13_g,    5,  0.125,  0,  1  },
    { "ex61 - Fox and Goodwin, 4. Example 1 (13)",                2,   0.0,   _pi_div_2,  -_4_div_pi,  ex61_f,    ex61_K,   ex61_g,    5,  0.125,  0,  1  },
    { "ex62 - Fox and Goodwin, 5. Example 2 (17) d = +1",         2,  -1.0,   1.0,         1.0 / _pi,  ex62_f1,   ex62_K,   ex62_g,    6,  0.250,  0,  1  },
    { "ex62 - Fox and Goodwin, 5. Example 2 (17) d = -1",         2,  -1.0,   1.0,        -1.0 / _pi,  ex62_f2,   ex62_K,   ex62_g,    6,  0.250,  0,  1  },
    { "ex63 - Fox and Goodwin, 11. Fredholm's Equation... (49)",  1,   0.0,   1.0,         1.0,        ex63_f,    ex63_K,   ex63_g,    1,  0.250,  0,  1  },
    { "ex71 - R. E. Scraton, Example 1, (10) lambda = 2",         2,  -1.0,   1.0,        -2.0,        ex71_f1,   ex71_K,   ex71_g,    6,  0.250,  0,  1  },
    { "ex71 - R. E. Scraton, Example 1. (10) lambda = 1/2",       2,  -1.0,   1.0,        -0.5,        ex71_f2,   ex71_K,   ex71_g,    6,  0.250,  0,  1  },
    { "ex72 - R. E. Scraton, Example 1. (13)",                    2,  -1.0,   1.0,         0.5,        ex72_f,    ex72_K,   ex72_g,    6,  0.250,  0,  1  },
    { "ex81 - Baker et al, Example 2 (39)",                       1,   0.0,   1.0,         1.0,        ex81_f,    ex81_K,   ex81_g,    6,  0.125,  0,  1  },
    { "ex91 - Elliott, Two Examples (54)",                        2,  -1.0,   1.0,        -1.0,        ex91_f,    ex91_K,   ex91_g,   11,  0.200,  0,  1  },
    { "ex101 - Rahbar and Hashemizadeh, Examples (13)",           2,  -1.0,   1.0,        -1.0,        ex101_f,   ex101_K,  ex101_g,   2,  0.250,  0,  1  },
    { "ex102 - Rahbar and Hashemizadeh, Examples (14)",           2,   0.0,   1.0,        -1.0,        ex102_f,   ex102_K,  ex102_g,   6,  0.125,  0,  1  },
    { "ex103 - Rahbar and Hashemizadeh, Examples (16)",           2,   0.0,   1.0,         1.0,        ex103_f,   ex103_K,  ex103_g,   6,  0.125,  0,  1  },
    { "ex104 - Rahbar and Hashemizadeh, Examples (17)",           2,  -1.0,   1.0,        -0.5,        ex104_f,   ex104_K,  ex104_g,   8,  0.250,  0,  1  },
    { "ex111 - Pilkington, Nystrom's method for solving... (II)", 2,  -1.0,   1.0,         1.0 / _pi,  ex111_f,   ex111_K,  ex111_g,   6,  0.250,  0,  1  },
    { "ex121 - Richardson, Introductory Example",                 2,   0.0,   1.0,        -0.5,        ex121_f,   ex121_K,  ex121_g,   4,  0.125,  1,  1  },
};
#define NTESTS  (sizeof(examples)/sizeof(examples[0]))


/*
 * main - run the examples
 */

int main( int argc, char *argv[] ) {
    unsigned    i   = 0;

    for ( ; i < NTESTS; ++i ) {
        unsigned    n       = examples[i].degree;
        unsigned    N       = n + 1;
        double      lower   = examples[i].lowerBound;
        double      upper   = examples[i].upperBound;
        vector_t    *b      = NULL;
        vector_t    *bf     = NULL;

        if ( ! examples[i].execute )
            continue;

        printf( "%s\n\n", examples[i].name );
            
        b   = fredholm( n, lower, upper, examples[i].g, examples[i].K, examples[i].lambda, examples[i].FredholmType );
        bf  = ChebyshevCoeff( N, lower, upper, examples[i].f );
            
        if ( b )
            vector_print( stdout, "x\n", b,  NULL );
        else
            printf( "x - no solution\n\n" );

        vector_print( stdout, "f(x)\n", bf, NULL );

        if ( b && examples[i].output ) {
            double  x   = lower;
            double  dx  = examples[i].dx;

            printf( "\tx\tf(x)\t\tTn(x)\t\tX(x)\t\tf(x)-X(x)\tTn(x)-X(x)\n\n" );
        
            for ( ; x <= upper; x += dx ) {
                double  y   = examples[i].f( x );
                double  yf  = ChebyshevEval( x, bf, N, lower, upper );
                double  yx  = ChebyshevEval( x, b,  N, lower, upper );

                if ( examples[i].K == ex62_K )
                    printf( "\t%+5.3f\t%+.4f\t\t%+.4f\t\t%+.4f\t\t%+.4f\t\t%+.4f\n", x, y, yf, yx, ( y - yx ), ( yf - yx ) );
                else if ( examples[i].K == ex91_K )
                    printf( "\t%+5.3f\t%+.6f\t%+.6f\t%+.6f\t%+.6f\t%+.6f\n", x, y, yf, yx, ( y - yx ), ( yf - yx ) );
                else
                    printf( "\t%+5.3f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\n", x, y, yf, yx, ( y - yx ), ( yf - yx ) );
            }

            putchar( '\n' );
        }

        putchar( '\n' );

        if ( b )
            vector_free( b );
        
        vector_free( bf );
    }
    
    return 0;
}
#endif
