/*
% cc -DTEST=1 -I ../libm -o inteq int.c ../libm/matrix.o ../libm/chebyshev.o ../libm/combination.o -lm
 *
 * INTEGRAL - generate an approximate solution to a Fredholm integral
 *		equation of the first or second kind using Chebyshev polynominals.
 *
 * written: rwmarejka@mac.com
 */

// Feature Test Macrox

#define	_POSIX_SOURCE

// Include Files

#include "matrix.h"

#include <assert.h>
#include <stdio.h>

// Constants and Macros

// External References

extern  vector_t    *fredholm( unsigned, double, double, double (*)( double ), double (*)( double, double ), double, unsigned );


// External Declarations

static const char *number_fmt   = "%+.5f";


/*
 * ChebyshevTMatrix - compute the t matrix.
 *
 *  t[i][j] = integral Ti(x) Tj(x) dx
 */

static matrix_t *ChebyshevTMatrix( unsigned n ) {
    int         i;
    matrix_t    *t  = matrix_alloc( n, n );

    for ( i = 0; i < n; ++i ) {
        int j = ( i % 2 )? 1 : 0;

        for ( ; j < n; j += 2 ) {
            int    diff = i - j;
            int    sum  = i + j;

            t->m[i][j]    = -( 1.0 / ( sum * sum - 1 ) + 1.0 / ( diff * diff - 1 ) );
        }
    }

    return t;
}

/*
 * Fredholm - solve a Fredholm integral equation of the first or second kind.
 */

vector_t    *fredholm( unsigned degree, double lower, double upper, double (*g)( double ), double (*K)( double, double ), double lambda, unsigned fredholmType ) {
    vector_t    *b;         // g(x)
    vector_t    *x;         // f(x)
    matrix_t    *a;         // K(x,y)
    matrix_t    *t;         // t matrix
    matrix_t    *A;         // A = a x t + I
    const char  *AMatrixTitle;

    assert( fredholmType == 1 || fredholmType == 2 );

    b       = ChebyshevCoeff( degree + 1, lower, upper, g, REGULAR, ZERO ); // g(x) vector
    vector_set( b, 0, vector_get( b, 0 ) * 0.5 );                           // half the first term of g(x)

    vector_print( stdout, "g(x)\n", b, number_fmt );

	a   = Chebyshev2Coeff( degree + 1, lower, upper, K );                   // K(x,y) matrix

    {                                                                       // apply the summation primes to K(x,y)
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

    matrix_print( stdout, "K(x,y)\n", a, number_fmt );

    t   = ChebyshevTMatrix( degree + 1 );                                   // compute t matrix
#if defined(DEBUG)
    matrix_print( stdout, "t matrix\n", t, number_fmt );
#endif

    A   = matrix_mul( a, t );                                               // A = K x t
    
    if ( ( upper - lower ) != 2.0 )                                         // apply the integration differential coefficient, (b-a)/2
        matrix_mul_constant( A, ( upper - lower ) * 0.5 );

    if ( lambda != 1.0 )                                                    // apply the integral coefficient, lanbda
        matrix_mul_constant( A, lambda );

    if ( fredholmType == 2 ) {                                              // for Fredholm equaitons of the 2nd kind add f(x), A += I
        matrix_add_identity( A );
        AMatrixTitle    = "A = K x t + I\n";
    } else if ( fredholmType == 1 )
        AMatrixTitle    = "A = K x t\n";

    matrix_print( stdout, AMatrixTitle, A, number_fmt );

    x   = matrix_solve( A, b );                                             // f(x) vector = solve ( K x t + I ) X = b

    vector_free( b );                                                       // release the working space
    matrix_free( a );
    matrix_free( t );
    matrix_free( A );

    if ( x )
        vector_set( x, 0, vector_get( x, 0 ) * 2.0 );                       // double the first term of f(x)
    
	return x;
}

#if defined(TEST)

#include <math.h>

// Type Definitions

/*
 * Test cases are described by an IntegralEqation_t.
 */

typedef struct {
    char        *name;                          // name of test case
    int         FredholmType;                   // Fredholm type { 1, 2 }
    double      lowerBound;                     // lower bound of range
    double      upperBound;                     // upper bound of range
    double      lambda;                         // coefficient of integral
    double      (*f)( double );                 // f(x) - solution
    double      (*K)( double, double );         // K(x,y) - kernel function
    double      (*g)( double );                 // g(x) -
    unsigned    degree;                         // degree of solution
    double      dx;                             // increment for range report
    unsigned    execute;                        // execute test case { 0, 1 }
    unsigned    output;                         // output range report { 0, 1 }
} IntegralEquation_t;

// Externla References

// External Declarations

static  const double _pi       = 3.14159265358979323846L;  /* pi           */
static  const double _pi_div_2 = 1.57079632679489661923L;  /* pi / 2       */
static  const double _4_div_pi = 1.27323954473516268615L;  /* 4 / pi       */
static  const double _2_div_pi = 0.63661977236758134308L;  /* 2 / pi       */
static  const double _e        = 2.71828182845904523536L;  /* e            */

/*
 * Test case functions K(x,y), g(x) and f(x) are defined here.
 */

// Test 0 - f(x) = 1, K(x,y) = 1, g(x) = 3, [-1,+1]

double  ex00_f( double x           ) {   return 1.0;   }
double  ex00_g( double x           ) {   return 3.0;   }
double  ex00_K( double x, double y ) {   return 1.0;   }

// Test 1 - f(x) = cos(x), K(x,y) = cos(x) x cos(y), g(x) = B x cos(x), [0,pi]

double  ex01_f( double x           ) { return cos( x );                          }
double  ex01_g( double x           ) { return ( 1.0 + _pi_div_2 ) * cos( x );    }
double  ex01_K( double x, double y ) { return cos( x ) * cos( y );               }

//  Test 2 - f(x) = B x exp(x), K(x,y) = g(x) = exp(x), [-1,+1]

double  ex04_f( double x           ) {  return ( _e * exp( x ) ) / ( _e * _e +_e - 1.0 ); }
double  ex04_g( double x           ) {  return exp( x );                                  }
double  ex04_K( double x, double y ) {  return exp( x );                                  }



// ex03 - Fox & Parker, Chapter 6 (6)

double  ex03_f( double x           ) {  return 4.0 / 3.0 + 2.0 * x - x * x; }
double  ex03_g( double x           ) {  return -( x * x );                  }
double  ex03_K( double x, double y ) {  return  ( x + y );                  }



// ex61 - Fox and Goodwin, Example 1 (13)
// ex105 - Rahbar and Hashemizadeh, Examples (11)

double  ex61_f( double x           ) {  return sin( x );                }
double  ex61_g( double x           ) {  return -_2_div_pi * cos( x );   }
double  ex61_K( double x, double y ) {  return cos( x - y );            }

// ex62 - Fox and Goodwin, Example 2 (17)

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

// ex63 - Fox and Goodwin, Example 1 (49)

double  ex63_f1( double x           ) {  return 1.0;                    }
double  ex63_f2( double x           ) {  return x;                      }
double  ex63_g1( double x           ) {  return 0.5 + x;                }
double  ex63_g2( double x           ) {  return 1.0 / 3.0 + 0.5 * x;    }
double  ex63_K ( double x, double y ) {  return x + y;                  }

// ex64 - Fox and Goodwin, Example 1 (51)
// ex81 - Baker et al, Example 2 (39)

double  ex64_f( double x           ) { return x;                        }
double  ex64_g( double x           ) {
    double  x2  = x * x;
    double  t   = 1.0 + x2;
    
    return ( sqrt( t * t * t ) - x2 * x ) / 3.0;
}

double  ex64_K( double x, double y ) { return sqrt( x * x + y * y );    }



// ex91 - Elliott, Two Examples (54)

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



// ex82 - Baker et al, Example 2 (17)

double  ex82_f( double x           ) {  return 4.0 - 6.0 * x;   }
double  ex82_g( double x           ) {  return x;               }
double  ex82_K( double x, double y ) {  return x + y;           }

// ex83 - Baker et al, Example 2 (42)
// ex13 - Numerical Soution of Intergal Equations, pg 104

double  ex83_f( double x           ) {  return exp( x );                                }
double  ex83_g( double x           ) {  return ( exp( x + 1.0 ) - 1.0 ) / ( x + 1.0 );  }
double  ex83_K( double x, double y ) {  return exp( x * y );                            }



// ex71 - R. E. Scraton, Example 1. (10)

double  ex71_f1( double x           ) { return cosh( x ) / ( 0.5 * 2.0 * sinh( 2.0 ) + ( 2.0 - 1.0 ) );     }
double  ex71_f2( double x           ) { return cosh( x ) / ( 0.5 * 0.5 * sinh( 2.0 ) + ( 0.5 - 1.0 ) );     }
double  ex71_g ( double x           ) { return -cosh( x );                                                  }
double  ex71_K ( double x, double y ) { return cosh( x + y );                                               }

// ex72 - R. E. Scraton, Example 3. (13)

double  ex72_f( double x           ) {  return cosh( x ) / cosh( 1.0 );                                 }
double  ex72_g( double x           ) {  return 1.0;                                                     }
double  ex72_K( double x, double y ) {  return ( y <= x )? 1.0 + y - x - x * y : 1.0 - y + x - x * y;   }



// ex12 - Numerical Soution of Intergal Equations, pg 98

double  ex12_f( double x           ) {  return sin( x );            }
double  ex12_g( double x           ) {  return sin( x ) - 0.25 * x; }
double  ex12_K( double x, double y ) {  return x * y;               }

// ex14 - Numerical Solution of Integral Equations, pg 102

double  ex14_f( double x           ) {  return ( 7.0 * _e ) / ( _e * _e - 1.0 ) * ( exp( x ) - exp( -x ) ) - 6.0 * x;   }
double  ex14_g( double x           ) {  return x * x * x;                                                               }
double  ex14_K( double x, double y ) {  return ( x <= y ) ? x * ( 1.0 - y ) : y * ( 1.0 - x );                          }



// ex111 - Pilkington, Nystrom's method for solving... (II)

double  ex111_f( double x           ) { return 1.0;                                                 }
double  ex111_g( double x           ) { return 1.0 + ( atan( 1.0 + x ) + atan( 1.0 - x ) ) / _pi;   }
double  ex111_K( double x, double y ) { return 1.0 / ( 1 + ( x - y ) * ( x - y ) );                 }



// ex121 - Richardson, Introductory Example

double  ex121_f( double x           ) { return exp( -x );                                   }
double  ex121_g( double x           ) { return exp( -x ) - 0.5 + 0.5 * exp( -( x + 1 ) );   }
double  ex121_K( double x, double y ) { return ( x + 1.0 ) * exp( -x * y );                 }



// ex106 - Rahbar and Hashemizadeh, Examples (12)

double  ex106_f( double x           ) { return ( 2.0 / 3.0 ) * sqrt( x );   }
double  ex106_g( double x           ) { return sqrt( x );                   }
double  ex106_K( double x, double y ) { return sqrt( x * y );               }

// ex101 - Rahbar and Hashemizadeh, Examples (13)

double  ex101_f( double x           ) { return 3.0 * x ;    }
double  ex101_g( double x           ) { return x;           }
double  ex101_K( double x, double y ) { return x * y;       }

// ex102 - Rahbar and Hashemizadeh, Examples (14), (15)

double  ex102_f( double x           ) { return sin( x ) / cos( 1.0 );   }
double  ex102_g( double x           ) { return x;                       }
double  ex102_K( double x, double y ) { return ( x < y )? x : y;        }

// ex103 - Rahbar and Hashemizadeh, Examples (16)

double  ex103_f( double x           ) { return exp( x );            }
double  ex103_g( double x           ) { return _e * x + 1.0;        }
double  ex103_K( double x, double y ) { return ( x < y )? x : y;    }

// ex104 - Rahbar and Hashemizadeh, Examples (17)

double  ex104_f( double x           ) {
    static const double  e2 = _e * _e;
    static const double  c2 = 0.125 * ( e2 * e2 + 6.0 * e2 + 1 ) / ( e2 + 1.0 );
    static const double  c1 = c2 + 1.0 / ( e2 + 1.0 );
    double               ex = exp( x );

    return 0.5 * x * ex + c1 * ex + c2 * exp( -x );
}
double  ex104_g( double x           ) { return exp( x );        }
double  ex104_K( double x, double y ) { return fabs( x - y );   }


// examples - define the test cases as array of IntegralEquation_t

IntegralEquation_t  examples[] = {
    { "Test 0 - f(x)=K(x,y)=1, g(x)=3, [-1,+1]",        2,  -1.0,   1.0,         1.0,        ex00_f,    ex00_K,   ex00_g,    2,  0.500,        0,  1  },
    { "Test 1",                                         2,   0.0,   _pi,         1.0,        ex01_f,    ex01_K,   ex01_g,    9,  0.125 * _pi,  0,  1  },
    { "Test 2",                                         2,  -1.0,   1.0,         1.0,        ex04_f,    ex04_K,   ex04_g,    8,  0.500,        0,  1  },

    { "ex03 - Fox & Parker, Chapter 6 (6)",             2,  -1.0,   1.0,        -1.0,        ex03_f,    ex03_K,   ex03_g,    3,  0.500,        0,  0  },

    { "ex61 - Fox and Goodwin, (13)",                   2,   0.0,   _pi_div_2,  -_4_div_pi,  ex61_f,    ex61_K,   ex61_g,    8,  0.125 * _pi,  0,  1  },
    { "ex62 - Fox and Goodwin, (17) d = +1",            2,  -1.0,   1.0,         1.0 / _pi,  ex62_f1,   ex62_K,   ex62_g,    8,  0.250,        0,  1  },
    { "ex62 - Fox and Goodwin, (17) d = -1",            2,  -1.0,   1.0,        -1.0 / _pi,  ex62_f2,   ex62_K,   ex62_g,    8,  0.250,        0,  1  },
    { "ex63 - Fox and Goodwin, (49) case 1",            1,   0.0,   1.0,         1.0,        ex63_f1,   ex63_K,   ex63_g1,   1,  0.250,        0,  0  },
    { "ex63 - Fox and Goodwin, (49) case 2",            1,   0.0,   1.0,         1.0,        ex63_f2,   ex63_K,   ex63_g2,   1,  0.250,        0,  0  },
    { "ex64 - Fox and Goodwin, (51)",                   1,   0.0,   1.0,         1.0,        ex64_f,    ex64_K,   ex64_g,    5,  0.125,        0,  1  },

    { "ex91 - Elliott, (54)",                           2,  -1.0,   1.0,        -1.0,        ex91_f,    ex91_K,   ex91_g,   11,  0.200,        0,  1  },

    { "ex82 - Baker et al, (17)",                       1,   0.0,   1.0,         1.0,        ex82_f,    ex82_K,   ex82_g,    1,  0.125,        0,  0  },
    { "ex83 - Baker et al, (42)",                       1,   0.0,   1.0,         1.0,        ex83_f,    ex83_K,   ex83_g,    5,  0.0625,       1,  1  },

    { "ex71 - R. E. Scraton, (10) lambda = 2",          2,  -1.0,   1.0,        -2.0,        ex71_f1,   ex71_K,   ex71_g,    8,  0.250,        0,  1  },
    { "ex71 - R. E. Scraton, (10) lambda = 1/2",        2,  -1.0,   1.0,        -0.5,        ex71_f2,   ex71_K,   ex71_g,    8,  0.250,        0,  1  },
    { "ex72 - R. E. Scraton, (13)",                     2,  -1.0,   1.0,         0.5,        ex72_f,    ex72_K,   ex72_g,   10,  0.125,        0,  1  },

    { "ex12 - Delves and Walsh, pg 98",                 2,   0.0,   _pi_div_2,  -0.25,       ex12_f,    ex12_K,   ex12_g,    6,  0.125 * _pi,  0,  1  },
    { "ex14 - Delves and Walsh, pg 102",                2,   0.0,   1.0,         1.0,        ex14_f,    ex14_K,   ex14_g,   10,  0.125,        0,  1  },
        
    { "ex111 - Pilkington, (II)",                       2,  -1.0,   1.0,         1.0 / _pi,  ex111_f,   ex111_K,  ex111_g,  10,  0.250,        0,  1  },

    { "ex121 - Richardson, Introductory Example",       2,   0.0,   1.0,        -0.5,        ex121_f,   ex121_K,  ex121_g,   5,  0.125,        0,  1  },
    
    { "ex106 - Rahbar and Hashemizadeh, (12)",          2,   0.0,   1.0,         1.0,        ex106_f,   ex106_K,  ex106_g,  10,  0.125,        0,  1  },
    { "ex101 - Rahbar and Hashemizadeh, (13)",          2,  -1.0,   1.0,        -1.0,        ex101_f,   ex101_K,  ex101_g,   2,  0.250,        0,  0  },
    { "ex102 - Rahbar and Hashemizadeh, (14), (15)",    2,   0.0,   1.0,        -1.0,        ex102_f,   ex102_K,  ex102_g,  10,  0.125,        0,  1  },
    { "ex103 - Rahbar and Hashemizadeh, (16)",          2,   0.0,   1.0,         1.0,        ex103_f,   ex103_K,  ex103_g,  10,  0.125,        0,  1  },
    { "ex104 - Rahbar and Hashemizadeh, (17)",          2,  -1.0,   1.0,        -0.5,        ex104_f,   ex104_K,  ex104_g,  10,  0.250,        0,  1  },
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
        bf  = ChebyshevCoeff( N, lower, upper, examples[i].f, REGULAR, ZERO );

        if ( b )
            vector_print( stdout, "x\n", b,  number_fmt );
        else
            printf( "x - no solution\n\n" );

        vector_print( stdout, "f(x)\n", bf, number_fmt );

        if ( b && examples[i].output ) {
            double  x   = lower;
            double  dx  = examples[i].dx;

            printf( "\tx\tf(x)\t\tTn(x)\t\tX(x)\t\tf(x)-X(x)\tTn(x)-X(x)\n\n" );
        
            for ( ; x <= upper; x += dx ) {
                double      y   = examples[i].f( x );
                double      yf  = ChebyshevEval( x, bf, N, lower, upper, REGULAR );
                double      yx  = ChebyshevEval( x, b,  N, lower, upper, REGULAR );
                const char  *fmt;

                if ( examples[i].K == ex62_K )
                    fmt =  "\t%+5.3f\t%+.4f\t\t%+.4f\t\t%+.4f\t\t%+.4f\t\t%+.4f\n";
                else if ( examples[i].K == ex91_K )
                    fmt =  "\t%+5.3f\t%+.6f\t%+.6f\t%+.6f\t%+.6f\t%+.6f\n";
                else
                    fmt =  "\t%+5.3f\t%+.6f\t%+.6f\t%+.6f\t%+.6f\t%+.6f\n";

                printf( fmt, x, y, yf, yx, ( y - yx ), ( yf - yx ) );
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
