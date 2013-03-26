c
c INTEGRAL - generate an approximate solution to a Fredholm integral equation
c		of the second kind.
c
c		 |* b
c		 |
c	f(x) - l | K(x,y) * f(y) dy = g(x),  a <= x <= b
c		 |
c		*| a
c
c written: Richard.Marejka@Canada.Sun.COM
c
c ident	: @(#) integral.f 1.5 94/02/08
c
	external func
	external fern
	implicit real*8(a-h,o-z)
	logical flag
	dimension d(100,100), b(100), x(100), a(100,100), value(100), dot(100)
	pi	= 4.d0 * datan( 1.d0 )
	ns	= 1
	ambda	= 1.d0
	mm	= 100
	elb	= 0.d0
	ub	= pi
	num	= 11
	write(06,11)
	read(05,12) mm
	call cheff( b, mm, ns, elb, ub, func )
	b(1) = b(1) / 2.d0
	write(06,1)(i,b(i),i=1,mm)
	call cheff2( mm, a, ns, elb, ub, fern )
	if ( mm .gt. 10 ) go to 10
	write(06,2)((a(m,n),n=1,mm),m=1,mm)
10	continue
	call matrix( mm, a, d, ambda, elb, ub )
	if ( mm .gt. 10 ) go to 20
	write(06, 3)((d(m,n),n=1,mm),m=1,mm)
20	continue
	call solve( mm, d, b, x, flag )
	if ( .not. flag ) go to 7
	write(06,4)(i,x(i),i=1,mm)
	x(1) = 2.d0 * x(1)
	call points( elb, ub, value, dot, x, ns, mm, num )
	write(06,5) num
	write(06,6)(dot(i), value(i), i=1,num)
	go to 9
7	write(06,8)
9	continue
1	format('Chebyshev coefficients of g(x)',//,10(2x,i2,3x,f22.16,/),/)
2	format('Chebyshev coefficients of K(x,y)',//,10(10(2x,f10.6),/),/)
3	format('matrix to be solved',//,10(10(2x,f10.6),/),/)
4	format('Chebyshev coefficients of f(x)',//,10(2x,i2,3x,f22.16,/),//)
5	format(2x,'value of f(x) at ',i3,' equally spaced points in the interval',//)
6	format(2x,'abscissa',10x,'ordinate',//,100(2x,f10.7,5x,f22.16,/),/)
8	format('there is not a unique solution to this matrix',/,/)
11	format('Enter the order of the solution')
12	format(i3)
	stop
	end
c g(x)	
	double precision function func(x)
	implicit real*8(a-h,o-z)
	func	= ( 1.d0 + 2.d0 * datan( 1.d0 ) ) * cos( x )
	return
	end
c K(x,y)
	double precision function fern( x, y )
	implicit real*8(a-h,o-z)
	fern	= dcos( x ) * dcos( y )
	return
	end
c
c CHEFF - generate the Chebyshev coefficients to f(x)
c
c	real*8	a(mm)	resulting coefficients		out		
c	integer mm	dimension of a			in
c	integer neve	even/odd/none switch		in
c	real*8	elb	lower bound of range		in
c	real*8	ub	upper bound of range		in
c
	subroutine cheff( a, mm, neve, elb, ub, f )
	implicit real*8(a-h,o-z)
	dimension a(mm), b(100)
	re(x)	= 0.5 * ( ub + elb + ( ub-elb ) * x )
	adam	= neve
	eve	= neve
	if ( neve .lt. 3 ) go to 1
	adam	= 2.d0
	eve	= 1.d0
1	xmm	= mm
	pi	= 4.d0 * datan( 1.d0 )
	um	= pi / ( adam * xmm - eve )
	fo	= 2.d0 * adam / ( adam * xmm - eve )
	do 4 j = 1, mm
4		b(j)	= f( re( dcos( ( j - 1 ) * um ) ) )
	do 3 i = 1, mm
		xi	= i
		sum	= 0.d0
		bit	= 0.5
		ui	= ( adam * xi - eve ) * um
		do 2 j = 1, mm
			add	= bit * b(j) * dcos( ( j - 1 ) * ui )
			sum	= sum + add
2			bit	= 1.d0
		if ( neve .eq. 3 ) go to 3
		sum	= sum - 0.5 * add
3		a(i)	= fo * sum
	a(mm)	= 0.5 * a(mm)
	return
	end
c
c CHESUM - evaluate a Chebyshev series at a point
c
c	real*8	xr	point at which to evaluate	in
c	real*8	p	Chebyshev coefficients		in
c	integer	mm	dimension of p			in
c	integer	ns	even/odd/none switch		in
c	real*8	elb	lower bound of range		in
c	real*8	ub	upper bound of range		in
c
	double precision function chesum( xr, p, mm, ns, elb, ub )
	implicit real*8(a-h,o-z)
	dimension p(mm)
	x	= ( 2.d0 * xr - ( ub + elb ) ) / ( ub - elb )
	y	= x
	go to ( 2, 1, 1 ), ns
1	y	= 2.d0 * x * x - 1.d0
2	c0	= 0.d0
	c1	= 0.d0
	c2	= 0.d0
	emul	= 2.d0 * y
	do 3 i = 1, mm
		ii	= mm + 1 - i
		c2	= c1
		c1	= c0
3		c0	= emul * c1 - c2 + p(ii)
	go to ( 4, 4, 5 ), ns
4	chesum	= 0.5d0 * ( c0 - c2 )
	go to 6
5	chesum	= x * ( c0 - c1 )
6	continue
	return 
	end
c
c CHEFF2 - generate the Chebyshev coefficients to K(x,y)
c
c	integer	mm	dimension of a matrix (square)	in
c	real*8	a	Cheyshev coefficient matrix	out
c	integer	ns	even/odd/none switch		in
c	real*8	elb	lower bound of range		in
c	real*8	ub	upper bound of range		in
c
	subroutine cheff2( mm, a, ns, elb, ub, f )
	implicit real*8(a-h,o-z)
	dimension a(100,100), b(100), v(100), xyk(100,100)
	value(i)	= dcos( ( i - 1 ) * pi / ( mm - 1 ) )
	re(x)		= 0.5d0 * ( ub + elb + ( ub - elb ) * x )
	pi		= 4.d0 * datan( 1.d0 )
	nn		= mm
	do 60 m = 1,mm
		do 70 n = 1,mm
			xyk(m,n)	= f( re ( value( m ) ), re( value( n ) ) )
70		continue
60	continue

	do 20 m = 1,mm
 		xm	= re( value( m ) )
		do 20 n = 1,mm
			yn	= re( value( n ) )
			do 30 k = 1,mm
				do 40 j=1,mm
					v(j)	= xyk(k,j)
40				continue
				v(mm)	= 0.5d0 * v(mm)
				b(k)	= chesum( yn, v, nn, ns, elb, ub )
30			continue
			b(mm)	= 0.5d0 * b(mm)
			a(m,n)	= 4.d0 * chesum( xm, b, nn, ns, elb, ub, ) / ( ( mm - 1 ) **2 )
20	continue
	do 50 j = 1, mm
		a(1,j)	= a(1,j)  / 2.d0
		a(j,1)	= a(j,1)  / 2.d0
		a(mm,j)	= a(mm,j) / 2.d0
		a(j,mm)	= a(j,mm) / 2.d0
50	continue
	return
	end
c
c MATRIX - integrate the right-hand side and add the identity matrix
c
c	integer mm	dimension of a					in
c	real*8	b	Chebyshev coefficient matrix for K(x,y)		in
c	real*8	d	Resultant matrix				out
c	real*8	ambda	lambda coefficient				in
c	real*8	elb	lower bound of range				in
c	real*8	ub	upper bound of range				in
c
	subroutine matrix( mm, b, d, ambda, elb, ub )
	implicit real*8(a-h,o-z)
	dimension c(100,100), b(100,100), d(100,100)
	data c/10000*0.d0/
	fact	= (ub - elb ) / ( 2.d0 * ambda )
	do 10 m = 1, mm
		if ( mod( m, 2 ) .eq. 0 ) go to 2
		k	= 1
		go to 3
2		k	= 2
3		continue
		do 10 n = k, mm, 2
			c(m,n)	= 0.5 * ( 1/(m+n-1.d0) - 1/(m+n-3.d0) + 1/(m-n+1.d0) - 1/(m-n-1.d0) )
10	continue
	do 20 m = 1, mm
		do 20 l = 1, mm
			sum	= 0.d0
			do 40 n = 1, mm
				sum	= sum + c(l,n) * b(m,n)
40			continue
			d(m,l)	= fact * sum
20	continue
	do 30 i = 1, mm
		d(i,i)	= d(i,i) + 1.d0
30	continue
	return
	end
c
c SOLVE - Gauss-Jordan solution to a system of linear equations.
c
c		[A][x] = [b]
c
c	integer	mm	dimension of matrix and vectors		in
c	real*8	a	sqaure matrix to solve			in/out
c	real*8	b	particular solution			in/out
c	real*8	x	result vector				out
c	logical	flag	status of solution			out
c
	subroutine solve( mm, a, b, x, flag )
	implicit real*8(a-h,o-z)
	logical flag
	dimension a(100,100), b(100), x(100)
	do 10 k = 1, mm-1
		do 45 i = k+1, mm
			fact	= a(i,k) / a(k,k)
			do 50 j = k+1, mm
				a(i,j)	= a(i,j) - a(k,j) * fact
50			continue
			b(i)	= b(i) - b(k) * fact
45		continue
10	continue
	do 55 i = 1, mm
		k	= mm + 1 - i
		j	= k + 1
		sum	= 0.d0
65		if ( j .gt. mm ) go to 60
		sum	= sum + a(k,j) * x(j)
		j	= j + 1
		go to 65
60		continue
		x(k)	= ( b(k) - sum ) / a(k,k)
55	continue
	flag	= .true.
	go to 70
99	flag	= .false.
70	continue
	return
	end
c
c POINTS - evaluate a Chebyshev series at a number of uniform points
c
c	real*8	elb		lower bound of range		in
c	real*8	ub		upper bound of range		in
c	real*8	value(num)	Y coordinates over range	out
c	real*8	dot(num)	X coordinates over range	out
c	real*8	x(mm)		Chebysehv coefficients		in
c	integer ns		even/odd/none switch		in
c	integer mm		dimension of x			in
c	integer	num		number of points to generate	in
c
	subroutine points( elb, ub, value, dot, x, ns, mm, num )
	implicit real*8(a-h,o-z)
	dimension value(num), x(100), dot(num)
	dx	= ( ub - elb ) / ( num - 1 )
	do 10 i = 1, num
		xr		= elb + ( i - 1 ) * dx
		value(i)	= chesum( xr, x, mm, ns, elb, ub )
		dot(i)		= xr
10	continue
	return 
	end
