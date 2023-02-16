      subroutine jmap(n,a,no)
c...
c...  subroutine jmap prints the jacobian matrix of an nth-order system of
c...  algebraic equations or first-order ordinary differential equations
c...
c...  argument list
c...
c...     n       number of algebraic or first-order ordinary differen-
c...             tial equations (odes) for which the jacobian matrix is
c...             to be mapped, i.e., the order of the algebraic or ode
c...             system (input)
c...
c...     a       two-dimensional array containing the jacobian matrix
c...             of the nth-order algebraic or ode system (output)
c...
c...     no      logical unit number for output from jmap.
c...
c...  from this point on, the discussion is in terms of first-order
c...                                                  -
c...  subroutine jmap prints a two-dimensional map of j with the numbers
c...  0 to 9 indicating the relative order-of-magnitude of the individ-
c...  ual elements of the matrix.  the values of the row subscript i are
c...  printed down the left side of the map and the values of the column
c...  subscript j are printed across the top.  the map is printed in
c...  sections 100 columns wide.  thus if the differential equation
c...  system is greater than 100th-order, successive sections of the map
c...  will be printed vertically.  these can then be joined together to
c...  make up the complete map.
c...
      implicit none
      integer n, no, n1, n10, n100, i, nsect, nsec, jl, ju, jup1
      integer jcol, jd, irow, js
      real a, scale
      character sym, symbol
c...  variable dimensions for the arrays passed to subroutine jmap
      dimension a(n,n)
c...
c...  absolute dimensions for the arrays used in subroutine jmap
      dimension sym(11),symbol(100),n1(100),n10(100),n100(100)
c...
c...  define the symbols used in printing the jacobian map
      data sym/" ","e","1","2","3","4","5","6","7","8","9"/
c...
c...  set the scale factor.
      data scale /1./
c...
c...  set the arrays used to print the units (n1), tens (n10) and
c...  hundreds (n100) places in the column subscripts of the jacobian
c...  map
      data n1(1),n1(2),n1(3),n1(4),n1(5),n1(6),n1(7),n1(8),n1(9),n1(10)/
     1     0,1,2,3,4,5,6,7,8,9/
      data n10/10*0,10*1,10*2,10*3,10*4,10*5,10*6,10*7,10*8,10*9/
      do 1 i=11,100
         n1(i)=n1(i-10)
1     continue
c...                                   -
c...  compute the number of sections of the printed jacobian map
      nsect=n/100+1
c...
c...  step through the calculation and printing of the individual sec-
c...  tions of the printed jacobian map
      do 4 nsec=1,nsect
c...
c...  set the lower and upper limits of the column subscripts of the
c...  current section to be printed
         jl=(nsec-1)*100
         ju=jl+99
c...
c...  if the upper limit exceeds the order of the system (v. equations
c...  (1)), set the upper limit to the order of the system
         if(n.lt.ju)ju=n
c...
c...  if the printed map contains more than one section, branch to 11
         if(nsec.ne.1)go to 11
c...
c...  the map contains only one section
         jl=1
         ju=99
         if(n.lt.ju)ju=n
         jup1=ju+1
c...
c...  begin the printing of the one-section map
         write(no,908)
c...
c...  if the system is less than 11th-order branch to 5 to skip printing
c...  the tens place in the column subscript
         if(jup1.lt.11)go to 5
c...
c...  print the tens and ones places in the column subscripts
         write(no,901)(n10(jcol),jcol=11,jup1)
5        write(no,902)( n1(jcol),jcol= 2,jup1)
         go to 8
c...
c...  for a map with more than one section, the column subscript
c...  increases from the lower column limit to the upper limit for
c...  a given section
11       jd=ju-jl+1
         do 12 jcol=1,jd
c...
c...  store the hundreds place for the current section of the map
         n100(jcol)=nsec-1
12       continue
c...
c...  begin the printing of the current section of the map
         if(nsec.eq.1)write(no,908)
c...
c...  print the hundreds, tens, unit places in the column subscripts of
c...  the current section of the map
         write(no,904)(n100(jcol),jcol=1,jd)
         write(no,904)( n10(jcol),jcol=1,jd)
         write(no,904)(  n1(jcol),jcol=1,jd)
c...
c...  determine which symbol is to be printed in each position of the
c...  current section of the jacobian map
8     do 9 irow=1,n
         do 10 jcol=jl,ju
            js=jcol-jl+1
            if(abs(a(irow,jcol)).le.1.e-100*scale)symbol(js)=sym(1)
            if(abs(a(irow,jcol)).gt.1.e-100*scale)symbol(js)=sym(2)
            if(abs(a(irow,jcol)).gt.0.00010*scale)symbol(js)=sym(3)
            if(abs(a(irow,jcol)).gt.0.00100*scale)symbol(js)=sym(4)
            if(abs(a(irow,jcol)).gt.0.01000*scale)symbol(js)=sym(5)
            if(abs(a(irow,jcol)).gt.0.10000*scale)symbol(js)=sym(6)
            if(abs(a(irow,jcol)).gt.1.00000*scale)symbol(js)=sym(7)
            if(abs(a(irow,jcol)).gt.10.0000*scale)symbol(js)=sym(8)
            if(abs(a(irow,jcol)).gt.100.000*scale)symbol(js)=sym(9)
            if(abs(a(irow,jcol)).gt.1000.00*scale)symbol(js)=sym(10)
            if(abs(a(irow,jcol)).gt.10000.0*scale)symbol(js)=sym(11)
10       continue
c...
c...  print the current section of the jacobian matrix, one row at a
c...  time
         jd=ju-jl+1
         if(nsec.eq.1)write(no,906)irow,(symbol(js),js=1,jd)
         if(nsec.gt.1)write(no,903)irow,(symbol(js),js=1,jd)
9     continue
c...
c...  return to 4 and repeat the calculation and printing of the next
c...  section of the jacobian map if required
4     continue
      return
901   format(1h ,/,20x,90i1)
902   format(11x,99i1)
903   format(i8,2x,100a1)
904   format(10x,100i1)
906   format(i8,3x,100a1)
908   format(1h ,/,
     1 58h dependent variable column index j (for yj) is printed hor,
     2 58hizontally                                                 ,//,
     3 58h derivative row index i (for dyi/dt = fi(y1,y2,...,yj,...,,
     4 58hyn) is printed vertically                                 ,//,
     5 58h jacobian matrix element in the map with indices i,j is fo,
     6 58hr pfi/pyj where p denotes a partial derivative            ,//)
      end
