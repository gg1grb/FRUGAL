      subroutine sslugm(n, b, x, nelt, ia, ja, a, isym, nsave,
     $     itol, tol, itmax, iter, err, ierr, iunit, rwork, lenw,
     $     iwork, leniw )
      implicit none



      integer   locrb, locib
      parameter (locrb=1, locib=11)

      integer  n, nelt, ia(nelt), ja(nelt), isym, nsave, itol
      integer  leniw,itmax, iter, ierr, iunit, lenw, iwork(leniw)
      integer   lociw, locl, locnc, locju, locnr, myitol
      integer   locw, locrgw, locdin, locu, locjl, lociu, locil
      integer   j, locigw, jbgn, jend, icol, nl, nu
      real     b(n), x(n), a(nelt), tol, err, rwork(lenw)
      external ssmv, sslui

c
c         change the slap input matrix ia, ja, a to slap-column format.
c***first executable statement  sslugm
      ierr = 0
      err  = 0.0
      if( nsave.le.1 ) then
         ierr = 3
         return
      endif
c      write(*,*) 'sslugm:  changing data structure'
      call ss2y( n, nelt, ia, ja, a, isym )
c
c         count number of non-zero elements preconditioner ilu matrix.
c         then set up the work arrays.  we assume maxl=kmp=nsave.
c      write(*,*) 'sslugm:  counting non-zeros'
      nl = 0
      nu = 0
      do 20 icol = 1, n
c         don't count diagonal.
         jbgn = ja(icol)+1
         jend = ja(icol+1)-1
         if( jbgn.le.jend ) then

            do 10 j = jbgn, jend
               if( ia(j).gt.icol ) then
                  nl = nl + 1
                  if( isym.ne.0 ) nu = nu + 1
               else
                  nu = nu + 1
               endif
 10         continue
         endif
 20   continue
c         
c      write(*,*) 'nl,nu =', nl,nu
      locigw = locib
      locil = locigw + 20
      locjl = locil + n+1
      lociu = locjl + nl
      locju = lociu + nu
      locnr = locju + n+1
      locnc = locnr + n
      lociw = locnc + n
c
      locl = locrb
      locdin = locl + nl
      locu = locdin + n
      locrgw = locu + nu
      locw = locrgw + 1+n*(nsave+6)+nsave*(nsave+3)
c
c         check the workspace allocations.
c      write(*,*) 'sslugm:  checking the workspace allocations.'
      call schkw( 'sslugm', lociw, leniw, locw, lenw, ierr, 
     &		iter, err )
      if( ierr.ne.0 ) return
c
      iwork(1) = locil
      iwork(2) = locjl
      iwork(3) = lociu
      iwork(4) = locju
      iwork(5) = locl
      iwork(6) = locdin
      iwork(7) = locu
      iwork(9) = lociw
      iwork(10) = locw
c
c         compute the incomplete lu decomposition.
c      write(*,*) 'sslugm:  computing ilu factors'
      call ssilus( n, nelt, ia, ja, a, isym, nl, iwork(locil),
     $     iwork(locjl), rwork(locl), rwork(locdin), nu, 
     $	   iwork(lociu),
     $     iwork(locju), rwork(locu), iwork(locnr), iwork(locnc) )
c         
c         perform the incomplet lu preconditioned generalized minimum
c         residual iteration algorithm.  the following sgmres 
c         defaults are used maxl = kmp = nsave, jscal = 0,
c         jpre = -1, nrmax = itmax/nsave
      iwork(locigw  ) = nsave
      iwork(locigw+1) = nsave
      iwork(locigw+2) = 0
      iwork(locigw+3) = -1
      iwork(locigw+4) = itmax/nsave
      myitol = 0
c      
c      write(*,*) 'sslugm:  starting gmres iteration, itmax',itmax

      call sgmres( n, b, x, nelt, ia, ja, a, isym, ssmv, sslui,
     $     myitol, tol, itmax, iter, err, ierr, iunit, rwork, 
     $	   rwork,
     $     rwork(locrgw), lenw-locrgw, iwork(locigw), 20,
     $     rwork, iwork )
c
c	write(*,*) ' iter, itmax = ', iter, itmax

      if( iter.gt.itmax ) ierr = 2


      return
c------------- last line of sslugm follows ----------------------------
      end


      subroutine ss2y(n, nelt, ia, ja, a, isym )
      implicit none


c***begin prologue  ss2y
c***date written   871119   (yymmdd)
c***revision date  881213   (yymmdd)
c***category no.  d2a4, d2b4
c***keywords  library=slatec(slap),
c             type=single precision(ss2y-s),
c             linear system, slap sparse
c***author  seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-300
c             livermore, ca 94550 (415) 423-3141
c             seager@lll-crg.llnl.gov
c***purpose  slap triad to slap column format converter.
c            routine to convert from the slap triad to slap column
c            format.
c***description
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym
c     real    a(nelt)
c
c     call ss2y( n, nelt, ia, ja, a, isym )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :inout    integer ia(nelt).
c ja     :inout    integer ja(nelt).
c a      :inout    real a(nelt).
c         these arrays should hold the matrix a in either the slap
c         triad format or the slap column format.  see "long 
c         description", below.  if the slap triad format is used
c         this format is translated to the slap column format by
c         this routine.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all nonzero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the lower
c         triangle of the matrix is stored.
c
c *precision:           single precision
c
c       the sparse linear algebra package (slap) utilizes two matrix
c       data structures: 1) the  slap triad  format or  2)  the slap
c       column format.  the user can hand this routine either of the
c       of these data structures.  if the slap triad format is give
c       as input then this routine transforms it into slap column
c       format.  the way this routine tells which format is given as
c       input is to look at ja(n+1).  if ja(n+1) = nelt+1 then we
c       have the slap column format.  if that equality does not hold
c       then it is assumed that the ia, ja, a arrays contain the 
c       slap triad format.
c
c***references  (none)
c***routines called  qs2i1r
c***end prologue  ss2y
      integer n, nelt, ia(nelt), ja(nelt), isym
      integer i, itemp, ibgn, iend, icol, j
      real    a(nelt), temp
c
c         check to see if the (ia,ja,a) arrays are in slap column 
c         format.  if it's not then transform from slap triad.
c***first executable statement  ss2lt
      if( ja(n+1).eq.nelt+1 ) return
c
c         sort into ascending order by column (on the ja array).
c         this will line up the columns.
c
      call qs2i1r( ja, ia, a, nelt, 1 )
c         
c         loop over each column to see where the column indicies change 
c         in the column index array ja.  this marks the beginning of the 
c         next column.
c         
cvd$r novector
      ja(1) = 1
      do 20 icol = 1, n-1
         do 10 j = ja(icol)+1, nelt
            if( ja(j).ne.icol ) then
               ja(icol+1) = j
               goto 20
            endif
 10      continue
 20   continue
      ja(n+1) = nelt+1
c         
c         mark the n+2 element so that future calls to a slap routine 
c         utilizing the ysmp-column storage format will be able to tell.
c         
      ja(n+2) = 0
c
c         now loop thru the ia(i) array making sure that the diagonal
c         matrix element appears first in the column.  then sort the
c         rest of the column in ascending order.
c
      do 70 icol = 1, n
         ibgn = ja(icol)
         iend = ja(icol+1)-1
         do 30 i = ibgn, iend
            if( ia(i).eq.icol ) then
c         swap the diag element with the first element in the column.
               itemp = ia(i)
               ia(i) = ia(ibgn)
               ia(ibgn) = itemp
               temp = a(i)
               a(i) = a(ibgn)
               a(ibgn) = temp
               goto 40
            endif
 30      continue
 40      ibgn = ibgn + 1
         if( ibgn.lt.iend ) then
            do 60 i = ibgn, iend
               do 50 j = i+1, iend
                  if( ia(i).gt.ia(j) ) then
                     itemp = ia(i)
                     ia(i) = ia(j)
                     ia(j) = itemp
                     temp = a(i)
                     a(i) = a(j)
                     a(j) = temp
                  endif
 50            continue
 60         continue
         endif
 70   continue

      return
c------------- last line of ss2y follows ----------------------------
      end



      subroutine schkw( name, lociw, leniw, locw, lenw,
     $     ierr, iter, err )
      implicit none

c***begin prologue  schkw
c***date written   880225   (yymmdd)
c***revision date  881213   (yymmdd)
c***category no.  r2
c***keywords  library=slatec(slap),
c             type=single precision(schkw-s),
c             slap, error checking, workspace checking
c***author  seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-300
c             livermore, ca 94550 (415) 423-3141
c             seager@lll-crg.llnl.gov
c***purpose  slap work/iwork array bounds checker.
c            this routine checks the work array lengths  and  inter-
c            faces to the slatec  error  handler  if  a  problem  is 
c            found.
c***description
c *usage:
c     character*(*) name
c     integer lociw, leniw, locw, lenw, ierr, iter
c     real    err
c
c     call schkw( name, lociw, leniw, locw, lenw, ierr, iter, err )
c
c *arguments:
c name   :in       character*(*).
c         name of the calling routine.  this is used in the output
c         message, if an error is detected.
c lociw  :in       integer.
c         location of the first free element in the integer workspace
c         array.
c leniw  :in       integer.
c         length of the integer workspace array.
c locw   :in       integer.
c         location of the first free element in the real workspace
c         array.
c lenrw  :in       integer.
c         length of the real workspace array.
c ierr   :out      integer.
c         return error flag.
c               ierr = 0 => all went well.
c               ierr = 1 => insufficient storage allocated for 
c                           work or iwork.
c iter   :out      integer.
c         set to 0 if an error is detected.
c err    :out      real.
c         set to a very large number if an error is detected.
c
c *precision:           single precision
c
c***references  (none)
c***routines called  r1mach, xerrwv
c***end prologue  schkw
      character*(*) name
      character*72 mesg
      integer lociw, leniw, locw, lenw, ierr, iter
      real    err, r1mach
      external r1mach, xerrwv
c
c         check the integer workspace situation.
c***first executable statement  schkw
      ierr = 0
      if( lociw.gt.leniw ) then
         ierr = 1
         iter = 0
         err = r1mach(2)
         mesg = name // ': integer work array too short. '//
     $        ' iwork needs i1: have allocated i2.'
         call xerrwv( mesg, len(mesg), 1, 1, 2, lociw, leniw,
     $        0, 0.0, 0.0 )
      endif
c
c         check the real workspace situation.
      if( locw.gt.lenw ) then
         ierr = 1
         iter = 0
         err = r1mach(2)
         mesg = name // ': real work array too short. '//
     $        ' rwork needs i1: have allocated i2.'
         call xerrwv( mesg, len(mesg), 1, 1, 2, locw, lenw,
     $        0, 0.0, 0.0 )
      endif

      return
c------------- last line of schkw follows ----------------------------
      end


      subroutine sslui(n, b, x, nelt, ia, ja, a, isym, rwork, iwork)
      implicit none


c***begin prologue  sslui
c***date written   871119   (yymmdd)
c***revision date  881213   (yymmdd)
c***category no.  d2a4
c***keywords  library=slatec(slap),
c             type=single precision(sslui-s),
c             non-symmetric linear system solve, sparse, 
c             iterative precondition
c***author  greenbaum, anne, courant institute
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-300
c             livermore, ca 94550 (415) 423-3141
c             seager@lll-crg.llnl.gov
c***purpose  slap msolve for ldu factorization.
c            this routine  acts as an  interface between  the   slap
c            generic mslove calling convention and the routine  that 
c            actually computes:     -1
c                              (ldu)  b = x.
c***description
c       see the "description" of sslui2 for the gory details.
c***routines called  sslui2
c***end prologue  sslui
      integer n, nelt, ia(nelt), ja(nelt), isym, iwork(*)
      integer locl, locdin, locu, locju, locil, locjl, lociu
      real    b(n), x(n), a(nelt), rwork(*)
c
c         pull out the locations of the arrays holding the ilu
c         factorization.
c***first executable statement  sslui
      locil = iwork(1)
      locjl = iwork(2)
      lociu = iwork(3)
      locju = iwork(4)
      locl = iwork(5)
      locdin = iwork(6)
      locu = iwork(7)
c
c         solve the system lux = b
      call sslui2(n, b, x, iwork(locil), iwork(locjl), rwork(locl),
     $     rwork(locdin), iwork(lociu), iwork(locju), rwork(locu) )
c         

      return
c------------- last line of sslui follows ----------------------------
      end


      subroutine sslui2(n, b, x, il, jl, l, dinv, iu, ju, u )
      implicit none


c***begin prologue  sslui2
c***date written   871119   (yymmdd)
c***revision date  881213   (yymmdd)
c***category no.  d2a4
c***keywords  library=slatec(slap),
c             type=single precision(sslui2-s),
c             non-symmetric linear system solve, sparse, 
c             iterative precondition
c***author  greenbaum, anne, courant institute
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-300
c             livermore, ca 94550 (415) 423-3141
c             seager@lll-crg.llnl.gov
c***purpose  slap back solve for ldu factorization.
c            routine  to  solve a system of the form  l*d*u x  =  b,
c            where l is a unit  lower  triangular  matrix,  d  is  a 
c            diagonal matrix, and u is a unit upper triangular matrix.
c***description
c *usage:
c     integer n, il(n+1), jl(nl), iu(nu), ju(n+1)
c     real    b(n), x(n), l(nl), dinv(n), u(nu)
c
c     call sslui2( n, b, x, il, jl, l, dinv, iu, ju, u )
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c b      :in       real b(n).
c         right hand side.
c x      :out      real x(n).
c         solution of l*d*u x = b.
c nel    :in       integer.
c         number of non-zeros in the el array.
c il     :in       integer il(n+1).
c jl     :in       integer jl(nl).
c  l     :in       real     l(nl).
c         il, jl, l contain the unit  lower triangular factor of the
c         incomplete decomposition of some matrix stored in slap row
c         format.  the diagonal of ones *is* stored.  this structure
c         can   be   set up  by   the  ssilus routine.   see 
c         "description", below  for more   details about   the  slap
c         format.
c dinv   :in       real dinv(n).
c         inverse of the diagonal matrix d.
c nu     :in       integer.
c         number of non-zeros in the u array.     
c iu     :in       integer iu(n+1).
c ju     :in       integer ju(nu).
c u      :in       real     u(nu).
c         iu, ju, u contain the unit upper triangular factor  of the
c         incomplete decomposition  of  some  matrix stored in  slap
c         column format.   the diagonal of ones  *is* stored.   this
c         structure can be set up  by the ssilus routine.  see
c         "description", below   for  more   details about  the slap
c         format.
c
c *description:
c       this routine is supplied with  the slap package as a routine
c       to  perform  the  msolve operation  in   the  sir and   sbcg
c       iteration routines for  the  drivers ssilur and sslubc.   it
c       must  be called  via   the  slap  msolve  calling   sequence
c       convention interface routine sslui.
c         **** this routine itself does not conform to the ****
c               **** slap msolve calling convention ****
c
c       il, jl, l should contain the unit lower triangular factor of
c       the incomplete decomposition of the a matrix  stored in slap
c       row format.  iu, ju, u should contain  the unit upper factor
c       of the  incomplete decomposition of  the a matrix  stored in
c       slap column format this ilu factorization can be computed by
c       the ssilus routine.  the diagonals (which is all one's) are
c       stored.
c
c       =================== s l a p column format ==================
c       this routine  requires that  the matrix a  be stored in  the
c       slap column format.  in this format the non-zeros are stored
c       counting down columns (except  for the diagonal entry, which
c       must appear first in each  "column") and  are stored  in the
c       real array a.  in other words, for each column in the matrix
c       put the diagonal entry in a.  then put in the other non-zero
c       elements going down   the  column (except  the diagonal)  in
c       order.  the ia array holds the row  index for each non-zero.
c       the ja array holds the offsets into the ia, a arrays for the
c       beginning of   each    column.    that  is,    ia(ja(icol)),
c       a(ja(icol)) points to the beginning of the icol-th column in
c       ia and  a.  ia(ja(icol+1)-1),  a(ja(icol+1)-1) points to the
c       end  of   the icol-th  column.  note   that  we  always have
c       ja(n+1) = nelt+1, where  n  is the number of columns in  the
c       matrix and  nelt   is the number of non-zeros in the matrix.
c       
c       here is an example of the  slap column  storage format for a
c       5x5 matrix (in the a and ia arrays '|'  denotes the end of a 
c       column):
c
c           5x5 matrix      slap column format for 5x5 matrix on left.
c                              1  2  3    4  5    6  7    8    9 10 11
c       |11 12  0  0 15|   a: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
c       |21 22  0  0  0|  ia:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
c       | 0  0 33  0 35|  ja:  1  4  6    8  9   12
c       | 0  0  0 44  0|
c       |51  0 53  0 55|
c
c       ==================== s l a p row format ====================
c       this routine requires  that the matrix a  be  stored  in the
c       slap  row format.   in this format  the non-zeros are stored
c       counting across  rows (except for the diagonal  entry, which
c       must appear first in each "row") and  are stored in the real
c       array a.  in other words, for each row in the matrix put the
c       diagonal entry in  a.   then   put  in the   other  non-zero
c       elements   going  across the  row (except   the diagonal) in
c       order.   the  ja array  holds   the column   index for  each
c       non-zero.   the ia  array holds the  offsets into  the ja, a
c       arrays  for   the   beginning  of   each  row.   that    is,
c       ja(ia(irow)),  a(ia(irow)) points  to  the beginning  of the
c       irow-th row in ja and a.   ja(ia(irow+1)-1), a(ia(irow+1)-1)
c       points to the  end of the  irow-th row.  note that we always
c       have ia(n+1) =  nelt+1, where  n  is  the number of rows  in
c       the matrix  and nelt  is the  number   of  non-zeros in  the
c       matrix.
c       
c       here is an example of the slap row storage format for a  5x5
c       matrix (in the a and ja arrays '|' denotes the end of a row):
c
c           5x5 matrix         slap row format for 5x5 matrix on left.
c                              1  2  3    4  5    6  7    8    9 10 11
c       |11 12  0  0 15|   a: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
c       |21 22  0  0  0|  ja:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
c       | 0  0 33  0 35|  ia:  1  4  6    8  9   12
c       | 0  0  0 44  0|  
c       |51  0 53  0 55|  
c
c       with  the slap  format  the "inner  loops" of  this  routine
c       should vectorize   on machines with   hardware  support  for
c       vector gather/scatter operations.  your compiler may require
c       a  compiler directive  to  convince   it that there  are  no
c       implicit vector  dependencies.  compiler directives  for the
c       alliant fx/fortran and cri cft/cft77 compilers  are supplied
c       with the standard slap distribution.
c
c *precision:           single precision
c *see also:
c       ssilus
c***references  (none)
c***routines called  (none)
c***end prologue  sslui2
      integer n, il(*), jl(*), iu(*), ju(*)
      integer i, icol, irow, j, jbgn, jend
      real    b(n), x(n), l(*), dinv(n), u(*)
c         
c         solve  l*y = b,  storing result in x, l stored by rows.
c***first executable statement  sslui2
      do 10 i = 1, n
         x(i) = b(i)
 10   continue
      do 30 irow = 2, n
         jbgn = il(irow)
         jend = il(irow+1)-1
         if( jbgn.le.jend ) then
clll. option assert (nohazard)
cdir$ ivdep
cvd$ assoc
cvd$ nodepchk
            do 20 j = jbgn, jend
               x(irow) = x(irow) - l(j)*x(jl(j))
 20         continue
         endif
 30   continue
c         
c         solve  d*z = y,  storing result in x.
      do 40 i=1,n
         x(i) = x(i)*dinv(i)
 40   continue
c         
c         solve  u*x = z, u stored by columns.
      do 60 icol = n, 2, -1
         jbgn = ju(icol)
         jend = ju(icol+1)-1
         if( jbgn.le.jend ) then
clll. option assert (nohazard)
cdir$ ivdep
cvd$ nodepchk
            do 50 j = jbgn, jend
               x(iu(j)) = x(iu(j)) - u(j)*x(icol)
 50         continue
         endif
 60   continue
c         

      return
c------------- last line of sslui2 follows ----------------------------
      end


      subroutine ssilus(n, nelt, ia, ja, a, isym, nl, il, jl,
     $     l, dinv, nu, iu, ju, u, nrow, ncol)
      implicit none



c***begin prologue  ssilus
c***date written   871119   (yymmdd)
c***revision date  881213   (yymmdd)
c***category no.  d2a4, d2b4
c***keywords  library=slatec(slap),
c             type=single precision(ssilus-s),
c             non-symmetric linear system, sparse, 
c             iterative precondition, incomplete lu factorization
c***author  greenbaum, anne, courant institute
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-300
c             livermore, ca 94550 (415) 423-3141
c             seager@lll-crg.llnl.gov
c***purpose  incomplete lu decomposition preconditioner slap set up.
c            routine to generate the incomplete ldu decomposition of a 
c            matrix.  the  unit lower triangular factor l is stored by 
c            rows and the  unit upper triangular factor u is stored by 
c            columns.  the inverse of the diagonal matrix d is stored.
c            no fill in is allowed.
c***description
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym
c     integer nl, il(n+1), jl(nl), nu, iu(n+1), ju(nu)
c     integer nrow(n), ncol(n)
c     real    a(nelt), l(nl), u(nu), dinv(n)
c
c     call ssilus( n, nelt, ia, ja, a, isym, nl, il, jl, l, 
c    $    dinv, nu, iu, ju, u, nrow, ncol )
c
c *arguments:
c n      :in       integer

c         order of the matrix.
c nelt   :in       integer.
c         number of elements in arrays ia, ja, and a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       real a(nelt).
c         these arrays should hold the matrix a in the slap column
c         format.  see "description", below. 
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all nonzero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the lower 
c         triangle of the matrix is stored.
c nl     :out      integer.
c         number of non-zeros in the el array.
c il     :out      integer il(n+1).
c jl     :out      integer jl(nl).
c l      :out      real     l(nl).
c         il, jl, l  contain the unit ower  triangular factor of  the
c         incomplete decomposition  of some  matrix stored  in   slap
c         row format.     the   diagonal  of ones  *is*  stored.  see
c         "description", below for more details about the slap format.
c nu     :out      integer.
c         number of non-zeros in the u array.     
c iu     :out      integer iu(n+1).
c ju     :out      integer ju(nu).
c u      :out      real     u(nu).
c         iu, ju, u contain   the unit upper triangular factor of the
c         incomplete  decomposition    of some matrix  stored in slap
c         column  format.   the diagonal of ones   *is*  stored.  see 
c         "description", below  for  more  details  about  the   slap 
c         format.
c nrow   :work     integer nrow(n).
c         nrow(i) is the number of non-zero elements in the i-th row
c         of l.
c ncol   :work     integer ncol(n).
c         ncol(i) is the number of non-zero elements in the i-th 
c         column of u.
c
c *description
c       il, jl, l should contain the unit  lower triangular factor of
c       the incomplete decomposition of the a matrix  stored in slap
c       row format.  iu, ju, u should contain  the unit upper factor
c       of the  incomplete decomposition of  the a matrix  stored in
c       slap column format this ilu factorization can be computed by
c       the ssilus routine.  the diagonals (which is all one's) are
c       stored.
c
c *precision:           single precision
c *see also:
c       silur
c***references  1. gene golub & charles van loan, "matrix computations",
c                 john hopkins university press; 3 (1983) ibsn 
c                 0-8018-3010-9.
c***routines called  (none)
c***end prologue  ssilus
      integer n, nelt, ia(nelt), ja(nelt), isym, nl, il(nl), jl(nl)
      integer nu, iu(nu), ju(nu), nrow(n), ncol(n)
      integer i, ibgn, iend, j, jbgn, jend, icol, irow
      integer itemp , jtemp, indx, indx1, indx2, indxr1, indxr2
      integer indxc2, indxc1, k, kc, kr
      real    a(nelt), l(nl), dinv(n), u(nu), temp

c         
c         count number of elements in each row of the lower triangle.
c***first executable statement  ssilus
      !write(*,*) 'ssilus:  counting'
      do 10 i=1,n
         nrow(i) = 0
         ncol(i) = 0
 10   continue
cvd$r noconcur
cvd$r novector
      do 30 icol = 1, n
         jbgn = ja(icol)+1
         jend = ja(icol+1)-1
         if( jbgn.le.jend ) then
            do 20 j = jbgn, jend
               if( ia(j).lt.icol ) then
                  ncol(icol) = ncol(icol) + 1
               else
                  nrow(ia(j)) = nrow(ia(j)) + 1
                  if( isym.ne.0 ) ncol(ia(j)) = ncol(ia(j)) + 1
               endif
 20         continue
         endif
 30   continue
      ju(1) = 1
      il(1) = 1
      do 40 icol = 1, n
         il(icol+1) = il(icol) + nrow(icol)
         ju(icol+1) = ju(icol) + ncol(icol)
         nrow(icol) = il(icol)
         ncol(icol) = ju(icol)
 40   continue
c         
c         copy the matrix a into the l and u structures.
      !write(*,*) 'ssilus:  copying'
      do 60 icol = 1, n
         dinv(icol) = a(ja(icol))
         jbgn = ja(icol)+1
         jend = ja(icol+1)-1
         if( jbgn.le.jend ) then
            do 50 j = jbgn, jend
               irow = ia(j)
               if( irow.lt.icol ) then
c         part of the upper triangle.
                  iu(ncol(icol)) = irow
                  u(ncol(icol)) = a(j)
                  ncol(icol) = ncol(icol) + 1
               else
c         part of the lower triangle (stored by row).
                  jl(nrow(irow)) = icol
                  l(nrow(irow)) = a(j)
                  nrow(irow) = nrow(irow) + 1
                  if( isym.ne.0 ) then
c         symmetric...copy lower triangle into upper triangle as well.
                     iu(ncol(irow)) = icol
                     u(ncol(irow)) = a(j)
                     ncol(irow) = ncol(irow) + 1
                  endif
               endif
 50         continue
         endif
 60   continue
c
c         sort the rows of l and the columns of u.
      !write(*,*) 'ssilus:  sorting'
      do 110 k = 2, n
         jbgn = ju(k)
         jend = ju(k+1)-1
         if( jbgn.lt.jend ) then
            do 80 j = jbgn, jend-1
               do 70 i = j+1, jend
                  if( iu(j).gt.iu(i) ) then
                     itemp = iu(j)
                     iu(j) = iu(i)
                     iu(i) = itemp
                     temp = u(j)
                     u(j) = u(i)
                     u(i) = temp
                  endif
 70            continue
 80         continue
         endif
         ibgn = il(k)
         iend = il(k+1)-1
         if( ibgn.lt.iend ) then
            do 100 i = ibgn, iend-1
               do 90 j = i+1, iend
                  if( jl(i).gt.jl(j) ) then
                     jtemp = ju(i)
                     ju(i) = ju(j)
                     ju(j) = jtemp
                     temp = l(i)
                     l(i) = l(j)
                     l(j) = temp
                  endif
 90            continue
 100        continue
         endif
 110  continue
c
c         perform the incomplete ldu decomposition.
      !write(*,*) 'ssilus:  factoring'
      do 300 i=2,n
c         
c           i-th row of l
         indx1 = il(i)
         indx2 = il(i+1) - 1
         if(indx1 .gt. indx2) go to 200
         do 190 indx=indx1,indx2
            if(indx .eq. indx1) go to 180
            indxr1 = indx1
            indxr2 = indx - 1
            indxc1 = ju(jl(indx))
            indxc2 = ju(jl(indx)+1) - 1
            if(indxc1 .gt. indxc2) go to 180
 160        kr = jl(indxr1)
 170        kc = iu(indxc1)
            if(kr .gt. kc) then
               indxc1 = indxc1 + 1
               if(indxc1 .le. indxc2) go to 170
            elseif(kr .lt. kc) then
               indxr1 = indxr1 + 1
               if(indxr1 .le. indxr2) go to 160
            elseif(kr .eq. kc) then
               l(indx) = l(indx) - l(indxr1)*dinv(kc)*u(indxc1)
               indxr1 = indxr1 + 1
               indxc1 = indxc1 + 1
               if(indxr1 .le. indxr2 .and. indxc1 .le. indxc2) go to 160
            endif
 180        l(indx) = l(indx)/dinv(jl(indx))
 190     continue
c         
c         ith column of u
 200     indx1 = ju(i)
         indx2 = ju(i+1) - 1
         if(indx1 .gt. indx2) go to 260
         do 250 indx=indx1,indx2
            if(indx .eq. indx1) go to 240
            indxc1 = indx1
            indxc2 = indx - 1
            indxr1 = il(iu(indx))
            indxr2 = il(iu(indx)+1) - 1
            if(indxr1 .gt. indxr2) go to 240
 210        kr = jl(indxr1)
 220        kc = iu(indxc1)
            if(kr .gt. kc) then
               indxc1 = indxc1 + 1
               if(indxc1 .le. indxc2) go to 220
            elseif(kr .lt. kc) then
               indxr1 = indxr1 + 1
               if(indxr1 .le. indxr2) go to 210
            elseif(kr .eq. kc) then
               u(indx) = u(indx) - l(indxr1)*dinv(kc)*u(indxc1)
               indxr1 = indxr1 + 1
               indxc1 = indxc1 + 1
               if(indxr1 .le. indxr2 .and. indxc1 .le. indxc2) go to 210
            endif
 240        u(indx) = u(indx)/dinv(iu(indx))
 250     continue
c         
c         ith diagonal element
 260     indxr1 = il(i)
         indxr2 = il(i+1) - 1
         if(indxr1 .gt. indxr2) go to 300
         indxc1 = ju(i)
         indxc2 = ju(i+1) - 1
         if(indxc1 .gt. indxc2) go to 300
 270     kr = jl(indxr1)
 280     kc = iu(indxc1)
         if(kr .gt. kc) then
            indxc1 = indxc1 + 1
            if(indxc1 .le. indxc2) go to 280
         elseif(kr .lt. kc) then
            indxr1 = indxr1 + 1
            if(indxr1 .le. indxr2) go to 270
         elseif(kr .eq. kc) then
            dinv(i) = dinv(i) - l(indxr1)*dinv(kc)*u(indxc1)
            indxr1 = indxr1 + 1
            indxc1 = indxc1 + 1
            if(indxr1 .le. indxr2 .and. indxc1 .le. indxc2) go to 270
         endif
c         
 300  continue
c         
c         replace diagonal lts by their inverses.
cvd$ vector
      do 430 i=1,n
         dinv(i) = 1./dinv(i)
 430  continue
c         

      return
c------------- last line of ssilus follows ----------------------------
      end


      subroutine sgmres(n, b, x, nelt, ia, ja, a, isym, matvec,
     $	   msolve,
     $     itol, tol, itmax, iter, err, ierr, iunit, sb, sx, 
     $     rgwk, lrgw, igwk, ligw, rwork, iwork )
      implicit none


      integer  n, nelt, ia(nelt), ja(nelt), isym, itol, itmax, iter
      integer  iunit, lrgw, ligw, igwk(ligw), iwork(*)
      real     b(n), x(n), tol, err, sb(n), sx(n), a(nelt) 
      real     rgwk(lrgw), rwork(*)
      integer  jpre, kmp, maxl, nms, maxlp1, nmsl, nrsts, nrmax
      integer  i, iflag, lr, ldl, lhes, lgmr, lq, lv, lw
      integer  lz, ierr, lzm1, lxl, jscal
      real     bnrm, rhol, sum, snrm2, r1mach
      external matvec, msolve

c
c***first executable statement  sgmres
      ierr = 0
c   ---------------------------------------------------------------
c         load method parameters with user values or defaults.
c   --------------------------------------------------------------
      maxl = igwk(1)
      if (maxl .eq. 0) maxl = 10 
      if (maxl .gt. n) maxl = n
      kmp = igwk(2)
      if (kmp .eq. 0) kmp = maxl
      if (kmp .gt. maxl) kmp = maxl
      jscal = igwk(3)
      jpre = igwk(4)
c         check for consistent values of itol and jpre.
      if( itol.eq.1 .and. jpre.lt.0 ) goto 650
      if( itol.eq.2 .and. jpre.ge.0 ) goto 650
      nrmax = igwk(5)
      if( nrmax.eq.0 ) nrmax = 10
c         if nrmax .eq. -1, then set nrmax = 0 to 
c		turn off restarting.
      if( nrmax.eq.-1 ) nrmax = 0
c         if input value of tol is zero, set it to its 
c		default value.
      if( tol.eq.0.0e0 ) tol = 500.0*r1mach(3)
c
c         initialize counters.
      iter = 0
      nms = 0
      nrsts = 0
c   ---------------------------------------------------------------
c         form work array segment pointers.
c   --------------------------------------------------------------
      maxlp1 = maxl + 1
      lv = 1
      lr = lv + n*maxlp1
      lhes = lr + n + 1
      lq = lhes + maxl*maxlp1
      ldl = lq + 2*maxl
      lw = ldl + n
      lxl = lw + n
      lz = lxl + n
c
c     load igwk(6) with required minimum length of the rgwk array.
      igwk(6) = lz + n - 1
      if( lz+n-1.gt.lrgw ) goto 640
c   --------------------------------------------------------------
c         calculate scaled-preconditioned norm of rhs vector b.
c   --------------------------------------------------------------
      if (jpre .lt. 0) then
         call msolve(n, b, rgwk(lr), nelt, ia, ja, a, isym,
     $        rwork, iwork)
         nms = nms + 1
      else
         call scopy(n, b, 1, rgwk(lr), 1)
      endif
      if( jscal.eq.2 .or. jscal.eq.3 ) then
         sum = 0.e0
         do 10 i = 1,n
            sum = sum + (rgwk(lr-1+i)*sb(i))**2
 10      continue
         bnrm = sqrt(sum)
      else
         bnrm = snrm2(n,rgwk(lr),1)
      !write(*,*) 'sgmres: jpre, jscal =', jpre, jscal
      !write(*,*) 'sgmres: bnrm, norm(b) =', bnrm, snrm2(n,b,1)
      endif
c   -------------------------------------------------------------
c         calculate initial residual.
c   -------------------------------------------------------------
      call matvec(n, x, rgwk(lr), nelt, ia, ja, a, isym)
      do 50 i = 1,n
         rgwk(lr-1+i) = b(i) - rgwk(lr-1+i)
 50   continue
c   -------------------------------------------------------------
c         if performing restarting, then load the residual into the 
c         correct location in the rgwk array.
c   --------------------------------------------------------------
 100  continue
      if( nrsts.gt.nrmax ) goto 610
      if( nrsts.gt.0 ) then
c    copy the curr residual to different loc in the rgwk array.
         call scopy(n, rgwk(ldl), 1, rgwk(lr), 1)
      endif
c   ------------------------------------------------------------
c         use the spigmr algorithm to solve the linear system a*z = r.
c   -------------------------------------------------------------
      call spigmr(n, rgwk(lr), sb, sx, jscal, maxl, maxlp1, kmp,
     $       nrsts, jpre, matvec, msolve, nmsl, rgwk(lz), rgwk(lv),
     $       rgwk(lhes), rgwk(lq), lgmr, rwork, iwork, rgwk(lw),
     $       rgwk(ldl), rhol, nrmax, b, bnrm, x, rgwk(lxl), itol,
     $       tol, nelt, ia, ja, a, isym, iunit, iflag, err)
      !write(*,*) 'sgmres: rhol =', rhol
      iter = iter + lgmr
      nms = nms + nmsl
c
c         increment x by the current approximate solution z of a*z = r.
c
      lzm1 = lz - 1
      do 110 i = 1,n
         x(i) = x(i) + rgwk(lzm1+i)
 110  continue
      if( iflag.eq.0 ) goto 600
      if( iflag.eq.1 ) then
         nrsts = nrsts + 1
         goto 100
      endif
      if( iflag.eq.2 ) goto 620
c   ------------------------------------------------------------------
c         all returns are made through this section.
c   ------------------------------------------------------------------
c         the iteration has converged.
c
 600  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 0
      return
c
c         max number((nrmax+1)*maxl) of linear iterations performed.
 610  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 1
      return
c
c         gmres failed to reduce last residual in maxl iterations.
c         the iteration has stalled.
 620  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 2
      return
c         error return.  insufficient length for rgwk array.
 640  continue
      err = tol
      ierr = -1
      return
c         error return.  inconsistent itol and jpre values.
 650  continue
      err = tol
      ierr = -2

      return
c------------- last line of sgmres follows ----------------------------
      end


      subroutine scopy(n,sx,incx,sy,incy)
      implicit none

c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(*),sy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue

      return
      end


      subroutine qs2i1r( ia, ja, a, n, kflag )
      implicit none


c***begin prologue  qs2i1r
c***date written   761118   (yymmdd)
c***revision date  890125   (yymmdd)
c***category no.  n6a2a
c***keywords  library=slatec(slap),
c             type=integer(qs2i1r-i),
c             quicksort,singleton quicksort,sort,sorting
c***author  jones, r. e., (snla)
c           kahaner, d. k., (nbs)
c           seager, m. k., (llnl) seager@lll-crg.llnl.gov
c           wisniewski, j. a., (snla)
c***purpose  sort an integer array also moving an integer and real array
c            this routine sorts the integer array ia and makes the same
c            interchanges in the integer array ja and the real array a.
c            the array ia may be sorted in increasing order or decreas-
c            ing order. a slightly modified quicksort algorithm is used.
c
c***description
c     written by rondall e jones
c     modified by john a. wisniewski to use the singleton quicksort
c     algorithm. date 18 november 1976.
c
c     further modified by david k. kahaner
c     national bureau of standards
c     august, 1981
c
c     even further modification made to bring the code up to the 
c     fortran 77 level and make it more readable and to carry
c     along one integer array and one real array during the sort by
c     mark k. seager
c     lawrence livermore national laboratory
c     november, 1987
c     this routine was adapted from the isort routine.
c
c     abstract
c         this routine sorts an integer array ia and makes the same
c         interchanges in the integer array ja and the real array a.  
c         the array a may be sorted in increasing order or decreasing 
c         order.  a slightly modified quicksort algorithm is used.
c
c     description of parameters
c        ia - integer array of values to be sorted.
c        ja - integer array to be carried along.
c         a - real array to be carried along.
c         n - number of values in integer array ia to be sorted.
c     kflag - control parameter
c           = 1 means sort ia in increasing order.
c           =-1 means sort ia in decreasing order.
c
c***references
c     singleton, r. c., algorithm 347, "an efficient algorithm for 
c     sorting with minimal storage", cacm, vol. 12, no. 3, 1969, 
c     pp. 185-187.
c***routines called  xerror
c***end prologue  qs2i1r
cvd$r novector
cvd$r noconcur
      dimension il(21),iu(21)
      integer   n,ia(n),ja(n),it,iit,jt,jjt, k
      integer   kflag, ij, iu, il, l, j, nn, m, i, kk
      real      a(n), ta, tta, r
c
c***first executable statement  qs2i1r
      nn = n
      if (nn.lt.1) then
         call xerror ( 'qs2i1r- the number of values to be sorted was no
     $t positive.',59,1,1)
         return
      endif
      if( n.eq.1 ) return
      kk = iabs(kflag)
      if ( kk.ne.1 ) then
         call xerror ( 'qs2i1r- the sort control parameter, k, was not 1
     $ or -1.',55,2,1)
         return
      endif
c
c     alter array ia to get decreasing order if needed.
c
      if( kflag.lt.1 ) then
         do 20 i=1,nn
            ia(i) = -ia(i)
 20      continue
      endif
c
c     sort ia and carry ja and a along.
c     and now...just a little black magic...
      m = 1
      i = 1
      j = nn
      r = .375
 210  if( r.le.0.5898437 ) then
         r = r + 3.90625e-2
      else
         r = r-.21875
      endif
 225  k = i
c
c     select a central element of the array and save it in location 
c     it, jt, at.
c
      ij = i + ifix (float (j-i) *r)
      it = ia(ij)
      jt = ja(ij)
      ta = a(ij)
c
c     if first element of array is greater than it, interchange with it.
c
      if( ia(i).gt.it ) then
         ia(ij) = ia(i)
         ia(i)  = it
         it     = ia(ij)
         ja(ij) = ja(i)
         ja(i)  = jt
         jt     = ja(ij)
         a(ij)  = a(i)
         a(i)   = ta
         ta     = a(ij)
      endif
      l=j
c                           
c     if last element of array is less than it, swap with it.
c
      if( ia(j).lt.it ) then
         ia(ij) = ia(j)
         ia(j)  = it
         it     = ia(ij)
         ja(ij) = ja(j)
         ja(j)  = jt
         jt     = ja(ij)
         a(ij)  = a(j)
         a(j)   = ta
         ta     = a(ij)
c
c     if first element of array is greater than it, swap with it.
c
         if ( ia(i).gt.it ) then
            ia(ij) = ia(i)
            ia(i)  = it
            it     = ia(ij)
            ja(ij) = ja(i)
            ja(i)  = jt
            jt     = ja(ij)
            a(ij)  = a(i)
            a(i)   = ta
            ta     = a(ij)
         endif
      endif
c
c     find an element in the second half of the array which is 
c     smaller than it.
c
  240 l=l-1
      if( ia(l).gt.it ) go to 240
c
c     find an element in the first half of the array which is 
c     greater than it.
c
  245 k=k+1
      if( ia(k).lt.it ) go to 245
c
c     interchange these elements.
c
      if( k.le.l ) then
         iit   = ia(l)
         ia(l) = ia(k)
         ia(k) = iit
         jjt   = ja(l)
         ja(l) = ja(k)
         ja(k) = jjt
         tta   = a(l)
         a(l)  = a(k)
         a(k)  = tta
         goto 240
      endif
c
c     save upper and lower subscripts of the array yet to be sorted.
c
      if( l-i.gt.j-k ) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 260
c
c     begin again on another portion of the unsorted array.
c                                  
  255 m = m-1
      if( m.eq.0 ) go to 300
      i = il(m)
      j = iu(m)
  260 if( j-i.ge.1 ) go to 225
      if( i.eq.j ) go to 255
      if( i.eq.1 ) go to 210
      i = i-1
  265 i = i+1
      if( i.eq.j ) go to 255
      it = ia(i+1)
      jt = ja(i+1)
      ta =  a(i+1)
      if( ia(i).le.it ) go to 265
      k=i
  270 ia(k+1) = ia(k)
      ja(k+1) = ja(k)
      a(k+1)  =  a(k)
      k = k-1
      if( it.lt.ia(k) ) go to 270
      ia(k+1) = it
      ja(k+1) = jt
      a(k+1)  = ta
      go to 265
c
c     clean up, if necessary.
c
  300 if( kflag.lt.1 ) then
         do 310 i=1,nn
            ia(i) = -ia(i)
 310     continue
      endif

      return
c------------- last line of qs2i1r follows ----------------------------
      end

      subroutine xerrwv(messg,nmessg,nerr,level,ni,i1,i2,nr,r1,r2)
      implicit none


c***begin prologue  xerrwv
c***date written   800319   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  process an error message allowing 2 integer and 2 real
c            values to be included in the message.
c***description
c
c     abstract
c        xerrwv processes a diagnostic message, in a manner
c        determined by the value of level and the current value
c        of the library error control flag, kontrl.
c        (see subroutine xsetf for details.)
c        in addition, up to two integer values and two real
c        values may be printed along with the message.
c
c     description of parameters
c      --input--
c        messg - the hollerith message to be processed.
c        nmessg- the actual number of characters in messg.
c        nerr  - the error number associated with this message.
c                nerr must not be zero.
c        level - error category.
c                =2 means this is an unconditionally fatal error.
c                =1 means this is a recoverable error.  (i.e., it is
c                   non-fatal if xsetf has been appropriately called.)
c                =0 means this is a warning message only.
c                =-1 means this is a warning message which is to be
c                   printed at most once, regardless of how many
c                   times this call is executed.
c        ni    - number of integer values to be printed. (0 to 2)
c        i1    - first integer value.
c        i2    - second integer value.
c        nr    - number of real values to be printed. (0 to 2)
c        r1    - first real value.
c        r2    - second real value.
c
c     examples
c        call xerrwv('smooth -- num (=i1) was zero.',29,1,2,
c    1   1,num,0,0,0.,0.)
c        call xerrwv('quadxy -- requested error (r1) less than minimum (
c    1r2).,54,77,1,0,0,0,2,errreq,errmin)
c
c     latest revision ---  1 august 1985
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  fdump,i1mach,j4save,xerabt,xerctl,xerprt,xersav,
c                    xgetua
c***end prologue  xerrwv
      character*(*) messg
      character*20 lfirst
      character*37 form
      integer      lun(5), nunit, i1mach, isizei, lerr, llevel
      integer	   isizef, ifatal, i, kunit, iunit, lmessg
      integer	   i1, i2, nr, ni, nmessg, nerr, level
      integer      lkntrl, maxmes, j4save, junk, mkntrl
      integer	   kount, kdummy
      real	   r1, r2
c     get flags
c***first executable statement  xerrwv
      lkntrl = j4save(2,0,.false.)
      maxmes = j4save(4,0,.false.)
c     check for valid input
      if ((nmessg.gt.0).and.(nerr.ne.0).and.
     1    (level.ge.(-1)).and.(level.le.2)) go to 10
         if (lkntrl.gt.0) call xerprt('fatal error in...',17)
         call xerprt('xerror -- invalid input',23)
c        if (lkntrl.gt.0) call fdump
         if (lkntrl.gt.0) call xerprt('job abort due to fatal error.',
     1  29)
         if (lkntrl.gt.0) call xersav(' ',0,0,0,kdummy)
         call xerabt('xerror -- invalid input',23)
         return
   10 continue
c     record message
      junk = j4save(1,nerr,.true.)
      call xersav(messg,nmessg,nerr,level,kount)
c     let user override
      lfirst = messg
      lmessg = nmessg
      lerr = nerr
      llevel = level
      call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
c     reset to original values
      lmessg = nmessg
      lerr = nerr
      llevel = level
      lkntrl = max0(-2,min0(2,lkntrl))
      mkntrl = iabs(lkntrl)
c     decide whether to print message
      if ((llevel.lt.2).and.(lkntrl.eq.0)) go to 100
      if (((llevel.eq.(-1)).and.(kount.gt.min0(1,maxmes)))
     1.or.((llevel.eq.0)   .and.(kount.gt.maxmes))
     2.or.((llevel.eq.1)   .and.(kount.gt.maxmes).and.(mkntrl.eq.1))
     3.or.((llevel.eq.2)   .and.(kount.gt.max0(1,maxmes)))) go to 100
         if (lkntrl.le.0) go to 20
            call xerprt(' ',1)
c           introduction
            if (llevel.eq.(-1)) call xerprt
     1('warning message...this message will only be printed once.',57)
            if (llevel.eq.0) call xerprt('warning in...',13)
            if (llevel.eq.1) call xerprt
     1      ('recoverable error in...',23)
            if (llevel.eq.2) call xerprt('fatal error in...',17)
   20    continue
c        message
         call xerprt(messg,lmessg)
         call xgetua(lun,nunit)
         isizei = log10(float(i1mach(9))) + 1.0
         isizef = log10(float(i1mach(10))**i1mach(11)) + 1.0
         do 50 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
            do 22 i=1,min(ni,2)
               write (form,21) i,isizei
   21          format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
               if (i.eq.1) write (iunit,form) i1
               if (i.eq.2) write (iunit,form) i2
   22       continue
            do 24 i=1,min(nr,2)
               write (form,23) i,isizef+10,isizef
   23          format ('(11x,21hin above message, r',i1,'=,e',
     1         i2,'.',i2,')')
               if (i.eq.1) write (iunit,form) r1
               if (i.eq.2) write (iunit,form) r2
   24       continue
            if (lkntrl.le.0) go to 40
c              error number
               write (iunit,30) lerr
   30          format (15h error number =,i10)
   40       continue
   50    continue
c        trace-back
c        if (lkntrl.gt.0) call fdump
  100 continue
      ifatal = 0
      if ((llevel.eq.2).or.((llevel.eq.1).and.(mkntrl.eq.2)))
     1ifatal = 1
c     quit here if message is not fatal
      if (ifatal.le.0) return
      if ((lkntrl.le.0).or.(kount.gt.max0(1,maxmes))) go to 120
c        print reason for abort
         if (llevel.eq.1) call xerprt
     1   ('job abort due to unrecovered error.',35)
         if (llevel.eq.2) call xerprt
     1   ('job abort due to fatal error.',29)
c        print error summary
         call xersav(' ',-1,0,0,kdummy)
  120 continue
c     abort
      if ((llevel.eq.2).and.(kount.gt.max0(1,maxmes))) lmessg = 0
      call xerabt(messg,lmessg)

      return
      end

      integer function i1mach(i)

c
c  i/o unit numbers.
c
c    i1mach( 1) = the standard input unit.
c
c    i1mach( 2) = the standard output unit.
c
c    i1mach( 3) = the standard punch unit.
c
c    i1mach( 4) = the standard error message unit.
c
c  words.
c
c    i1mach( 5) = the number of bits per integer storage unit.
c
c    i1mach( 6) = the number of characters per integer storage unit.
c
c  integers.
c
c    assume integers are represented in the s-digit, base-a form
c
c               sign ( x(s-1)*a**(s-1) + ... + x(1)*a + x(0) )
c
c               where 0 .le. x(i) .lt. a for i=0,...,s-1.
c
c    i1mach( 7) = a, the base.
c
c    i1mach( 8) = s, the number of base-a digits.
c
c    i1mach( 9) = a**s - 1, the largest magnitude.
c
c  floating-point numbers.
c
c    assume floating-point numbers are represented in the t-digit,
c    base-b form
c
c               sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c               where 0 .le. x(i) .lt. b for i=1,...,t,
c               0 .lt. x(1), and emin .le. e .le. emax.
c
c    i1mach(10) = b, the base.
c
c  single-precision
c
c    i1mach(11) = t, the number of base-b digits.
c
c    i1mach(12) = emin, the smallest exponent e.
c
c    i1mach(13) = emax, the largest exponent e.
c
c  double-precision
c
c    i1mach(14) = t, the number of base-b digits.
c
c    i1mach(15) = emin, the smallest exponent e.
c
c    i1mach(16) = emax, the largest exponent e.
c

      integer imach(16),output,i
c
      equivalence (imach(4),output)
c

c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9 and the sel systems 85/86.
c
      data imach( 1) /   5 /
      data imach( 2) /   6 /
      data imach( 3) /   7 /
      data imach( 4) /   6 /
      data imach( 5) /  32 /
      data imach( 6) /   4 /
      data imach( 7) /   2 /
      data imach( 8) /  31 /
      data imach( 9) / 1000000 /
      data imach(10) /  16 /
      data imach(11) /   6 /
      data imach(12) / -64 /
      data imach(13) /  63 /
      data imach(14) /  14 /
      data imach(15) / -64 /
      data imach(16) /  63 /

      if (i .lt. 1  .or.  i .gt. 16) go to 999
      i1mach=imach(i)

      return
  999 write(output,1999) i
 1999 format(' i1mach - i out of bounds',i10)
      stop
      end


      real function r1mach(i)
c
c  single-precision machine constants
c
c  r1mach(1) = b**(emin-1), the smallest positive magnitude.
c
c  r1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c  r1mach(3) = b**(-t), the smallest relative spacing.
c
c  r1mach(4) = b**(1-t), the largest relative spacing.
c
c  r1mach(5) = log10(b)
c
c  to alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  removing the c from column 1.
c  on rare machines a static statement may need to be added.
c  (but probably more systems prohibit it than require it.)
c
c  where possible, octal or hexadecimal constants have been used
c  to specify the constants exactly which has in some cases
c  required the use of equivalent integer arrays.
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
      integer i	         ! put this in myself
      integer i1mach
c
      real rmach(5)
c
      equivalence (rmach(1),small(1))
      equivalence (rmach(2),large(1))
      equivalence (rmach(3),right(1))
      equivalence (rmach(4),diver(1))
      equivalence (rmach(5),log10(1))

c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9 and the sel systems 85/86.
c
      data rmach(1) / 100000 /
      data rmach(2) /  100000 /
      data rmach(3) / 100000 /
      data rmach(4) / 100000 /
      data rmach(5) /  100000 /

      if (i .lt. 1  .or.  i .gt. 5) goto 999
      r1mach = rmach(i)
      return
  999 write(i1mach(2),1999) i
 1999 format(' r1mach - i out of bounds',i10)
      stop
      end


      function j4save(iwhich,ivalue,iset)


c***begin prologue  j4save
c***refer to  xerror
c***routines called  (none)
c***description
c
c     abstract
c        j4save saves and recalls several global variables needed
c        by the library error handling routines.
c
c     description of parameters
c      --input--
c        iwhich - index of item desired.
c                = 1 refers to current error number.
c                = 2 refers to current error control flag.
c                 = 3 refers to current unit number to which error
c                    messages are to be sent.  (0 means use standard.)
c                 = 4 refers to the maximum number of times any
c                     message is to be printed (as set by xermax).
c                 = 5 refers to the total number of units to which
c                     each error message is to be written.
c                 = 6 refers to the 2nd unit for error messages
c                 = 7 refers to the 3rd unit for error messages
c                 = 8 refers to the 4th unit for error messages
c                 = 9 refers to the 5th unit for error messages
c        ivalue - the value to be set for the iwhich-th parameter,
c                 if iset is .true. .
c        iset   - if iset=.true., the iwhich-th parameter will be
c                 given the value, ivalue.  if iset=.false., the
c                 iwhich-th parameter will be unchanged, and ivalue
c                 is a dummy parameter.
c      --output--
c        the (old) value of the iwhich-th parameter will be returned
c        in the function value, j4save.
c
c     written by ron jones, with slatec common math library subcommittee
c    adapted from bell laboratories port library error handler
c     latest revision ---  1 august 1985
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***end prologue  j4save
      logical iset
      integer iparam(9), iwhich, ivalue, j4save
      save iparam
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,2,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/
c***first executable statement  j4save
      j4save = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue

      return
      end


*deck xerclr
      subroutine xerclr
      implicit none


c***begin prologue  xerclr
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  reset current error number to zero.
c***description
c
c     abstract
c        this routine simply resets the current error number to zero.
c        this may be necessary to do in order to determine that
c        a certain error has occurred again since the last time
c        numxer was referenced.
c
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xerclr
c***first executable statement  xerclr
      integer      junk, j4save

      junk = j4save(1,0,.true.)

      return
      end


      subroutine spigmr(n, r0, sr, sz, jscal, maxl, maxlp1, kmp, 
     $     nrsts, jpre, matvec, msolve, nmsl, z, v, hes, q, lgmr,
     $     rpar, ipar, wk, dl, rhol, nrmax, b, bnrm, x, xl,
     $     itol, tol, nelt, ia, ja, a, isym, iunit, iflag, err)
      implicit none


      integer  n,maxl,maxlp1,kmp,jpre,nmsl,lgmr,iflag,jscal,nrsts
      integer  nrmax,itol,nelt,isym
      real     rhol,bnrm,tol
      real     r0(*), sr(*), sz(*), z(*), v(n,*), hes(maxlp1,*)
      real     q(*), rpar(*), wk(*), dl(*)
      real     a(*), b(*), x(*), xl(*)
      integer  ipar(*), ia(*), ja(*)
      real     snrm2
c
c         local variables.
c
      integer    i, info, ip1, i2, j, k, ll, llp1
      real       r0nrm,c,dlnrm,prod,rho,s,snormw,tem, err
      integer	 iter, iunit, itmax, issgmr
      external   matvec, msolve

c
c         zero out the z array.
c***first executable statement  spigmr
      do 5 i = 1,n
         z(i) = 0.0e0
 5    continue
c
      iflag = 0
      lgmr = 0
      nmsl = 0
c         load itmax, the maximum number of iterations.
      itmax =(nrmax+1)*maxl
c   -------------------------------------------------------------------
c         the initial residual is the vector r0.
c         apply left precon. if jpre < 0 and this is not a restart.
c         apply scaling to r0 if jscal = 2 or 3.
c   -------------------------------------------------------------------
      if ((jpre .lt. 0) .and.(nrsts .eq. 0)) then
         call scopy(n, r0, 1, wk, 1)
         call msolve(n, wk, r0, nelt, ia, ja, a, isym, rpar, ipar)
         nmsl = nmsl + 1
      endif
      if (((jscal.eq.2) .or.(jscal.eq.3)) .and.(nrsts.eq.0)) then
         do 10 i = 1,n
            v(i,1) = r0(i)*sr(i)
 10      continue
      else
         do 20 i = 1,n
            v(i,1) = r0(i)
 20      continue
      endif
      r0nrm = snrm2(n, v, 1)
      iter = nrsts*maxl
c
c         call stopping routine issgmr.
c
      if (issgmr(n, b, x, xl, nelt, ia, ja, a, isym, msolve,
     $    nmsl, itol, tol, itmax, iter, err, iunit, v(1,1), z, wk,
     $    rpar, ipar, r0nrm, bnrm, sr, sz, jscal,
     $    kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
     $    hes, jpre) .ne. 0) return
      tem = 1.0e0/r0nrm
      call sscal(n, tem, v(1,1), 1)
c
c         zero out the hes array.
c
      do 50 j = 1,maxl
         do 40 i = 1,maxlp1
            hes(i,j) = 0.0e0
 40      continue
 50   continue
c   -------------------------------------------------------------------
c         main loop to compute the vectors v(*,2) to v(*,maxl).
c         the running product prod is needed for the convergence test.
c   -------------------------------------------------------------------
      prod = 1.0e0
      do 90 ll = 1,maxl
         lgmr = ll
c   -------------------------------------------------------------------
c        unscale  the  current v(ll)  and store  in wk.  call routine
c        msolve    to   compute(m-inverse)*wk,   where    m   is  the
c        preconditioner matrix.  save the answer in z.   call routine
c        matvec to compute  vnew  = a*z,  where  a is  the the system
c        matrix.  save the answer in  v(ll+1).  scale v(ll+1).   call
c        routine sorth  to  orthogonalize the    new vector vnew   =
c        v(*,ll+1).  call routine sheqr to update the factors of hes.
c   -------------------------------------------------------------------
	!write(*,*) 'jscal =', jscal
        if ((jscal .eq. 1) .or.(jscal .eq. 3)) then
           do 60 i = 1,n
              wk(i) = v(i,ll)/sz(i)
 60        continue
        else
           call scopy(n, v(1,ll), 1, wk, 1)
	   !write(*,*) 'v(:,ll)=', (v(i,ll),i=1,20)
        endif
        if (jpre .gt. 0) then
           call msolve(n, wk, z, nelt, ia, ja, a, isym, rpar, ipar)
	   !write(*,*) '1:  z(:) =', (z(i),i=1,20)
           nmsl = nmsl + 1
           call matvec(n, z, v(1,ll+1), nelt, ia, ja, a, isym)
	   !write(*,*) '1: v(:,ll+1)=', (v(i+1,ll+1),i=1,20)
        else
           call matvec(n, wk, v(1,ll+1), nelt, ia, ja, a, isym)
	   !write(*,*) '2: v(:,ll+1)=', (v(i+1,ll+1),i=1,20)
        endif
        if (jpre .lt. 0) then
           call scopy(n, v(1,ll+1), 1, wk, 1)
           call msolve(n,wk,v(1,ll+1),nelt,ia,ja,a,isym,rpar,ipar)
	   !write(*,*) '3: v(:,ll+1)=', (v(i+1,ll+1),i=1,20)
           nmsl = nmsl + 1
        endif
        if ((jscal .eq. 2) .or.(jscal .eq. 3)) then
           do 65 i = 1,n
              v(i,ll+1) = v(i,ll+1)*sr(i)
 65        continue
	   !write(*,*) '4: v(:,ll+1)=', (v(i+1,ll+1),i=1,20)
        endif
	!write(*,*) '1: maxlp1, ll, (hes(i,ll),i=1,maxlp1) :'
	!write(*,*) maxlp1, ll, (hes(i,ll),i=1,maxlp1)
        call sorth(v(1,ll+1), v, hes, n, ll, maxlp1, kmp, snormw)
        hes(ll+1,ll) = snormw
	!write(*,*) '2: maxlp1, ll, (hes(i,ll),i=1,maxlp1) :'
	!write(*,*) maxlp1, ll, (hes(i,ll),i=1,maxlp1)
        call sheqr(hes, maxlp1, ll, q, info, ll)
	!write(*,*) 'info =', info
        if (info .eq. ll) go to 120
c   -------------------------------------------------------------------
c         update rho, the estimate of the norm of the residual r0-a*zl.
c         if kmp <  maxl, then the vectors v(*,1),...,v(*,ll+1) are not
c         necessarily orthogonal for ll > kmp.  the vector dl must then
c         be computed, and its norm used in the calculation of rho.
c   -------------------------------------------------------------------
	!write(*,*) 'prod 1 =', prod
	!write(*,*) '2*maxl, ll, q(2*ll) =', 2*maxl, ll, q(2*ll)
        prod = prod*q(2*ll)
	!write(*,*) 'prod 2 =', prod
        rho = abs(prod*r0nrm)
	!write(*,*) 'rho 1 =', rho
        if ((ll.gt.kmp) .and.(kmp.lt.maxl)) then
           if (ll .eq. kmp+1) then
              call scopy(n, v(1,1), 1, dl, 1)
              do 75 i = 1,kmp
                 ip1 = i + 1
                 i2 = i*2
                 s = q(i2)
                 c = q(i2-1)
                 do 70 k = 1,n
                    dl(k) = s*dl(k) + c*v(k,ip1)
 70              continue
 75           continue
           endif
           s = q(2*ll)
           c = q(2*ll-1)/snormw
           llp1 = ll + 1
           do 80 k = 1,n
              dl(k) = s*dl(k) + c*v(k,llp1)
 80        continue
           dlnrm = snrm2(n, dl, 1)
           rho = rho*dlnrm
	!write(*,*) 'rho 2 =', rho
        endif
        rhol = rho
c   -------------------------------------------------------------------
c         test for convergence.  if passed, compute approximation zl.
c         if failed and ll < maxl, then continue iterating.
c   -------------------------------------------------------------------
        iter = nrsts*maxl + lgmr
	!write(*,*) 'spigmr calling issgmr: rhol =', rhol
        if (issgmr(n, b, x, xl, nelt, ia, ja, a, isym, msolve,
     $      nmsl, itol, tol, itmax, iter, err, iunit, dl, z, wk,
     $      rpar, ipar, rhol, bnrm, sr, sz, jscal,
     $      kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
     $      hes, jpre) .ne. 0) go to 200
        if (ll .eq. maxl) go to 100
c   -------------------------------------------------------------------
c         rescale so that the norm of v(1,ll+1) is one.
c   -------------------------------------------------------------------
        tem = 1.0e0/snormw
        call sscal(n, tem, v(1,ll+1), 1)
 90   continue
 100  continue
      if (rho .lt. r0nrm) go to 150
 120  continue
      iflag = 2
c
c         load approximate solution with zero.
c
      do 130 i = 1,n
         z(i) = 0.e0
 130  continue
      return
 150  iflag = 1
c
c         tolerance not met, but residual norm reduced.
c
      if (nrmax .gt. 0) then
c
c        if performing restarting (nrmax > 0)  calculate the residual
c        vector rl and  store it in the dl  array.  if the incomplete
c        version is being used (kmp < maxl) then dl has  already been
c        calculated up to a scaling factor.   use srlcal to calculate
c        the scaled residual vector.
c
         call srlcal(n, kmp, maxl, maxl, v, q, dl, snormw, prod,
     $        r0nrm)
      endif
c   -------------------------------------------------------------------
c         compute the approximation zl to the solution.  since the 
c         vector z was used as work space, and the initial guess
c         of the linear iteration is zero, z must be reset to zero.
c   -------------------------------------------------------------------
 200  continue
      ll = lgmr
      llp1 = ll + 1
      do 210 k = 1,llp1
         r0(k) = 0.0e0
 210  continue
      r0(1) = r0nrm
      call shels(hes, maxlp1, ll, q, r0)
      do 220 k = 1,n
         z(k) = 0.0e0
 220  continue
      do 230 i = 1,ll
         call saxpy(n, r0(i), v(1,i), 1, z, 1)
 230  continue
      if ((jscal .eq. 1) .or.(jscal .eq. 3)) then
         do 240 i = 1,n
            z(i) = z(i)/sz(i)
 240     continue
      endif
      if (jpre .gt. 0) then
         call scopy(n, z, 1, wk, 1)
         call msolve(n, wk, z, nelt, ia, ja, a, isym, rpar, ipar)
         nmsl = nmsl + 1
      endif

      return
c------------- last line of spigmr follows ----------------------------
      end

      subroutine xerror(messg,nmessg,nerr,level)
      implicit none


c***begin prologue  xerror
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  process an error (diagnostic) message.
c***description
c
c     abstract
c        xerror processes a diagnostic message, in a manner
c        determined by the value of level and the current value
c        of the library error control flag, kontrl.
c        (see subroutine xsetf for details.)
c
c     description of parameters
c      --input--
c        messg - the hollerith message to be processed, containing
c                no more than 72 characters.
c        nmessg- the actual number of characters in messg.
c        nerr  - the error number associated with this message.
c                nerr must not be zero.
c        level - error category.
c                =2 means this is an unconditionally fatal error.
c                =1 means this is a recoverable error.  (i.e., it is
c                   non-fatal if xsetf has been appropriately called.)
c                =0 means this is a warning message only.
c                =-1 means this is a warning message which is to be
c                   printed at most once, regardless of how many
c                   times this call is executed.
c
c     examples
c        call xerror('smooth -- num was zero.',23,1,2)
c        call xerror('integ  -- less than full accuracy achieved.',
c    1                43,2,1)
c        call xerror('rooter -- actual zero of f found before interval f
c    1ully collapsed.',65,3,0)
c        call xerror('exp    -- underflows being set to zero.',39,1,-1)
c
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  xerrwv
c***end prologue  xerror
      character*(*) messg
      integer       level, nerr, nmessg
c***first executable statement  xerror
      call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.,0.)

      return
      end

      subroutine xerabt(messg,nmessg)
      implicit none


c***begin prologue  xerabt
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  abort program execution and print error message.
c***description
c
c     abstract
c        ***note*** machine dependent routine
c        xerabt aborts the execution of the program.
c        the error message causing the abort is given in the calling
c        sequence, in case one needs it for printing on a dayfile,
c        for example.
c
c     description of parameters
c        messg and nmessg are as in xerror, except that nmessg may
c        be zero, in which case no message is being supplied.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  1 august 1982
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  (none)
c***end prologue  xerabt
      character*72 messg
      integer      nmessg
c      dimension messg(nmessg)
c***first executable statement  xerabt
      nmessg = 72
      stop 1
      end


      subroutine xerctl(messg1,nmessg,nerr,level,kontrl)
      implicit none
c***begin prologue  xerctl
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  allow user control over handling of errors.
c***description
c
c     abstract
c        allows user control over handling of individual errors.
c        just after each message is recorded, but before it is
c        processed any further (i.e., before it is printed or
c        a decision to abort is made), a call is made to xerctl.
c        if the user has provided his own version of xerctl, he
c        can then override the value of kontrol used in processing
c        this message by redefining its value.
c        kontrl may be set to any value from -2 to 2.
c        the meanings for kontrl are the same as in xsetf, except
c        that the value of kontrl changes only for this message.
c        if kontrl is set to a value outside the range from -2 to 2,
c        it will be moved back into that range.
c
c     description of parameters
c
c      --input--
c        messg1 - the first word (only) of the error message.
c        nmessg - same as in the call to xerror or xerrwv.
c        nerr   - same as in the call to xerror or xerrwv.
c        level  - same as in the call to xerror or xerrwv.
c        kontrl - the current value of the control flag as set
c                 by a call to xsetf.
c
c      --output--
c        kontrl - the new value of kontrl.  if kontrl is not
c                 defined, it will remain at its original value.
c                 this changed value of control affects only
c                 the current occurrence of the current message.
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  (none)
c***end prologue  xerctl
      character*20 messg1
      integer      nmessg,nerr, level, kontrl
c***first executable statement  xerctl
      return
      end

      subroutine xerprt(messg,nmessg)
      implicit none
c***begin prologue  xerprt
c***date written   790801   (yymmdd)
c***revision date  851213   (yymmdd)
c***category no.  r3
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  print error messages.
c***description
c
c     abstract
c        print the hollerith message in messg, of length nmessg,
c        on each file indicated by xgetua.
c     latest revision ---  1 august 1985
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  i1mach,xgetua
c***end prologue  xerprt
      integer lun(5), i1mach, iunit, ichar, last
      integer nunit, nmessg, kunit, lenmssg, lenmes
      character*(*) messg
c     obtain unit numbers and write line to each unit
c***first executable statement  xerprt
      call xgetua(lun,nunit)
      lenmes = len(messg)
      do 20 kunit=1,nunit
         iunit = lun(kunit)
         if (iunit.eq.0) iunit = i1mach(4)
         do 10 ichar=1,lenmes,72
            last = min0(ichar+71 , lenmes)
            write (iunit,'(1x,a)') messg(ichar:last)
   10    continue
   20 continue

      return
      end


      subroutine xersav(messg,nmessg,nerr,level,icount)
      implicit none


c***begin prologue  xersav
c***date written   800319   (yymmdd)
c***revision date  851213   (yymmdd)
c***category no.  r3
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  record that an error has occurred.
c***description
c
c     abstract
c        record that this error occurred.
c
c     description of parameters
c     --input--
c       messg, nmessg, nerr, level are as in xerror,
c       except that when nmessg=0 the tables will be
c       dumped and cleared, and when nmessg is less than zero the
c       tables will be dumped and not cleared.
c     --output--
c       icount will be the number of times this message has
c       been seen, or zero if the table has overflowed and
c       does not contain this message specifically.
c       when nmessg=0, icount will not be altered.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  1 august 1985
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  i1mach,xgetua
c***end prologue  xersav
      integer lun(5), i1mach, iunit, kunit, ii, i, level, nerr
      integer nmessg, nunit, icount, kountx
      character*(*) messg
      character*20 mestab(10),mes
      integer      nertab(10),levtab(10),kount(10)
      save mestab,nertab,levtab,kount,kountx
c     next two data statements are necessary to provide a blank
c     error table initially
      data kount(1),kount(2),kount(3),kount(4),kount(5),
     1     kount(6),kount(7),kount(8),kount(9),kount(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      data kountx/0/
c***first executable statement  xersav
      if (nmessg.gt.0) go to 80
c     dump the table
         if (kount(1).eq.0) return
c        print to each unit
         call xgetua(lun,nunit)
         do 60 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
c           print table header
            write (iunit,10)
   10       format (32h0          error message summary/
     1      51h message start             nerr     level     count)
c           print body of table
            do 20 i=1,10
               if (kount(i).eq.0) go to 30
               write (iunit,15) mestab(i),nertab(i),levtab(i),kount(i)
   15          format (1x,a20,3i10)
   20       continue
   30       continue
c           print number of other errors
            if (kountx.ne.0) write (iunit,40) kountx
   40       format (41h0other errors not individually tabulated=,i10)
            write (iunit,50)
   50       format (1x)
   60    continue
         if (nmessg.lt.0) return
c        clear the error tables
         do 70 i=1,10
   70       kount(i) = 0
         kountx = 0
         return
   80 continue
c     process a message...
c     search for this messg, or else an empty slot for this messg,
c     or else determine that the error table is full.
      mes = messg
      do 90 i=1,10
         ii = i
         if (kount(i).eq.0) go to 110
         if (mes.ne.mestab(i)) go to 90
         if (nerr.ne.nertab(i)) go to 90
         if (level.ne.levtab(i)) go to 90
         go to 100
   90 continue
c     three possible cases...
c     table is full
         kountx = kountx+1
         icount = 1
         return
c     message found in table
  100    kount(ii) = kount(ii) + 1
         icount = kount(ii)
         return
c     empty slot found for new message
  110    mestab(ii) = mes
         nertab(ii) = nerr
         levtab(ii) = level
         kount(ii)  = 1
         icount = 1

         return
      end


      subroutine xgetua(iunita,n)
      implicit none


c***begin prologue  xgetua
c***date written   790801   (yymmdd)
c***revision date  851111   (yymmdd)
c***category no.  r3c
c***keywords  error,xerror package
c***author  jones, r. e., (snla)
c***purpose  return unit number(s) to which error messages are being
c            sent.
c***description
c
c     abstract
c        xgetua may be called to determine the unit number or numbers
c        to which error messages are being sent.
c        these unit numbers may have been set by a call to xsetun,
c        or a call to xsetua, or may be a default value.
c
c     description of parameters
c      --output--
c        iunit - an array of one to five unit numbers, depending
c                on the value of n.  a value of zero refers to the
c                default unit, as defined by the i1mach machine
c                constant routine.  only iunit(1),...,iunit(n) are
c                defined by xgetua.  the values of iunit(n+1),...,
c                iunit(5) are not defined (for n .lt. 5) or altered
c                in any way by xgetua.
c        n     - the number of units to which copies of the
c                error messages are being sent.  n will be in the
c                range from 1 to 5.
c
c     latest revision ---  19 mar 1980
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xgetua
      integer     iunita(5), n, j4save, index, i 
c***first executable statement  xgetua
      n = j4save(5,0,.false.)
      do 30 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         iunita(i) = j4save(index,0,.false.)
   30 continue

      return
      end


      subroutine sscal(n,sa,sx,incx)
      implicit none

c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increments, 9/29/88.
c
      real sa,sx(*)
      integer i,ix,incx,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      do 10 i = 1,n 
        sx(ix) = sa*sx(ix)
        ix = ix + incx 
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue

      return
      end


      subroutine sorth(vnew, v, hes, n, ll, ldhes, kmp, snormw)
      implicit none


      integer n, ll, ldhes, kmp
      real vnew, v, hes, snormw
      dimension vnew(*), v(n,*), hes(ldhes,*)
      real	snrm2, sdot
c
c         internal variables.
c
      integer i, i0
      real arg, sumdsq, tem, vnrm
c
c         get norm of unaltered vnew for later use.
c***first executable statement  sorth
      vnrm = snrm2(n, vnew, 1)
      !write(*,*) 'sorth:  vnrm=', vnrm
      !write(*,*) 'vnew(1:20) =', (vnew(i),i=1,20)
c   -------------------------------------------------------------------
c         perform the modified gram-schmidt procedure on vnew =a*v(ll).
c         scaled inner products give new column of hes.
c         projections of earlier vectors are subtracted from vnew.
c   -------------------------------------------------------------------
      i0 = max0(1,ll-kmp+1)
      !write(*,*) 'sorth:  i0, ll =',  i0, ll
      do 10 i = i0,ll
         hes(i,ll) = sdot(n, v(1,i), 1, vnew, 1)
	 !write(*,*) 'i =', i
	 !write(*,*) 'v(:,i) =', (v(k,i),k=1,n)
	 !write(*,*) 'sorth:  hes(i,ll) =', hes(i,ll)
	 !write(*,*) 'sorth: stop!'
         tem = -hes(i,ll)
         call saxpy(n, tem, v(1,i), 1, vnew, 1)
 10   continue
c   -------------------------------------------------------------------
c         compute snormw = norm of vnew.  if vnew is small compared 
c         to its input value (in norm), then reorthogonalize vnew to 
c         v(*,1) through v(*,ll).  correct if relative correction 
c         exceeds 1000*(unit roundoff).  finally, correct snormw using 
c         the dot products involved.
c   -------------------------------------------------------------------
      snormw = snrm2(n, vnew, 1)
      if (vnrm + 0.001e0*snormw .ne. vnrm) return
      sumdsq = 0.0e0
      do 30 i = i0,ll
         tem = -sdot(n, v(1,i), 1, vnew, 1)
         if (hes(i,ll) + 0.001e0*tem .eq. hes(i,ll)) go to 30
         hes(i,ll) = hes(i,ll) - tem
         call saxpy(n, tem, v(1,i), 1, vnew, 1)
         sumdsq = sumdsq + tem**2
 30   continue
      if (sumdsq .eq. 0.0e0) return
      arg = amax1(0.0e0,snormw**2 - sumdsq)
      snormw = sqrt(arg)
c

      return
c------------- last line of sorth follows ----------------------------
      end


      subroutine sheqr(a, lda, n, q, info, ijob)
      implicit none


      integer lda, n, info, ijob
      real a(lda,*), q(*)
c
c         local variables.
c
      integer i, iq, j, k, km1, kp1, nm1
      real c, s, t, t1, t2
c
c***first executable statement  sheqr
      if (ijob .gt. 1) go to 70
c   -------------------------------------------------------------------
c         a new facorization is desired.
c   -------------------------------------------------------------------
c         qr decomposition without pivoting.
c
      info = 0
      !write(*,*) 'sheqr: n =', n
      !write(*,*) '1: (q(i),i=1,20) =',  (q(i),i=1,20)
      do 60 k = 1, n
         km1 = k - 1
         kp1 = k + 1
c
c           compute k-th column of r.
c           first, multiply the k-th column of a by the previous
c           k-1 givens rotations.
c
         if (km1 .lt. 1) go to 20
         do 10 j = 1, km1
            i = 2*(j-1) + 1
            t1 = a(j,k)
            t2 = a(j+1,k)
            c = q(i)
            s = q(i+1)
            a(j,k) = c*t1 - s*t2
            a(j+1,k) = s*t1 + c*t2
 10      continue
c
c         compute givens components c and s.
c
 20      continue
         iq = 2*km1 + 1
         t1 = a(k,k)
         t2 = a(kp1,k)
	 !write(*,*) 'k =', k
	 !write(*,*) 't1, t2 =', t1, t2
         if( t2.eq.0.0e0 ) then
            c = 1.0e0
            s = 0.0e0
         elseif( abs(t2).ge.abs(t1) ) then
            t = t1/t2
            s = -1.0e0/sqrt(1.0e0+t*t)
            c = -s*t
         else
            t = t2/t1
            c = 1.0e0/sqrt(1.0e0+t*t)
            s = -c*t
         endif
         q(iq) = c
         q(iq+1) = s
	 !write(*,*) '1: c,s =', c, s
	 !write(*,*) 'iq+1, q(iq+1) =',  iq+1, q(iq+1)
         a(k,k) = c*t1 - s*t2
         if( a(k,k).eq.0.0e0 ) info = k
 60   continue
      !write(*,*) '2: (q(i),i=1,20) =',  (q(i),i=1,20)
      return
c   -------------------------------------------------------------------
c         the old factorization of a will be updated.  a row and a 
c         column has been added to the matrix a.  n by n-1 is now 
c         the old size of the matrix.
c   -------------------------------------------------------------------
 70   continue
      !write(*,*) 'sheqr 2: n =', n
      nm1 = n - 1
c   -------------------------------------------------------------------
c         multiply the new column by the n previous givens rotations.
c   -------------------------------------------------------------------
      do 100 k = 1,nm1
         i = 2*(k-1) + 1
         t1 = a(k,n)
         t2 = a(k+1,n)
         c = q(i)
         s = q(i+1)
	 !write(*,*) '2: c,s =', c, s
	 !write(*,*) 'iq+1, q(iq+1) =',  iq+1, q(iq+1)
         a(k,n) = c*t1 - s*t2
         a(k+1,n) = s*t1 + c*t2
 100  continue
c   -------------------------------------------------------------------
c         complete update of decomposition by forming last givens 
c         rotation, and multiplying it times the column 
c         vector(a(n,n),a(np1,n)).
c   -------------------------------------------------------------------
      info = 0
      t1 = a(n,n)
      t2 = a(n+1,n)
      if ( t2.eq.0.0e0 ) then
         c = 1.0e0
         s = 0.0e0
      elseif( abs(t2).ge.abs(t1) ) then
         t = t1/t2
         s = -1.0e0/sqrt(1.0e0+t*t)
         c = -s*t
      else
         t = t2/t1
         c = 1.0e0/sqrt(1.0e0+t*t)
         s = -c*t
      endif
      iq = 2*n - 1
      q(iq) = c
      q(iq+1) = s
      a(n,n) = c*t1 - s*t2
	 !write(*,*) '3: c,s =', c, s
	 !write(*,*) 'iq+1, q(iq+1) =',  iq+1, q(iq+1)
      !write(*,*) '3: (q(i),i=1,20) =',  (q(i),i=1,20)
      if (a(n,n) .eq. 0.0e0) info = n

      return
c------------- last line of sheqr follows ----------------------------
      end


       subroutine srlcal(n, kmp, ll, maxl, v, q, rl, snormw, prod,
     $     r0nrm)
      implicit none


      integer n, kmp, ll, maxl
      real v, q, rl, snormw
      dimension v(n,*), q(*), rl(n)
c
c         internal variables.
c
      integer i, ip1, i2, k, llp1, llm1
      real    c, s, tem, prod, r0nrm
c
c***first executable statement  srlcal
      if (kmp .eq. maxl) then
c
c         calculate rl.  start by copying v(*,1) into rl.
c
         call scopy(n, v(1,1), 1, rl, 1)
         llm1 = ll - 1
         do 20 i = 1,llm1
            ip1 = i + 1
            i2 = i*2
            s = q(i2)
            c = q(i2-1)
            do 10 k = 1,n
               rl(k) = s*rl(k) + c*v(k,ip1)
 10         continue
 20      continue
         s = q(2*ll)
         c = q(2*ll-1)/snormw
         llp1 = ll + 1
         do 30 k = 1,n
            rl(k) = s*rl(k) + c*v(k,llp1)
 30      continue
      endif
c
c         when kmp < maxl, rl vector already partially calculated. 
c         scale rl by r0nrm*prod to obtain the residual rl.
c
      tem = r0nrm*prod
      call sscal(n, tem, rl, 1)

      return
c------------- last line of srlcal follows ----------------------------
      end


      subroutine shels(a, lda, n, q, b)
      implicit none


      integer lda, n
      real a(lda,*), b(*), q(*)
c
c         local variables.
c
      integer iq, k, kb, kp1
      real c, s, t, t1, t2
c
c         minimize(b-a*x,b-a*x).  first form q*b.
c
      do 20 k = 1, n
         kp1 = k + 1
         iq = 2*(k-1) + 1
         c = q(iq)
         s = q(iq+1)
         t1 = b(k)
         t2 = b(kp1)
         b(k) = c*t1 - s*t2
         b(kp1) = s*t1 + c*t2
 20   continue
c
c         now solve  r*x = q*b.
c
      do 40 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call saxpy(k-1, t, a(1,k), 1, b(1), 1)
 40   continue

      return
c------------- last line of shels follows ----------------------------
      end


      subroutine saxpy(n,sa,sx,incx,sy,incy)
      implicit none

c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(*),sy(*),sa
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue

      return
      end

      subroutine ssmv( n, x, y, nelt, ia, ja, a, isym )
      implicit none



c***begin prologue  ssmv
c***date written   871119   (yymmdd)
c***revision date  881213   (yymmdd)
c***category no.  d2a4, d2b4
c***keywords  library=slatec(slap),
c             type=single precision(ssmv-s),
c             matrix vector multiply, sparse
c***author  greenbaum, anne, courant institute
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-300
c             livermore, ca 94550 (415) 423-3141
c             seager@lll-crg.llnl.gov
c***purpose  slap column format sparse matrix vector product.
c            routine to calculate the sparse matrix vector product:
c            y = a*x.
c***description
c *usage:
c     integer  n, nelt, ia(nelt), ja(n+1), isym
c     real     x(n), y(n), a(nelt)
c
c     call ssmv(n, x, y, nelt, ia, ja, a, isym )
c         
c *arguments:
c n      :in       integer.
c         order of the matrix.
c x      :in       real x(n).
c         the vector that should be multiplied by the matrix.
c y      :out      real y(n).
c         the product of the matrix and the vector.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(n+1).
c a      :in       integer a(nelt).
c         these arrays should hold the matrix a in the slap column
c         format.  see "description", below. 
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all nonzero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
c
c *precision:           single precision
c *cautions:
c     this   routine   assumes  that  the matrix a is stored in slap 
c     column format.  it does not check  for  this (for  speed)  and 
c     evil, ugly, ornery and nasty things  will happen if the matrix 
c     data  structure  is,  in fact, not slap column.  beware of the 
c     wrong data structure!!!
c
c *see also:
c       ssmtv
c***references  (none)
c***routines called  (none)
c***end prologue  ssmv
      integer n, nelt, ia(nelt), ja(nelt), isym
      integer i, ibgn, iend, j, jbgn, jend, irow, icol
      real    a(nelt), x(n), y(n)
c
c         zero out the result vector.
c***first executable statement  ssmv
      do 10 i = 1, n
         y(i) = 0.0
 10   continue
c
c         multiply by a.
c
cvd$r noconcur
      do 30 icol = 1, n
         ibgn = ja(icol)
         iend = ja(icol+1)-1
clll. option assert (nohazard)
cdir$ ivdep
cvd$ nodepchk
         do 20 i = ibgn, iend
            y(ia(i)) = y(ia(i)) + a(i)*x(icol)
 20      continue
 30   continue
c
      if( isym.eq.1 ) then
c
c         the matrix is non-symmetric.  need to get the other half in...
c         this loops assumes that the diagonal is the first entry in
c         each column.
c
         do 50 irow = 1, n
            jbgn = ja(irow)+1
            jend = ja(irow+1)-1
            if( jbgn.gt.jend ) goto 50
            do 40 j = jbgn, jend
               y(irow) = y(irow) + a(j)*x(ia(j))
 40         continue
 50      continue
      endif

      return
c------------- last line of ssmv follows ----------------------------
      end


      real function sdot(n,sx,incx,sy,incy)

c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp

      return
      end


      real function snrm2 ( n, sx, incx)


      integer          next, i, j, n, nn, incx
      real   sx(*),  cutlo, cuthi, hitest, sum, xmax, zero, one
      data   zero, one /0.0e0, 1.0e0/
c
c     euclidean norm of the n-vector stored in sx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 4.441e-16,  1.304e19 /
c
      if(n .gt. 0) go to 10
         snrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( abs(sx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( sx(i) .eq. zero) go to 200
      if( abs(sx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / sx(i)) / sx(i)
  105 xmax = abs(sx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( abs(sx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( abs(sx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / sx(i))**2
         xmax = abs(sx(i))
         go to 200
c
  115 sum = sum + (sx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(abs(sx(j)) .ge. hitest) go to 100
   95    sum = sum + sx(j)**2
      snrm2 = sqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      snrm2 = xmax * sqrt(sum)
  300 continue

      return
      end


       function issgmr(n, b, x, xl, nelt, ia, ja, a, isym, msolve,
     $     nmsl, itol, tol, itmax, iter, err, iunit, r, z, dz,
     $     rwork, iwork, rnrm, bnrm, sb, sx, jscal,
     $     kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
     $     hes, jpre)



      integer   kmp, lgmr, maxl, maxlp1, jpre, nmsl, n, nelt
      real      dxnrm, rnrm, r0nrm, snormw, solnrm, prod
      real      snrm2, r1mach, sb(*), sx(*), q(*), v(n,*)
      real      b(*), x(*), r(*), z(*), dz(*), rwork(*)
      real      hes(maxlp1,*), xl(*), a(nelt)
      integer   ia(nelt), ja(nelt), iwork(*)

c
c         local variables.
c
      integer    i, itol, iter, itmax, issgmr, iunit, ielmax
      integer    jscal, isym
      real       tol, err, bnrm, fuzz, ratmax, rat
      real       c,dlnrm,rho,s,tem

c      common /solblk/ soln(1)
c      save solnrm

      external   msolve

c         
c***first executable statement issgmr
      issgmr = 0
      if ( itol.eq.0 ) then
c
c       use input from spigmr to determine if stop conditions are met.
c
         err = rnrm/bnrm
      endif
      if ( (itol.gt.0) .and. (itol.le.3) ) then
c
c       use srlcal to calculate the scaled residual vector. 
c       store answer in r.
c
         if ( lgmr.ne.0 ) call srlcal(n, kmp, lgmr, maxl, v, q, r,
     $                                snormw, prod, r0nrm)
         if ( itol.le.2 ) then
c         err = ||residual||/||righthandside||(2-norms).
            err = snrm2(n, r, 1)/bnrm
c
c         unscale r by r0nrm*prod when kmp < maxl.
c
            if ( (kmp.lt.maxl) .and. (lgmr.ne.0) ) then
               tem = 1.0e0/(r0nrm*prod)
               call sscal(n, tem, r, 1)
            endif
         elseif ( itol.eq.3 ) then
c         err = max |(minv*residual)(i)/x(i)|
c         when jpre .lt. 0, r already contains minv*residual.
            if ( jpre.gt.0 ) then
               call msolve(n, r, dz, nelt, ia, ja, a, isym, rwork,
     $              iwork)
               nmsl = nmsl + 1
            endif
c
c         unscale r by r0nrm*prod when kmp < maxl.
c
            if ( (kmp.lt.maxl) .and. (lgmr.ne.0) ) then
               tem = 1.0e0/(r0nrm*prod)
               call sscal(n, tem, r, 1)
            endif
c
            fuzz = r1mach(1)
            ielmax = 1
            ratmax = abs(dz(1))/amax1(abs(x(1)),fuzz)
            do 25 i = 2, n
               rat = abs(dz(i))/amax1(abs(x(i)),fuzz)
               if( rat.gt.ratmax ) then
                  ielmax = i
                  ratmax = rat
               endif
 25         continue
            err = ratmax
            if( ratmax.le.tol ) issgmr = 1
            if( iunit.gt.0 ) write(iunit,1020) iter, ielmax, ratmax
            return
         endif
      endif
c      if ( itol.eq.11 ) then
c
c       use sxlcal to calculate the approximate solution xl.
c
c         if ( (lgmr.ne.0) .and. (iter.gt.0) ) then
c            call sxlcal(n, lgmr, x, xl, xl, hes, maxlp1, q, v,
c     $		 r0nrm,
c     $           dz, sx, jscal, jpre, msolve, nmsl, rwork, iwork,
c     $           nelt, ia, ja, a, isym)
c         elseif ( iter.eq.0 ) then
c         copy x to xl to check if initial guess is good enough.
c            call scopy(n, x, 1, xl, 1)
c         else
c         return since this is the first call to spigmr on a restart.
c            return
c         endif
c
c         if ((jscal .eq. 0) .or.(jscal .eq. 2)) then
c         err = ||x-truesolution||/||truesolution||(2-norms).
c            if ( iter.eq.0 ) solnrm = snrm2(n, soln, 1)
c            do 30 i = 1, n
c               dz(i) = xl(i) - soln(i)
c 30         continue
c            err = snrm2(n, dz, 1)/solnrm
c         else
c            if (iter .eq. 0) then
c               solnrm = 0.e0
c               do 40 i = 1,n
c                  solnrm = solnrm + (sx(i)*soln(i))**2
c 40            continue
c               solnrm = sqrt(solnrm)
c            endif
c            dxnrm = 0.e0
c            do 50 i = 1,n
c               dxnrm = dxnrm + (sx(i)*(xl(i)-soln(i)))**2
c 50         continue
c            dxnrm = sqrt(dxnrm)
c            err = dxnrm/solnrm
c         endif
c      endif
c         
      if( iunit.ne.0 ) then
         if( iter.eq.0 ) then
            write(iunit,1000) maxl, kmp, n, itol
         endif
c	 write(6,*) 'itol, iter =', itol, iter
c	 write(6,*) 'rnrm, bnrm =', rnrm, bnrm
         write(iunit,1010) iter, rnrm/bnrm, err
      endif
      if ( err.le.tol ) issgmr = 1
c         
 1000 format(' generalized minimum residual(',i3,i3,') for ',
     $     'n, itol = ',i6, i6,
     $     /' iter','   natral err est','   error estimate')
 1010 format(1x,i4,1x,e16.7,1x,e16.7)
 1020 format(1x,' iter = ',i5, ' ielmax = ',i5,
     $     ' |r(ielmax)/x(ielmax)| = ',e12.5)

c------------- last line of issgmr follows ----------------------

c
      return
      end
