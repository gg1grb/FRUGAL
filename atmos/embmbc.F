      subroutine embmbc (data)
      implicit none
c
c=======================================================================
c     set boundary values
c     based on code by a.fanning
c
c     input:
c       data = array to be set
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c

#include "param.h"
c
      real data(imta,jmta)
c
c     set cyclic east west
c
      do j=1,jmta
        data(1,j) = data(imtam1,j)
        data(imta,j) = data(2,j)
      enddo
c
c     set closed North South
c
      do i=1,imta
        data(i,1) = data(i,2)
        data(i,jmta) = data(i,jmtam1)
      enddo

c
      return
      end

      subroutine setbca (a, imt, jmtorkm)
c
c=======================================================================
c     set zonal boundary condition on the first index of array "a" for
c     every second index. the first index corresponds to the "x"
c     or longitude direction.
c
c     input:
c      a = array in need of setting the zonal b.c.
c     output
c      a = array with zonal b.c. set
c=======================================================================
c
      integer jmtorkm
      dimension a(imt,jmtorkm)
      do k=1,jmtorkm
#ifdef cyclic
        a(1,k)   = a(imt-1,k)
        a(imt,k) = a(2,k)
#else
# ifdef obc_west
        a(1,k)   = a(2,k)
# else
        a(1,k)   = 0.0
# endif
# ifdef obc_east
        a(imt,k) = a(imt-1,k)
# else
        a(imt,k) = 0.0
# endif
#endif
      enddo
      return
      end


