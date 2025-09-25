c ---------------------------------------------------------------------

      subroutine rsphhead(line1,lmax,idp1,idp2,ndp)
      parameter(MAXL=100)
      parameter(MAXP=24)
      dimension lask(0:MAXL),mask(MAXP)

c     for strgrep
      parameter(MAXSUB=50)
      character*50 strout(MAXSUB)
      dimension il(MAXSUB),icol(MAXSUB)

      character*(*) line1

      call strgrep(line1,strout,il,nsub)
      read(strout(1),*) lmaxrd
      read(strout(2),10)(lask(i),i=0,lmaxrd)
10    format(100i1)
      read(strout(3),*) ndep
      read(strout(4),20)(mask(i),i=1,ndep)
20    format(24i1)

      do i=0,lmaxrd
       if(lask(i).eq.1) lmax=i
      enddo

      idp1=1
      do i=1,ndep
       if(i.ne.1) then
        if(mask(i).eq.1.and.mask(i-1).eq.0) idp1=i
       endif
       if(mask(i).eq.1.or.mask(i).eq.2) idp2=i
      enddo

      ndp=idp2-idp1+1

      end
