c ---------------------------------------------------------------------

      subroutine wsphhead(lmax,idp1,idp2,strt,lstrt)
      parameter(MAXL=100)
      parameter(MAXP=24)
      dimension lask(0:MAXL),mask(MAXP)

      character*120 str1,str2,strt

      if(lmax.gt.MAXL) stop 'wsphhead >>> lmax too large'
      if(idp2.gt.MAXP) stop 'wsphhead >>> idp2 too large'
      do i=0,MAXL
       if(i.le.lmax) then
        lask(i)=1
       else
        lask(i)=0
       endif
      enddo

      do i=1,24
       mask(i)=0
      enddo
      do i=idp1,idp2
       mask(i)=1
      enddo

      write(str1,10) lmax,(lask(i),i=0,lmax)
10    format(10x,i4,1x,100i1)
      write(str2,20) 24,(mask(i),i=1,24)
20    format(i4,1x,24i1)

      len1=istlen(str1)
      len2=istlen(str2)
      if((len1+len2).gt.120) stop 'wsphhead: increase dimension of strt'
      strt=str1(1:len1)//str2(1:len2)
      lstrt=len1+len2+1
      end
