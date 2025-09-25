      subroutine wspthead(nfil,lmax,mp1,iparsw,ipardps)
     
      parameter(MXPARTYP=7)
      dimension iparsw(MXPARTYP),ipardps(2,MXPARTYP)
      
      if(mp1.gt.10) stop 'wspthead >>> number of parameters too large'
      write(nfil,10) lmax,mp1,(iparsw(i),i=1,mp1)
10    format(i2,1x,i2,1x,10i1)
      do i=1,mp1
       if(iparsw(i).eq.1) then
        write(nfil,'(2i4)') ipardps(1,i),ipardps(2,i)
       endif
      enddo
       
      end
