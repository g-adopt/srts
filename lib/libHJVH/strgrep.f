c ------------------------------------------------------------
      subroutine strgrep(strin,strout,ilsub,nsub)

      parameter(MAXSUB=50)
      character*50 strout(MAXSUB)
      character*(*) strin
      dimension ilsub(MAXSUB)

      do i=1,MAXSUB
       strout(i)=''
      enddo

      l=istlen(strin)

      iwrd=0
      nw=0
      do i=1,l
       if(strin(i:i).eq.char(32).or.strin(i:i).eq.char(0).or.
     1    strin(i:i).eq.char(6)) then
        if(iwrd.eq.1) then
         iend=i-1
         ils=iend-ist+1
         ilsub(nw)=ils
         strout(nw)(1:ils)=strin(ist:iend)
         iwrd=0
        endif
       else
        if(iwrd.ne.1) then
         nw=nw+1
         iend=0
         ist=i
         iwrd=1
        endif
       endif
      enddo

      if(iwrd.eq.1) then
         iend=l
         ils=iend-ist+1
         ilsub(nw)=ils
         strout(nw)=strin(ist:iend)
      endif
      nsub=nw

      end
