      subroutine ylm02xlm (lmax,d0,d1)

      dimension d0(*),d1(*)

c    construct vector d1 with values of Xlm:
      ind=1
      d1(ind)=d0(ind)
      do l=1,lmax
       do m=0,l
        if(m.eq.0) then
         ind=ind+1
         d1(ind)=d0(ind)
        else
         ind=ind+1
         d1(ind)=d0(ind)
         ind=ind+1
         d1(ind)=d0(ind-1)
        endif
       enddo
      enddo

      end
