c------------------------------------------------------------

      Subroutine convpack(pack,n,unp)

c     This subroutine assumes that the physical dimension of pack
c     is n(n+1)/2 and the physical dimension of unp is n**2

      dimension pack(*),unp(*)

      indexp=1
      Do i=1,n
       lenn=n-i+1
       lent=n-i
       indexu=(i-1)*n+i
       indexut=(i*n)+i
       call scopy(lenn,pack(indexp),1,unp(indexu),1)
       call scopy(lent,pack(indexp+1),1,unp(indexut),n)
       indexp=indexp+lenn
      Enddo
 
      End
