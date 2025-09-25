      subroutine saxpy(n,s,a,inca,d,incd)
      dimension a(inca,*),d(incd,*)
      do i=1,n
        d(1,i)=d(1,i)+s*a(1,i)
      enddo
      return
      end
