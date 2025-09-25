      function sdot(n,a,inca,b,incb)
      dimension a(inca,*),b(incb,*)
      sum=0.
      do i=1,n
        sum=sum+a(1,i)*b(1,i)
      enddo
      sdot=sum
      return
      end

