c ---------------------------------------------------------------

      subroutine integ(x,nmi,nma,dx,rint)
      dimension x(*)

      f1=3.*x(nmi)/8.+7.*x(nmi+1)/6.+23.*x(nmi+2)/24.
      f2=3.*x(nma)/8.+7.*x(nma-1)/6.+23.*x(nmi-2)/24.
      f3=0.
      do i=nmi+3,nma-3
       f3=f3+x(i)
      enddo

      rint=dx*(f1+f2+f3)

      end
