      subroutine ylm(xlat,xlon,lmax,y,wk1,wk2,wk3)
c
      complex temp,fac,dfac
      dimension wk1(1),wk2(1),wk3(1),y(1)
c
c     wk1,wk2,wk3 should be dimensioned at least (lmax+1)*4
c
      data radian/57.2957795/    ! 360./2pi
c
c     transform to spherical coordinates
      theta=(90.-xlat)/radian
      phi=xlon/radian
c
c    loop over l values
      ind=0
      lm1=lmax+1
      do 10 il1=1,lm1
      l=il1-1
      call legndr(theta,l,l,wk1,wk2,wk3)
c
      fac=(1.,0.)
      dfac=cexp(cmplx(0.,phi))
c
c    loop over m values
      do 20 im=1,il1
      temp=fac*cmplx(wk1(im),0.)
      ind=ind+1
      y(ind)=real(temp)
      if(im.eq.1) goto 20
      ind=ind+1
      y(ind)=aimag(temp)
   20 fac=fac*dfac   ! calculates exp(im phi)
c
   10 continue
      return
      end
