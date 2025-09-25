      subroutine normylm(lmx,wtsnorm)

c     ..... This routines normalises the ylm to 1/2.......

c     we devide by sqrt(2) because ylm is normalised
c     to 1 for m=0 and 1/2 for m.ne.0. 
c     Therefore, surface integrals of Y_l0, as defined
c     in subroutines ylm, times some function will be sqrt(2) larger
c     than when the spherical harmonics were fully normalised.
c     In a typical surface inversion with model norm damping, derivatives
c     with respect to the m=0 ylms will be larger and therefore, for a
c     unit model coefficient, the model size will actually be twice
c     as large for Yl0 as for other Ylm (n.ne.0) 

      dimension wtsnorm(*)

      k=0
      do l=0,lmx
        do m=0,l
          k=k+1
          wtsnorm(k)=sqrt(2.)
          if(m.ne.0) then
            k=k+1
            wtsnorm(k)=sqrt(2.) ! correct for normalisation in ylm
          else
            wtsnorm(k)=1.   
          endif
        enddo
      enddo

      end
