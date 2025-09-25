c ---------------------------------------------------------------

      subroutine modtrf(lmx,ndim,xin,xout,evcs)

c     this subroutine transforms a model vector in sph style
c     defined with respect to the eigenvector basis
c     back to a the normal sph representation.
c     run inittrf before using this routine!

      dimension evcs(*)
      double precision alph,xin(*),xout(*)

      numatd=(lmx+1)**2

      do i=1,ndim
       idin=(i-1)*numatd
       idev=(i-1)*ndim
c      write(6,*) i,idin,idev
       do j=1,ndim
        idout=(j-1)*numatd
        alph=dble(evcs(idev+j))
c       write(6,*) '     ',idout,alph
c       write(6,*) (xin(idin+k),k=1,numatd)
        call daxpy(numatd,alph,xin(idin+1),1,xout(idout+1),1)
       enddo
      enddo

      end
