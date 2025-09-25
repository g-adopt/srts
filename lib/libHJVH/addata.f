c ------------------------------------------------------------

      subroutine addata(atadbl,ndim,ttfilenm,phwt,ttwt)

c     this subroutine adds a .ata matrix for travel time data to
c     a .ata matrix for phase velocity data

      parameter(MAXNRD=21)
      parameter(MAXLINV=20)
      parameter(MXLENY=(MAXLINV+1)**2)
      parameter(MXATD=MAXNRD*MXLENY)

      double precision atadbl(*)
      dimension dum(MXATD)
      character*80 ttfilenm

      if(ndim.gt.MXATD) stop 'addata: ndim.gt.MXATD'

      open(11,file=ttfilenm,status='old',form='unformatted')
c     read(11) nata,dtd
      read(11) lmax,nsplmn,nsplmx,nata

      if(nata.ne.ndim) stop 'addata: incompatible travel time ata file'

      ind=0
      do j=1,nata
        read(11) (dum(i),i=j,nata)
        do i=j,nata
         ind=ind+1
         atadbl(ind)=dble(phwt)*atadbl(ind)+dble(ttwt*dum(i))
        enddo
      enddo

      close(11)

      end

