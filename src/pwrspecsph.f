      program pwrspecsph

cp    March 2014
cp    power spectrum according to Bernard's equation
cp    Schuberth et al 2009 eq 1 and 2 - or Dahlen and Tromp B.95 B.104
cp    computed for an .sph model file slice

      parameter(MXL=40)
      parameter(MXLENY=(MXL+1)**2)

      character*256 rawfile
      integer lmax,ncoef,l,m,i
      double precision dum(0:MXLENY)
      double precision x(0:MXL,0:MXL*2)
      double precision rawsum,powerdeg

      pi=acos(-1.)
      write(6,*) 'Enter .raw file'
      read(5,'(a)') rawfile

cp    read in model rawfile
      open(45,file=rawfile)
      read(45,*) lmax
      write(6,*) 'rawfile=',rawfile
      write(6,*) 'model lmax =',lmax
      if(lmax.gt.MXL) stop 'stop >> lmax.gt.MXL'
      ncoef=(lmax+1)**2
      read(45,*) (dum(i),i=1,ncoef)
      close(45)

cp    initialize arrays
      do l=0,lmax
        do m=0,2*l
          x(l,m)=0.d0
        enddo
      enddo

cp    attribute dum(ind) to x(i,j)
      ind=0
      do l=0,lmax
        do m=0,2*l
          ind=ind+1
          x(l,m)=dum(ind)
cp          write(6,*) l,m,x(l,m)
        enddo
      enddo

cp    calculate power spectra for the cst function functions for each degree
      open(85,file='powerdeg.dat')
      open(95,file='power.dat')

      power=0.d0

      do l=0,lmax,1

      rawsum=0.d0
      powerdeg=0.d0

      do m=0,2*l
        rawsum=rawsum+x(l,m)**2
      enddo

      powerdeg=sqrt(1./(2.*l+1.)*rawsum)
      if(l.gt.0) power=power+powerdeg*powerdeg

      write(6,*) 'Spectrum amplitude at l= ',l,'is ',powerdeg
      write(85,*) l,powerdeg

      enddo

      power=sqrt(1./sqrt(4.*pi)*power)

      write(6,*) 'Total power is ',power
      write(95,*) power

      close(85)
      close(95)

      end program
  


