      program correlatorsph

cp    correlation C12 (equation 14 Smith and Masters 1989)
cp    and treatise on geophysics Forte 2007 Vol 1 Ch 23 Eq 4 
cp    adapted for sph parameterisation
cp    March 2014

      parameter(MXL=40)
      parameter(MXLENY=(MXL+1)**2)

      character*256 rawfile1,rawfile2
      integer lmax1,lmax2,ncoef1,ncoef2,l,m,i,lmax

      double precision dum1(0:MXLENY),dum2(0:MXLENY)
      double precision x1(0:MXL,0:MXL*2),x2(0:MXL,0:MXL*2)

      double precision x1sum,x2sum,x1sq,x2sq,cross,corr
      double precision x1sumdeg,x2sumdeg,x1sqdeg,x2sqdeg,crossdeg,corrdeg

      pi=acos(-1.)
      write(6,*) 'Enter .raw file 1'
      read(5,'(a)') rawfile1
      write(6,*) 'Enter .raw file 2'
      read(5,'(a)') rawfile2

      write(6,*) 'reading in .raw file 1'
cp    read in model rawfile 1
      open(45,file=rawfile1)
      read(45,*) lmax1
      write(6,*) 'rawfile 1=',rawfile1
      write(6,*) 'model lmax 1 =',lmax1
      if(lmax1.gt.MXL) stop 'stop >> lmax1.gt.MXL'
      ncoef1=(lmax1+1)**2
      read(45,*) (dum1(i),i=1,ncoef1)
      close(45)

      write(6,*) 'reading in .raw file 2'
cp    read in model rawfile 2
      open(55,file=rawfile2)
      read(55,*) lmax2
      write(6,*) 'rawfile 2=',rawfile2
      write(6,*) 'model lmax 2 =',lmax2
      if(lmax2.gt.MXL) stop 'stop >> lmax2.gt.MXL'
      ncoef2=(lmax2+1)**2
      read(55,*) (dum2(i),i=1,ncoef2)
      close(55)

      if(lmax1.gt.lmax2) lmax=lmax2
      if(lmax1.le.lmax2) lmax=lmax1
      write(6,*) 'using lmax is ',lmax

      write(6,*) 'copying rawfiles to model'
cp    initialize arrays
      do l=0,lmax
        do m=0,2*l
          x1(l,m)=0.d0
          x2(l,m)=0.d0
        enddo
      enddo

cp    attribute dum(ind) to x(i,j)
      ind=0
      do l=0,lmax
        do m=0,2*l
          ind=ind+1
          x1(l,m)=dum1(ind)
          x2(l,m)=dum2(ind)
cp          write(6,*) l,m,x1(l,m),x2(l,m)
        enddo
      enddo

      write(6,*) 'calculating corrrelation'
cp    calculate correlation for each degree
      open(85,file='corrdeg.dat')
      open(95,file='corr.dat')

      x1sum=0.d0
      x2sum=0.d0
      x1sq=0.d0
      x2sq=0.d0
      cross=0.d0
      corr=0.d0

      do l=0,lmax,1

      x1sumdeg=0.d0
      x2sumdeg=0.d0
      x1sqdeg=0.d0
      x2sqdeg=0.d0
      crossdeg=0.d0
      corrdeg=0.d0

      do m=0,2*l
        x1sumdeg=x1sumdeg+x1(l,m)**2
        x2sumdeg=x2sumdeg+x2(l,m)**2
        crossdeg=crossdeg+x1(l,m)*x2(l,m)
      enddo

      x1sqdeg=sqrt(x1sumdeg)
      x2sqdeg=sqrt(x2sumdeg)
      corrdeg=crossdeg/(x1sqdeg*x2sqdeg)

      write(6,*) 'Correlation at l= ',l,'is ',corrdeg
      write(85,*) l,corrdeg

      if(l.ne.0) x1sum=x1sum+x1sumdeg
      if(l.ne.0) x2sum=x2sum+x2sumdeg
      if(l.ne.0) cross=cross+crossdeg

      enddo

      x1sq=sqrt(x1sum)
      x2sq=sqrt(x2sum)
      corr=cross/(x1sq*x2sq)

      write(6,*) 'Total correlation is ',corr
      write(95,*) corr

      close(85)
      close(95)

      end program
      









