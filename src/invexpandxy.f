      program invexpandxy

c     least squares expansion in spherical harmonics of data
c     on a grid
cp    PK obtained from JR, 2013
cp    removed multiplication factor of 0.01 in writing out the .xyz file
cp    Paula Koelemeijer, February 2015

      parameter (MXL=40)
      parameter (MXLENY=(MXL+1)**2)
      dimension atd(MXLENY),y(MXLENY)
      double precision evc(MXLENY),x(MXLENY),eigv(MXLENY)
      double precision f1,sum,w,damp
      dimension wnorm(MXLENY)

      character*80 afl,evcfl,getunx,grdfl,ofl

c     call chekcl('|  :r:2:Input xyz file, output .raw file'
c    1          //'|-a:r:1:input .a file'
c    1          //'|-evc:r:1:evc file'
c    1          //'|-damp:o:1:damping'
c    1          //'|')

      read(5,*) grdfl
      read(5,*) ofl
      read(5,*) afl
      read(5,*) evcfl
      read(5,*) dum
      damp=dble(dum)

      open(21,file=afl,status='old',form='unformatted')
      open(22,file=evcfl,status='old',form='unformatted')
      open(23,file=grdfl,status='old')

      read(22) lmax
      leny=(lmax+1)**2

12    read(23,*,end=100) xlon,xlat,grdv
       read(21) xlon2,xlat2,(y(k),k=1,leny)
       if(xlon.ne.xlon2) stop 'inconsistent .a and .xyz file'
       if(xlat.ne.xlat2) stop 'inconsistent .a and .xyz file'
       do k=1,leny
        atd(k)=atd(k)+grdv*y(k)
       enddo
      goto 12

100   continue

      do i=1,leny
        read(22) eigv(i),(evc(k),k=1,leny)

        if(eigv(i).gt.1.d-7*eigv(1)) then
         sum=0.
         do j=1,leny
          sum=sum+dble(atd(j))*evc(j)
         enddo

         f1=1.d0/(eigv(i)+damp)
         w=sum*f1
         do j=1,leny
          x(j)=x(j)+w*evc(j)
         enddo
        endif
      enddo

c     normaliseer harmonics
      call normylm(lmax,wnorm)
      do i=1,leny
       x(i)=x(i)*dble(wnorm(i))
      enddo

      open(24,file=ofl,status='unknown')
      write(24,'(i3)') lmax
c-- Hendrik multiplied by 0.01 .... WHY?
      write(24,'(5e16.8)') (x(i)*.01,i=1,leny)
c--   write(24,'(5e16.8)') (x(i),i=1,leny)
      close(24)

      end







