      program mkexpmatxy

c     construct matrices to expand data on a nxn degree grid
c     to spherical harmonics
c     Jeroen Ritsema
cp    PK obtained from JR, 2013

      parameter(LMX=40)
      parameter(MXLENY=(LMX+1)**2)
      parameter(MXWORK=(LMX+1)*4)
      parameter(MXMAT=(MXLENY*(MXLENY+1))/2)

      dimension y(MXLENY),wk1(MXWORK),wk2(MXWORK),wk3(MXWORK)
      dimension wnorm(MXLENY)
      double precision ata(MXMAT),eigv(MXLENY)
      double precision d1(MXLENY),d2(MXLENY),d3(MXLENY),d4(MXLENY)
      double precision d5(MXLENY),d6(MXLENY),d7(MXLENY),d8(MXLENY)

      character*80 afl,evcfl,xyfl,getunx

c     call chekcl('| :r:1:input x,y(,z) file (order= lon, lat!!!)'
c    1          //'|-a:r:1:output .a file'
c    1          //'|-evc:r:1:output .evc file'
c    1          //'|-lmax:r:1:maximum degree and order for expansion'
c    1          //'|')

      read(5,*) xyfl
      read(5,*) afl
      read(5,*) evcfl
      read(5,*) lmax

      if(lmax.gt.LMX) stop 'lmx.gt.LMX'
      leny=(lmax+1)**2

      do i=1,MXMAT
       ata(i)=0. 
      enddo

      open(19,file=xyfl,status='old')
      open(21,file=afl,status='unknown',form='unformatted')
      open(22,file=evcfl,status='unknown',form='unformatted')

c     get normalisation
      call normylm(lmax,wnorm)

      nread=0
10    read(19,*,end=100) xlon,xlat
       nread=nread+1
       if(mod(nread,1000).eq.0) write(6,'(i8,'' points read'')') nread
       call ylm(xlat,xlon,lmax,y,wk1,wk2,wk3)
       do k=1,leny
        y(k)=y(k)*wnorm(k)
       enddo

       ind=0
       do jj=1,leny
        do ii=jj,leny
         ind=ind+1
         ata(ind)=ata(ind)+dble(y(ii)*y(jj))
        enddo
       enddo

       write(21) xlon,xlat,(y(k),k=1,leny)
      goto 10

100   continue

      write(6,*) 'decomposing......'
      write(22) lmax
      call ahouse2(leny,22,ata,d1,d2,d3,d4,d5,d6,d7,d8,eigv)
 
      end
