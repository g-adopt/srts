      program sphadd

c-- adds up the coefficients of two SPH files
cp  PK obtained from JR, 2013

      parameter (MXLH=40)
      parameter (MXLENY=(MXLH+1)**2)
      parameter (MXDEP=21)
      parameter (MXSPL=MXDEP+3)
      parameter (MXMSZ=MXSPL*MXLENY)
      parameter (MAXP=1024)
      dimension d0(MXLENY)
      dimension maskl1(MXLH),masks1(MXSPL)
      dimension maskl2(MXLH),masks2(MXSPL)
      dimension spl(MXDEP)
      dimension x(MXMSZ),y(MXMSZ),summ(MXMSZ)
      dimension spln(MAXP,MXDEP),rad(MAXP),wl(MAXP)
      character*80 getunx,mfl1,mfl2,outfl
      character*120 strt

c     call chekcl('|  :r:2:Input model .sph files '
c    1          //'|-o:r:1:output file (1-2)'
c    1          //'|')

111   format(a80)
      read(5,111) mfl1
      read(5,111) mfl2
      read(5,111) outfl

      pi4=16.*atan(1.)
      fc=sqrt(1./pi4)

c    read sph model1
c     write(6,*) 'file 1'
      call rmod2(mfl1,x,lmx1,nsmx1,maskl1,masks1,ndmx1,ndmn1)
      if(lmx1.gt.MXLH) stop'lmx1.gt.MXLH'
      if(ndmx1.gt.MXSPL) stop'ndmx1.gt.MXSPL'
      natd1=(lmx1+1)**2
c     write(6,*) 'lmx1,ndmx1', lmx1,ndmx1
c     write(6,*) 'maskl1', maskl1

c    read sph model2
c     write(6,*) 'file 2'
      call rmod2(mfl2,y,lmx2,nsmx2,maskl2,masks2,ndmx2,ndmn2)
      if(lmx2.gt.MXLH) stop'lmx2.gt.MXLH'
      if(ndmx2.gt.MXSPL) stop'ndmx2.gt.MXSPL'
      natd2=(lmx2+1)**2
c     write(6,*) 'lmx1,ndmx1', lmx2,ndmx2
c     write(6,*) 'maskl1', maskl2

      if(lmx1.eq.lmx2.and.nsmx1.eq.nsmx2.and.ndmx1.eq.ndmx2.and.ndmn1.eq.ndmn2) then
       do i=1,lmx1
        if(maskl1(i).ne.maskl2(i)) stop 'cannot deal with maskl1.ne.maskl2'
       enddo
       do i=1,nsmx1
        if(masks1(i).ne.masks2(i)) stop 'cannot deal with masks1.ne.masks2'
       enddo

       nl=ndmx1-ndmn1+1
       np=nl*natd1
       ind=(ndmn1-1)*natd1
       do j=1,np
         summ(ind+j)=x(ind+j)+y(ind+j)
       enddo


       open(21,file=outfl,status='unknown')
       call wsphhead(lmx1,ndmn1,ndmx1,strt,lstrt)
c      write(6,*) 'start= ',strt(1:lstrt)
       write(21,*) strt(1:lstrt)
       ind=(ndmn1-1)*natd1+1
       do i=ndmn1,ndmx1
        do j=0,lmx1
         ind1=ind+2*j
         write(21,'(11e12.4)')(summ(k),k=ind,ind1)
         ind=ind1+1
        enddo
       enddo
       close(21)
      else
       stop 'cannot deal with inequalities'
      endif

      end
