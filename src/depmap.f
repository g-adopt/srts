      program depmaphj_pk

c-- Program to calculate horizontal slice through a .sph model
c--
c-- Program from AD/JR - March 2013
cp  adapted to read from input line - not stdin

      parameter (MXLH=40)
      parameter (MXLENY=(MXLH+1)**2)
      parameter (MXDEP=21)
      parameter (MXNDP=200)
      parameter (MXSPL=MXDEP+3)
      parameter (MXMSZ=MXSPL*MXLENY)

      dimension yraw(MXLENY)
      dimension maskl(MXLH),masks(MXSPL)
      dimension x(MXMSZ)
      dimension spln(MXDEP)
      dimension dep(MXNDP),ldep(MXNDP)
      character*10 depstr(MXNDP)
      character*200 getunx,mfl,mofl,outfl,dumstr
      character*120 line1

33    format(a200)
c-- Get the modelname
      read(5,33) mfl
      lm = lnblnk(mfl)
      write(6,*) 'model=', mfl
      write(6,*) 'lm= ', lm
      if(mfl(lm-2:lm).eq.'sph') then 
        imod=1
        iwhole=1
        idptp=1
        ncrtopo=3
      else
        stop 'unknown model format'
      endif

c-- Read the depth from standard input
      i=0
10    i=i+1
      read(5,*,end=20) dep(i)
      idum=int(dep(i))
      if(idum.ge.1000) write(dumstr,'(i4)') idum
      if(idum.ge.100.and.idum.lt.1000) write(dumstr,'(i3)') idum
      if(idum.ge.10.and.idum.lt.100) write(dumstr,'(i2)') idum
      ldep(i) = lnblnk(dumstr)
      if(ldep(i).gt.10) stop 'depth too long'
      if(ldep(i).eq.0) goto 20
      if(ldep(i).eq.2) then
       depstr(i)='00'//dumstr(1:ldep(i))
       ldep(i)=ldep(i)+2
      else if(ldep(i).eq.3) then
       depstr(i)='0'//dumstr(1:ldep(i))
       ldep(i)=ldep(i)+1
      else
       depstr(i)=dumstr(1:ldep(i))
      endif
15    goto 10
20    continue
      ndep=i-1

      do i=1,ndep
         write(6,*) 'n= ',i,' dep= ',dep(i)
      enddo

c      call chekcl('|  :r:1:Input model file '
c     1          //'|-d:r:1,*:depth'
c     1          //'|')

c      i=0
c10    i=i+1
c      dumstr=getunx('-d',i,ldep(i))
c      if(ldep(i).gt.10) stop 'depth too long'
c      if(ldep(i).eq.0) goto 20
c      if(ldep(i).eq.2) then
c       depstr(i)='00'//dumstr(1:ldep(i))
c       ldep(i)=ldep(i)+2
c      else if(ldep(i).eq.3) then
c       depstr(i)='0'//dumstr(1:ldep(i))
c       ldep(i)=ldep(i)+1
c      else
c       depstr(i)=dumstr(1:ldep(i))
c      endif
c      dep(i)=reunx('-d',i,ll)
c15    goto 10
c20    continue
c      ndep=i-1

c      mfl=getunx(' ',1,lm)

      if(mfl(lm-2:lm).eq.'sph') then 
       iwhole=1
       idptp=1
       ncrtopo=3
       call splhsetup()
ca      else if(mfl(lm-2:lm).eq.'spe') then
ca       iwhole=1
ca       idptp=4
ca       ncrtopo=3
ca       call spesetup
      else if(mfl(lm-2:lm).eq.'c13') then
       iwhole=1
       idptp=2
       ncrtopo=2
      else if(mfl(lm-2:lm).eq.'p57') then
       iwhole=0507
       idptp=3
       ncrtopo=2
      else if(mfl(lm-2:lm).eq.'t75') then
       iwhole=0705
       idptp=2
       ncrtopo=2
      else
       stop 'unknown model format'
      endif

c    take directory out of model name
      ip1=lm
      do while(ip1.gt.0.and.mfl(ip1:ip1).ne.'/')
        ip1=ip1-1
      enddo
      mofl=mfl(ip1+1:lm-3)
      lo=istlen(mofl)

c    read model
      open(11,file=mfl,status='old')
      read(11,'(a)') line1
      call rsphhead(line1,lmx,ndmn,ndmx,ndp)
      write(6,*) lmx,ndmn,ndmx,ndp
      natd=(lmx+1)**2
      ind=(ndmn-1)*natd+1
      do i=ndmn,ndmx
       do j=0,lmx
        ind1=ind+2*j
        read(11,'(11e12.4)')(x(k),k=ind,ind1)
c       write(6,*) (x(k),k=ind,ind1)
        ind=ind1+1
       enddo
      enddo
      close(11)

      if(lmx.gt.MXLH) stop'lmx.gt.MXLH'
      if(ndmx.gt.MXSPL) stop'ndmx.gt.MXSPL'
      nbeg=max(ncrtopo+1,ndmn)
      ntot=ndmx-ncrtopo

      write(6,*) lmx,nbeg,ndmx,ntot

c    Calculate the spline basis functions at a regular grid
      do ii=1,ndep
       do ip=1,ntot
        call getfdp(iwhole,idptp,dep(ii),ip,fp)
        write(6,*) 'ii,dep(ii),ip,fp= ',ii,dep(ii),ip,fp
        spln(ip)=fp
       enddo

       do i=1,natd
        yraw(i)=0.
       enddo

c      calculate map
       do i=1,ntot
        ind=(i+ncrtopo-1)*natd+1
        write(6,*) i,ind
        call saxpy(natd,spln(i),x(ind),1,yraw,1)
       enddo

       outfl=mofl(1:lo)//depstr(ii)(1:ldep(ii))//'.raw'
       open(14,file=outfl,status='unknown')
       write(14,'(i3,x,i6,x,f11.3)') lmx
       write(14,'(5e16.8)') (yraw(i),i=1,natd)
       close(14)
      enddo


100   continue

      end

c -----------------------------------------------------

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

c --------------------------------------------------------------------
      SUBROUTINE LEGNDR(THETA,L,M,X,XP,XCOSEC)
      DIMENSION X(*),XP(*),XCOSEC(*)
      DOUBLE PRECISION SMALL,SUM,COMPAR,CT,ST,FCT,COT,FPI,X1,X2,X3,
     1F1,F2,XM,TH,DFLOAT
      DATA FPI/12.56637062D0/
      DFLOAT(I)=FLOAT(I)
      SUM=0.D0
      LP1=L+1
      TH=THETA
      CT=DCOS(TH)
      ST=DSIN(TH)
      MP1=M+1
      FCT=DSQRT(DFLOAT(2*L+1)/FPI)
      SFL3=SQRT(FLOAT(L*(L+1)))
      COMPAR=DFLOAT(2*L+1)/FPI
      DSFL3=SFL3
      SMALL=1.D-16*COMPAR
      DO 1 I=1,MP1
      X(I)=0.
      XCOSEC(I)=0.
    1 XP(I)=0.
      IF(L.GT.1.AND.ABS(THETA).GT.1.E-5) GO TO 3
      X(1)=FCT
      IF(L.EQ.0) RETURN
      X(1)=CT*FCT
      X(2)=-ST*FCT/DSFL3
      XP(1)=-ST*FCT
      XP(2)=-.5D0*CT*FCT*DSFL3
      IF(ABS(THETA).LT.1.E-5) XCOSEC(2)=XP(2)
      IF(ABS(THETA).GE.1.E-5) XCOSEC(2)=X(2)/ST
      RETURN
    3 X1=1.D0
      X2=CT
      DO 4 I=2,L
      X3=(DFLOAT(2*I-1)*CT*X2-DFLOAT(I-1)*X1)/DFLOAT(I)
      X1=X2
    4 X2=X3
      COT=CT/ST
      COSEC=1./ST
      X3=X2*FCT
      X2=DFLOAT(L)*(X1-CT*X2)*FCT/ST
      X(1)=X3
      X(2)=X2
      SUM=X3*X3
      XP(1)=-X2
      XP(2)=DFLOAT(L*(L+1))*X3-COT*X2
      X(2)=-X(2)/SFL3
      XCOSEC(2)=X(2)*COSEC
      XP(2)=-XP(2)/SFL3
      SUM=SUM+2.D0*X(2)*X(2)
      IF(SUM-COMPAR.GT.SMALL) RETURN
      X1=X3
      X2=-X2/DSQRT(DFLOAT(L*(L+1)))
      DO 5 I=3,MP1
      K=I-1
      F1=DSQRT(DFLOAT((L+I-1)*(L-I+2)))
      F2=DSQRT(DFLOAT((L+I-2)*(L-I+3)))
      XM=K
      X3=-(2.D0*COT*(XM-1.D0)*X2+F2*X1)/F1
      SUM=SUM+2.D0*X3*X3
      IF(SUM-COMPAR.GT.SMALL.AND.I.NE.LP1) RETURN
      X(I)=X3
      XCOSEC(I)=X(I)*COSEC
      X1=X2
      XP(I)=-(F1*X2+XM*COT*X3)
    5 X2=X3
      RETURN
      END


c ------------------------------------------------------------------------------

      subroutine wint2ch(int,ch,lrch)

      character*2 ch

      if(int.lt.10) then
       lrch=1
       write(ch,'(i1)') int
      else if(int.ge.10.and.int.lt.100) then
       lrch=2
       write(ch,'(i2)') int
      else
       stop 'wint2ch; int.gt.99'
      endif

      end

c ------------------------------------------------------------------------
