      program mk3d_res_ESEP

cp    Program to multiply an .sph file with the tomographic resolution matrix
cp    according to Ritsema et al. 2007
cp    PK obtained from JR, 2013

cp    This version has been adapted to work on either 1 parameter or 2 parameter
cp    input files, e.g. both S-wave and P-wave velocity structure. 
cp    Paula Koelemeijer - Feb 2014

      parameter (MXLH=12)
      parameter (MXLENY=(MXLH+1)**2)
      parameter (MXDEP=21*2)
      parameter (MXMSZ=MXDEP*MXLENY)
      parameter (MXSPHR=MXLENY)
      double precision evc(MXMSZ),eigv(MXMSZ),x(MXMSZ),xout(MXMSZ)
      double precision sum,eta,w,f1,tr
      character*80 getunx

      dimension twts(MXMSZ),twtsinv(MXMSZ)

c     these declarations are related to the smoothness constraints
      parameter (MXSTRP=30)
      parameter (MXEVC=MXSTRP*MXSTRP)
      character*130 wtsfl,evcfl,mifl,mofl
      character*130 line

c     declarations for getmulkern
      parameter(MXPARTYP=7)
      dimension iparsw(MXPARTYP),ipardps(2,MXPARTYP)
      dimension iparsw2(MXPARTYP),ipardps2(2,MXPARTYP)
      dimension parwts(MXPARTYP)

c     call chekcl('|   :r:3:Input evc, model input, and model output '
c    1          //'|-c:r:1: cut off factor., f*biggest ev'
c    1          //'|-wf:o:1: wts file'
c    1          //'|')

111   format(a130)
      read(5,111) evcfl
      read(5,111) mifl
      read(5,111) mofl
      imofl=istlen(mofl)
      write(6,*) 'mofl= ', mofl
      write(6,*) 'imofl= ', imofl

      read(5,*) damp
      write(6,*) 'damp= ',damp 
      read(5,111) wtsfl

c--   read evc file
      open(32,file=evcfl,form='unformatted',status='old')
      read(32) lmaxh,numatd2,ndep,icrust,idensi1,idum1,ismth
      read(32) mp1,(iparsw(j),j=1,mp1),(parwts(j),j=1,mp1),
     1          ((ipardps(i,j),i=1,2),j=1,mp1)
      write(6,*) 'lmaxh=', lmaxh
      if(lmaxh.gt.MXLH) stop 'lmaxh.gt.MXLH'
      if(ndep.gt.MXDEP) stop 'ndep.gt.MXDEP'

c--   read input model
      open(29,file=mifl,status='old')
      call rspthead(29,lmx,mp2,iparsw2,ipardps2,ndp)
 
c      call rsphhead(line,lmx,ndmn,ndmx,ndp)
      write(6,*) lmx,ndp
      write(6,*) lmaxh,ndep
      if(lmx.gt.MXLH.or.ndp.gt.MXDEP) stop 'lmx.gt.MXLH.or.ndp.gt.MXDEP'

      natd   = (lmx+1)**2
      lenatd = natd*ndp

      do i=1,lenatd
       x(i)=0.
       xout(i)=0.
      enddo

c--- Hendrik had it as follows:
c     ind=(ndmn-1)*natd+1
c     do i=ndmn,ndmx
c      do j=0,lmx
c       ind1=ind+2*j
c       read(29,'(11e12.4)')(x(k),k=ind,ind1)
c       ind=ind1+1
c      enddo
c     enddo
c--- Jeroen:
      ind=1
      do i=1,ndp
       do j=0,lmx
        ind1=ind+2*j
        read(29,'(11e12.4)')(x(k),k=ind,ind1)
        ind=ind1+1
       enddo
      enddo
      close(29)

      if(lmaxh.ne.lmx)    stop 'STOP >>>> lmaxh.ne.lmx'
      if(ndep.ne.ndp)  then
         write(6,*) 'ndep= ', ndep, ' ndp= ', ndp
         stop 'STOP >>>> ndep.ne.ndp'
      endif
      if(natd.ne.numatd2) stop 'inconsistent numatd'

      if(ismth.eq.1) then
        open(22,file=wtsfl,form='unformatted',status='old')
        read(22) lmaxw,nsmn,nsmx,ndepw,etaz,etah,etai,iderh,iderv
        if(etaz.ne.0) stop' vertical smoothing is not supported!'
        if(lmaxw.ne.lmx.or.ndepw.ne.ndp) then
          stop 'mk3d_res; inconsitent weights file'
        endif
        do i=1,ndp
          id=(i-1)*natd
          read(22) (twts(id+j),j=1,natd)
        enddo
        close(22)
      else
        do i=1,lenatd
          twts(i)=1.
        enddo
      endif

      do i=1,lenatd
       twtsinv(i)=1./twts(i)
      enddo

      i=0
      inxt=1
      tr=0.

c     do i=1,lenatd
      do while(i.lt.lenatd.and.inxt.eq.1)
        i=i+1
        read(32,end=100) eigv(i),(evc(k),k=1,lenatd)
        if(eigv(i).lt.stpfact) then
          inxt=0
          write(6,*) 'STOP BUILDING MODEL'
          write(6,*) 'LAST EIGENVECTOR USED IS NR ',i
          ilstev=i
        endif
        if(i.eq.1) then
         eta=eigv(1)*damp
         write(6,*) eigv(1),eta
         stpfact=eta/5000.
        endif

        sum=0.
        do j=1,lenatd
          sum=sum+twtsinv(j)*x(j)*evc(j)
        enddo
c       f1=1./(eigv(i)+eta)
        f1=eigv(i)/(eigv(i)+eta)
        w=sum*f1
        do j=1,lenatd
          xout(j)=xout(j)+w*evc(j)
        enddo
      enddo
      close(32)

100   continue
      if(i.ne.lenatd.and.i.ne.ilstev) then
        write(6,*) 'REDUCED EVC FILE  !!!???'
        write(6,*) (i-1),'  eigenvectors used'
      endif

      write(6,*) 'eikel'
      if(ismth.eq.1) then
       write(6,*) 'applying a priori model weights'
       do i=1,lenatd
        xout(i)=xout(i)*twts(i)
       enddo
      endif

c    reweight the crustal thickness
      if(icrust.eq.1) then
       do i=1,natd
        xout(i)=xout(i)*1000.
       enddo
      endif

      ind=1
      open(26,file=mofl(1:imofl)//'.spt',status='unknown')
      call wspthead(26,lmx,mp1,iparsw,ipardps)
      Do i=1,ndp
       Do j=0,lmx
        ind1=ind+2*j
        write(26,'(11e12.4)')(xout(k),k=ind,ind1)
        ind=ind1+1
       Enddo
      Enddo
      close(26)

      end

c------------------------------------------------------
      subroutine wspthead(nfil,lmax,mp1,iparsw,ipardps)

      parameter(MXPARTYP=7)
      dimension iparsw(MXPARTYP),ipardps(2,MXPARTYP)

      if(mp1.gt.10) stop 'wspthead >>> number of parameters too large for format statement'
      write(nfil,10) lmax,mp1,(iparsw(i),i=1,mp1)
10    format(i2,x,i2,x,10i1)
      do i=1,mp1
       if(iparsw(i).eq.1) then
        write(nfil,'(2i4)') ipardps(1,i),ipardps(2,i)
       endif
      enddo

      end


c ------------------------------------------------------

      subroutine rspthead(nfil,lmx,mp2,iparsw2,ipardps2,ndp)
 
      parameter(MXPARTYP=7)
      dimension iparsw2(MXPARTYP),ipardps2(2,MXPARTYP)
 
      read(nfil,10) lmx,mp2,(iparsw2(i),i=1,mp2)
      if(mp2.gt.MXPARTYP) stop'cannot store more than MXPARTYP parameter types'
10    format(i2,x,i2,x,10i1)
      ndp=0
      do i=1,mp2
       if(iparsw2(i).eq.1) then
        read(nfil,*) ipardps2(1,i),ipardps2(2,i)
        ndp=ndp+ipardps2(2,i)-ipardps2(1,i)+1
       endif
      enddo
 
      end

