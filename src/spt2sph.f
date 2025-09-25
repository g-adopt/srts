      program spt2sph

cp  program to convert .spt file (all parts of model, e.g. ES, EP) to separate .sph files
cp  PK obtained from JR, 2013

      parameter (MXLH=40)
      parameter (MXLENY=(MXLH+1)**2)
      parameter (MXDEP=42)
      parameter (MXMSZ=MXDEP*MXLENY)
      double precision x(MXMSZ)

c     declarations for getmulkern
      parameter(MXPARTYP=7)
      dimension iparsw(MXPARTYP),ipardps(2,MXPARTYP)
      character*3 cparm
      character*80 sptfl,getunx,ofil(MXPARTYP)
      character*120 strt

c     call chekcl('|   :r:1:Input spt model file '
c    1          //'|')
111   format(a80)
      read(5,111) sptfl
      lspt = istlen(sptfl)
      write(6,*) 'sptfl= ', sptfl, ' lspt= ', lspt

      ofil(1)=sptfl(1:lspt-3)//'.SP.sph'
      ofil(2)=sptfl(1:lspt-3)//'.ES.sph'
      ofil(3)=sptfl(1:lspt-3)//'.EP.sph'
      ofil(4)=sptfl(1:lspt-3)//'.DE.sph'
      ofil(5)=sptfl(1:lspt-3)//'.ER.sph'
      ofil(6)=sptfl(1:lspt-3)//'.ZS.sph'
      ofil(7)=sptfl(1:lspt-3)//'.ZP.sph'

      open(10,file=sptfl,status='old')

cp      read(10,'(a3)') cparm
cp      write(6,*) 'cparm= ', cparm
      read(10,10) lmax,mp1,(iparsw(i),i=1,mp1)
10    format(i2,x,i2,x,10i1)
      if(lmax.gt.MXLH) stop 'spt2sph >>> lmax.gt.MXLH'
      if(mp1.gt.MXPARTYP) stop 'spt2sph >>> mp1.gt.MXPARTYP'
      write(6,*) lmax,mp1
      npartot=0
      do i=1,mp1
       if(iparsw(i).eq.1) then
        read(10,'(2i4)') ipardps(1,i),ipardps(2,i)
        npartot=npartot+ipardps(2,i)-ipardps(1,i)+1
       endif
      enddo

      ind=1
      do i=1,npartot
       do j=0,lmax
        ind1=ind+2*j
        read(10,'(11e12.4)')(x(k),k=ind,ind1)
        ind=ind1+1
       enddo
      enddo

      ind=1
      do i=1,mp1
       if(iparsw(i).eq.1) then
        write(6,*) 'writing to ',ofil(i)
        open(11,file=ofil(i),status='unknown')
        idp1=ipardps(1,i)
        idp2=ipardps(2,i)
        call wsphhead(lmax,idp1,idp2,strt,lstrt)
        write(11,*) strt(1:lstrt)
        do j=idp1,idp2
         do k=0,lmax
          ind1=ind+2*k
          write(11,'(11e12.4)')(x(l),l=ind,ind1)
          ind=ind1+1
         enddo
        enddo
        close(11)
       endif
      enddo


      end

c ---------------------------------------------------------------------
 
      subroutine wsphhead(lmax,idp1,idp2,strt,lstrt)
      parameter(MAXL=100)
      parameter(MAXP=24)
      dimension lask(0:MAXL),mask(MAXP)
 
      character*120 str1,str2,strt
 
      if(lmax.gt.MAXL) stop 'wsphhead >>> lmax too big for format statement'
      if(idp2.gt.MAXP) stop 'wsphhead >>> idp2 too large'
      do i=0,MAXL
       if(i.le.lmax) then
        lask(i)=1
       else
        lask(i)=0
       endif
      enddo
 
      do i=1,24
       mask(i)=0
      enddo
      do i=idp1,idp2
       mask(i)=1
      enddo
 
      write(str1,10) lmax,(lask(i),i=0,lmax)
10    format(10x,i4,x,100i1)
      write(str2,20) 24,(mask(i),i=1,24)
20    format(i4,x,24i1)
 
      len1=istlen(str1)
      len2=istlen(str2)
      if((len1+len2).gt.120) stop 'wsphhead >>> increase dimension of strt'
      strt=str1(1:len1)//str2(1:len2)
      lstrt=len1+len2+1
      end

