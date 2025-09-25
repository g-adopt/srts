c ----------------------------------------------------------------------

      subroutine rmod(infl,x,lmx,nsmx,mskl,msks,ndmx,ndmn)

      character*80 infl
      character*200 line,dum1
      character*2 c1
      dimension x(*)
      dimension mskl(*), msks(*)

      open(10,file=infl,status='old')
      read(10,*) lmx,dum,nsmx
      rewind(10)
      write(6,*) 'RMOD: lmx= ',lmx,' nsmx= ', nsmx

      if(lmx.gt.100.or.nsmx.gt.100) stop'lmx.gt.100.or.nsmx.gt.100'

      read(10,'(a)') line
      ll=len(line)

c     find out where lmx is on the line
      do i=1,ll-1
       c1=line(i:i+1)
       read(c1,'(i2)') ltst
       if(ltst.eq.lmx) imb=i+3
      enddo

      dum1=line(imb:imb+lmx-1)
      read(dum1,11) (mskl(i),i=1,lmx)
11    format(100i1)
      write(6,12) (mskl(i),i=1,lmx)
12    format(' RMOD: mskl = ',100i1)

c     find out where lmx is on the line
      do i=imb+lmx,ll-1
       c1=line(i:i+1)
       read(c1,'(i2)') ltst
       if(ltst.eq.nsmx) imb=i+3
      enddo

      dum1=line(imb:imb+nsmx-1)
      read(dum1,11) (msks(i),i=1,nsmx)
      write(6,13) (msks(i),i=1,nsmx)
13    format(' RMOD: msks = ',100i1)

      ndmn=1
      do while(msks(ndmn).ne.1.and.msks(ndmn).ne.2)
       ndmn=ndmn+1
      enddo

      ndmx=ndmn
      do while((msks(ndmx).eq.1.or.msks(ndmx).eq.2).and.ndmx.le.nsmx)
       ndmx=ndmx+1
      enddo
      ndmx=ndmx-1

      write(6,*) 'RMOD: ndmn= ',ndmn,' ndmx= ', ndmx

      if(mskl(lmx).ne.1) stop 'do not know how to deal with this '

      natd=(lmx+1)**2
      ind=(ndmn-1)*natd+1
      do i=ndmn,ndmx
       do j=0,lmx
        ind1=ind+2*j
        read(10,'(11e12.4)',end=100)(x(k),k=ind,ind1)
        ind=ind1+1
       Enddo
      Enddo

      goto 200

 100  stop 'incompatible model header'

 200  continue
      
      end
