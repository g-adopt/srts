c -----------------------------------------------------

      subroutine inittrf(nsmn,nsmx,dir,ider,eval,evcs)

c     this subroutine initialises a transformation matrix defined in terms
c     of eigenvectors.

      dimension eval(*),evcs(*)

      character*80 dir,rid
      character*200 fl
      character*2 chmn,chmx,chdr

      call wint2ch(nsmn,chmn,lcmn)
      call wint2ch(nsmx,chmx,lcmx)
      call wint2ch(ider,chdr,lcdr)
      rid=chmn(1:lcmn)//'_'//chmx(1:lcmx)//'.der'//chdr(1:lcdr)
      lrid=istlen(rid)

      ld=istlen(dir)
      fl=dir(1:ld)//'/evmat.'//rid(1:lrid)
      lfl=istlen(fl)
      write(6,*) 'OPENING EVC FILE :',fl(1:lfl)

      open (15,file=fl,form='unformatted',status='old')
      read(15) ndim
      if (ndim.ne.(nsmx-nsmn+1)) stop 'inittrf, unexpected ndim'
      do i=1,ndim
       id=(i-1)*ndim
       read(15) eval(i),(evcs(id+j),j=1,ndim)
      enddo
      close(15)

      end
