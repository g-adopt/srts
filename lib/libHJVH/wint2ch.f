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
