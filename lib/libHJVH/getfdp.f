      subroutine getfdp(iw,idptp,dep,ip,fp)

c     iw: 1 = whole mantle anything else is split
c     idptp: 1 = sph, 2 = chebyshev, 3 = legendre
c     ip: parameter number numbering starts at 1

      rcmb=3480.
      r670=5701.
      rmoho=6346.619
      rearth=6371.

      r=rearth-dep

      nlow=mod(iw,100)+1
      nup=((iw-nlow+1)/100)+1
      if(ip.le.nup) then
       iup=1
      else
       iup=0
      endif

c     whole mantle parameterisation or split
      if(iw.eq.1) then
       xd=-1.+2.*(r-rcmb)/(rmoho-rcmb)
       ipc=ip
      else
       if(r.lt.r670) then
        if(iup.ne.0) goto 1000
        xd=-1.+2.*(r-rcmb)/(r670-rcmb)
        ipc=ip-nup
       else
        if(iup.ne.1) goto 1000
        xd=-1.+2.*(r-r670)/(rmoho-r670)
        ipc=ip
       endif
      endif


      if(xd.gt.1.or.xd.lt.-1) stop 'getfd: bad depth!'

c     spline, chebyshev, or legendre
      if(idptp.eq.1) then
       fp=splh(ipc-1,xd)

c-- Jeroen Ritsema: only sph format for now ...
c     else if(idptp.eq.2) then
c      fp=tn(ipc-1,xd)
c     else if(idptp.eq.3) then
c      fp=pn(ipc-1,xd)

      else
       stop 'getfdp: unknown idptp'
      endif

      return

1000  continue
      fp=0.
      return

      end
