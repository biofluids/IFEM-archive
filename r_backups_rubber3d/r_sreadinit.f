      subroutine r_sreadinit
      implicit real*8 (a-h,o-z)
      include 'r_common'
c     assign nonlinear initial conditions
      ninp=0
 1000 read(1,*) ndum,x1,y1,z1
      ndumtest=ndum-1
      if (ndumtest .ge. 0) then
         ninp=ninp+1
         if (initdir .eq. 1) then
            dis(1,ndum)=x1
            dis(2,ndum)=y1
		  dis(3,ndum)=z1
c
            xindis(1,ndum)=x1
            xindis(2,ndum)=y1
		  xindis(3,ndum)=z1
         endif
         goto 1000
      endif
      return
      end	