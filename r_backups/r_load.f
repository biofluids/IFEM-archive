c     
c     load assignment
c     
      subroutine r_load
      implicit real*8 (a-h,o-z)
      include 'r_common'
      if (ntprint .eq. 1) then
         write(4,100) iti,tfun(2)
 100     format(1x, 'tfun(',i4,')=',e23.7,';'/)
         write(4,101) iti,iti*dt_solid
 101     format(1x, 'time(',i4,')=',e23.7,';'/)
      endif

      xtime=tfun(ntfun)
c
      do k=1,numeb
         boupo(k,1)=boup(k,1)*xtime
         boupo(k,2)=boup(k,2)*xtime
	   boupo(k,3)=boup(k,3)*xtime
	enddo
c     concentrated load  
      do i=1,numfn
         fnodo(nodefn(i),ndirfn(i))=fnod(nodefn(i),ndirfn(i))*
     $        xtime
	enddo
      return
      end






