!     
!     load assignment
!     
subroutine r_load
  use r_common
  implicit none

  integer :: i
  real*8 :: xtime
!      if (ntprint .eq. 1) then
!         write(4,100) iti,tfun(2)
! 100     format(1x, 'tfun(',i4,')=',e23.7,';'/)
!         write(4,101) iti,iti*dt
! 101     format(1x, 'time(',i4,')=',e23.7,';'/)
!      endif

  xtime=tfun(ntfun)

  do i=1,numeb
     boupo(i,1)=boup(i,1)*xtime
     boupo(i,2)=boup(i,2)*xtime
     boupo(i,3)=boup(i,3)*xtime
  enddo

 !...concentrated load  
  do i=1,numfn
     fnodo(nodefn(i),ndirfn(i))=fnod(nodefn(i),ndirfn(i))*xtime
  enddo
  
  return
end subroutine r_load






