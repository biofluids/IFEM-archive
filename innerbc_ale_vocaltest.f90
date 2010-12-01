subroutine innerbc_ale_vocaltest(disp,ids,node_alebc,nn_alebc,nsd,nn,xold)
! Test vocal case moving boundary, rigid roatation for a rectangular
 use run_variables, only: its,tt,dt
 implicit none
integer nsd
integer nn
real(8) disp(nsd,nn)
integer ids(nsd,nn)
integer nn_alebc
integer node_alebc(nn_alebc)
real(8) xold(nsd,nn)
integer i
real(8) theta1
real(8) theta2
real(8) x1
real(8) y1
real(8) dx
real(8) dy
real(8) tmp
theta2 = 3.14/2-3.14*tt
theta1 = 3.14/2-3.14*(tt-dt)


do i=1,nn_alebc
	x1=xold(1,node_alebc(i))
	y1=xold(2,node_alebc(i))
if (y1 .lt. 28) then
	dx=y1/sin(theta1)*(cos(theta2)-cos(theta1))
	dy=y1/sin(theta1)*(sin(theta2)-sin(theta1))
        ids(:,node_alebc(i))=0
        disp(1,node_alebc(i))=dx
        disp(2,node_alebc(i))=dy   ! Test the displacement only in x and is a constant
else
	tmp=56.0-y1
	dx=tmp/sin(theta1)*(cos(theta2)-cos(theta1))
        dy=tmp/sin(theta1)*(sin(theta1)-sin(theta2))
        ids(:,node_alebc(i))=0
        disp(1,node_alebc(i))=dx
        disp(2,node_alebc(i))=dy   ! Test the displacement only in x and is a constant
end if
end do
!==========================================
! Test non-moving
!	ids(:,node_alebc(i))=0
!        disp(1,node_alebc(i))=0.0
!        disp(2,node_alebc(i))=0.0   ! Test the displacement only in x and is a constant
return

end subroutine innerbc_ale_vocaltest
