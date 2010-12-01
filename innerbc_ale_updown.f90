subroutine innerbc_ale_updown(disp,ids,node_alebc,nn_alebc,nsd,nn,xold)
! Two rigid vocal folds move up and down at fixed frequency
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
real(8) uplimit
real(8) downlimit
real(8) T ! movement frequency
real(8) vel ! movement velocity
real(8) dy ! y direction displacement
real(8) A
real(8) omiga
real(8) dgap
integer node

T=1.0d0
uplimit=0.48
downlimit=0.4
dgap=(1-2*downlimit)-(1-2*uplimit)
A=0.5*dgap*3.14/T
omiga=2*3.14/T
vel=A*(sin(omiga*tt))
dy=vel*dt

do i=1,nn_alebc
	node=node_alebc(i)
	ids(:,node)=0
	disp(1,node)=0.0
	if ( xold(2,node) .lt. 0.5 ) then
		disp(2,node)=dy
		if (xold(2,node) .eq. 0.0d0) disp(2,node)=0.0d0
	else
		disp(2,node)=-dy
		if (xold(2,node) .eq. 1.0d0) disp(2,node)=0.0d0
	end if
end do




return 
end subroutine
