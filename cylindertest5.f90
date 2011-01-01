subroutine cylindertest5(disp,ids,node_alebc,nn_alebc,nsd,nn,accel,initvelocity)
 use run_variables, only: its,tt,dt
 implicit none
integer nsd
integer nn
real(8) disp(nsd,nn)
integer ids(nsd,nn)
integer nn_alebc
integer node_alebc(nn_alebc)
integer i
real(8) f
real(8) A
real(8) tmp
real(8) accel(nsd,nn)
real(8) initvelocity(nsd,nn)

f=5.0
A=5.0
disp(:,:)=0.0d0
accel(:,:)=0.0d0
initvelocity(:,:)=0.0d0

do i=1,nn_alebc
	ids(:,node_alebc(i))=0
	disp(1,node_alebc(i))=0  
	disp(2,node_alebc(i))=(5*sin(2*3.14*5*tt))-(5*sin(2*3.14*5*(tt-dt)))
	initvelocity(1,node_alebc(i))=0
	initvelocity(2,node_alebc(i))=2*3.14*f*A
	accel(1,node_alebc(i))=0
	accel(2,node_alebc(i))=-(2*3.14*f)*(2*3.14*f)*A*sin(2*3.14*f*tt)
end do

return

end subroutine cylindertest5 
