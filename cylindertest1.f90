subroutine cylindertest1(disp,ids,node_alebc,nn_alebc,nsd,nn)
 use run_variables, only: its,tt,dt
 implicit none
integer nsd
integer nn
real(8) disp(nsd,nn)
integer ids(nsd,nn)
integer nn_alebc
integer node_alebc(nn_alebc)
integer i
real(8) tmp

tmp=1.0*(sin(2*3.14*tt)-sin(2*3.14*(tt-dt)))
disp(:,:)=0.0d0

do i=1,nn_alebc
	ids(:,node_alebc(i))=0
	disp(1,node_alebc(i))=0  
	!disp(2,node_alebc(i))=(.5*sin(2*3.14*tt))-(.5*sin(2*3.14*(tt-dt)))
	disp(2,node_alebc(i))=0  
end do

return

end subroutine cylindertest1
