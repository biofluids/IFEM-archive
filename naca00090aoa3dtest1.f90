subroutine naca00090aoa3dtest1(disp,ids,node_alebc,nn_alebc,nsd,nn,xold)
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
real(8) deg
real(8) predeg
real(8) change 
real(8) xx
real(8) yy
real(8) dist
real(8) angle
real(8) newang
real(8) newxx
real(8) newyy
integer node 
real(8) xold(nsd,nn)

disp(:,:)=0.0d0
deg=5*sin(2*3.14*tt)
predeg=5*sin(2*3.14*(tt-dt))
change=(deg-predeg)*(3.14/180)

do i=1,nn_alebc
	node=node_alebc(i)
	ids(:,node_alebc(i))=0
	xx=xold(1,node)-3.4925
	yy=xold(2,node)
	!disp(1,node_alebc(i))=((xx*cos(change))-(yy*sin(change)))-xx  
	!disp(2,node_alebc(i))=((xx*sin(change))+(yy*cos(change)))-yy
	disp(1,node_alebc(i))=0
	disp(2,node_alebc(i))=(.5*sin(2*3.14*tt))-(.5*sin(2*3.14*(tt-dt))) 
	disp(3,node_alebc(i))=0 
end do

return

end subroutine  naca00090aoa3dtest1
