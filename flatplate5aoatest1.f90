subroutine flatplate5aoatest1(disp,ids,node_alebc,nn_alebc,nsd,nn,xold)
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
real(8) f ! movement frequency
integer node
real(8) length
real(8) norm 
real(8) amp 
real(8) alpha 

length=3.669362      
!norm=minval(xold(1,node_alebc(:)))
norm=0.681048  
f=10.0   
amp=0.001   
alpha=5.0  

do i=1,nn_alebc
	node=node_alebc(i)
	ids(:,node)=0.0
	disp(1,node)=amp*sin(3.14*((xold(1,node)-norm)/length))*cos(f*2*3.14*tt)*f*2*3.14*dt*sin(alpha*(3.14/180))    
	disp(2,node)=amp*sin(3.14*((xold(1,node)-norm)/length))*cos(f*2*3.14*tt)*f*2*3.14*dt*cos(alpha*(3.14/180))
end do


return 
end subroutine
