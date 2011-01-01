subroutine flatplate0aoatest2(disp,ids,node_alebc,nn_alebc,nsd,nn,xold)
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
real(8) omega
real(8) k
real(8) x

length=.6
k=(2*3.14)/(2*length)      
norm=minval(xold(1,node_alebc(:)))
f=50.0 
omega=2*3.14*50  
amp=0.02  
alpha=0.0  

do i=1,nn_alebc
	node=node_alebc(i)
	ids(:,node)=0.0
	x=(xold(1,node)-norm)/(2*length)   
        disp(1,node)=0.0
        disp(2,node)=-1*abs(amp*sin(k*x)*sin(omega*tt))  
end do


return 
end subroutine
