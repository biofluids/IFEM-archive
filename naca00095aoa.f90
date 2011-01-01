subroutine naca00095aoa(disp,ids,node_alebc,nn_alebc,nsd,nn,accel,initvelocity,xold)
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
real(8) accel(nsd,nn)
real(8) initvelocity(nsd,nn)

length=2.475135-1.873513 
k=(2*3.14)/(2*length)      
norm=1.873513 
f=55.0 
omega=2*3.14*f  
amp=0.01  
alpha=5.0
accel(:,:)=0.0
initvelocity(:,:)=0.0  

do i=1,nn_alebc
	node=node_alebc(i)
	ids(:,node)=0.0
        x=xold(1,node)-norm  
	disp(1,node)=sin((alpha*3.14)/180)*((-1*abs(amp*sin(k*x)*sin(omega*tt)))-(-1*abs(amp*sin(k*x)*sin(omega*(tt-dt))))) 
        disp(2,node)=cos((alpha*3.14)/180)*((-1*abs(amp*sin(k*x)*sin(omega*tt)))-(-1*abs(amp*sin(k*x)*sin(omega*(tt-dt))))) 
	!initvelocity(1,node)=sin((alpha*3.14)/180)*((-1*abs(amp*sin(k*x)*omega*cos(omega*0)))) 
	!initvelocity(2,node)=cos((alpha*3.14)/180)*((-1*abs(amp*sin(k*x)*omega*cos(omega*0))))
	accel(1,node)=sin((alpha*3.14)/180)*(-1)*amp*omega*omega*sin(k*x)*sin(omega*tt) 
	accel(2,node)=cos((alpha*3.14)/180)*(-1)*amp*omega*omega*sin(k*x)*sin(omega*tt) 
end do


return 
end subroutine
