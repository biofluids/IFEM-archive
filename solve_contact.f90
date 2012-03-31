

subroutine solve_contact(x_fix,x_wall,norm_fix,norm_wall,x_con,nn_con)

  use fluid_variables,only:nsd
  use interface_variables, only:maxmatrix
  use mpi_variables
  real(8) x_fix(2),x_wall(2),norm_fix(2),norm_wall(2)
  real(8) x_con(2,maxmatrix)
  integer nn_con

  real(8) x1,y1,nx,ny,x2,y2,mx,my,a,b
  real(8) thelta2, thelta1
  real(8) C1,C2,C3,C4,C5,C6,r2,root1,root2

  real(8) x3,y3
  integer i,j
  real(8) pi,thelta,eps
  integer nop !num of points to regenerate
  nop=10
!  nn_con=0
  pi=3.14159
  eps=0.001
  x1=x_fix(1)
  y1=x_fix(2)
  y2=x_wall(2)
  nx=norm_fix(1)
  ny=norm_fix(2)
  mx=norm_wall(1)
  my=norm_wall(2)
  C1=sqrt(mx**2*nx**2 + mx**2*ny**2 + my**2*nx**2 + my**2*ny**2)
  C2=mx*ny*x1 - mx*nx*y1 + mx*nx*y2
  C3=my*ny*x1 - my*nx*y1 + mx*ny*y2
  C4=mx*ny*x1 + my*nx*x1 - mx*nx*y1 + mx*nx*y2 + my*ny*y1 - my*ny*y2
  C5=mx*ny + my*nx
  C6=mx*ny - my*nx

  if((abs(C5).lt.eps) .or. (abs(C6).lt.eps))then
   if(abs(mx).lt.eps) then
    do i=1,30
       nn_con=nn_con+1
       x_con(1,nn_con)=x1
       x_con(2,nn_con)=y1+(y2-y1)/10.0*(i-1)
    end do
   else
     x2=x1-(y2-y1)*my/mx
     do i=1,30
        nn_con=nn_con+1
	x_con(1,nn_con)=x1+(x2-x1)/10.0*(i-1)
	x_con(2,nn_con)=y1+(y2-y1)/10.0*(i-1)
     end do
   end if
    goto 100
  end if


  x2=(y1*C1-y2*C1+C4)/C5

!  if(x2*x1.gt.0) then
  
    a=(C2-(my*nx*(y1*C1-y2*C1+C4))/C5)/C6
    b=(C3-(my*ny*(y1*C1-y2*C1+C4))/C5)/C6
!  else
  if( ((x1-a)*nx+(y1-b)*ny)*((x2-a)*mx+(y2-b)*my).lt.0) then
    x2=(y2*C1-y1*C1+C4)/C5
    a=(C2-(my*nx*(y2*C1-y1*C1+C4))/C5)/C6
    b=(C3-(my*ny*(y2*C1-y1*C1+C4))/C5)/C6
  end if



  r2=sqrt((x1-a)**2+(y1-b)**2)
!  if(myid==0) then
!    write(*,*)'cal=',(x1-a)/r2,(y1-b)/r2
!    write(*,*)'norm=',nx,ny
!  end if
!if(myid==0)write(*,*)'curv for contact=',1.0/r2

  thelta1=atan((y1-b)/(x1-a))
  if((x1-a).gt.0) then
    if( (y1-b).gt.0) then
       thelta1=thelta1
    else
       thelta1=2.0*pi+thelta1
    end if
  else
    thelta1=thelta1+pi
  end if

  thelta2=atan((y2-b)/(x2-a))
  if((x2-a).gt.0) then
     if( (y2-b).gt.0) then
         thelta2=thelta2
     else
          thelta2=2.0*pi+thelta2
     end if
  else
    thelta2=thelta2+pi
  end if

!  if(myid==0)write(*,*)'x1=',x1,a+r2*cos(thelta1),'y1=',y1,r2*sin(thelta1)+b
!  if(myid==0)write(*,*)'x2=',x2,a+r2*cos(thelta2),'y2=',y2,r2*sin(thelta2)+b
!  if(myid==0)write(*,*)'thelta1=',thelta1,'thelta2=',thelta2
  y3=b+r2*sin(thelta1+1.0/180.0*pi)
!  if(abs(y3).gt.abs(y1)) then
  if(y3.lt.y1) then
   if(thelta2.lt.thelta1) thelta2=thelta2+2.0*pi

    do i=1,nop
       thelta=(thelta2-thelta1)/real(nop)*(i-1)+thelta1
       x3=a+r2*cos(thelta)
       y3=b+r2*sin(thelta)
       nn_con=nn_con+1
       x_con(1,nn_con)=x3
       x_con(2,nn_con)=y3
!       if(myid==0)write(*,*)'x3=',x3,'y3=',y3,'thelta=',thelta
    end do
  else
    if(thelta2.gt.thelta1) thelta2=thelta2-2.0*pi

    do i=1,nop
       thelta=(thelta2-thelta1)/real(nop)*(i-1)+thelta1
       x3=a+r2*cos(thelta)
       y3=b+r2*sin(thelta)
       nn_con=nn_con+1
       x_con(1,nn_con)=x3
       x_con(2,nn_con)=y3
!       if(myid==0)write(*,*)'x3=',x3,'y3=',y3,'thelta=',thelta
    end do
  end if

100 continue     
end subroutine solve_contact










