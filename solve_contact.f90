

subroutine solve_contact(x_fix,x_wall,norm_fix,norm_wall,x_con,nn_con)

  use fluid_variables,only:nsd
  use interface_variables, only:maxmatrix
  use mpi_variables
  real(8) x_fix(nsd),x_wall(nsd),norm_fix(nsd),norm_wall(nsd)
  real(8) x_con(nsd,maxmatrix)
  integer nn_con

  real(8) x1,y1,nx,ny,x2,y2,mx,my,a,b
  real(8) thelta2, thelta1
  real(8) C1,C2,C3,C4,C5,C6,r2,root1,root2

  real(8) x3,y3
  integer i,j
  real(8) pi,thelta,eps
  pi=3.14159
  eps=0.0001
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
if(myid==0)write(*,*)'curv for contact=',1.0/r2

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
!if(myid==0)write(*,*)'t1=',thelta1/pi*180,'t2=',thelta2/pi*180,'b=',b
  y3=b+r2*sin(thelta1+1.0/180.0*pi)
!if(myid==0)write(*,*)'y3=',y3,'y2=',y2
  if(abs(y3).gt.abs(y1)) then
    if(thelta2.lt.thelta1) thelta2=thelta2+2.0*pi

    do i=1,30
       thelta=(thelta2-thelta1)/30.0*(i-1)+thelta1
       x3=a+r2*cos(thelta)
       y3=b+r2*sin(thelta)
       nn_con=nn_con+1
       x_con(1,nn_con)=x3
       x_con(2,nn_con)=y3
    end do
  else
    if(thelta2.gt.thelta1) thelta2=thelta2-2.0*pi

    do i=1,30
       thelta=(thelta2-thelta1)/30.0*(i-1)+thelta1
       x3=a+r2*cos(thelta)
       y3=b+r2*sin(thelta)
       nn_con=nn_con+1
       x_con(1,nn_con)=x3
       x_con(2,nn_con)=y3
    end do
  end if

!  do i=1,51
!     x3=(x2-x1)/50.0*(i-1)+x1
!     root1=b+sqrt(r2-(x3-a)**2)
!     root2=b-sqrt(r2-(x3-a)**2)
!     if((root1.lt.y1).and.(root1.gt.-0.5)) then
!       nn_con=nn_con+1
!       x_con(1,nn_con)=x3
!       x_con(2,nn_con)=root1
!     else if((root2.lt.y1).and.(root2.gt.-0.5)) then
!       nn_con=nn_con+1
!       x_con(1,nn_con)=x3
!       x_con(2,nn_con)=root2
!     end if
!  end do
100 continue     
end subroutine solve_contact










