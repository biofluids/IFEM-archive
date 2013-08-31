
subroutine indicator_initialize(I_fluid,I_fluid_center,x,x_center,nn,nn_center,nsd)

  real(8) I_fluid(nn),I_fluid_center(nn_center)
  real(8) x(nsd,nn),x_center(nsd,nn_center)
  integer nn,nn_center,nsd

  integer i,j


  I_fluid_center(:)=0.0
  I_fluid(:)=0.0
!  do j=1,nn_center
!     temp=sqrt((x_center(1,j)-3.5)**2+(x_center(2,j)-1.0)**2)
!     if(temp.le.0.3) I_fluid_center(j)=1.0
!  end do
!  do j=1,nn
!     temp=sqrt((x(1,j)-3.5)**2+(x(2,j)-1.0)**2)
!     if(temp.le.0.3)I_fluid(j)=1.0
!  end do

  do j=1,nn_center
     if(x_center(2,j).gt.1.0)I_fluid_center(j)=1.0
  end do
  do j=1,nn
     if(x(2,j).gt.1.0)I_fluid(j)=1.0
  end do


!  do j=1,nn_center
!     temp=sqrt((x_center(1,j)+1.6)**2+(x_center(3,j)-0.4)**2)
!     if((temp.le.0.1).and.(x_center(1,j).lt.-1.5).and.(x_center(3,j).lt.0.5)) I_fluid_center(j)=1.0
!     if((x_center(1,j).lt.-1.5).and.(x_center(3,j).lt.0.4))I_fluid_center(j)=1.0
!     if((x_center(1,j).lt.-1.6).and.(x_center(3,j).lt.0.5))I_fluid_center(j)=1.0
!  end do
!  do j=1,nn
!     temp=sqrt((x(1,j)+1.6)**2+(x(3,j)-0.4)**2)
!     if((temp.le.0.1).and.(x(1,j).lt.-1.5).and.(x(3,j).lt.0.5))I_fluid(j)=1.0
!     if((x(1,j).lt.-1.5).and.(x(3,j).lt.0.4))I_fluid(j)=1.0
!     if((x(1,j).lt.-1.6).and.(x(3,j).lt.0.5))I_fluid(j)=1.0
!  end do


!  do j=1,nn_center
!     temp=sqrt((x_center(1,j)+0.05)**2+(x_center(2,j)+0.05)**2)
!     if((temp.le.0.05).and.(x_center(1,j).lt.0.0).and.(x_center(2,j).lt.0.0))I_fluid_center(j)=1.0
!     if((x_center(1,j).lt.-0.05).and.(x_center(2,j).lt.0.0))I_fluid_center(j)=1.0
!     if((x_center(1,j).lt.0.0).and.(x_center(2,j).lt.-0.05))I_fluid_center(j)=1.0
!  end do

!  do j=1,nn
!     temp=sqrt((x(1,j)+0.05)**2+(x(2,j)+0.05)**2)
!     if((temp.le.0.05).and.(x(1,j).lt.0.0).and.(x(2,j).lt.0.0))I_fluid(j)=1.0
!     if((x(1,j).lt.-0.05).and.(x(2,j).lt.0.0))I_fluid(j)=1.0
!     if((x(1,j).lt.0.0).and.(x(2,j).lt.-0.05))I_fluid(j)=1.0
!  end do


!  do j=1,nn_center
!     temp=sqrt((x_center(1,j)+1.0)**2+(x_center(2,j)+0.1)**2+(x_center(3,j)-1.5)**2)
!     if(temp.lt.0.4-max_hg) I_fluid_center(j)=1.0
!     if(abs(temp-0.4).le.max_hg) I_fluid_center(j)=-0.5*(temp-0.4)/max_hg+0.5
!  end do 
!  do j=1,nn
!     temp=sqrt((x(1,j)+1.0)**2+(x(2,j)+0.1)**2+(x(3,j)-1.5)**2)
!     if(temp.le.0.4-max_hg)I_fluid(j)=1.0
!     if(abs(temp-0.4).le.max_hg)I_fluid(j)=-0.5*(temp-0.4)/max_hg+0.5
!  end do


end subroutine indicator_initialize



