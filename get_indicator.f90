!==============================================
!get indicator considering weight summation
!===============================================

subroutine get_indicator(xp,xg,hg,infdomain,I_fluid_center,x_center,corr_Ip,I_c)
			 !Sp_sum_g,Sp_sum_0d)

  use interface_variables
  use fluid_variables, only:nn,ne,nsd
!  use mpi_variables

  real(8) xp(nsd,maxmatrix), xg(nsd,nn)    ! coor of interface points and fluid nodes
  real(8) Ip(maxmatrix),Ig(nn)             !indicator of interface points and fluid nodes
  integer infdomain(maxmatrix)             ! the elements that each inter points belong to
  real(8) hg(ne)                           !spacing
  real(8) I_fluid_center(ne)
  real(8) x_center(nsd,ne)
  real(8) corr_Ip(maxmatrix)
  real(8) I_c
  real(8) Sp_sum_g(2,nn), Sp_sum_0d(2,maxmatrix) !weight summation

  integer i,j
  real(8) dx(nsd),hs,Sp,length

  real(8) A(nn_inter,nn_inter)
  real(8) B(nn_inter)   !A*x=B
  integer IPIV(nn_inter), INFO
  Ip(1:nn_inter)=0.0
  Ig(1:nn)=0.0
  B(:)=0.0
  I_c=0.5


  do i=1,nn_inter
     Sp_sum_0d(1:2,i)=0.0
     hs=hg(infdomain(i))
     do j=1,nn_inter
	dx(:)=abs(xp(:,i)-xp(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	A(i,j)=Sp
	Sp_sum_0d(2,i)=Sp_sum_0d(2,i)+Sp
     end do   !interpolate from interfacial points

     do j=1,ne
	dx(:)=abs(xp(:,i)-x_center(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	Ip(i)=Ip(i)+I_fluid_center(j)*Sp
	Sp_sum_0d(1,i)=Sp_sum_0d(1,i)+Sp
     end do   !interpolate from center points

!Sp_sum_0d is the summation of the weights for both the inter and center points
     
!      A(i,:)=A(i,:)/Sp_sum_0d(2,i)
     Ip(i)=Ip(i)/Sp_sum_0d(1,i)   ! consider weighting summation

  end do

  B(1:nn_inter)=I_c-Ip(1:nn_inter)
  
  call DGESV(nn_inter,1,A,nn_inter,IPIV,B,nn_inter,INFO)

  corr_Ip(1:nn_inter)=B(1:nn_inter)  !calculate delta Ip
!================================================================
!  do i=1,nn
!     Sp_sum_g(1:2,i)=0.0
!     do j=1,nn_inter
!	hs=hg(infdomain(j))
!	dx(:)=abs(xg(:,i)-xp(:,j))
!	call B_Spline(dx,hs,nsd,Sp)
!	Sp_sum_g(2,i)=Sp_sum_g(2,i)+Sp
!        Ig(i)=Ig(i)+corr_Ip(j)*Sp
!     end do !interpolate from the inter points
!	Ig(i)=Ig(i)/Sp_sum_g(2,i)    
!     if(Sp_sum_g(2,i).lt.1.0e-3) then
!	Ig(i)=0.0
!     end if
!     do j=1,ne
!	hs=hg(j)
!	dx(:)=abs(xg(:,i)-x_center(:,j))
!	call B_Spline(dx,hs,nsd,Sp)
!	Sp_sum_g(1,i)=Sp_sum_g(1,i)+Sp
!!	Ig(i)=Ig(i)+I_fluid_center(j)*Sp
!     end do !interpolate from the center points
!     do j=1,ne
!	hs=hg(j)
!	dx(:)=abs(xg(:,i)-x_center(:,j))
!	call B_Spline(dx,hs,nsd,Sp)
!	Ig(i)=Ig(i)+I_fluid_center(j)*Sp/Sp_sum_g(1,i)
!    end do
!!     Ig(i)=Ig(i)/Sp_sum_g(i) !devided by weight summation
!
!  end do

end subroutine get_indicator













    
  














        














