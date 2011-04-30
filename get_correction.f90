!=====================================
! get correction term
!=====================================

subroutine get_correction(x_inter,x_center,hg,infdomain_inter,corr_Ip,I_fluid_center)

  use interface_variables
  use fluid_variables, only:nn,ne,nsd
  use allocate_variables, only:center_domain,nn_center_domain

  real(8) x_inter(nsd,maxmatrix)
  real(8) x_center(nsd,ne)
  real(8) hg(ne)
  integer infdomain_inter(maxmatrix)
  real(8) corr_Ip(maxmatrix)
  real(8) I_fluid_center(ne)
!  real(8) x(nsd,nn),I_fluid(nn)

  integer i,j,ie
  real(8) dx(nsd),hs,Sp
  real(8) A(nn_inter,nn_inter)
  real(8) B(nn_inter)  !Ax=B
  integer IPIV(nn_inter), INFO

  A(:,:)=0.0
  B(:)=0.0

  do i=1,nn_inter
     hs=hg(infdomain_inter(i))
     do j=1,nn_inter
	dx(:)=abs(x_inter(:,i)-x_inter(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	A(i,j)=Sp
     end do

     do j=1,nn_center_domain
	ie=center_domain(j)
	dx(:)=abs(x_inter(:,i)-x_center(:,ie))
	call B_Spline(dx,hs,nsd,Sp)
	B(i)=B(i)+I_fluid_center(ie)*Sp
     end do
  end do

     B(:)=0.5-B(:)

     call DGESV(nn_inter,1,A,nn_inter,IPIV,B,nn_inter,INFO)

     corr_Ip(1:nn_inter)=B(1:nn_inter) !calculate delta Ip
!  I_fluid(1:nn)=0.0
!  do i=1,nn
!     do j=1,nn_inter
!        dx(:)=abs(x(:,i)-x_inter(:,j))
!        call B_Spline(dx,hs,nsd,Sp)
!        I_fluid(i)=I_fluid(i)+corr_Ip(j)*Sp
!     end do
!     do j=1,ne
!        dx(:)=abs(x(:,i)-x_center(:,j))
!	call B_Spline(dx,hs,nsd,Sp)
!	I_fluid(i)=I_fluid(i)+I_fluid_center(j)*Sp
!     end do
!  end do

end subroutine get_correction
