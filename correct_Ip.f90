!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!get the correction of Ip for interfacial points!!!!
!!!!get the correction of Ig for fluid nodes!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine correct_Ip(xp,xg,Ip,Ig,hg,I_c,infdomain_cI,corr_Ip)

  use interface_variables
  use fluid_variables, only:nsd,ne,nn

  real(8) xp(nsd,maxmatrix)       !coordinates of interfacial points
  real(8) xg(nsd,nn)              !coordinates of fluid nodes
  real(8) Ip(maxmatrix)           !indicator of interfacial points
  real(8) Ig(nn)                  !indicator of fluid nodes
  real(8) Ig_n
  real(8) hg(ne)
  integer infdomain_cI(maxmatrix)
  real(8) corr_Ip(maxmatrix)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8) I_c                    
!!!constant interface indicator.Initially equals to 0.0.
!!!Should be only calculated in the first time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer i,j
  real(8) A(nn_inter,nn_inter)
  real(8) dIp(nn_inter)
  real(8) B(nn_inter)     ! A*dIp=B
  real(8) dx(nsd)
  real(8) hs
  real(8) Sp
  integer IPIV(nn_inter)
  integer INFO
  A(:,:) = 0.0
  if (I_c.lt.1.0e-6) then
     do i=1,nn_inter
	I_c=I_c+Ip(i)
     end do
     I_c=I_c/nn_inter
  end if 
! if Ic=0.0, then calculate the initial Ic
! Ic is a constant over time
     
  do i=1,nn_inter
     hs=hg(infdomain_cI(i))
     do j=1,nn_inter
	dx(:)=abs(xp(:,i)-xp(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	A(i,j) = Sp
     end do
  end do ! construct A matrix

  do i=1,nn_inter
     B(i)=I_c-Ip(i)
  end do ! construct B matrix

  call DGESV(nn_inter,1,A,nn_inter,IPIV,B,nn_inter,INFO)

  corr_Ip(1:nn_inter) = B(1:nn_inter)

  do i=1,nn
     Ig_n = 0.0
     do j=1,nn_inter
	hs=hg(infdomain_cI(j))
	dx(:)=abs(xg(:,i)-xp(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	Ig_n=Ig_n+B(j)*Sp
     end do
     Ig(i)=Ig(i)+Ig_n
  end do

write(*,*)'corr_Ip=',corr_Ip(1:nn_inter)
stop

end subroutine correct_Ip

















