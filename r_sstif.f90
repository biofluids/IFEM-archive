!
!     element tangent stiffness matrix assemblage
!
subroutine r_sstif(ocpp,ocuu,ocup,ne,w,toxj,lx,ly,lz)
  use r_common
  use solid_variables
  implicit none

  real*8 :: ocpp
  real*8 :: ocuu(6,6),ocup(6)
  integer,intent(in) :: ne,lx,ly,lz
  real*8 :: w
  real*8 :: toxj(3,3)
  
  real*8 :: sdensit
  integer :: ni,i,j,k,nu1,nv1,nw1,mu1,nk,ntem
  real*8 :: xac(3),xve(3)
  real*8 :: fu,fkup,fp,fkpp,totalh



  sdensit=density_solid

!ccccccccccc
!     i-u
!ccccccccccc

  do ni=1,nis
     do i=1,3
        nu1=(i-1)*nn_solid+solid_fem_con(ne,ni)
        mu1=(i-1)*nis+ni
!ccccccccccc
!     fu
!ccccccccccc
        call r_scalfu(fu,i,ni)
        !call r_scalfu_curr(fu,i,ni,ne,lx,ly,lz)
        predrf(nu1)=predrf(nu1)-fu*w
        do nk=1,nump
           call r_scalkup(fkup,ocup,i,nk,ni)
           xkup(mu1,nk,ne)=xkup(mu1,nk,ne)+fkup*w
        enddo
     enddo
  enddo

  do i=1,nump
     call r_scalfp(fp,ocpp,i)
     xfp(i,ne)=xfp(i,ne)+fp*w
     do j=1,nump
        call r_scalkpp(fkpp,ocpp,i,j)
        xkpp(i,j,ne)=xkpp(i,j,ne)+fkpp*w
     enddo
  enddo
!cccccccccccccccccccccccccc
!     inertia forces
!cccccccccccccccccccccccccc
  do i=1,3
     xac(i)=0.0d0
     xve(i)=0.0d0
     do k=1,nis
        ntem=solid_fem_con(ne,k)
        xac(i)=xac(i)+h(k)*solid_accel(i,ntem)
        xve(i)=xve(i)+h(k)*solid_vel(i,ntem)
     enddo
  enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  xac(1)=xac(1)+xmg(1)
  xac(2)=xac(2)+xmg(2)
  xac(3)=xac(3)+xmg(3)
  totalh=0
  do ni=1,nis
     nu1=             solid_fem_con(ne,ni)
     nv1=  nn_solid + solid_fem_con(ne,ni)
     nw1=2*nn_solid + solid_fem_con(ne,ni)

     predrf(nu1) = predrf(nu1) - w*sdensit*h(ni)*xac(1) - w*xviss*h(ni)*xve(1)
     predrf(nv1) = predrf(nv1) - w*sdensit*h(ni)*xac(2) - w*xviss*h(ni)*xve(2)
     predrf(nw1) = predrf(nw1) - w*sdensit*h(ni)*xac(3) - w*xviss*h(ni)*xve(3)
!	totalh=totalh+h(ni)

  enddo

  return
end subroutine r_sstif



