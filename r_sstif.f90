!
!     element tangent stiffness matrix assemblage
!
subroutine r_sstif(ocpp,ocuu,ocup,xkup,xkpp,xfp,ne,w,toxj,vel,acc,solid_fem_con,cstr_element)
  use r_common
  use solid_variables
  implicit none

  integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con    !...connectivity for solid FEM mesh
  real(8) :: cstr_element(6)  !...Cauchy stress

  real(8) :: ocpp
  real(8) :: ocuu(6,6),ocup(6)
  real(8) :: xkup(3*nen_solid,nup,ne_solid)
  real(8) :: xkpp(nup,nup,ne_solid)
  real(8) :: xfp(nup,ne_solid)
  integer,intent(in) :: ne
  real(8) :: vel(nsd_solid,nen_solid)  !...element nodes velocity
  real(8) :: acc(nsd_solid,nen_solid)  !...element nodes acceleration
  real(8) :: w
  real(8) :: toxj(3,3)
  
  real(8) :: sdensit
  integer :: ni,isd,ip,jp,k,nu1,nv1,nw1,mu1,nk
  real(8) :: xac(1:nsd_solid),xve(1:nsd_solid)
  real(8) :: fu,fkup,fp,fkpp,totalh



  sdensit=density_solid

!ccccccccccc
!     i-u
!ccccccccccc

  do ni=1,nen_solid
     do isd=1,3
        nu1 = (isd-1)*nn_solid+solid_fem_con(ne,ni)
        mu1 = (isd-1)*nen_solid+ni
!ccccccccccc
!     fu
!ccccccccccc
        call r_scalfu(fu,isd,ni)
        !call r_scalfu_curr(fu,isd,ni,cstr_element)
        predrf(nu1) = predrf(nu1) + fu*w
        do nk=1,nump
           call r_scalkup(fkup,ocup,isd,nk,ni)
           xkup(mu1,nk,ne) = xkup(mu1,nk,ne) + fkup*w
        enddo
     enddo
  enddo

  do ip=1,nump
     call r_scalfp(fp,ocpp,ip)
     xfp(ip,ne) = xfp(ip,ne) + fp*w
     do jp=1,nump
        call r_scalkpp(fkpp,ocpp,ip,jp)
        xkpp(ip,jp,ne) = xkpp(ip,jp,ne) + fkpp*w
     enddo
  enddo
!cccccccccccccccccccccccccc
!     inertia forces
!cccccccccccccccccccccccccc
  do isd=1,3
     xac(isd)=0.0d0
     xve(isd)=0.0d0
     do k=1,nen_solid
        xac(isd) = xac(isd)+h(k)*acc(isd,k)
        xve(isd) = xve(isd)+h(k)*vel(isd,k)
     enddo
  enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  xac(1) = xac(1) + xmg(1)
  xac(2) = xac(2) + xmg(2)
  xac(3) = xac(3) + xmg(3)
  totalh=0
  do ni=1,nen_solid
     nu1=             solid_fem_con(ne,ni)
     nv1=  nn_solid + solid_fem_con(ne,ni)
     nw1=2*nn_solid + solid_fem_con(ne,ni)

     predrf(nu1) = predrf(nu1) - w*sdensit*h(ni)*xac(1) - w*xviss*h(ni)*xve(1)
     predrf(nv1) = predrf(nv1) - w*sdensit*h(ni)*xac(2) - w*xviss*h(ni)*xve(2)
     predrf(nw1) = predrf(nw1) - w*sdensit*h(ni)*xac(3) - w*xviss*h(ni)*xve(3)
!   totalh=totalh+h(ni)

  enddo

  return
end subroutine r_sstif



