!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module: adaptive_meshing
!  
!  Axel Gerstenberger, October 2003
!

module adaptive_meshing
  use fluid_variables, only: nsd
  implicit none
  save

  integer,parameter :: adapt_mesh = 1

  integer,parameter :: bgmesh_filename_length = 15
  character(len = bgmesh_filename_length),parameter :: bgmesh_filename = "ifem_bgmesh.msh"

  integer :: nn_bg,ne_bg,ne_surf_bg
  integer,parameter :: nen_surf_bg = 3

  integer,parameter :: nen_bg = 4

contains

subroutine write_background_mesh(nn,x,ne,nen,ien,ne_surf,nen_surf,ien_surf,d)
  use global_constants
  use delta_nonuniform
  use meshgen_gmsh
  use output_gmsh
  implicit none

  integer,intent(in) :: nn
  real(8),intent(in) :: x(1:nsd,1:nn)
  integer,intent(in) :: ne,nen
  integer,intent(in) :: ien(1:nen,1:ne)
  integer,intent(in) :: ne_surf,nen_surf
  integer,intent(in) :: ien_surf(1:nen_surf,1:ne_surf)
  real(8),intent(in) :: d(1:3,1:nn)

  real(8) :: x_bg(nsd,nn_bg)
  integer :: ien_bg(nen_bg,ne_bg)
  integer :: ien_surf_bg(nen_surf_bg,ne_surf_bg)
  integer :: bid_bg(nn_bg)
  real(8) :: dvolume(nn)

  real(8) :: shrknode(maxconn,nn_bg)             !...shape function for each node, contains the weights
  integer :: cnn(maxconn,nn_bg),ncnn(nn_bg)  !...connectivity arrays for domain of incluence for each solid node

  real(8) :: scal(1:nn)          !...element size calculated on original mesh
  real(8) :: scal_bg(1:nn_bg)    !...element size projected on background mesh
  real(8) :: dgrad(1:nn)

  integer,parameter :: one = 1,file = 21
  real(8) :: wave_length

  !...read backgroundmesh from input file
  call read_gmsh(x_bg,ien_bg,ien_surf_bg,bid_bg,nn_bg,ne_bg,ne_surf_bg,nen_bg,nen_surf_bg,bgmesh_filename_length,bgmesh_filename)

  !...generate error information
  !wave_length = maxval( x_orig(1,:) ) - minval( x_orig(1,:) )
  !scal_orig(1,:) = 0.25 + 0.15*sin(4*pi*x_orig(1,:)/wave_length)*sin(6*pi*x_orig(3,:)/wave_length)
  !scal_orig(1,:) = 1.0d0
  scal(:) = sqrt(d(1,:)**2 + d(2,:)**2 + d(3,:)**2)

  call velocity_gradient(d(1:3,1:nn),x,dgrad,ien)

  call calculate_element_size(nn,dgrad)



  !...initialize RKPM shape functions for transfer between calculation mesh and background mesh
  !call delta_initialize(nn_bg,x_bg,nn,ne,x,ien,dvolume,cnn,ncnn,shrknode)
  !write(*,*) "Initialize mapping between original and background mesh"

  !...transfer mesh density to backgroundmesh
  !call delta_exchange(one,scal_bg,nn_bg,dgrad,nn,ndelta,dvolume,cnn,ncnn,shrknode,delta_exchange_fluid_to_solid)
  !write(*,*) "Transfer data from original to background mesh"




  call write_2d_scalar(nn,ne_surf,nen_surf,nsd,ien_surf,x,dgrad,18,"Velocity Gradient",28,"output_velocity_gradient.pos")

  call write_3d_scalar(nn,ne,nen,nsd,ien,x,dgrad,21,"Characteristic Length",22,"output_bgmesh_orig.pos")
  !call write_3d_scalar(nn_bg,ne_bg,nen_bg,nsd,ien_bg,x_bg,dgrad,21,"Characteristic Length",20,"output_bgmesh_bg.pos")


end subroutine write_background_mesh


subroutine velocity_gradient(dloc,xloc,dgradloc,ien)
  use fluid_variables, only: nn,ne,nen,nsd,nsdpad,nenpad,nquad,sq
  implicit none

  real(8),intent(in)  :: dloc(1:nsd,1:nn)
  real(8),intent(in)  :: xloc(1:nsd,1:nn)
  real(8),intent(out) :: dgradloc(1:nn)
  integer,intent(in)  :: ien(1:nen,1:ne)

  real(8) :: x(nsdpad,nenpad)

  real(8) :: vel_mag(nn)
  real(8) :: sh(0:nsdpad,nenpad)
  real(8) :: xr(nsdpad,nsdpad), cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)
  real(8) :: vel_grad_ele(nsd)
  real(8) :: det,vel_grad_mag
  integer :: inl, ine, isd, iq, ntem, in, inen


  !...calculate velocity magnitude
  vel_mag(:) = sqrt(dloc(1,:)**2 + dloc(2,:)**2 + dloc(3,:)**2)

  !...calculate d(Na)/d(x)

  dgradloc(:) = 0.0d0

  do ine=1,ne

     do inl=1,nen
        do isd=1,nsd
           x(isd,inl)=xloc(isd,ien(inl,ine))
        enddo
     enddo

     do iq=1,nquad
        
        if (nen.eq.4) then
           include "sh3d4n.h"
        else if (nen.eq.8) then
           include "sh3d8n.h"
        end if

        !...calculate the first derivative, d(v^abs)/d(x)
        vel_grad_ele(:) = 0.0d0
        do inl=1,nen
           vel_grad_ele(1) = vel_grad_ele(1) + sh(1,inl)*vel_mag(ien(inl,ine))
           vel_grad_ele(2) = vel_grad_ele(2) + sh(2,inl)*vel_mag(ien(inl,ine))
           vel_grad_ele(3) = vel_grad_ele(3) + sh(3,inl)*vel_mag(ien(inl,ine))  
        enddo

        vel_grad_mag = sqrt(vel_grad_ele(1)**2 + vel_grad_ele(2)**2 + vel_grad_ele(3)**2)

        do inl=1,nen
           dgradloc(ien(inl,ine)) = dgradloc(ien(inl,ine)) + sh(0,inl)*vel_grad_mag
        enddo

     enddo
  enddo


  do in=1,nn

     ntem=0

     do ine = 1,ne
        do inen=1,nen
           if (ien(inen,ine) == in) then
              ntem = ntem + 1
              goto 541
           endif
        enddo
 541 enddo

     dgradloc(in) = dgradloc(in)/ntem

  enddo

end subroutine velocity_gradient


subroutine calculate_element_size(nn,dgrad)
  implicit none

  integer,intent(in) :: nn
  real(8),intent(inout) :: dgrad(nn)

  integer :: nperc
  real(8) :: mingrad,maxgrad,range
  real(8) :: dgradtest(nn)
  integer :: i,in

  real(8),parameter :: esizemin = 0.17, esizemax = 0.55


  dgrad(:) = 1/dgrad(:)  !...large values become small element lengths
  dgrad(:) = log(1 + dgrad(:))






 !...95 percentile

  nperc = nint(0.005*nn)
  write(*,*) nperc

 !...cut maximum peaks
  maxgrad = maxval(dgrad(:))
  mingrad = minval(dgrad(:))
  write(*,*) mingrad,maxgrad
  dgradtest(:) = dgrad(:)
  do i = 1,nperc
     write(*,*) maxloc(dgradtest(:))
     dgradtest(maxloc(dgradtest(:))) = mingrad - 100.0d0
  enddo
  maxgrad = maxval(dgradtest(:))
  where(dgradtest(:) < mingrad)
     dgrad(:) = maxgrad
  endwhere


  maxgrad = maxval(dgrad(:))
  mingrad = minval(dgrad(:))
  write(*,*) mingrad,maxgrad
  dgradtest(:) = dgrad(:)
  do i = 1,nperc
     write(*,*) minloc(dgradtest(:))
     dgradtest(minloc(dgradtest(:))) = maxgrad + 100.0d0
  enddo
  mingrad = minval(dgradtest(:))
  where(dgradtest(:) > maxgrad)
     dgrad(:) = mingrad
  endwhere

  maxgrad = maxval(dgrad(:))
  mingrad = minval(dgrad(:))
  write(*,*) mingrad,maxgrad



  range = maxgrad - mingrad

  do in = 1,nn
     dgrad(in) = esizemin + (esizemax - esizemin) * (dgrad(in) - mingrad) / range
  enddo




  !where(x(3,:) < -1.5)
  !   dgrad(:) = 0.3
  !endwhere
  !where(x(3,:) > 1.5)
  !   dgrad(:) = 0.3
  !endwhere
  !where(x(1,:) > 7.0)
  !   dgrad(:) = 0.3
  !endwhere

  !where(scal_bg(1,:) <= 0.2)
  !   scal_bg(1,:) = 0.4
  !elsewhere
  !   scal_bg(1,:) = 0.1
  !endwhere

end subroutine calculate_element_size


end module adaptive_meshing