module adaptive_meshing_2
  implicit none
  save



contains

subroutine exchange_mesh
  use restart_lib
  !use fluid_variables, only: nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface
  implicit none

  integer :: nn_old,ne_old,nen_old,ne_surf_old,nen_surf_old,etype_old,nnface_old,neface_old

 !...check for restart files
  call restart_file_check
 !...read header information from old mesh
  call restart_read_header(nn_old,ne_old,nen_old,ne_surf_old,nen_surf_old,etype_old,nnface_old,neface_old)
 !...exchange data between the two meshes
  call exchange(nn_old,ne_old,nen_old,ne_surf_old,nen_surf_old,etype_old,nnface_old,neface_old)

end subroutine exchange_mesh


subroutine exchange(nn_old,ne_old,nen_old,ne_surf_old,nen_surf_old,etype_old,nnface_old,neface_old)
  use solid_variables, only: nn_solid,ne_solid,nen_solid,nn_solid_1,ne_solid_1,n_solid,nsurface,nsd_solid
  use fluid_variables, only: nsd,ndf,nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface
  use restart_lib
  use delta_nonuniform
  use meshgen_gmsh
  use output_gmsh
  implicit none

  integer,intent(inout) :: nn_old,ne_old,nen_old,ne_surf_old,nen_surf_old,etype_old,nnface_old,neface_old

  integer :: nts_start

  real(8) :: tt
  integer :: klok

  integer :: ien(nen,ne)
  integer :: ien_surf(nen_surf,ne_surf)
  integer :: rng(neface,ne)
  real(8) :: x(nsd,nn),xref(nsd,nn)
  real(8) :: hg(ne)
  real(8) :: d(ndf,nn)
  integer :: bid(nn)

  integer :: ien_old(nen_old,ne_old)
  integer :: ien_surf_old(nen_surf_old,ne_surf_old)
  integer :: rng_old(neface_old,ne_old)
  real(8) :: x_old(nsd,nn_old),xref_old(nsd,nn_old)
  real(8) :: hg_old(ne_old)
  real(8) :: d_old(ndf,nn_old)
  integer :: bid_old(nn_old)

  integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con
  integer,dimension(1:ne_solid,1:nsurface)  :: solid_surface
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel


  real(8) :: shrknode(maxconn,nn)             !...shape function for each node, contains the weights
  integer :: cnn(maxconn,nn),ncnn(nn)  !...connectivity arrays for domain of incluence for each solid node

  real(8) :: dvolume(nn_old)





 !...read old mesh information
  call restart_read_data(nts_start,nn_old,ne_old,nen_old,ne_surf_old,nen_surf_old,etype_old,nnface_old,neface_old,             &
                         nn_solid,ne_solid,nen_solid,nn_solid_1,ne_solid_1,n_solid,nsurface,   &
                         nsd,ndf,nsd_solid,                                                    &
                         tt,klok,                                                              &
                         ien_old,ien_surf_old,rng_old,x_old,xref_old,hg_old,d_old,bid_old,     &
                         solid_fem_con,solid_surface,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel)

 !...read new mesh information
  call read_gmsh(x,ien,ien_surf,bid,nn,ne,ne_surf,nen,nen_surf,15,"input_fluid.msh")
    write(*,*) "read_gmsh done"
 !...prepare remaining variables depending on new mesh information
  xref(:,:) = x(:,:)
  call lenght(x,ien,hg)


 !...initialize RKPM shape functions for transfer between calculation mesh and background mesh
  call delta_initialize(nn,x,nn_old,ne_old,x_old,ien_old,dvolume,cnn,ncnn,shrknode)
  write(*,*) "Initialize mapping between original and background mesh"

 !...transfer unknown vector to new mesh
  call delta_exchange(ndf,d,nn,d_old,nn_old,ndelta,dvolume,cnn,ncnn,shrknode,delta_exchange_fluid_to_solid)
  write(*,*) "Transfer data from original to background mesh"

  call write_2d_scalar(nn_old,ne_surf_old,nen_surf_old,nsd,ien_surf_old,x_old,sqrt(d_old(1,:)**2 + d_old(2,:)**2 + d_old(3,:)**2),18,"Velocity Magnitude",33,"output_velocity_magnitude_old.pos")
  call write_2d_scalar(nn,ne_surf,nen_surf,nsd,ien_surf,x,sqrt(d(1,:)**2 + d(2,:)**2 + d(3,:)**2),18,"Velocity Magnitude",33,"output_velocity_magnitude_new.pos")

  call restart_write_data(nts_start,nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface,             &
                          nn_solid,ne_solid,nen_solid,nn_solid_1,ne_solid_1,n_solid,nsurface,   &
                          nsd,ndf,nsd_solid,                                                    &
                          tt,klok,                                                              &
                          ien,ien_surf,rng,x,xref,hg,d,bid,                                     &
                          solid_fem_con,solid_surface,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel)




end subroutine exchange

end module adaptive_meshing_2