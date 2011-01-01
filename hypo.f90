!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   *.fi files are used to shorten hypo.f (keeping the overview)
!   the include command reads these files and replaces the include line
!   with the content of these files

subroutine hypo
  use global_simulation_parameter
  use global_constants
  use run_variables
  use delta_nonuniform
  use solid_variables
  use fluid_variables
  use r_common, only: ninit
  use meshgen_fluid
  use meshgen_solid
  use form
  use ensight_output
  use ale_variables
  use mpi_variables ! call mpi variable module
  implicit none
  include 'mpif.h'
!==============================	  
! Definition of variables
  integer :: klok,j
  integer ie, inen
!  integer tmp_index(ne)
!============================
! Variables for boudary equations
  integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
  real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration 
  real(8) time
  real(8) time1_begin
  real(8) time1_end
!==========================================
! Variables define for ALE feature
integer node_alebc(nn_alebc)
real(8) kp(nsd,nn)
real(8) kw(nsd,nn)
real(8) e_ale(nsd,nn)
real(8) disp_ale(nsd,nn)
real(8) kdg(nsd,nn)
!========================================
include "hypo_declaration_ale.fi"
!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
!===================================
! Define varibales on each processor
  include "hypo_declaration_part.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information

  include "hypo_restart_file_check.fi"
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"
!===================================
! Prepare for MPI
  include "hypo_prepare_part.fi"
  include "hypo_prepare_com_node.fi"
!=============================
! define the influence domain matrix
 ! integer infdomain(nn_solid)
      call mpi_barrier(mpi_comm_world,ierror)
  write(*,*) 'myid', myid, 'nn_local', nn_local, 'ne_local', ne_local !id for debuger
!=============================
! call the subroutine to set up ginflow and linflow


if (edge_inflow .ne. 0) then
call edgeele(edge_inflow,rng,neface,ne,bc4el,ne_inflow)
end if  
  
  if (restart == 0) then
	if (myid == 0) then
     	include 'hypo_write_output.fi'
	end if
  else
     include "hypo_restart_read.fi"
  endif

!================================
	if (myid == 0) write(*,*) 'apply initial jac'
jac(:,:)=1.0d0
refvel(:,:,:)=0.0d0
call read_alebc(node_alebc,nn_alebc) ! read in node index on the inner boundary
call formid_ale(kid,rng,ien) ! fix outer boundary
disp_ale(:,:)=0.0d0


!=================================================================
!						 time loop	
!=================================================================
  time_loop: do its = nts_start,nts !.....count from 1 or restart-timestep to number of timesteps
      call mpi_barrier(mpi_comm_world,ierror)
	if (myid == 0) then
     	write (6,*) ' '
     	write (6,*) 'TIME STEP = ', its
     	write (6,*) ' '
     	write (7,*) ' '
     	write (7,*) 'TIME STEP = ', its
     	write (7,*) ' '

!=================================================================
! Write restart information in binary file

     	include "hypo_restart_write.fi"
	end if
     tt = tt + dt    !....update real time
     klok = klok + 1 !....update counter for output
	if (myid == 0) then
     	write (6,'("  physical time = ",f7.3," s")') tt
     	write (7,'("  physical time = ",f7.3," s")') tt
!=================================================================
! Adding the ALE feature for testing purpose
	write(*,*) 'begin ALE'
	end if

time=mpi_wtime()

!----------------------------------------------------
xold(:,:)=x(:,:) ! store previous step moving mesh
dold(:,:)=d(:,:) ! store previous step velcotiy and pressure
refvelo(:,:,:)=refvel(:,:,:) ! store previous step referential velocity
jaco(:,:)=jac(:,:) ! store previous step determinant of deformation 
!disp_ale(:,:)=0.0d0
!----------------------------------------------------
kinner=inner ! using same # of inner space to solve both 2 gmres
kouter=outer ! using same # of outer space to solve both 2 gmres
!kinner=100 ! Need inputs
!kouter=outer ! Need inputs

!call read_alebc(node_alebc,nn_alebc) ! read in node index on the inner boundary
!call formid_ale(kid,rng,ien) ! fix outer boundary


!call innerbc_ale(disp_ale,kid,node_alebc,nn_alebc,nsd,nn) ! apply displacement on inner boundary
!call innerbc_ale_vocaltest(disp_ale,kid,node_alebc,nn_alebc,nsd,nn,xold) ! Apply boundary movement for test vocal case
!call innerbc_ale_updown(disp_ale,kid,node_alebc,nn_alebc,nsd,nn,xold) ! Apply boundary movement for rigid vocal folds 
! move up and down
!call innerbc_ale_updown_larger(disp_ale,kid,node_alebc,nn_alebc,nsd,nn,xold) ! Apply boundary movement for rigid vocal folds move up and down using the geometry similar as Krane's 2006 paper
call innerbc_ale_re(disp_ale,kid,node_alebc,nn_alebc,nsd,nn,xold)


!do iit=1,2
kp(:,:)=0.0d0
kw(:,:)=0.0d0 ! clear the matrix
e_ale(:,:)=0.0d0 ! no extern force for this case

time1_begin=mpi_wtime()

call blockm(xold,e_ale,disp_ale,kw,kp,ien,jac,ne_intlocal,ien_intlocal)
time1_end=mpi_wtime()
if (myid == 0) write(*,*) 'Time evaluate block', time1_end-time1_begin

!=====================================================================
! Communicate residual and diagonal preconditioner
!        call communicate_res(global_com,nn_global_com,local_com,nn_local_com,kp,nsd,nn)
!        call communicate_res(global_com,nn_global_com,local_com,nn_local_com,kw,nsd,nn)
      call mpi_barrier(mpi_comm_world,ierror)
        call communicate_res_ad(kp,nsd,nn,send_address,ad_length)
      call mpi_barrier(mpi_comm_world,ierror)
        call communicate_res_ad(kw,nsd,nn,send_address,ad_length)
if (myid == 0) write(*,*) 'First communication in hypo'
!=======================================================================
!call setid(kp,kid,nsd)
call setid_pa(kp,nsd,nn,kid,node_local,nn_local)

call getnorm_pa(kp,nsd,nn,node_local,nn_local,res_l)
res_l=sqrt(res_l)
kdg(:,:)=0.0d0 ! clear matrix
if (myid ==0) write(*,*) 'res_l', res_l, 'land-over-mu', landa_over_mu

 call gmres_new(xold,kw,kp,kdg,ien,kid,jac, &
                        ne_intlocal,ien_intlocal,node_local,nn_local, &
                        global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length)
! call gmresm(xref,kid,kw,kp,kdg,ien,kz,kv,kzg,kavg,ksm,kavloc,kh,ky,kcc,kss)
call getnorm_pa(kdg,nsd,nn,node_local,nn_local,del_l)
 del_l = sqrt(del_l)
 call mpi_barrier(mpi_comm_world,ierror)
 call update(kp,disp_ale,kdg,nsd)

time=mpi_wtime()-time
if (myid == 0) write(*,*) 'Time for Mesh updating equation', time

!write(*,*) 'ALE iteration', 'ALE error norm', del_l
!end do

 
if (myid == 0) write(*,*) '**********   Mesh Updated   **********'

time=mpi_wtime()

call add(x,xold,disp_ale,nsd)
!call defgrad(x,xref,ien,f,jac,finv,jacinv) ! old deformation gradient subroutine only for 3-D
call defgrad_new(x,xref,ien,f,jac,finv,jacinv) ! works for both 2-D and 3-D

!call defgrad_new(x,xref,ien,f_new,jac_new,finv_new,jacinv_new)
!write(*,*) 'diff between old and new in f', maxval(abs(jac_new(:,:)*jacinv_new(:,:)))
call velocity(x,xold,meshvel,refvel,finv,d,ien)
!??????????????????
!??? do not know what 'af' is 
!??????????????????

call formid(id,rng,ien)  ! apply essential B.C 
call formd(d,rng,ien)
! Not sure about the moving inner ALE boundary
! Moving inner bloundary condition
call moving_bc(d,id,node_alebc,nn_alebc,meshvel)

! ??????????????????????????????????????????
f_fluids(:,:)=0.0d0
!===========================================================================
if (edge_inflow .ne. 0) then
! Apply nature boundary condition
if (myid == 0) write (*,*) 'Apply nature BC'
        call nature_pre(x,d,ien,rng,bc4el,ne_inflow,edge_inflow,pin,res_bc)
        f_fluids(1,:) =  res_bc(1,:)
end if

time=mpi_wtime()-time
if (myid == 0) write(*,*) 'Time for Updating mesh_vel and JAC', time

time=mpi_wtime()

iit=1
del_l=1.0d0

do 100, while((iit .le. nit) .and. (del_l .ge. 10.0e-6 ))    
p(:,:)=0.0d0
w(:,:)=0.0d0
!call block_ale(xref,d,dold,af,p,w,hg,ien,f,finv,jac,jaco,refvel,refvelo) ! get residuals
call block_ale(x,d,dold,p,w,hg,ien,f_fluids,rng,f_stress,meshvel,ne_intlocal,ien_intlocal,node_local,nn_local)

!call communicate_res(global_com,nn_global_com,local_com,nn_local_com,p,ndf,nn)
!call communicate_res(global_com,nn_global_com,local_com,nn_local_com,w,ndf,nn)
      call mpi_barrier(mpi_comm_world,ierror)
        call communicate_res_ad(p,ndf,nn,send_address,ad_length)
      call mpi_barrier(mpi_comm_world,ierror)
        call communicate_res_ad(w,ndf,nn,send_address,ad_length)
if(myid==0) write(*,*) 'fluid solver first communication'
        call setid_pa(p,ndf,nn,id,node_local,nn_local)


!call setid(p,id,ndf) ! set residual be zero at essential B.C DOF
!write(*,*) 'p', p(:,:)
call getnorm_pa(p,ndf,nn,node_local,nn_local,res_l)
res_l =sqrt(res_l)

if(myid==0) write(*,*) 'res_l', res_l

dg(:,:)=0.0d0

!call gmres_ale(xref,d,dold,id,af,w,p,dg,hg,ien,z,v,zg,avg,sm, &
!	   vloc,avloc,h_gmres,y_gmres,cc,ss,finv,jac,jaco,refvel,refvelo)

call gmres_ale(x,d,dold,w,p,dg,hg,ien,f_fluids,id,meshvel, &
		ne_intlocal,ien_intlocal,node_local,nn_local, &
                global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length)

call getnorm_pa(dg,ndf,nn,node_local,nn_local,del_l)
del_l = sqrt(del_l)
call update(p,d,dg,ndf)

if (myid == 0) write(*,10)  iit, res_l, del_l
10 format('Block RES and Nowton Iteration Error', I3, e15.7, e15.7)

iit=iit+1
100 continue ! iteration loop

time=mpi_wtime()-time
if (myid == 0) write(*,*) 'Time for Solving N-S', time

time=mpi_wtime()
if (myid == 0) then
     include "hypo_write_output.fi"
end if
time=mpi_wtime()-time
if (myid == 0) write(*,*) 'Time for Outputting data', time

  enddo time_loop



end subroutine hypo
