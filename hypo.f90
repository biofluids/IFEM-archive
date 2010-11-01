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
  implicit none

!==============================	  
! Definition of variables
  integer :: klok,j

  integer infdomain(nn_solid)
  real(8) mass_center(2)

!==========================================
! Variables define for ALE feature
integer node_alebc(80)
integer nn_alebc
real(8) kp(nsd,nn)
real(8) kw(nsd,nn)
real(8) e_ale(nsd,nn)
real(8) disp_ale(nsd,nn)
real(8) kdg(nsd,nn)
include "hypo_declaration_ale.fi"
!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information

  include "hypo_restart_file_check.fi"
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"
!=============================
 ! integer infdomain(nn_solid)
  
  
  if (restart == 0) then
     include 'hypo_write_output.fi'
  else
     include "hypo_restart_read.fi"
  endif

!================================
write(*,*) 'apply initial jac'
jac(:,:)=1.0d0
refvel(:,:,:)=0.0d0



!=================================================================
!						 time loop	
!=================================================================
  time_loop: do its = nts_start,nts !.....count from 1 or restart-timestep to number of timesteps

     write (6,*) ' '
     write (6,*) 'TIME STEP = ', its
     write (6,*) ' '
     write (7,*) ' '
     write (7,*) 'TIME STEP = ', its
     write (7,*) ' '

!=================================================================
! Write restart information in binary file

     include "hypo_restart_write.fi"

     tt = tt + dt    !....update real time
     klok = klok + 1 !....update counter for output

     write (6,'("  physical time = ",f7.3," s")') tt
     write (7,'("  physical time = ",f7.3," s")') tt
!=================================================================
! Adding the ALE feature for testing purpose
write(*,*) 'begin ALE'
!----------------------------------------------------
xold(:,:)=x(:,:) ! store previous step moving mesh
dold(:,:)=d(:,:) ! store previous step velcotiy and pressure
refvelo(:,:,:)=refvel(:,:,:) ! store previous step referential velocity
jaco(:,:)=jac(:,:) ! store previous step determinant of deformation 
!----------------------------------------------------
nn_alebc=80
kinner=inner ! using same # of inner space to solve both 2 gmres
kouter=outer ! using same # of outer space to solve both 2 gmres
call read_alebc(node_alebc,nn_alebc) ! read in node index on the inner boundary
call formid_ale(kid,rng,ien) ! fix outer boundary
call innerbc_ale(disp_ale,kid,node_alebc,nn_alebc,nsd,nn) ! apply displacement on inner boundary
!do iit=1,2
kp(:,:)=0.0d0
kw(:,:)=0.0d0 ! clear the matrix
e_ale(:,:)=0.0d0 ! no extern force for this case
call blockm(xold,e_ale,disp_ale,kw,kp,ien,jac)
call setid(kp,kid,nsd)
call getnorm(kp,kp,nsd*nn,res_l)
res_l=sqrt(res_l)

kdg(:,:)=0.0d0 ! clear matrix
 call gmres_new(xold,kw,kp,kdg,ien,kid,jac)
! call gmresm(xref,kid,kw,kp,kdg,ien,kz,kv,kzg,kavg,ksm,kavloc,kh,ky,kcc,kss)
call getnorm(kdg,kdg,nsd*nn,del_l)
 del_l = sqrt(del_l)
 call update(kp,disp_ale,kdg,nsd)
!write(*,*) 'ALE iteration', 'ALE error norm', del_l
!end do
 
write(*,*) '**********   Mesh Updated   **********'

call add(x,xold,disp_ale,nsd)
!call defgrad(x,xref,ien,f,jac,finv,jacinv) ! old deformation gradient subroutine only for 3-D
call defgrad_new(x,xref,ien,f,jac,finv,jacinv) ! works for both 2-D and 3-D

!call defgrad_new(x,xref,ien,f_new,jac_new,finv_new,jacinv_new)
!write(*,*) 'diff between old and new in f', maxval(abs(jac_new(:,:)*jacinv_new(:,:)))
call velocity(x,xold,meshvel,refvel,finv,d,ien)

!??????????????????
!??? do not know what 'af' is 
!??????????????????
af(:)=0.0d0

call formid(id,rng,ien)  ! apply essential B.C 
call formd(d,rng,ien)
! Not sure about the moving inner ALE boundary
! Moving inner bloundary condition
call moving_bc(d,id,node_alebc,nn_alebc,meshvel)

! ??????????????????????????????????????????

do iit=1,nit
p(:,:)=0.0d0
w(:,:)=0.0d0
!call block_ale(xref,d,dold,af,p,w,hg,ien,f,finv,jac,jaco,refvel,refvelo) ! get residuals

call block_newale(x,d,dold,p,w,hg,ien,f_fluids,rng,f_stress,meshvel)

call setid(p,id,ndf) ! set residual be zero at essential B.C DOF
!write(*,*) 'p', p(:,:)
call getnorm(p,p,ndf*nn,res_l)
res_l =sqrt(res_l/nq)

dg(:,:)=0.0d0

!call gmres_ale(xref,d,dold,id,af,w,p,dg,hg,ien,z,v,zg,avg,sm, &
!	   vloc,avloc,h_gmres,y_gmres,cc,ss,finv,jac,jaco,refvel,refvelo)

call gmres_newale(x,d,dold,id,w,p,dg,hg,ien,    &
                z,v,zg,avg,sm,vloc,avloc,h_gmres,y_gmres,cc,ss,f_fluids,meshvel)

call getnorm(dg,dg,ndf*nn,del_l)
del_l = sqrt(del_l/nq)
call update(p,d,dg,ndf)

write(*,*) 'del_l', del_l
end do

write(*,*) 'After ALE'

f_fluids(:,:)=0.0

     include "hypo_write_output.fi"

  enddo time_loop



end subroutine hypo
