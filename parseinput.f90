module parseinput
  implicit none
  save

  !Read_Real,Read_Int are stored in read.f
  ! If choise has been added, to make the program have two options of data exchanging 1---rkpm 2---FEM
  ! Febn12nd 2008 Xingshi
contains

subroutine parseinput_solid

  use global_simulation_parameter
  use solid_variables
  use meshgen_solid
  use delta_nonuniform, only: ndelta
  use r_common
  use mpi_variables ! call mpi variable module
  implicit none
   include 'mpif.h'    
  real(8) :: ftemp
  integer :: i,error_id
  integer :: nbe(nel)
  integer,parameter :: one = 1
  integer :: file_in,echo_out
  common /filename/file_in,echo_out
  real(8) :: shift1(3)

  file_in=8
  open(file_in,file='input_solid.in',status='old',action='read')
if (myid == 0) then
  write(*,*) 'reading input_solid.in'
end if
  ! Initial displacement
  CALL Read_Int(ninit,1)
  CALL Read_Int(initdir,1)
if (myid == 0) then
  if (ninit.eq.1)   write(*,*) 'apply initial condition' 
  if (initdir.eq.1) write(*,*) 'apply initial displacement'
  if (initdir.eq.2) write(*,*) 'apply initial velocity'
end if
 !...Gauss integration point
  CALL Read_Int(iquad_solid,1)
  CALL Read_Int(nen_solid,1)
if (myid ==0) then
  write(*,*) 'gauss integration type     : ',iquad_solid
  write(*,*) 'number of nodes per element: ',nen_solid
end if
 !...number of dimensions in the solid domain
  CALL Read_Int(nsd_solid,1)
if (myid == 0) then
  write(*,*) 'number of space dimensions in solid domain: ',nsd_solid
end if
 !...rigid body or not
  CALL Read_Int(nrigid,1)
if (myid ==0) then
  if (nrigid.eq.1) then
     write(*,*) 'treating the solid as a RIGID BODY'
  else
     write(*,*) 'treating the solid as a DEFORMABLE BODY'
  endif
end if
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccc from the original main_dat
  CALL Read_Int(n_tec_ens,1)
  CALL Read_Int(ndelta,1)
if (myid ==0) then
  if (ndelta .eq. 1) then
     write(*,*) ' RKPK data exchange method', ndelta

  else if(ndelta .eq. 2) then
     write(*,*) 'FEM data exchage method',ndelta
  else
     write(*,*) 'Currently we only offer these two options...'
  endif
end if
!cccccccccccccccccccccccccccccccccccccccccccccc
!
!     use ibm information
!
  CALL Read_Int(nn_solid_1,1)
  CALL Read_Int(ne_solid_1,1)
  CALL Read_Int(nump,1)
  
  CALL Read_Real(solid_scale,nsd_solid)
  CALL Read_Int(n_solid,1)

  nn_solid = nn_solid_1 * n_solid
  ne_solid = ne_solid_1 * n_solid
if (myid == 0) then
  if (mno .lt. nn_solid) then !...see r_common.f90
     write(*,*) 'boost mno in r_common.f'
     stop
  endif

  write(*,*) 'number of nodes    (for one object)  = ',nn_solid_1
  write(*,*) 'number of elements (for one object)  = ',ne_solid_1
  write(*,*) 'number of multiple defined solids    = ',n_solid
  write(*,*) 'number of nodes    (all objects)     = ',nn_solid
  write(*,*) 'number of elements (all objects)     = ',ne_solid
end if
  allocate(shift(nsd_solid,n_solid),stat=error_id)
  do i=1,n_solid
    CALL Read_Real(shift1,nsd_solid)
    shift(1:nsd_solid,i)=shift1(1:nsd_solid)
  enddo

 !...time step
  CALL Read_Int(ntfun,1)

 !...pressure force
  CALL Read_Int(numfn,1)
  CALL Read_Int(numeb,1)
  do i=1,numeb
     read(8,*) nbe(i),nface(i),boup(i,1:nsd_solid)
  enddo

 !...concentrated force
  do i=1,numfn
     read(8,*) nodefn(i),ndirfn(i),ftemp
     fnod(nodefn(i),ndirfn(i)) = ftemp*1.0d5
  enddo

  CALL Read_Int(material_type,1)	!1=hyperelastic material, 2=linear elastic material  
  ! Read in 2 young's modules for 2 solid parts=====================
  CALL Read_Real(young_mod,1)		! young's modulus (elastic material)
  !============================================================
  CALL Read_Real(Poisson,1)			! Poisson ratio (elastic material_
  CALL Read_Real(rc1,1)				! constants C1 (hyperelastic material)
  CALL Read_Real(rc2,1)				! constants C2 (hyperelastic material)
  CALL Read_Real(rk,1)				! constants Ck (hyperelastic material)
  CALL Read_Real(density_solid,1)	! Density solid-fluid
if (myid ==0) then
  if (material_type==1) write(*,*) 'The solid is HYPERELASTIC material'
  if (material_type==2) write(*,*) 'The solid is LINEAR ELASTIC material'
  write(*,*) 'Youngs modulus =', young_mod
  write(*,*) 'Poisson ratio=', Poisson
  write(*,*) 'C1  = ', rc1
  write(*,*) 'C2  = ', rc2
  write(*,*) 'kappa   = ',rk
  write(*,*) 'density = ',density_solid
end if
  CALL Read_Real(xviss,1)
if (myid ==0) then
  write(*,*) 'xviss  = ',xviss
end if
  vnorm = 0.0d0

  CALL Read_Real(vnorm,1)
  CALL Read_Real(fnorm,1)

 !...gravity acceleration
  CALL Read_Real(xmg,nsd_solid)
if (myid ==0) then
  write(*,*) 'gravity acceleration'
  write(*,*) ' xmg   =',xmg(1:nsd_solid)
end if

  close(8)

  prec(1:nump*ne_solid)=0.0d0

  return
end subroutine parseinput_solid


!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   sets parameters from standard input or control file                  c
!   ---------------------------------------------------------------------c
!   written by Lucy Zhang, Axel Gerstenberger                            c
!   ---------------------------------------------------------------------c
!   Northwestern University                                              c
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine parseinput_fluid
  use run_variables
  use fluid_variables
  use ale_variables
  use mpi_variables ! call mpi variable module
  implicit none

  integer :: restart_onoff, steady_onoff,hg_vol_onoff, taudt_onoff
  integer :: static_onoff, conserve_onoff, stokes_onoff
  real(8) :: fix(5),read_delta(2)
  integer :: idelta,i,ibc,idf,isd
  integer :: file_in,echo_out
  common /filename/file_in,echo_out
  integer,parameter :: one = 1

  file_in=5
  echo_out=7
  OPEN(file_in,file='input_fluid.in',STATUS='old',ACTION='read')
  OPEN(echo_out,file='summary.dat',STATUS="unknown")
if (myid == 0) then
  write(*,*) 'reading input_fluid.in'
end if
  CALL Read_Int(restart_onoff,1)
  if (restart_onoff.eq.1) then
  !   restart=.TRUE.
     restart = 1
  elseif (restart_onoff.ne.0) then
  !   restart=.FALSE.
     restart = restart_onoff
  endif
  CALL Read_Int(restart_freq,1)

  CALL Read_Int(steady_onoff,1)
  if (steady_onoff.eq.1) then
     steady=.TRUE.
  elseif (steady_onoff.eq.0) then
     steady=.FALSE.
  endif

  CALL Read_Int(stokes_onoff,1)
  if (stokes_onoff.eq.1) then
     stokes=.TRUE.
  elseif (stokes_onoff.eq.0) then
     stokes=.FALSE.
  endif

  CALL Read_Int(hg_vol_onoff,1)
  if (hg_vol_onoff.eq.1) then
     hg_vol=.TRUE.
  elseif (hg_vol_onoff.eq.0) then
     hg_vol=.FALSE.
  endif

  CALL Read_Int(taudt_onoff,1)
  if (taudt_onoff.eq.1) then
     taudt=.TRUE.
  elseif (taudt_onoff.eq.0) then
     taudt=.FALSE.
  endif

  CALL Read_Int(static_onoff,1)
  if (static_onoff.eq.1) then
     static=.TRUE.
  elseif (static_onoff.eq.0) then
     static=.FALSE.
  endif

  CALL Read_Int(conserve_onoff,1)
  if (conserve_onoff.eq.1) then
     conserve=.TRUE.
  elseif (conserve_onoff.eq.0) then
     conserve=.FALSE.
  endif

  CALL Read_Int(nn,1)
  CALL Read_Int(ne,1)
  CALL Read_Int(nrng,1)
  CALL Read_Int(nen,1)
  CALL Read_Int(ndf,1)
  CALL Read_Int(nsd,1)
  if (nsd==2) then
    udf=1
	vdf=2
	pdf=3
  elseif (nsd==3) then
    udf=1
	vdf=2
	wdf=3
	pdf=4
  endif

!write(*,*) 'nn', nn
!write(*,*) 'ne', ne
!write(*,*) 'nrng', nrng
!write(*,*) 'ndf', ndf
!write(*,*) 'nsd', nsd

  CALL Read_Int(iquad,1)
  CALL Read_Int(nit,1)
  CALL Read_Int(nts,1)
  CALL Read_Int(ntsbout,1)
  CALL Read_Int(idisk,1)
  CALL Read_Int(inner,1)
  CALL Read_Int(outer,1)
  CALL Read_Int(iscaling,1)
  CALL Read_Int(maxconn,1)

  CALL Read_Real(dt,1)
  CALL Read_Real(t_start,1)
  CALL Read_Real(alpha,1)

  do i=1,nrng
                do idf=1,ndf+1
        CALL Read_Real(fix(idf),1)
                end do
!        write(*,*) 'fix', fix(:)
     ibc=int(fix(1))
     do idf=1,ndf
        bv(idf,ibc) = fix(idf+1)
        if(abs(bv(idf,ibc)+999.0).gt.1.0e-6) bc(idf,ibc) = 1
     enddo
  enddo

! Read in the nature boundary condition
  call Read_Int(edge_inflow,1)
  call Read_Int(ne_inflow,1)
  call Read_Real(pin,1)
!======================================
! Read in number of nodes on ALE moving boundary
 Call Read_Int(nn_alebc,1)

  CALL Read_Real(landa_over_mu,1)
  CALL Read_Real(ic,ndf)
  CALL Read_Real(gravity,nsd)
  CALL Read_Real(interface,nsd)
  CALL Read_Real(vis_liq,1)
  CALL Read_Real(den_liq,1)
  CALL Read_Real(ref_lgt,1)
  CALL Read_Real(ref_vel,1)
  CALL Read_Real(ref_den,1)

  CALL Read_Real(read_delta,2)
  idelta=int(read_delta(1))
  delta(idelta)=read_delta(2)

  CALL Read_Real(read_delta,2)
  idelta=int(read_delta(1))
  delta(idelta)=read_delta(2)

  CALL Read_Real(read_delta,2)
  idelta=int(read_delta(1))
  delta(idelta)=read_delta(2)

  CALL Read_Real(read_delta,2)
  idelta=int(read_delta(1))
  delta(idelta)=read_delta(2)

  CALL Read_Real(read_delta,2)
  idelta=int(read_delta(1))
  delta(idelta)=read_delta(2)

  CALL Read_Real(read_delta,2)
  idelta=int(read_delta(1))
  delta(idelta)=read_delta(2)

  CALL Read_Real(read_delta,2)
  idelta=int(read_delta(1))
  delta(idelta)=read_delta(2)

  CALL Read_Real(turb_kappa,1)

!       further defaults
  if (ntsbout.eq.0) ntsbout = nts + 1
  if (steady) alpha = 1.0

  if      (nen.eq.3) then
     etype = tri
     neface = 3
     nnface = 2
  else if ((nen.eq.4).and.(nsd.eq.3)) then
     etype = tet
     neface = 4
     nnface = 3
  else if ((nen.eq.4).and.(nsd.eq.2)) then
     etype = qud
     neface = 4
     nnface = 2
  else if (nen.eq.8) then
     etype = hex
     neface = 6
     nnface = 4
  end if

  if ( nq.eq.0) nq  = ndf * nn
  twod = .false.
  if((nn.gt.2*ne).and.(nen.eq.8)) then
     twod = .true.
     hg_vol = .true.
  endif

  CLOSE(5)

  if(nsd.eq.3) then
	call shape
  elseif(nsd.eq.2) then
	call shape2d
  endif

  return
end subroutine parseinput_fluid


end module parseinput
