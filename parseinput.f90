module parseinput
  implicit none
  save

  !Read_Real,Read_Int are stored in read.f

contains

subroutine parseinput_solid
  use global_simulation_parameter
  use solid_variables
  use meshgen_solid
  use delta_nonuniform, only: ndelta
  use r_common
  implicit none
     
  
  real(8) :: ftemp
  real(8) :: fbacc(3)

  integer :: n_ibmfem,i,error_id
  integer :: nbe(nel)
  integer,parameter :: fileunit = 21

  nsurface = 4

  open(fileunit,file='input_solid.in',status='old',readonly)
  write(*,*) 'reading input_solid.in'

  read(fileunit,*) ninit,initdir
  if (ninit.eq.1)   write(*,*) 'apply initial condition' 
  if (initdir.eq.1) write(*,*) 'apply initial displacement'
  if (initdir.eq.2) write(*,*) 'apply initial velocity'

 !...Gauss integration point
  read(fileunit,*) iquad_solid,nen_solid
  write(*,*) 'gauss integration type     : ',iquad_solid
  write(*,*) 'number of nodes per element: ',nen_solid

 !...number of dimensions in the solid domain
  read(fileunit,*) nsd_solid
  write(*,*) 'number of space dimensions in solid domain: ',nsd_solid

 !...rigid body or not
  read(fileunit,*) nrigid
  if (nrigid.eq.1) then
     write(*,*) 'treating the solid as a RIGID BODY'
  else
     write(*,*) 'treating the solid as a DEFORMABLE BODY'
  endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccc from the original main_dat
  read(fileunit,*) n_ibmfem, n_tec_ens


  read(fileunit,*) ndelta
  if (ndelta /= 1) then
     write(*,*) " currently, only 1 deltafunction type is implemented ..."
     stop
  else
     write(*,*) 'ndelta=',ndelta
  endif

!cccccccccccccccccccccccccccccccccccccccccccccc
!
!     use ibm information
!
  read(fileunit,*) nn_solid_1,ne_solid_1
  read(fileunit,*) nump
  
  read(fileunit,*) solid_scale(1),solid_scale(2),solid_scale(3)
  
  read(fileunit,*) n_solid

  !if (n_solid .gt. n_solid_max) then
  !   write(*,*) 'BOOST n_solid_max in r_common'
  !   stop
  !endif

  nn_solid = nn_solid_1 * n_solid
  ne_solid = ne_solid_1 * n_solid

  if (mno .lt. nn_solid) then !...see r_common.f90
     write(*,*) 'boost maxnn_solids in hypo.f and delta_nonuniform'
     stop
  endif

  write(*,*) 'number of nodes    (for one object)  = ',nn_solid_1
  write(*,*) 'number of elements (for one object)  = ',ne_solid_1
  write(*,*) 'number of multiple defined solids    = ',n_solid
  write(*,*) 'number of nodes    (all objects)     = ',nn_solid
  write(*,*) 'number of elements (all objects)     = ',ne_solid

  allocate(shift(nsd_solid,n_solid),stat=error_id)
  do i=1,n_solid
     read(fileunit,*) shift(1,i),shift(2,i),shift(3,i)
  enddo

 !...time step
  read(fileunit,*) ntfun

 !...pressure force
  read(fileunit,*) numfn,numeb      
  do i=1,numeb
     read(fileunit,*) nbe(i),nface(i),boup(i,1),boup(i,2),boup(i,3)
  enddo

 !...concentrated force
  do i=1,numfn
     read(fileunit,*) nodefn(i),ndirfn(i),ftemp
     fnod(nodefn(i),ndirfn(i)) = ftemp*1.0d5
  enddo



  read(fileunit,*) rc1,rc2,rk,density_solid
  write(*,*) 'C1      = ',rc1
  write(*,*) 'C2      = ',rc2
  write(*,*) 'kappa   = ',rk
  write(*,*) 'density = ',density_solid

  read(fileunit,*) xviss
  write(*,*) 'xviss  = ',xviss

  vnorm = 0.0d0

  read(fileunit,*) vnorm,fnorm,vtol,ftol
  read(fileunit,*) alpha_solid,beta_solid

 !...gravity acceleration
  read(fileunit,*) xmg(1),xmg(2),xmg(3)  !gravity
  write(*,*) 'gravity acceleration'
  write(*,*) ' xmg(1)   =',xmg(1)
  write(*,*) ' xmg(2)   =',xmg(2)
  write(*,*) ' xmg(3)   =',xmg(3)

 !...body force
  read(fileunit,*) fbacc(1),fbacc(2),fbacc(3)
  write(*,*) 'body force acceleration - not activated'
  write(*,*) ' fbacc(1) =',fbacc(1)
  write(*,*) ' fbacc(2) =',fbacc(2)
  write(*,*) ' fbacc(3) =',fbacc(3)

  close(fileunit)

  prec(1:nump*ne_solid)=0.0d0

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
  use delta_nonuniform, only: maxconn
  use meshgen_gmsh
  use adaptive_meshing
  use restart_lib
  implicit none

  !character*32 key, keyaux
  !character*8 onoff
  !logical fctrl, getkey, isatty
  !logical enough
  !data enough /.false./
  integer :: restart_onoff, steady_onoff,hg_vol_onoff, taudt_onoff
  integer :: static_onoff, conserve_onoff, stokes_onoff
  real(8) :: fix(5),fixd(4),read_delta(2)
  integer :: idelta,i,ibc,idf,isd
  integer :: file_in,echo_out
  common /filename/file_in,echo_out

  integer,parameter :: one = 1

  nen_surf = 3

  file_in=5
  echo_out=7
  OPEN(file_in,file='input_fluid.in',STATUS='old',ACTION='read')
  OPEN(echo_out,file='summary.dat',STATUS="unknown")

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
  elseif (hg_vol_onoff.eq.0) then
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

  if (gmsh_input == 1) then 
     call read_gmsh_node_elem_numbers(nn,ne,ne_surf,fluid_mesh_filename_length,fluid_mesh_filename)  !...get node and element numbers from Gmsh input
     if (adapt_mesh == 1) then
        call read_gmsh_node_elem_numbers(nn_bg,ne_bg,ne_surf_bg,bgmesh_filename_length,bgmesh_filename)
     endif
  endif

  CALL Read_Int(nrng,1)
  CALL Read_Int(nen,1)
  CALL Read_Int(ndf,1)
  CALL Read_Int(nsd,1)
  CALL Read_Int(iquad,1)
  CALL Read_Int(nit,1)
  CALL Read_Int(nts,1)
  CALL Read_Int(ntsbout,1)
  CALL Read_Int(idisk,1)
  CALL Read_Int(inner,1)
  CALL Read_Int(outer,1)

  CALL Read_Int(kinner,1)
  CALL Read_Int(kouter,1)
  CALL Read_Int(iscaling,1)
  CALL Read_Int(maxconn,1)

  CALL Read_Real(dt,1)
  CALL Read_Real(t_start,1)
  CALL Read_Real(alpha,1)

  do i=1,nrng
     CALL Read_Real(fix,5)
     ibc=int(fix(1))
     do idf=1,ndf
        bv(idf,ibc) = fix(idf+1)
        if(abs(bv(idf,ibc)+999.0).gt.1.0e-6) bc(idf,ibc) = 1
     enddo
  enddo

  !amp = 0.4
  !spring = 11.302
  !damper = 2.0
  !mass = 35
  do i=1,nrng
     CALL Read_Real(fixd,4)
     ibc=int(fixd(1))
     do isd=1,nsd
        bvd(isd,ibc)=fixd(isd+1)
        if(abs(bvd(isd,ibc)+999.0).gt.1.0e-6) bcd(isd,ibc) = 1
     enddo
  enddo


  CALL Read_Real(landa_over_mu,1)

  CALL Read_Real(ic,ndf)

  CALL Read_Real(gravity,3)
  CALL Read_Real(interface,3)

  CALL Read_Int(surf,2)
  CALL Read_Int(hydro,1)

  CALL Read_Real(vis_liq,1)
  CALL Read_Real(vis_gas,1)
  CALL Read_Real(den_liq,1)
  CALL Read_Real(den_gas,1)
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

  if      (nen.eq.4) then
     etype = tet
     neface = 4
     nnface = 3
  else if (nen.eq.8) then
     etype = hex
     neface = 6
     nnface = 4
  end if

  !if (ndf.eq.0) ndf = 4
  !if (nen.eq.0) nen = 8
  if ( nq.eq.0) nq  = ndf * nn
  !if (nqf.eq.0) nqf = nn
  twod = .false.
  if((nn.gt.2*ne).and.(nen.eq.8)) then
     twod = .true.
     hg_vol = .true.
  endif



  close(file_in)
  close(echo_out)

  call shape
  !call shape2d

end subroutine parseinput_fluid


end module parseinput