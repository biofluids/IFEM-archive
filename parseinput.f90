module parseinput
  implicit none
  save

  !Read_Real,Read_Int are stored in read.f

contains



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
!  use delta_nonuniform, only: maxconn
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
  integer :: file_in,echo_out,maxconn
  common /filename/file_in,echo_out

  integer,parameter :: one = 1

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
     CALL Read_Real(fix,ndf+1)
     ibc=int(fix(1))
     do idf=1,ndf
        bv(idf,ibc) = fix(idf+1)
        if(abs(bv(idf,ibc)+999.0).gt.1.0e-6) bc(idf,ibc) = 1
     enddo
  enddo


  do i=1,nrng
     CALL Read_Real(fixd,nsd+1)
     ibc=int(fixd(1))
     do isd=1,nsd
        bvd(isd,ibc)=fixd(isd+1)
        if(abs(bvd(isd,ibc)+999.0).gt.1.0e-6) bcd(isd,ibc) = 1
     enddo
  enddo


  CALL Read_Real(landa_over_mu,1)
  CALL Read_Real(ic,ndf)

  CALL Read_Real(gravity,nsd)
  CALL Read_Real(interface,nsd)

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

  !if (ndf.eq.0) ndf = 4
  !if (nen.eq.0) nen = 8
  if ( nq.eq.0) nq  = ndf * nn
  !if (nqf.eq.0) nqf = nn
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