c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	parseinput.fcm                                                       c
c	---------------------------------------------------------------------c
c	sets parameters from standard input or control file                  c
c	---------------------------------------------------------------------c
c     written by Lucy Zhang - 12/4/02	                                                    c
c	---------------------------------------------------------------------c
c	Northwestern University                                              c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine parseinput_fluid
	use run_variables, only:ntsbout,nts,dt
	use fluid_variables
	use delta_nonuniform, only: maxconn
	use global_simulation_parameter, only:nrestart
	implicit none
	
!	character*32 key, keyaux
!	character*8 onoff
!	logical fctrl, getkey, isatty
!	logical enough
!	data enough /.false./
!	logical restart_onoff, 
	logical steady_onoff,hg_vol_onoff, taudt_onoff
	logical static_onoff, conserve_onoff, stokes_onoff
	real*8  fix(5),read_delta(2)
	integer idelta,i,ibc,idf
	integer file_in,echo_out
	common /filename/file_in,echo_out

	file_in=5
	echo_out=7
	OPEN(file_in,file='fluid.in',STATUS='old')
	OPEN(echo_out,file='summary.dat',STATUS="unknown")

c	CALL Read_Int(restart_onoff,1)
c	if (restart_onoff.eq.1) then
c		restart=.TRUE.
c	elseif (restart_onoff.eq.0) then
c		restart=.FALSE.
c	endif

	CALL Read_Int(nrestart,1)

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
      CALL Read_Int(iscaling,1)
      CALL Read_Int(maxconn,1)

	CALL Read_Real(dt,1)
      CALL Read_Real(t_start,1)
      CALL Read_Real(alpha,1)

	do i=1,nrng
		CALL Read_Real(fix,5)
		ibc=int(fix(1))
		do idf=1,ndf
			bv(idf,ibc)=fix(idf+1)
			if(abs(bv(idf,ibc)+999.0).gt.1.0e-6) bc(idf,ibc) = 1
		enddo
	enddo


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

c       further defaults
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

	if (ndf.eq.0) ndf = 4
	if (nen.eq.0) nen = 8
	if ( nq.eq.0) nq  = ndf * nn
	if (nqf.eq.0) nqf = nn
	twod = .false.
	if((nn.gt.2*ne).and.(nen.eq.8)) then
	   twod = .true.
	   hg_vol = .true.
	endif

      CLOSE(5)

	return
      END



	SUBROUTINE Read_Real(var,count)

      IMPLICIT REAL*8  (a-h,o-z)
      IMPLICIT INTEGER (i-n)

      REAL*8 var(*)
      INTEGER count,file_in,echo_out
      COMMON/filename/ file_in,echo_out
      CHARACTER line*132

*  Read lines in until reaching a non-blank line that is not commented out

      i = 0
 10   CONTINUE
      READ(file_in,'(a132)',err=911) line
      linelen = 133
 11   CONTINUE
      linelen = linelen - 1
      IF (line(linelen:linelen) .EQ. ' ' .and. 
     +       linelen .GT. 0) GOTO 11
      IF (echo_out .GT. 0)WRITE(echo_out,15)(line(k:k),k=1,linelen)
 15	FORMAT(132a1)
        
      IF (INDEX('!@#%/',line(1:1)) .NE. 0) GOTO 10
      IF (linelen .LE. 0) GOTO 10
        
*  Convert the string beyond the equals sign into the var(count) array

      j = INDEX(line,'=') + 1
 20   CONTINUE
      READ(line(j:linelen),*,err=40) var(i+1)
      i = i + 1
      IF (i .LT. count) THEN
 30		CONTINUE
          k = INDEX(line(j:linelen),' ')
          j = j + k
          IF (k .EQ. 1) GOTO 30
          IF (k .EQ. 0) GOTO 10
          GOTO 20
       ENDIF

 40   continue
      do while (i .lt. count)
        var(i+1) = var(i)
        i = i+1
      enddo


*  Now distribute the array var(count) to all nodes

      RETURN

*  Read Error routine

 911  CONTINUE
      print*,'Read_Real error - ',i,count
      STOP'Read_Real error'

      END


********************************************************************
      SUBROUTINE Read_Int(ivar,count)

      IMPLICIT REAL*8  (a-h,o-z)
      IMPLICIT INTEGER (i-n)

      INTEGER ivar(*),count,file_in,echo_out
      COMMON/filename/ file_in,echo_out
      CHARACTER line*132

* Read lines in until reaching a non-blank line that is not commented out

        i = 0
   10   CONTINUE
      READ(file_in,'(a132)',err=911) line
      linelen = 133
   11 CONTINUE
      linelen = linelen - 1
      IF (line(linelen:linelen) .EQ. ' ' .and.
     +     linelen .GT. 0) GOTO 11
      IF (echo_out .GT. 0)WRITE(echo_out,15)(line(k:k),k=1,linelen)
   15 FORMAT(132a1)
 
      IF (INDEX('!@#%/',line(1:1)) .NE. 0) GOTO 10
      IF (linelen .LE. 0) GOTO 10

* Convert the string beyond the equals sign into the ivar(count) array

      j = INDEX(line,'=') + 1
   20 CONTINUE
      READ(line(j:linelen),*,err=10) ivar(i+1)
      i = i + 1
      IF (i .LT. count) THEN
   30    CONTINUE
         k = INDEX(line(j:linelen),' ')
         j = j + k
         IF (k .EQ. 1) GOTO 30
         IF (k .EQ. 0) GOTO 10
         GOTO 20
      ENDIF
* Now distribute the array ivar(count) to all nodes


      RETURN

* Read Error routine

  911 CONTINUE
      print*,'Read_Int error - ',i,count
      STOP'Read_Int error'

      END