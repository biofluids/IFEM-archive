!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine echoinput
      use run_variables
      use fluid_variables
      implicit none

      !character*8 date
      character(len=3) :: yon
      integer :: io,idelta,i

!     if (myid.ne.0) return

      io = 7
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(io,'(/"Control information")')
      write(io,'( "-------------------"/)')
      write(io,'(" Nodes........................(nn) = ",i7)') nn
      write(io,'(" Elements.....................(ne) = ",i7)') ne
      write(io,'(" Total integration points..(nquad) = ",i7)') nquad
      write(io,'(" Time steps....................(nts) = ",i5)') nts
      write(io,'(" Iterations....................(nit) = ",i5)') nit
      write(io,'(" Time steps b/output.......(ntsbout) = ",i5)') ntsbout
      write(io,'(" Space dimensions..............(nsd) = ",i5)') nsd
      write(io,'(" Degrees of freedom............(ndf) = ",i5)') ndf
      WRITe(io,'(" Number of element nodes.......(nen) = ",i5)') nen
      write(io,'(" Integration points..........(iquad) = ",i5)') iquad
      write(io,'(" Hydrostatic surface.........(hydro) = ",i5)') hydro
      write(io,'(" Scaling..................(iscaling) = ",i5)') iscaling
      write(io,'(" Inner GMRES iterations......(inner) = ",i5)') inner 
      write(io,'(" Outer GMRES iterations......(outer) = ",i5)') outer 
      write(io,'(" Restart...................(restart) = ",a5)') restart
      write(io,'(" Stokes.....................(stokes) = ",a5)') yon(stokes)
      write(io,'(" Steady.....................(steady) = ",a5)') yon(steady)
      write(io,'(" Mass conservation........(conserve) = ",a5)') yon(conserve)
      write(io,'(" Static Problem.............(static) = ",a5)') yon(static)
      write(io,'(" Lenght based on volume.....(hg_vol) = ",a5)') yon(hg_vol)
      write(io,'(" Dt in tau...................(taudt) = ",a5)') yon(taudt)
      write(io,'(" 2D computation...............(twod) = ",a5)') yon(twod)
      write(io,'(" Time step.................= ",e15.8)') dt
      write(io,'(" Initial time..............= ",e15.8)') t_start
      write(io,'(" Alpha.....................= ",e15.8)') alpha   
      write(io,'(" Reference length..........= ",e15.8)') ref_lgt 
      write(io,'(" Reference velocity........= ",e15.8)') ref_vel 
      write(io,'(" Reference density.........= ",e15.8)') ref_den 
      write(io,'(" Density of liquid.........= ",e15.8)') den_liq 
      write(io,'(" Density of gas............= ",e15.8)') den_gas 
      write(io,'(" Viscosity of liquid.......= ",e15.8)') vis_liq 
      write(io,'(" Viscosity of gas..........= ",e15.8)') vis_gas 
      write(io,'(" Gravity in x..............= ",e15.8)') gravity(1)
      write(io,'(" Gravity in y..............= ",e15.8)') gravity(2)
      write(io,'(" Gravity in z..............= ",e15.8)') gravity(3)
      write(io,'(" Interface in x............= ",e15.8)') interface(1)
      write(io,'(" Interface in y............= ",e15.8)') interface(2)
      write(io,'(" Interface in z............= ",e15.8)') interface(3)

      do idelta=0,21
      if (delta(idelta).ne.0.0) then
          write(io,'(" Delta_",i1,".............. = ",e15.8)') idelta, delta(idelta)
      end if
      end do

      if (delta(4).gt.1.0e-6) then
         write(io,'(/"Smagorinsky/Hitachi turbulence model")')
         write(io,'( "------------------------------------"/)')
         write(io,'( " Kappa........(turb_kappa) = ",e15.8)') turb_kappa
      else
         write(io,'(/"No turbulence modeling")')
      end if

          write(io,'(/"Boundary condition information" )')
          write(io,'( "------------------------------"/)')
      do i=1,nrng
	      write(io,14) i,bc(1:ndf,i) !,bcf(i)
 14   format (' boundary no. ',i2,3x,5i9)
      enddo
      do i=1,nrng
	      write(io,15) i,bv(1:ndf,i) !,bvf(i)
 15   format (' Boundary no. ',i2,3x,5f9.3)
      enddo
          write(io,'(/"Initial condition information" )')
          write(io,'( "-----------------------------"/)')
      write(io,16) ic(1:4) !,icf
 16   format (' Initial value : ',5f9.3)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      end

      function yon(switch)

      character(len=3) :: yon
      logical switch

      if (switch) then
          yon = "yes"
      else
          yon = " no"
      end if

      return
      end
