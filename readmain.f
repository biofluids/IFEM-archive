      subroutine readmain
      implicit real*8 (a-h,o-z)

      include 'r_common'
      include 'main_common'

      ntf=9988  ! file="main_dat"
      read(ntf,*) n_step_run,iramp
      n_step_wr_pc_file=n_step_run
      read(ntf,*) n_step_wr_ib_user_files
      read(ntf,*) nrestart
      if (nrestart .eq. 1) then
         read(ntf,*) xfactor
      endif
      read(ntf,*) xshift,zshift
      read(ntf,*) nslay
      read(ntf,*) r_nu, r_rho, t_step
      read(ntf,*) r_len_fluid_exp1
      read(ntf,*) r_len_fluid_exp2
      read(ntf,*) r_len_fluid_exp3
      read(ntf,*) n_ibmfem, n_tec_ens, n_dispforce
      read(ntf,*) n_mass, n_periodicvel
      read(ntf,*) mnk, mxk
      read(ntf,*) n_ibmfv,applyv,applyf

      !unit_length=r_len_fluid_exp1/n_ce1
	!reset unit_length=1.0
	unit_length=1.0
      unit_mass=r_rho*(unit_length**3)
	!reset unit_time
      !unit_time=t_step/n_iter_step
	unit_time=1.0
      unit_divergence = 1.0d0/unit_time
      unit_force= unit_mass*unit_length/(unit_time**2)
      unit_frc_density = unit_force   / (unit_length ** 3)
      unit_pressure   = unit_force    / (unit_length ** 2)
      unit_resistance = unit_pressure*unit_time/(unit_length ** 3)
      unit_velocity   = unit_length/unit_time
      unit_viscosity  = r_nu/((unit_length ** 2)/unit_time)
      unit_density   = unit_mass/unit_length**3
      unit_flow = unit_length**3/unit_time
      unit_volume = unit_length**3

      write(*,201) unit_viscosity
 201  format(1x,'unit_viscosity=',e23.5)

      write(*,202) unit_length
 202  format(1x,'unit_length=',e23.5)

      write(*,203) unit_time
 203  format(1x,'unit_time=',e23.5)

      write(*,204) unit_velocity
 204  format(1x,'unit_velocity=',e23.5)

      write(*,206) unit_pressure
 206  format(1x,'unit_pressure=',e23.5)

      write(*,207) unit_force
 207  format(1x,'unit_force=',e23.5)

      write(*,208) unit_mass
 208  format(1x,'unit_mass=',e23.5)

      applyv=applyv/unit_velocity
      applyf=applyf/unit_force

      write(*,*) 'applyv=',applyv
      write(*,*) 'applyf=',applyf

      read(ntf,*) npmax
      write(*,*) 'npmax=',npmax
      if (npmax .ne. 0) then
         do 30 i=1,npmax
            read(ntf,*) npout(i)
 30      continue
      endif

      read(ntf,*) nteclay
      read(ntf,*) threshold
      write(*,*) 'nteclay=',nteclay
      write(*,*) 'threshold=',threshold
c
      read(ntf,*) nmove,nbou
      read(ntf,*) ndelta

      write(*,*) 'ndelta=',ndelta
c      cma=0.8d1
      read(ntf,*) cma

      cma=cma/unit_mass

      write(*,*) 'cma=',cma

      return
      end

