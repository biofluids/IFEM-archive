c     
c     adina input coordinate
c
      subroutine r_input
      implicit real*8 (a-h,o-z)
      include 'r_common'
      include 'main_common'      
c
      include 'ibd0_exchange_pars.fh'
      include 'ibg_variables_cfd.fh'
      include 'ibg_variables_cloud.fh'
      include 'ibg_variables_domain.fh'
      include 'ibg_variables_fluid.fh'
      include 'ibg_variables_marker.fh'
      include 'ibg_variables_point.fh'
      include 'ibg_variables_run.fh'
      include 'ibg_variable_equivalences.fh'
      include 'iba_application_parameters.fh'
      include 'ibg_change_me_ptcon_var_decl.fh'
      include 'ibg_change_me_ptcon_var_common_equiv.fh'
      include 'ibg_parameters_run.fh' 
c
c     print control
c     if ninit=1 initial condition
c     if initdir=1/2 initial disp/vel

c	read in coortable	

      read(1,*) ninit,initdir
      if (ninit.eq.1) write(*,*) 'apply initial condition' 
      if (initdir.eq.1) write(*,*) 'apply initial displacement'
      if (initdir.eq.2) write(*,*) 'apply initial velocity'
c     if npr=1/0 plot.m gives the coor/dispc     if ntprint=1 time.m gives the timefunstion

      read(1,*) npr,npdis,ntprint

      do 31 i=1,npr
         read(1,*) ndprint(i),nprint(i)
 31   continue

c     data input check and iteration control
c     nbouc -- boundary change influence the distributed forces
      read(1,*) nbouc,nchkread
c     Gauss integration point
      read(1,*) nint,nis
      write(*,*) 'number of gauss integration point=',nint
      write(*,*) 'number of nodes per element=',nis
c
c     use ibm information
c
      read(1,*) nnd,numel,nnda,numela
      read(1,*) nump
      write(*,*) 'number of nodes=',nnd
      write(*,*) 'number of elements=',numel
c
      do 42 i=1,numel
         read(1,*) ntx,(nea(i,k),k=1,nis),ntemps1,ntemps2,
     $        ntemps3,ntemps4,ntemps5
 42   continue
c
c for 2-D case only 
      do j=1,nnd
         read(1,*) ntemps1,xlay(j),(coor(j,i),i=1,2), ntemps2
	enddo
	coor(:,:)=coor(:,:)*xnmag

      read(1,*) numskew,numgb,intnum,numct
      write(*,*) 'numgb=',numgb

      ntether=0
      do 7 j=1,numgb
         read(1,*) ndirgb(j)
         read(1,*) numdir(j)
         do 77 k=1,numdir(j)
            read(1,*) nodegb(j,k),
     $           nxt(1,nodegb(j,k)),nxt(3,nodegb(j,k))
 77      continue
c
         if (ndirgb(j) .eq. 111111) then
            do 38 i=1,numdir(j)
               if (nxt(1,nodegb(j,i)) .ne. 0) then
                  ntether=ntether+1
               endif
               if (nxt(3,nodegb(j,i)) .ne. 0) then
                  ntether=ntether+1
               endif
 38         continue
         endif
 7    continue  
c     skew system
      write(*,*) 'number of tether points=',ntether
      do 9 i=1,intnum
         read(1,*) numint(i),ninsk(i)
 9    continue
      do 10 i=1,numskew
         read(1,*) xang(i)
 10   continue
c     slave and master nodes
      do 11 i=1,numct
 11   continue
c     time step
      read(1,*) ntfun
c     data print out step
      read(1,*) nina
c     pressure force
      read(1,*) numfn,numeb      
      do 12 i=1,numeb
         read(1,*) nbe(i),nface(i),boup(i,1),boup(i,2)
 12   continue
c     concentrated force
      do 13 i=1,numfn
         read(1,*) nodefn(i),ndirfn(i),ftem
         fnod(nodefn(i),ndirfn(i))=ftem*1.0d5/unit_force
 13   continue
c     body force
      read(1,*) fbacc(1),fbacc(2)
      write(*,*) 'body force=',fbacc(1),'i+',fbacc(2),'j'
c     time functions
      do 14 i=1,2
         read(1,*) nfuns(i)
         if (nfuns(i) .eq. 2) then
            read(1,*) xome(i)
         endif
 14   continue
c
      read(1,*) thic,nreact
      if (nreact .eq. 1) then
         read(1,*) nrtp
         do 700 i=1,nrtp
            read(1,*) ntt,mtt
            ndraf(i)=ntt
            nraf(i)=mtt
 700     continue
      endif
      read(1,*) rc1,rc2,rk,sdensi
      rc1=rc1/unit_pressure
      rc2=rc2/unit_pressure
      rk=rk/unit_pressure
      sdensi=sdensi/unit_density
      iflag=0
c
c     in units
c
      read(1,*) xk,xtedis,xstretch,xvisc,xviss

      xk=xk/(unit_pressure*unit_length)
      xtedis=xtedis/unit_length
      xvisc=xvisc/unit_viscosity
      xviss=xviss/unit_viscosity
      
      write(*,*) 'xk=',xk
      write(*,*) 'xtedis=',xtedis
      write(*,*) 'xvisc=',xvisc
      write(*,*) 'xviss=',xviss

      read(1,*) vnorm,fnorm,vtol,ftol
      read(1,*) alpha,beta

      read(1,*) xmg1,xmg2  !gravity

      xmg1=xmg1/unit_length*unit_time**2
      xmg2=xmg2/unit_length*unit_time**2

      write(*,*) 'xmg1=',xmg1
      write(*,*) 'xmg2=',xmg2

      return
      end