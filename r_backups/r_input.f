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
	real*8 coor_1(mno,3)
	integer nea_1(nel,9)
c
c     print control
c     if ninit=1 initial condition
c     if initdir=1/2 initial disp/vel

c	read in coortable.in	

      read(1,*) ninit,initdir
      if (ninit.eq.1) write(*,*) 'apply initial condition' 
      if (initdir.eq.1) write(*,*) 'apply initial displacement'
      if (initdir.eq.2) write(*,*) 'apply initial velocity'
c     if npr=1/0 plot.m gives the coor/dispc     if ntprint=1 time.m gives the timefunstion

      read(1,*) npr,npdis,ntprint

      do i=1,npr
         read(1,*) ndprint(i),nprint(i)
	enddo

c     data input check and iteration control
c     nbouc -- boundary change influence the distributed forces
      read(1,*) nbouc,nchkread
c     Gauss integration point
      read(1,*) nint,nis
      write(*,*) 'number of gauss integration point=',nint
      write(*,*) 'number of nodes per element=',nis

c	number of dimensions in the solid domain
	read(1,*) nsd_solid
	write(*,*) 'number of dimensions in the solid domain=',nsd_solid
c	rigid body or not
	read(1,*) nrigid
	if (nrigid.eq.1) then
	 write(*,*) 'treating the solid as a RIGID BODY'
	else
	 write(*,*) 'treating the solid as a NON-RIGID BODY'
	endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc from the original main_dat
      read(1,*) n_ibmfem, n_tec_ens, n_dispforce
      read(1,*) n_ibmfv,applyv,applyf

      write(*,*) 'applyv=',applyv
      write(*,*) 'applyf=',applyf

      read(1,*) npmax
      write(*,*) 'npmax=',npmax
      if (npmax .ne. 0) then
         do 30 i=1,npmax
            read(1,*) npout(i)
 30      continue
      endif

      read(1,*) threshold
      write(*,*) 'threshold=',threshold

      read(1,*) ndelta
      write(*,*) 'ndelta=',ndelta

ccccccccccccccccccccccccccccccccccccccccccccccc
c	
c     use ibm information
c
      read(1,*) nnd,numel,nnda,numela
      read(1,*) nump
	read(1,*) n_solid
	if (n_solid .gt. n_solid_max) then
		write(*,*) 'BOOST n_solid_max in r_common'
		stop
	endif

c	if (n_solid .gt.1) then
	  do i=1,n_solid
		read(1,*) shift(1,i),shift(2,i),shift(3,i)
	  enddo
c	endif

      write(*,*) 'number of nodes=',nnd
      write(*,*) 'number of elements=',numel
c
      do i=1,numel
         read(1,*) ntx,(nea_1(i,k),k=1,nis)
	enddo
c
c for 2-D case only 
      do j=1,nnd
         read(1,*) ntemps1,(coor_1(j,i),i=1,3), ntemps2
	enddo
	coor_1(:,1)=coor_1(:,1)*xnmag
	coor_1(:,2)=coor_1(:,2)*xnmag
	coor_1(:,3)=coor_1(:,3)*xnmag

	do i=1,n_solid
         coor(nnd*(i-1)+1:nnd*i,1)=coor_1(1:nnd,1)+shift(1,i)
	   coor(nnd*(i-1)+1:nnd*i,2)=coor_1(1:nnd,2)+shift(2,i)
	   coor(nnd*(i-1)+1:nnd*i,3)=coor_1(1:nnd,3)+shift(3,i)
	   nea(numel*(i-1)+1:numel*i,1:nis)=
     +		 nea_1(1:numel,1:nis)+(i-1)*nnd
	enddo

	nnd=nnd*n_solid
	numel=numel*n_solid


      read(1,*) numskew,numgb,intnum,numct
      write(*,*) 'numgb=',numgb

      ntether=0
      do j=1,numgb
         read(1,*) ndirgb(j)
         read(1,*) numdir(j)
         do k=1,numdir(j)
            read(1,*) nodegb(j,k),
     $       nxt(1,nodegb(j,k)),nxt(2,nodegb(j,k)),nxt(3,nodegb(j,k))
	   enddo
c
         if (ndirgb(j) .eq. 111111) then
            do i=1,numdir(j)
               if (nxt(1,nodegb(j,i)) .ne. 0) then
                  ntether=ntether+1
               endif
               if (nxt(2,nodegb(j,i)) .ne. 0) then
                  ntether=ntether+1
               endif
               if (nxt(3,nodegb(j,i)) .ne. 0) then
                  ntether=ntether+1
               endif
		  enddo
         endif
	enddo  
c     skew system
      write(*,*) 'number of tether points=',ntether
      do i=1,intnum
         read(1,*) numint(i),ninsk(i)
	enddo
      do i=1,numskew
         read(1,*) xang(i)
	enddo
c     slave and master nodes
      do i=1,numct
	enddo
c     time step
      read(1,*) ntfun
c     data print out step
      read(1,*) nina
c     pressure force
      read(1,*) numfn,numeb      
      do i=1,numeb
         read(1,*) nbe(i),nface(i),boup(i,1),boup(i,2),boup(i,3)
	enddo
c     concentrated force
      do i=1,numfn
         read(1,*) nodefn(i),ndirfn(i),ftem
         fnod(nodefn(i),ndirfn(i))=ftem*1.0d5
	enddo
c     body force
      read(1,*) fbacc(1),fbacc(2),fbacc(3)
      write(*,*) 'body force=',fbacc(1),'i+',fbacc(2),'j'
c     time functions
      do i=1,3
         read(1,*) nfuns(i)
         if (nfuns(i) .eq. 2) then
            read(1,*) xome(i)
         endif
	enddo
c
      read(1,*) thic,nreact
      if (nreact .eq. 1) then
         read(1,*) nrtp
         do i=1,nrtp
            read(1,*) ntt,mtt
            ndraf(i)=ntt
            nraf(i)=mtt
	   enddo
      endif
      read(1,*) rc1,rc2,rk,sdensi

      iflag=0
c
c     in units
c
      read(1,*) xk,xtedis,xstretch,xvisc,xviss
      
      write(*,*) 'xk=',xk
      write(*,*) 'xtedis=',xtedis
      write(*,*) 'xvisc=',xvisc
      write(*,*) 'xviss=',xviss

      read(1,*) vnorm,fnorm,vtol,ftol
      read(1,*) alpha,beta

      read(1,*) xmg1,xmg2,xmg3  !gravity

      write(*,*) 'xmg1=',xmg1
      write(*,*) 'xmg2=',xmg2
	write(*,*) 'xmg3=',xmg3

      return
      end