cSUBROUTINE DELTA
c Lucy Zhang
c 11/06/02
c Northwestern University

c This subroutine calculate the delta function which is used for both
c interpolation and distribution of the velocities and forces respectively
c between fluids and solids domain.

c There are 3 options for calculating
c 1. RKPM - cubic spline for non-uniform spacing
c 2. RKPM - cubic spline for uniform spacing
c 3. Original delta function for uniform spacing 

      subroutine diracdelta(data_solids,nn_solids, coord_pt, 
     +	data_fluids,nn_fluids,
     +	ndelta,shrknode,cnn,ncnn,maxconn,ibuf)

      integer ibuf,ndelta,maxconn
	parameter (maxnn_solids=1000)
	real* 8 shrknode(maxconn,maxnn_solids)

	!solids variables
      integer nn_solids
	integer ncnn(maxnn_solids),cnn(maxconn,maxnn_solids)
	real* 8 coord_pt(3,nn_solids),data_solids(3,nn_solids)

	!fluids variables
	integer nn_fluids
      real* 8 data_fluids(3,nn_fluids)

	!local variables
	integer i,j,k,inn,i0,j0,k0,icnn,pt
	real* 8 x(3)

c  ibuf=1 interpolation of velocity from the fluids domain to the solids domain
c  ibuf=2 distribution of forces from the solids domain to the fluids domain

ccccccccccccccccccccccccccccccccccccccccccccc
c Initialization of velocities and forces
ccccccccccccccccccccccccccccccccccccccccccccc
c      if (ibuf.eq.1) then  !velocity interpolation
c         write(*,*) 'velocity interpolation in process'
c      elseif (ibuf.eq.2) then !force distribution
c         write(*,*) 'force distribution in process'
c      endif

cccccccccccccccccccccccccccccccccccccccccccccc
c Calculate the delta functions
cccccccccccccccccccccccccccccccccccccccccccccc

c    If non-uniform grid 

      if (ndelta. eq. 1) then
         if (ibuf.eq.1) then    !velocity interpolation
		  data_solids(:,:)=0
            do inn=1,nn_solids
               do icnn=1,ncnn(inn)
                  pt=cnn(icnn,inn)
                  data_solids(1:3,inn) = data_solids(1:3,inn)
     +                 + data_fluids(1:3,pt) * shrknode(icnn,inn)
               enddo
            enddo
         elseif (ibuf.eq.2) then !force distribution
		  data_fluids(:,:)=0
            do inn=1,nn_solids
               do icnn=1,ncnn(inn)
                  pt=cnn(icnn,inn)
                  data_fluids(1:3,pt) = data_fluids(1:3,pt)
     +                 + data_solids(1:3,inn) * shrknode(icnn,inn)
               enddo
            enddo
         endif
      else

c    if uniform grid
c         do inn=1,nn_solids
c            if (ndelta .eq.2) then
c               call delta_rkpm_uniform(deltaeachaxis,inn,
c     +              dlptlocal_number,coord_pt)
c            elseif (ndelta .eq.3) then
c               call delta_original(deltaeachaxis,inn,
c     +              dlptlocal_number,coord_pt)
c            endif
c
c            do k=0,n_c_del_oneaxis-1
c               do j=0,n_c_del_oneaxis-1
c                  do i=0,n_c_del_oneaxis-1
c                     wt_delta(i,j,k) = deltaeachaxis(i,ix)
c     $                    * deltaeachaxis(j,iy)
c     $                    * deltaeachaxis(k,iz)
c                  enddo
c               enddo
c            enddo
c         
c            i0=modrealinto( coord_pt(ix,ipt), mn_ce1, mx_ce1)
c     $           - n_lo_gc_del_oneaxis
c            j0=modrealinto( coord_pt(iy,ipt), mn_ce2, mx_ce2)
c     $           - n_lo_gc_del_oneaxis
c            k0=modrealinto( coord_pt(iz,ipt), mn_ce3, mx_ce3)
c     $           - n_lo_gc_del_oneaxis
c         
c            if (ibuf .eq. 2) then !force distribution
c               do k = 0, n_c_del_oneaxis - 1 
c                  do j = 0, n_c_del_oneaxis - 1 
c                     do i = 0, n_c_del_oneaxis - 1       
c                        do idim = ix, iz
c                           data_fluids(i0+i,j0+j,k0+k,idim)
c     $                          = data_fluids(i0+i,j0+j,k0+k,idim)
c     $                          + data_solids(idim,ipt) * wt_delta(i,j,k)
c                        
c                        enddo
c                     enddo
c                  enddo
c               enddo
c            
c            elseif (ibuf .eq. 1) then !velocity interpolation
c               do k=0,n_c_del_oneaxis - 1 
c                  do j=0,n_c_del_oneaxis - 1 
c                     do i=0,n_c_del_oneaxis - 1 
c                        do idim = ix, iz
c                           data_solids(idim,ipt)
c     $                          = data_solids(idim,ipt)
c     $                          + data_fluids(i0+i,j0+j,k0+k, idim)
c     $                          * wt_delta(i,j,k)
c                        enddo
c                     enddo
c                  enddo
c               enddo
c            endif 
c
c         enddo ! end of loop in nn_solids
c
      endif ! end of option for delta function

      
      return
      end