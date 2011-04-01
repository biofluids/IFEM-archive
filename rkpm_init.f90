subroutine rkpm_init(x_solid,nn_solid,xna,nsd,nn,dwjp,sp_radius)
! calculate RKPM coefficient parallelly a
use mpi_variables
use delta_nonuniform, only: cnn, ncnn, shrknode
use fluid_variables, only: maxconn
implicit none
include 'mpif.h'

real(8) x_solid(nsd,nn_solid)
integer nsd
integer nn_solid
real(8) xna(nsd,nn)
integer nn
real(8) dwjp(nn)
real(8) sp_radius(nsd,nn)
!-------------------------------------
real(8) tmp_solid
real(8) tmp_ncpus
integer npart
integer i
integer j
real(8) x(nsd)
integer ninf
integer inf(maxconn)
integer :: maxinf,mininf,totinf
integer :: maxinf_r,mininf_r,totinf_r
real(8) :: b(4), bd(nsd,4)
real(8) shp
integer error_id
integer n
integer nnum
integer isd
real(8) y(nsd)
real(8) a(nsd)
real(8) avginf
integer cnn_r(maxconn,nn_solid)
integer ncnn_r(nn_solid)
real(8) shrknode_r(maxconn,nn_solid)
real(8) tmp(maxconn,nn_solid)


!------------------------------------

if (myid == 0) write(*,*) '*** Calculate RKPM coefficient ***'

  if (allocated(shrknode)) then
     deallocate(shrknode)
  end if
  if (allocated(cnn)) then
     deallocate(cnn)
  end if
  if (allocated(ncnn)) then
     deallocate(ncnn)
  end if

  allocate(shrknode(maxconn,nn_solid),stat=error_id)
  allocate(cnn(maxconn,nn_solid)     ,stat=error_id)
  allocate(ncnn(nn_solid)            ,stat=error_id)


  maxinf = 0
  mininf = 9999
  avginf = 0
  cnn(:,:)=0
  ncnn(:)=0
  shrknode(:,:)=0.0d0
  totinf=0

tmp_solid=dble(nn_solid)
tmp_ncpus=dble(ncpus)


npart=ceiling(tmp_solid/tmp_ncpus)


do j=1,npart
	i=(j-1)*ncpus+myid+1
	if (i .le. nn_solid) then
		     x(1:nsd)=x_solid(1:nsd,i) !get solids point coordinate

		     ninf=0
		     inf(:)=0
! get a list of influence nodes from the fluids grid
		     call getinf_rkpm(inf,ninf,x,xna,sp_radius,nn,nsd,maxconn)
		    cnn(1:ninf,i)=inf(1:ninf)
		     ncnn(i)=ninf
!write(*,*) 'myid', myid, 'node index', i
!write(*,*) 'cnn', cnn(1:ninf,i)
		    	 if (ninf > maxinf) then
               		         maxinf = ninf
 		  	  elseif (ninf < mininf) then
                	        mininf = ninf
		  	 endif
	   	  totinf = totinf + ninf
! calculate the correction function
    		 if (nsd==2) then
        	        call correct2d(b,bd,x,xna,sp_radius,dwjp,nn,1,inf,ninf,maxconn)
   		  elseif (nsd==3) then
        	        call correct3d(b,bd,x,xna,sp_radius,dwjp,nn,1,inf,ninf,maxconn)
  	 	 endif


		     do n = 1, ninf
		        nnum = inf(n)
		        do isd = 1,nsd
		           y(isd) = xna(isd,nnum)
		           a(isd) = sp_radius(isd,nnum)
		        enddo

	                 if (nsd==2) then
        	      	  call RKPMshape2d(shp,b,bd,x,y,a,dwjp(nnum))
                	 elseif (nsd==3) then
                          call RKPMshape3d(shp,b,bd,x,y,a,dwjp(nnum))
               		 endif
    		    shrknode(n,i)=shp
	            enddo
	end if
enddo
maxinf_r=0
mininf_r=0
totinf_r=0
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(maxinf,maxinf_r,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
call mpi_allreduce(mininf,mininf_r,1,mpi_integer,mpi_min,mpi_comm_world,ierror)
call mpi_allreduce(totinf,totinf_r,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)

  avginf = totinf_r/nn_solid
if (myid == 0) then
  write(6,'("  Maximum Influence Nodes = ",i7)') maxinf_r
  write(6,'("  Minimum Influence Nodes = ",i7)') mininf_r
  write(6,'("  Average Influence Nodes = ",f7.2)') avginf
end if

ncnn_r(:)=0
cnn_r(:,:)=0
shrknode_r(:,:)=0.0d0

!tmp(:,:)=1.0


call mpi_allreduce(ncnn(1),ncnn_r(1),nn_solid,mpi_integer,mpi_sum,mpi_comm_world,ierror)
call mpi_allreduce(cnn(1,1),cnn_r(1,1),nn_solid*maxconn,mpi_integer,mpi_sum,mpi_comm_world,ierror)
call mpi_allreduce(shrknode(1,1),shrknode_r(1,1),nn_solid*maxconn,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

!write(*,*) 'cnn', cnn(:,1)


ncnn(:)=ncnn_r(:)
cnn(:,:)=cnn_r(:,:)
shrknode(:,:)=shrknode_r(:,:)

!if (myid == 0) then
!open(unit=8466, file='rkpmsh_pa.txt', status='unknown')
!do isd=1,nn_solid
!write(8466,*) 'node', isd, 'rkcoeff',shrknode(:,isd)
!write(8466,*) 'ninf', cnn(:,isd)
!write(8466,*) 'should be 1', sum(shrknode(:,isd))
!end do
!close(8466)
!end if


return
end 
