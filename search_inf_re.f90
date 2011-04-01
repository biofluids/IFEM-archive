subroutine search_inf_re(xyz_solid, xyz_fluid, nn_fluid,nn_solid, nsd, ne_solid, nen_solid, ien_solid,&
			 flag_fnode,node_local,nn_local)
! Seach all the fluid nodes overlapping with solid domain and give flag_node =1 
use mpi_variables
implicit none
include 'mpif.h'
integer nn_fluid
integer nn_solid
integer nsd
integer ne_solid
integer nen_solid
integer ien_solid(ne_solid,nen_solid)
integer flag_fnode(nn_fluid)
real(8) xyz_solid(nsd,nn_solid)
real(8) xyz_fluid(nsd,nn_fluid)
integer node_local(nn_local)
integer nn_local
integer i
integer finf
real(8) x(nsd)
integer flag_fnode_s(nn_fluid)
integer inn
integer maxconn
integer ien_solid_tmp(nen_solid,ne_solid)

maxconn=30
flag_fnode(:)=0
flag_fnode_s(:)=0

do i=1,ne_solid
	ien_solid_tmp(:,i)=ien_solid(i,:)
end do

do i=1,nn_local
	inn=node_local(i)
	x=xyz_fluid(1:nsd,inn)
	finf=0
	call getinf_el_3d(finf,x,xyz_solid,nn_solid,nsd,ne_solid,nen_solid,ien_solid_tmp,maxconn)
	if (finf .ne. 0) then
		flag_fnode_s(inn)=1
	end if
end do


call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(flag_fnode_s(1),flag_fnode(1),nn_fluid,mpi_integer,mpi_sum,mpi_comm_world,ierror)
   
      return
end
