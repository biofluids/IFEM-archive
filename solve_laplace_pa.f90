subroutine solve_laplace_pa(source,nsd,nn,nn_solid,ien,ne,nen,x_fluid,node_sbc,nn_sbc,I_fluid,flag_fnode,&
		ne_local,ien_local,nn_local,node_local,send_address,ad_length,&
		global_com,nn_global_com,local_com,nn_local_com)
! 1 decide interior, outer, boundary fluid nodes
! 2 set id matrix for laplace equation
! 3 sovle for laplace equation
use delta_nonuniform, only: cnn, ncnn
use mpi_variables
implicit none
include 'mpif.h'
real(8) source(nsd,nn)
integer nsd
integer nn
integer nn_solid
integer ien(nen,ne)
integer ne
integer nen
integer node_sbc(nn_sbc)
integer nn_sbc
real(8) x_fluid(nsd,nn)
real(8) I_fluid(nn)
integer flag_fnode(nn)
!----------------------
integer flag_node(nn)
integer flag_el(ne)
integer i
integer j
integer inn
integer count_el
integer node
real(8) dg(nn)
integer lp_id(nn)
real(8) p_inter(nn)
real(8) w_inter(nn)
integer icount
real(8) time
!--------------------
integer ne_local
integer ien_local(ne_local)
integer nn_local
integer node_local(nn_local)
integer ad_length
integer send_address(ad_length,2)
integer nn_global_com
integer global_com(nn_global_com)  ! global node index for communication
integer nn_local_com
integer local_com(nn_local_com)  ! local index in the communication region on each processor

flag_node(:)=0 ! set out node
flag_el(:)=0
p_inter(:)=0.0
w_inter(:)=0.0
dg(:)=0.0

lp_id(:)=0
I_fluid(:)=0.0

do i=1,nn_solid
	do j=1,ncnn(i)
		node=cnn(j,i)
		flag_node(node)=1 ! set interior node
		I_fluid(node)=1.0
	end do
end do


do i=1,nn_sbc
	inn=node_sbc(i)
	do j=1,ncnn(inn)
		node=cnn(j,inn)
		flag_node(node)=2 ! set boundary node
		I_fluid(node)=0.0
		lp_id(node) = 1 ! set id vector for laplace eq
	end do
end do
!------------------------
! reset inside boundary based on fnode
!do i=1,nn
!        if (flag_fnode(i) == 1) then
!                I_fluid(i)=1.0
!                flag_node(i)=1
!		lp_id(i) = 0
!        end if
!end do

time = mpi_wtime()
call block_Laplace(x_fluid,I_fluid,p_inter,w_inter,ien,ien_local,ne_local,source)

! commute residual for each processor
call communicate_res_ad_sub(p_inter,1,nn,send_address,ad_length)
call communicate_res_ad_sub(w_inter,1,nn,send_address,ad_length)


call setid_pa(p_inter,1,nn,lp_id,node_local,nn_local)
time = mpi_wtime() - time
if (myid == 0) write(*,*) 'Laplace block time', time


call gmres_Laplace_pa(x_fluid,I_fluid,w_inter,p_inter,dg,ien,lp_id,&
	ne_local,ien_local,nn_local,node_local,send_address,ad_length,&
	global_com,nn_global_com,local_com,nn_local_com)

!do i=1,nn


if (myid == 0)	write(*,*) 'dg', maxval(dg(:)), minval(dg(:))
!end do

dg(:) = abs(dg(:)) / maxval(dg(:))

    I_fluid(:)=dg(:)+I_fluid(:)

return

end 
