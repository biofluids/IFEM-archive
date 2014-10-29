subroutine solve_laplace(source,nsd,nn,nn_solid,ien,ne,nen,x_fluid,node_sbc,nn_sbc,I_fluid,flag_fnode)
! 1 decide interior, outer, boundary fluid nodes
! 2 set id matrix for laplace equation
! 3 sovle for laplace equation
use delta_nonuniform, only: cnn, ncnn
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
integer, allocatable ::  lp_el(:)
integer node
real(8) dg(nn)
integer lp_id(nn)
real(8) p_inter(nn)
real(8) w_inter(nn)
integer icount

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
	end do
end do
!------------------------
! reset inside boundary based on fnode
!do i=1,nn
!        if (flag_fnode(i) == 1) then
!                I_fluid(i)=1.0
!                flag_node(i)=1
!        end if
!end do


count_el=0
do i=1,ne
	do j=1,nen
		node=ien(j,i)
		if (flag_node(node) .eq. 2) then
		lp_id(node)=1  ! set id vector for laplace eq
		flag_el(i)=1 ! elements need to be looped
		else
		w_inter(node)=1.0 ! prevent diagonal preconditioner to be zero
		end if
	end do
	if (flag_el(i) == 1) then
		count_el=count_el+1
	end if
end do

allocate(lp_el(count_el))
!do i=1,nn
!write(*,*) 'id flag_node', i, flag_node(i)
!end do 
count_el=0


do i=1,ne
	if (flag_el(i) == 1) then
		count_el=count_el+1
		lp_el(count_el)=i ! element index needs to be looped
	end if
end do


!open(unit=8406, file='lpsource.txt', status='new')
!do i=1,nn_solid
!j=ncnn(i)
!write(8406,*) 'cnn', i,j,'aaa', cnn(1:J,i)
!end do
!close(8406)
!source(:,:)=0.0

call block_Laplace(x_fluid,I_fluid,p_inter,w_inter,ien,lp_el,count_el,source)
continue
call setid(p_inter,lp_id,1)

!open(unit=8406, file='source.txt', status='unknown')
!do i=1,nn
!j=ncnn(i)
!write(8406,*) 'node', i ,'lp_id', lp_id(i), 'source', source(:,i)
!end do
!close(8406)

!write(*,*) 'w_inter', minval( w_inter(:)), maxval(w_inter(:))


call gmres_Laplace(x_fluid,I_fluid,w_inter,p_inter,dg,ien,lp_id,lp_el,count_el)

!do i=1,nn
	!write(*,*) 'dg', maxval(dg(:)), minval(dg(:))
!end do

dg(:) = abs(dg(:)) / maxval(dg(:))


    I_fluid(:)=dg(:)+I_fluid(:)
return

deallocate(lp_el)
end 
