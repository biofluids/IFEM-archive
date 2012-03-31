!!!!!!!!!!!!!!!Get interface elements!!!!!!!!!!!!!!


subroutine find_inter(inf_inter,ien,nn_inter)

  use interface_variables,only:maxmatrix
  use allocate_variables,only:ne_inter,inter_ele
  use fluid_variables, only:nn,ne,nen
  use mpi_variables
  include 'mpif.h'
  integer inf_inter(maxmatrix),ien(nen,ne),nn_inter
  integer i,j,icount,jcount,flag,ie,inl,node
  integer ncount
  integer temp_ele(maxmatrix)
  integer flag_node(nn)
  integer base,top
!  real(8) I_fluid_center(ne),I_fluid(nn)
!  integer its


! find the interface element that contains the interface points
  ncount = 0
  temp_ele(:) = 0
  do i=1,nn_inter
     do j=1,ncount
	if (inf_inter(i) == temp_ele(j)) then
	   go to 100
	end if
     end do

	ncount = ncount + 1
	temp_ele(ncount) = inf_inter(i)

100 continue
  end do
!=======================
  if(ncount.ge.maxmatrix) then
    write(*,*)'enlarge maxatrix in get_inter_ele'
    stop
  end if
!=========================
  ne_inter = ncount
  if(allocated(inter_ele)) then
     deallocate(inter_ele)
  end if
  allocate(inter_ele(ne_inter))
  inter_ele(1:ne_inter)=temp_ele(1:ne_inter)


if(myid==0) then
write(*,*)'ne_inter=',ne_inter
end if

call mpi_barrier(mpi_comm_world,ierror)


end subroutine find_inter









