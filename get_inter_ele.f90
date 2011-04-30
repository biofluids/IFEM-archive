!!!!!!!!!!!!!!!Get interface elements!!!!!!!!!!!!!!


subroutine get_inter_ele(inf_inter,inf_den,ien)!,&
!			I_fluid_center,I_fluid,its)

  use interface_variables,only:nn_inter,maxmatrix
  use allocate_variables
  use fluid_variables, only:nn,ne,nen
  use mpi_variables
  include 'mpif.h'
  integer inf_inter(maxmatrix),inf_den(maxmatrix),ien(nen,ne)
  integer i,j,icount,jcount,flag,ie,inl,node
  integer ncount
  integer temp_ele(maxmatrix)
  integer flag_node(nn)
  integer base,top
!  real(8) I_fluid_center(ne),I_fluid(nn)
!  integer its

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

  ncount=0
  temp_ele(:)=0
  do i=1,nn_inter
     do j=1,ncount
	if(inf_den(i)==temp_ele(j)) then
	   goto 200
	end if
     end do
     
     ncount=ncount+1
     temp_ele(ncount)=inf_den(i)
200 continue
   end do
!===========================
  if(ncount.gt.maxmatrix) then
    write(*,*)'enlarge maxamtrix in get_inter_ele'
    stop
   end if
!============================
   ne_inter_den=ncount
  if(allocated(inter_ele_den)) then
    deallocate(inter_ele_den)
  end if
  allocate(inter_ele_den(ne_inter_den))
  inter_ele_den(1:ne_inter_den)=temp_ele(1:ne_inter_den)
if(myid==0) then
write(*,*)'ne_inter=',ne_inter
write(*,*)'ne_inter_den=',ne_inter_den
end if

!=====================================
  flag_node(:)=0
  do icount=1,ne_inter
     ie=inter_ele(icount)
     do inl=1,nen
	node=ien(inl,ie)
	flag_node(node)=1
     end do
  end do
  ncount=0
  do ie=1,ne
     flag=0
     do inl=1,nen
	node=ien(inl,ie)
	if(flag_node(node)==1) then
	  flag=1
	  goto 987
	end if
     end do
987 continue
    if(flag==1) then
	ncount=ncount+1
	temp_ele(ncount)=ie
    end if
  end do
  ne_regen_ele=ncount
  if(allocated(regen_ele)) then
     deallocate(regen_ele)
  end if
  allocate(regen_ele(ne_regen_ele))
  regen_ele(1:ne_regen_ele)=temp_ele(1:ne_regen_ele)
if(myid==0) then
  write(*,*)'ne_regen_ele=',ne_regen_ele
end if

!========================================================
!allocate local regen element for each processor
!========================================================
if(ne_regen_ele.le.ncpus) then
  if(allocated(regen_ele_loc)) then
    deallocate(regen_ele_loc)
  end if
  allocate(regen_ele_loc(1))
  if(myid+1.le.ne_regen_ele) then
    ne_regen_ele_loc=1
    regen_ele_loc(1)=regen_ele(myid+1)
  else
    ne_regen_ele_loc=0 
    regen_ele_loc(1)=0
  end if
else
  base=floor(real(ne_regen_ele)/real(ncpus))
  top=ne_regen_ele-base*ncpus
  if(myid+1.le.top) then
     ne_regen_ele_loc=base+1
  else
     ne_regen_ele_loc=base
  end if
  if(allocated(regen_ele_loc)) then
    deallocate(regen_ele_loc)
  end if
  allocate(regen_ele_loc(ne_regen_ele_loc))
  do icount=1,ne_regen_ele_loc
     regen_ele_loc(icount)=regen_ele(myid+1+(icount-1)*ncpus)
  end do
end if

!=========
call mpi_barrier(mpi_comm_world,ierror)


end subroutine get_inter_ele









