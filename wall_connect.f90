
subroutine wall_connect(x_inter,I_fluid_center,nn_con_ele,con_ele,infdomain_inter)

  use fluid_variables,only:nsd,nen,nn,ne
  use interface_variables,only:nn_center,maxmatrix,nn_inter,ele_refine
  use allocate_variables,only:inter_ele,ne_inter
  use mpi_variables
  include 'mpif.h'

  real(8) x_inter(nsd,maxmatrix)
  real(8) I_fluid_center(nn_center)
  integer nn_con_ele,con_ele(nn_con_ele)
  integer flag_con(nn_con_ele),infdomain_inter(maxmatrix)
  integer nn_elelist, elelist(nn_con_ele)

  integer i,j,icount,jcount,isd,inl,node
  integer nn_inter_tmp,flag
  real(8) x_inter_tmp(nsd,maxmatrix)
  integer nump

  nump=ele_refine**nsd

  flag_con(:)=0
  nn_elelist=0

  do i=1,ne_inter
     do j=1,nn_con_ele
        if(inter_ele(i)==con_ele(j))flag_con(j)=1
     end do
  end do
! set the element flag = 1 for the element that has interfacial points in it

  do i=3,nn_con_ele-2
     if(flag_con(i-2)==1.and.flag_con(i+2)==1.and.flag_con(i)==1) then
 if(myid==0)write(*,*)'element to change=',con_ele(i)
        do inl=1,nump
           I_fluid_center(nn_center-nump*nn_con_ele+nump*(i-1)+inl)=1.0
I_fluid_center(nn_center-nump*nn_con_ele+nump*(i-2)+inl)=0.6
I_fluid_center(nn_center-nump*nn_con_ele+nump*(i)+inl)=0.6
        end do
     end if
     if(flag_con(i-1)==1.and.flag_con(i+1)==1.and.flag_con(i)==1) then
        nn_elelist=nn_elelist+1
        elelist(nn_elelist)=con_ele(i)
     end if
  end do

! find out the element id : the adjacent elements all contain interfacial
! elements
  nn_inter_tmp=0
  do i=1,nn_inter
     flag=1
     do j=1,nn_elelist
        if(infdomain_inter(i)==elelist(j))flag=0
     end do
     if(flag==1) then
        nn_inter_tmp=nn_inter_tmp+1
        x_inter_tmp(:,nn_inter_tmp)=x_inter(:,i)
     end if
  end do
  nn_inter=nn_inter_tmp
  x_inter(:,1:nn_inter)=x_inter_tmp(:,1:nn_inter)

end subroutine wall_connect



