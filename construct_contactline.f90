

subroutine construct_contactline(x,x_inter,x_center,I_fluid,I_fluid_center,hg,corr_Ip,ien,&
			nn_con_ele,con_ele,norm_con_ele,rngface,ne_intlocal,ien_intlocal,vol_nn,d,flag_contact)

  use fluid_variables, only:nsd,nn,ne,nen,neface,ndf
  use interface_variables
!  use allocate_variables,only:ne_inter
  use mpi_variables
  include 'mpif.h'


  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center),I_fluid(nn),I_fluid_center(nn_center),hg(ne)
  real(8) corr_Ip(maxmatrix),vol_nn(nn)
  integer ien(nen,ne)

  integer nn_con_ele,con_ele(nn_con_ele)
  real(8) norm_con_ele(nsd,nn_con_ele)
  integer rngface(neface,ne)

  integer ne_intlocal,ien_intlocal(ne_intlocal)

  integer nn_inter_regen
  real(8) x_inter_regen(nsd,maxmatrix)
  real(8) d(ndf,nn),vel_fluid(nsd,nn)
  integer infdomain_inter(maxmatrix)
 
  integer i,j,icount,jcount

  integer flag_contact,flag

  integer nn_inter_regen_temp
  real,dimension(:,:),allocatable :: x_inter_regen_temp

  real(8) temp
  vel_fluid(1:nsd,:)=d(1:nsd,:)

flag_contact=1
  call points_on_contact_sur(x,x_inter,x_center,I_fluid,I_fluid_center,hg,corr_Ip,ien,&
                      nn_con_ele,con_ele,norm_con_ele,rngface,nn_inter_regen,x_inter_regen)



  if(allocated(x_inter_regen_temp)) then
    deallocate(x_inter_regen_temp)
  end if
  allocate(x_inter_regen_temp(nsd,nn_inter_regen))

  nn_inter_regen_temp=nn_inter_regen
  x_inter_regen_temp(:,1:nn_inter_regen)=x_inter_regen(:,1:nn_inter_regen)

  nn_inter_regen=0

  do i=1,nn_inter_regen_temp
     flag=1
     do j=1,nn_inter_regen
        temp=0.0
        do icount=1,nsd
           temp=temp+(x_inter_regen_temp(icount,i)-x_inter_regen(icount,j))**2
        end do
        if(sqrt(temp).lt.max_hg/4.0) flag=0
     end do
     if(flag==1) then
       nn_inter_regen=nn_inter_regen+1
       x_inter_regen(:,nn_inter_regen)=x_inter_regen_temp(:,i)
     end if
  end do


!open(10,file='pointsoncontact.dat',status='unknown')
!if(myid==0) then
!do j=1,nn_inter_regen
!   write(10,*)x_inter_regen(:,j)
!end do
!end if
!close(10)


if(nn_inter_regen==0) flag_contact=0
if(flag_contact==0) goto 222
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,ne_intlocal,ien_intloal)

!  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)


  call generate_contactline(x_inter_regen,nn_inter_regen,x,x_inter,x_center,hg,I_fluid_center,&
			corr_Ip,ien,nn_con_ele,con_ele,norm_con_ele,vel_fluid,&
				vol_nn,infdomain_inter)

222 continue

endsubroutine construct_contactline
