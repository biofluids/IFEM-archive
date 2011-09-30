
subroutine blockgmres_correction(x_inter,RW,nn_inter,dv,avloc,nn_inter_loc,nsd)

  use mpi_variables
  use interface_variables, only:maxmatrix,hsp

  real(8) x_inter(nsd,maxmatrix)
  real(8) RW(nsd+1,nn_inter)
  integer nn_inter,nn_inter_loc,nsd
  real(8) dv(nn_inter),avloc(nn_inter)

  integer i,j,icount,jcount,INFO
  real(8) dx(nsd),Sp,temp
  real(8) vec(nsd+1)

  integer loc_index
  vec(1)=1.0
  do loc_index=1,nn_inter_loc
     i=myid+1+(loc_index-1)*ncpus
     do j=1,nn_inter
!	dx(:)=abs(x_inter(:,i)-x_inter(:,j))
	vec(2:nsd+1)=x_inter(:,i)-x_inter(:,j)
	call B_Spline1(abs(vec(2:nsd+1)),hsp,nsd,Sp,INFO)
	if(INFO==1)then
!	   vec(2:nsd+1)=x_inter(:,i)-x_inter(:,j)
	   temp=0.0
	   do icount=1,nsd+1
	      temp=temp+vec(icount)*RW(icount,i)
	   end do
	   avloc(i)=avloc(i)+Sp*temp*dv(j)
	end if
     end do
  end do
end subroutine blockgmres_correction
