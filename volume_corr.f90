!==============================!
! do volume correction         !
!==============================!

subroutine volume_corr(x_inter_ini,x_inter,arc_inter,norm_inter)

  use interface_variables  !include vol_corr and total_length
  use fluid_variables, only:nsd
  use mpi_variables
  implicit none

  real(8) x_inter_ini(nsd,maxmatrix)
  real(8) x_inter(nsd,maxmatrix)
  real(8) arc_inter(maxmatrix)
  real(8) norm_inter(nsd,maxmatrix)

  integer icount
  real(8) x0(nsd),x1(nsd),x(nsd) !corr for initial , after update, after correction
  real(8) cons     !volchange/total length
  real(8) norm(nsd)   !normal

  vol_corr=0.0
  do icount=1,nn_inter
       vol_corr=vol_corr+arc_inter(icount)*((x_inter(1,icount)-x_inter_ini(1,icount))*norm_inter(1,icount)+&
                                            (x_inter(2,icount)-x_inter_ini(2,icount))*norm_inter(2,icount))
  end do
  if(myid==0)write(*,*)'vol_corr=',vol_corr

  cons=vol_corr/total_length

  do icount=1,nn_inter
     x0(1:nsd)=x_inter_ini(1:nsd,icount)
     x1(1:nsd)=x_inter(1:nsd,icount)
     norm(1:nsd)=norm_inter(1:nsd,icount)
     if(nsd==2) then
	x(1)=-(cons*x0(1)-cons*x1(1)+norm(1)*x1(1)**2-norm(1)*x0(1)*x1(1)- &
		norm(2)*x1(1)*x0(2)+norm(2)*x1(1)*x1(2))/ &
		(norm(1)*x0(1)-norm(1)*x1(1)+norm(2)*x0(2)-norm(2)*x1(2))
	x(2)=-(cons*x0(2)-cons*x1(2)+norm(2)*x1(2)**2-norm(1)*x0(1)*x1(2)+ &
		norm(1)*x1(1)*x1(2)-norm(2)*x0(2)*x1(2))/ &
		(norm(1)*x0(1)-norm(1)*x1(1)+norm(2)*x0(2)-norm(2)*x1(2))
     end if
     x_inter(1:nsd,icount)=x(1:nsd)
  end do

end subroutine










































