!
!!             ifnter indicator coupling with correction
!!============================================

subroutine center_indicator_update(x_center,x_inter,I_fluid_center,corr_Ip,infdomain_inter,hg)

  use interface_variables
  use fluid_variables, only:ne,nsd
  use allocate_variables, only: ne_regen_ele,regen_ele
  use mpi_variables

  real(8) x_center(nsd,ne)
  real(8) x_inter(nsd,maxmatrix)
  real(8) I_fluid_center(ne)
  real(8) corr_Ip(maxmatrix)
  integer infdomain_inter(maxmatrix)
  real(8) hg(ne)

  integer i,j,ie,icount
  real(8) dx(nsd),hs,Sp
  real(8) I_fluid_center_temp(ne_regen_ele)
  real(8) err_center

  do icount=1,ne_regen_ele
     ie=regen_ele(icount)
     I_fluid_center_temp(icount)=0.0
     hs=hg(infdomain_inter(i))
     do j=1,nn_inter
        dx(:)=abs(x_center(:,ie)-x_inter(:,j))
        call B_Spline(dx,hs,nsd,Sp)
        I_fluid_center_temp(icount)=I_fluid_center_temp(icount)+corr_Ip(j)*Sp
     end do

     do j=1,ne
        dx(:)=abs(x_center(:,ie)-x_center(:,j))
        call B_Spline(dx,hs,nsd,Sp)
        I_fluid_center_temp(icount)=I_fluid_center_temp(icount)+I_fluid_center(j)*Sp
     end do
  end do
  err_center=0.0

  do icount=1,ne_regen_ele
     ie=regen_ele(icount)
     if(abs(I_fluid_center(ie)-I_fluid_center_temp(icount)).gt.err_center) then
        err_center=abs(I_fluid_center(ie)-I_fluid_center_temp(icount))
     end if
     I_fluid_center(ie)=I_fluid_center_temp(icount)
  end do
if(myid==0) write(*,*)'max err for center indicator=',err_center
end subroutine center_indicator_update
