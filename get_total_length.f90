!==========================!
! get total length         !
!==========================!

subroutine get_total_length(x_inter,infdomain,hg)

  use fluid_variables, only:nsd,ne
  use interface_variables
  implicit none
  real(8) x_inter(nsd,nn_inter)
  integer infdomain(maxmatrix)
  real(8) hg(ne)

  real(8) hs, Bsum,dx(nsd),temp
  integer i,j,isd
  real(8) Sp
  total_length=0.0
  do i=1,nn_inter
     hs=hg(infdomain(i))
     Bsum=0.0
     do j=1,nn_inter
        dx(:)=abs(x_inter(:,i)-x_inter(:,j))
        temp=0.0
        do isd=1,nsd
           temp=temp+dx(isd)**2
        end do
        temp=sqrt(temp)
        call B_Spline_0order(temp,hs,Sp)
        Bsum=Bsum+Sp
     end do
     total_length=total_length+hs/Bsum
  end do

end subroutine get_total_length
