

subroutine get_Ia_ls(x_inter,x_center,infdomain_inter,ne_inter,inter_ele,hg,I_fluid_center,Ic_inter)

  use interface_variables
  use fluid_variables,only:nsd,ne
  use mpi_variables

  real(8) x_inter(nsd,maxmatrix)
  real(8) x_center(nsd,ne)
  integer infdomain_inter(maxmatrix)
  integer ne_inter
  integer inter_ele(ne)
  real(8) hg(ne)
  real(8) I_fluid_center(ne)
  real(8) Ic_inter

  integer i,j,icount,jcount
  real(8) A(nn_inter,ne_inter)
  real(8) B(nn_inter)
  real(8) dx(nsd),hs,Sp
  integer INFO
  character(1) TRANS
  real(8) WORK(2*ne_inter)
  integer LWORK

  A(:,:)=0.0
  B(:)=0.0
  TRANS='N'
  LWORK=2*ne_inter

  do i=1,ne_inter
     I_fluid_center(inter_ele(i))=0.0
  end do

  do i=1,nn_inter
     do j=1,ne_inter
	hs=hg(inter_ele(j))
	dx(1:nsd)=abs(x_inter(1:nsd,i)-x_center(1:nsd,inter_ele(j)))
	call B_Spline(dx,hs,nsd,Sp)
	A(i,j)=Sp
     end do
     do j=1,ne
	hs=hg(j)
	dx(1:nsd)=abs(x_inter(1:nsd,i)-x_center(1:nsd,j))
	call B_Spline(dx,hs,nsd,Sp)
	B(i)=B(i)+I_fluid_center(j)*Sp
     end do
     B(i)=Ic_inter-B(i)
  end do

  call DGELS(TRANS,nn_inter,ne_inter,1,A,nn_inter,B,nn_inter,WORK,LWORK,INFO)

  do i=1,ne_inter
     I_fluid_center(inter_ele(i))=B(i)
if(myid==0) then
  write(*,*)'I_fluid_center',inter_ele(i),B(i)
end if
     if(B(i).gt.1.0) then
	I_fluid_center(inter_ele(i))=1.0
     else if(B(i).lt.0.0) then
	I_fluid_center(inter_ele(i))=0.0
!     else
!	I_fluid_center(inter_ele(i))=Ic_inter
     end if
  end do

end subroutine get_Ia_ls














