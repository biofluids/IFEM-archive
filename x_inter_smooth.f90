!========================================

!smooth x_inter using polynomial

!==========================================

subroutine x_smooth(x_inter,x_inter_smooth,infdomain_inter,hg)

  use interface_variables
  use fluid_variables,only:nsd,ne

  real(8) x_inter(nsd,maxmatrix),x_inter_smooth(nsd,maxmatrix)
  integer infdomain_inter(maxmatrix)
  real(8) hg(ne)

  integer i,j,icount,isd
  integer inter_index(nn_inter,50)
  integer inter_num(nn_inter)
  real(8) hd   !distance
  real(8) hs   !critirial
  real(8) A(50,4),B(50)
  integer row
  character(1) TRANS
  real(8) workls(400)
  integer INFO
  real(8) a0,a1,a2,a3,b2,b3,c
  real(8) err
!  B(:)=1.0
  TRANS='N'
  do i=1,nn_inter
     hs=hg(infdomain_inter(i))/2.0
     row=0
     B(:)=1.0
     do j=1,nn_inter
	hd=sqrt((x_inter(1,i)-x_inter(1,j))**2+(x_inter(2,i)-x-inter(2,j))**2)
        if(hd.le.hs) then
	  row=row+1
	  A(row,1)=x_inter(1,j)**3
	  A(row,2)=x_inter(1,j)**2*x_inter(2,j)
	  A(row,3)=x_inter(1,j)*x_inter(2,j)**2
	  A(row,4)=x_inter(2,j)**3
	end if
     end do

     call DGELS(TRANS,row,4,1,A(1:row,1:4),row,B(1:row),row,workls,400,INFO)

     a0=B(1)
     a1=B(2)*x_inter(2,i)
     a2=B(3)*x_inter(2,i)**2
     a3=B(4)*x_inter(2,i)**3-1.0
     b2=3.0*a0*a2-a1**2
     b3=27.0*(a0**2)*a3-9.0*a0*a1*a2+2.0*(a1**3)
     c=((-b3+(b3**2+4.0*b2**3)**0.5)/2.0)**(1.0/3.0)
     x_inter_smooth(1,i)=1.0/(3.0*a0)*(-a1+c-b2/c)
     x_inter_smooth(2,i)=x_inter(2,i)

     err=abs(x_inter(1,i)-x_inter_smooth(1,i))
    write(*,*)'err_smooth=',err
  end do

  end subroutine x_smooth



