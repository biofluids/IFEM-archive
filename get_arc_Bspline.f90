!=======================================================
!calculate the arclength for each interfacial points
!=======================================================

subroutine get_arc_Bspline(arc_total,x_inter,arc_inter,infdomain_inter,hg)

  use fluid_variables,only:nsd,nn,ne,nen
  use interface_variables
  use mpi_variables
  real(8) x_inter(nsd,maxmatrix)       !coordinates for interfacial points
  real(8) arc_inter(maxmatrix)         !arclength for interfacial points
  integer infdomain_inter(maxmatrix)
  real(8) hg(ne)
  real(8) arc_total

  integer i,j,flag
  real(8) hs,Sp,dx(nsd)
  real(8) Bsum                         !summation of Bspline function
  real(8) temp
  real(8) total_length

  real(8) A(nn_inter,nn_inter),B(nn_inter)
  integer IPIV(nn_inter),INFO
  flag=1
  A(:,:)=0.0
  if(flag==2) then
  do i=1,nn_inter
     B(i)=arc_inter(i)
     hs=arc_inter(i)
     do j=1,nn_inter
        dx(:)=abs(x_inter(:,j)-x_inter(:,i))
        temp=sqrt(dx(1)**2+dx(2)**2)
	call B_Spline_0order(temp,hs,Sp)
	A(i,j)=Sp
     end do
  end do

  call DGESV(nn_inter,1,A,nn_inter,IPIV,B,nn_inter,INFO)

  arc_inter(1:nn_inter)=B(1:nn_inter)
  total_length=0.0
  do i=1,nn_inter
     total_length=total_length+arc_inter(i)
  end do
  end if
  
  if(flag==1) then
  total_length=0.0
  if(nsd==2) then
    do i=1,nn_inter
       Bsum=0.0
!       hs=hg(infdomain_inter(i))/3
        hs=arc_inter(i)
!        hs=0.0005
!         num=0
       do j=1,nn_inter
	  dx(:)=abs(x_inter(:,j)-x_inter(:,i))
	  temp=sqrt(dx(1)**2+dx(2)**2)
	  if(temp.le.3*hs)then
!	    num=num+1
	  end if
	  call B_Spline_0order(temp,hs,Sp)
	  Bsum=Bsum+Sp
       end do
!       arc_exact(i)=6*hs/num
       arc_inter(i)=hs/Bsum
       total_length=total_length+arc_inter(i)
    end do
  else if(nsd==3) then
    write(*,*)'no 3D right now'
    stop
  end if
  end if
if(myid==0) then
write(*,*)'total_length=',total_length
end if
  arc_total=total_length
end subroutine get_arc_Bspline





