!========================================

!smooth x_inter using polynomial

!==========================================

subroutine x_smooth(x_inter,infdomain_inter,hg)

  use interface_variables
  use fluid_variables,only:nsd,ne

  real(8) x_inter(nsd,maxmatrix),x_inter_smooth(nsd,maxmatrix)
  integer infdomain_inter(maxmatrix)
  real(8) hg(ne)

  integer i,j,icount,isd,ncount,inter_count
  integer inter_index(nn_inter,50)
  integer inter_num(nn_inter)
  real(8) hd   !distance
  real(8) hs   !critirial
  real(8) hint
  real(8) A(50,4),B(50)
  integer row
  character(1) TRANS
  real(8) workls(400)
  integer INFO
  real(8) a0,a1,a2,a3,b2,b3,cc
  real(8) err
  integer ecount
!  B(:)=1.0
  TRANS='N'
  inter_count=0
  ecount=0
  do i=1,nn_inter
     hint=hg(infdomain_inter(i))/30.0
     hs=hint
100 continue
     ncount=0
     do j=1,nn_inter
        hd=sqrt((x_inter(1,i)-x_inter(1,j))**2+(x_inter(2,i)-x_inter(2,j))**2)
	if(hd.le.hs) then
	  ncount=ncount+1
	  inter_index(i,ncount)=j
	end if
     end do
     if(ncount.lt.7) then
	hs=hs+hint
        go to 100
     end if
    inter_num(i)=ncount
  end do
!write(*,*)'inter_num(:)',inter_num(:)



  do i=1,nn_inter
!  do i=17,17
     row=inter_num(i)
!     write(*,*)'i=',i,'inter_index=',inter_index(i,1:inter_num(i))
     do j=1,inter_num(i)
	  A(j,1)=x_inter(1,inter_index(i,j))**3
	  A(j,2)=x_inter(1,inter_index(i,j))**2
	  A(j,3)=x_inter(1,inter_index(i,j))
	  A(j,4)=1.0
	  B(j)=x_inter(2,inter_index(i,j))
     end do

     call DGELS(TRANS,row,4,1,A(1:row,1:4),row,B(1:row),row,workls,400,INFO)
     x_inter_smooth(1,i)=x_inter(1,i)
     x_inter_smooth(2,i)=B(1)*x_inter(1,i)**3+B(2)*x_inter(1,i)**2+B(3)*x_inter(1,i)+B(4)
     err=abs(x_inter_smooth(2,i)-x_inter(2,i))


!    write(*,*)'err_smooth=',err,'row=',row,INFO
    if(err.lt.1.0e-3)then
      inter_count=inter_count+1
      x_inter(1:nsd,inter_count)=x_inter_smooth(1:nsd,i)
!      x_inter(1:nsd,inter_count)=x_inter(1:nsd,i)
    else
       ecount=ecount+1
	write(*,*)'eleminate point',i
    end if
  end do
!stop
!  x_inter(1:nsd,:)=x_inter_smooth(1:nsd,:)
  write(*,*)'# of points eleminate',ecount 
  nn_inter=inter_count
  end subroutine x_smooth



