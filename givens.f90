subroutine givens(h,M,beta)
! Solve the least square problem M * x = beta 
! with Givens rotation and Back substitution 
! h is a hessenberg matrix
! M is dimension of h, h=dim(M+1,M)
! use beta as the out puts beta=x

integer M
real(8) h(M+1,M)
real(8) beta(M+1)
real(8) x(M+1)

integer i,j,k
real(8) tmp,tmp1,tmp2
real(8) s1,c1
real(8) tmph(2,M)

! Givens rotation to get up-triangle matrix
do i=1,M
	tmp=sqrt(h(i,i)**2 + h(i+1,i)**2)
	s1=h(i+1,i)/tmp
	c1=h(i,i)/tmp
	tmp1=c1*beta(i)
	tmp2=-s1*beta(i)
	beta(i)=tmp1
	beta(i+1)=tmp2;
	
	do j=i,M
		tmph(1,j)=c1*h(i,j)+s1*h(i+1,j)
		tmph(2,j)=-s1*h(i,j)+c1*h(i+1,j)
	end do
	h(i:i+1,:)=tmph(1:2,:)
	tmph(:,:)=0.0d0
end do
!do i=1,M+1
!write(*,*) ' h up-tri', h(i,:)
!end do

!write(*,*) 'beta', beta(:)
! Back substitution used to solve a up-triangle system 
x(M)=beta(M)/h(M,M)

do j=M-1,1,-1
	tmp=0.0d0
	do k=j+1,M
		tmp=tmp+x(k)*h(j,k)
	end do
	x(j)=(beta(j)-tmp)/h(j,j);
end do

beta(:)=x(:)

return
end