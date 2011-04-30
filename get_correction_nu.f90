!=====================================
! get correction term
!=====================================

subroutine get_correction_nu(x_inter,x_center,hg,infdomain_inter,corr_Ip,I_fluid_center)

  use interface_variables
  use fluid_variables, only:nn,ne,nsd
  use allocate_variables, only:center_domain,nn_center_domain
  use mpi_variables

  real(8) x_inter(nsd,maxmatrix)
  real(8) x_center(nsd,ne)
  real(8) hg(ne)
  integer infdomain_inter(maxmatrix)
  real(8) corr_Ip(maxmatrix)
  real(8) I_fluid_center(ne)
!  real(8) x(nsd,nn),I_fluid(nn)

  integer i,j,ie,icount,jcount
  real(8) dx(nsd),hs,Sp
  real(8) A(nn_inter,nn_inter)
  real(8) BB(nn_inter)  !Ax=B
  integer IPIV(nn_inter), INFO

  real(8) M(nsd+1,nsd+1),B(nsd+1,nn_inter),P(nsd+1) ! used for nonuniform
  real(8) vec(nsd+1),hsg !M*B=P from rkpm
  integer IP(nsd+1)
  real(8) temp

  real(8) debugA
  real(8) temph

  P(:)=0.0
  P(1)=1.0
  vec(1)=1.0
  M(:,:)=0.0
  B(:,:)=0.0
  do i=1,nn_inter
!     hs=hg(infdomain_inter(i))
     hs=hsp
     M(:,:)=0.0
     B(:,i)=0.0
     P(:)=0.0
     P(1)=1.0
     do j=1,nn_center_domain
	ie=center_domain(j)
	hsg=hg(ie)
        dx(:)=abs(x_inter(:,i)-x_center(:,ie))
	call B_Spline(dx,hs,nsd,Sp)
	vec(2:nsd+1)=x_inter(:,i)-x_center(:,ie)
	do icount=1,nsd+1
	   do jcount=1,nsd+1
	      M(icount,jcount)=M(icount,jcount)+vec(icount)*vec(jcount)*Sp/(hs**nsd)*(hsg**nsd)
	   end do
	end do
     end do
     call DGESV(nsd+1,1,M,nsd+1,IP,P,nsd+1,INFO)
     B(:,i)=P(:)
  end do
  A(:,:)=0.0
  BB(:)=0.0

  do i=1,nn_inter
!     hs=hg(infdomain_inter(i))
     hs=hsp
     do j=1,nn_inter
!      do j=1,2
	dx(:)=abs(x_inter(:,i)-x_inter(:,j))
	call B_Spline(dx,hs,nsd,Sp)
!=====
	vec(2:nsd+1)=x_inter(:,i)-x_inter(:,j)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount,i)
	end do
!=======
	A(i,j)=Sp*temp
     end do
  end do

  do i=1,nn_inter
!     hs=hg(infdomain_inter(i))
     hs=hsp
     do j=1,nn_center_domain
	ie=center_domain(j)
	hsg=hg(ie)
	dx(:)=abs(x_inter(:,i)-x_center(:,ie))
	call B_Spline(dx,hs,nsd,Sp)

	vec(2:nsd+1)=x_inter(:,i)-x_center(:,ie)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount,i)
	end do
	BB(i)=BB(i)+I_fluid_center(ie)*temp*Sp/(hs**nsd)*(hsg**nsd)
!	debugA=debugA+temp*Sp/(hs**nsd)*(hsg**nsd)
     end do
!if((myid==0).and.(i==1))write(*,*)'debugA=',debugA
  end do

     BB(:)=0.5-BB(:)
!if(myid==0)write(*,*)'BB=',BB(:)
!if(myid==0)write(*,*)'================='
     call DGESV(nn_inter,1,A,nn_inter,IPIV,BB,nn_inter,INFO)

     corr_Ip(1:nn_inter)=BB(1:nn_inter) !calculate delta Ip
!  I_fluid(1:nn)=0.0
!  do i=1,nn
!     do j=1,nn_inter
!        dx(:)=abs(x(:,i)-x_inter(:,j))
!        call B_Spline(dx,hs,nsd,Sp)
!        I_fluid(i)=I_fluid(i)+corr_Ip(j)*Sp
!     end do
!     do j=1,ne
!        dx(:)=abs(x(:,i)-x_center(:,j))
!	call B_Spline(dx,hs,nsd,Sp)
!	I_fluid(i)=I_fluid(i)+I_fluid_center(j)*Sp
!     end do
!  end do

end subroutine get_correction_nu