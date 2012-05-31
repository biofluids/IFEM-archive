

subroutine res_slipbc_precon(p,w,d,x,rngface,ien,spbcnode,spbcele,ne_local,ien_local,I_fluid)

  use fluid_variables
  use mpi_variables
  use interface_variables, only:vis_inter
  include 'mpif.h'
  
  real(8) p(ndf,nn),d(ndf,nn),x(nsd,nn),w(ndf,nn)
  integer rngface(neface,ne),ien(nen,ne)
  integer spbcnode(nn_spbc)  ! includes overlapping boundary nodes
  integer spbcele(ne_spbc)      ! includes overlapping boundary elements
  integer ne_local,ien_local(ne_local)
  real(8) I_fluid(nn)

  integer ieface,inface,inl,irng,ie,i,j,icount,jcount,iec,node(2),isd
  real(8) u(nsd,2),l,temp(nsd)
  integer flag
  real(8) sh(2,2),const,u1(nsd),u2(nsd)
  real(8) mu(2),lam(2)
!  lamda=10.0  ! slip length
if(lambda.lt.0.0)goto 2000
  const=0.577350269
!====================================!
!two guass quad points               !
!====================================!
  sh(1,1)=0.5*(1+const)
  sh(2,1)=0.5*(1-const)
  sh(1,2)=0.5*(1-const)
  sh(2,2)=0.5*(1+const)


  do ie=1,ne_spbc
     iec=spbcele(ie)
     flag=0
     do j=1,ne_local
        if(iec==ien_local(j)) flag=1
     end do
if(flag==1) then
     do ieface=1,neface
        irng=rngface(ieface,iec)
	if(irng.ne.0) then
	   do inface=1,nnface
	      inl=mapping(ieface,inface,etype)
	      node(inface)=ien(inl,iec)
	   end do
! right now all the calculation is based on 2D case, which the boundary edge has only two nodes
	   l= sqrt((x(1,node(1))-x(1,node(2)))**2+(x(2,node(1))-x(2,node(2)))**2)
!	   u(1:nsd)=0.5*(d(1:nsd,node(1))+d(1:nsd,node(2)))
!	   u(1:nsd)=0.0
	   mu(1)=vis_liq+(vis_inter-vis_liq)*I_fluid(node(1))
           u1(1:nsd)=d(1:nsd,node(1))
	   mu(2)=vis_liq+(vis_inter-vis_liq)*I_fluid(node(2))
           u2(1:nsd)=d(1:nsd,node(2))
           u(1:nsd,1)=(u2(1:nsd)-u1(1:nsd))*0.5*(-const)+(u1(1:nsd)+u2(1:nsd))*0.5
           u(1:nsd,2)=(u2(1:nsd)-u1(1:nsd))*0.5*(+const)+(u1(1:nsd)+u2(1:nsd))*0.5
	   lam(1)=((mu(2)-mu(1))*0.5*(-const)+(mu(1)+mu(2))*0.5)/lambda
	   lam(2)=((mu(2)-mu(1))*0.5*(+const)+(mu(1)+mu(2))*0.5)/lambda

	   do inface=1,2
	      temp(1:nsd)=sh(1,inface)*lam(1)*u(1:nsd,1)+sh(2,inface)*lam(2)*u(1:nsd,2)
	      p(1:nsd,node(inface))=p(1:nsd,node(inface))-l*temp(1:nsd)*0.5 
	      w(1:nsd,node(inface))=w(1:nsd,node(inface))+l*0.5*(lam(1)*sh(1,inface)*0.5*(1+const)+ &
								 lam(2)*sh(2,inface)*0.5*(1-const))     
	   end do
	end if
     end do
end if
  end do
!call mpi_barrier(mpi_comm_world,ierror)
!  do i=1,nn_spbc1
!     p(1:nsd,spbcnode1(i))=p(1:nsd,spbcnode1(i))+dp(1:nsd,spbcnode1(i))
!  end do
2000 continue
end subroutine res_slipbc_precon

