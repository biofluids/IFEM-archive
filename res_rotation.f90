

subroutine res_rotation(p,norm_node,nn_spbc,spbcnode)

  use fluid_variables, only:ndf,nn,nsd,ne,lambda
  use mpi_variables

  real(8) p(ndf,nn)
  real(8) norm_node(nsd,nn)
  integer nn_spbc,spbcnode(nn_spbc)

  real(8) norm(nsd),tang(nsd),tang1(nsd),tang2(nsd)
  real(8) temp(nsd),axis(nsd)
  integer i,j,node

  real(8) modmod

if(lambda.lt.0.0) goto 2000
if(nsd==2) then
  do i=1,nn_spbc
     node=spbcnode(i)
     norm(1:nsd)=norm_node(1:nsd,node)
     tang(1)=-norm(2)
     tang(2)=norm(1)
     temp(1)=norm(1)*p(1,node)+norm(2)*p(2,node)
     temp(2)=tang(1)*p(1,node)+tang(2)*p(2,node)
     p(1,node)=temp(1)
     p(2,node)=temp(2)
  end do

else if(nsd==3) then
  do i=1,nn_spbc
     node=spbcnode(i)
     norm(1:nsd)=norm_node(1:nsd,node)
     axis(:)=0.0
     if(abs(norm(1)).lt.0.99) then
	axis(1)=1.0
     else if(abs(norm(2)).lt.0.99) then
	axis(2)=1.0
     else
	axis(3)=1.0
     end if

     tang1(1)=norm(2)*axis(3)-norm(3)*axis(2)
     tang1(2)=norm(3)*axis(1)-norm(1)*axis(3)
     tang1(3)=norm(1)*axis(2)-norm(2)*axis(1)
  
     modmod=sqrt(tang1(1)**2+tang1(2)**2+tang1(3)**2)
     tang1(:)=tang1(:)/modmod     

     tang2(1)=tang1(2)*norm(3)-tang1(3)*norm(2)
     tang2(2)=tang1(3)*norm(1)-tang1(1)*norm(3)
     tang2(3)=tang1(1)*norm(2)-tang1(2)*norm(1)  

     modmod=sqrt(tang2(1)**2+tang2(2)**2+tang2(3)**2)
     tang2(:)=tang2(:)/modmod

     temp(1)=norm(1)*p(1,node)+norm(2)*p(2,node)+norm(3)*p(3,node)  
     temp(2)=tang1(1)*p(1,node)+tang1(2)*p(2,node)+tang1(3)*p(3,node)
     temp(3)=tang2(1)*p(1,node)+tang2(2)*p(2,node)+tang2(3)*p(3,node)

     p(1,node)=temp(1)
     p(2,node)=temp(2)
     p(3,node)=temp(3)
!if(myid==0)write(*,*)'p=',node,p(:,node)
  end do
end if
2000 continue
end subroutine res_rotation
