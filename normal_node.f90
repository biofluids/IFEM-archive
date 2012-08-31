!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
! used to calculate normal for every !
! node in the domian                 !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!

subroutine normal_node(norm_node,xloc,ien,spbcele,spbcnode,index_bcnode,pbnode)

  use global_constants
  use run_variables
  use fluid_variables
  use mpi_variables
  implicit none

  real(8) norm_node(nsd,nn),xloc(nsd,nn),x(nsd,nen)
  integer ien(nen,ne)
  integer spbcele(ne_spbc),spbcnode(nn_spbc),index_bcnode(2,ne_spbc)
  integer pbnode(2,nn_pb)
  
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)
  real(8) ph(1:nsd,nen) ! 1st derivative of shape function times weight
  integer i,j,icount,jcount,ie,inl,isd,iq,node
  real(8) temp
  norm_node(:,:)=0.0
  index_bcnode(:,:)=0
  do ie=1,ne
     do inl=1,nen
	x(1:nsd,inl)=xloc(1:nsd,ien(inl,ie))
     end do
  

do iq=1,nquad
  if(nsd==2) then
    if(nen==3) then
	include "sh2d3n.h"
    elseif(nen==4) then
	include "sh2d4n.h"
    end if
  elseif(nsd==3) then
    if(nen==4) then
	include "sh3d4n.h"
    elseif(nen==8) then
	include "sh3d8n.h"
    end if
  end if

  eft0=abs(det)*wq(iq)
  ph(1:nsd,1:nen)=sh(1:nsd,1:nen)*eft0
  do inl=1,nen
     node=ien(inl,ie)
     norm_node(1:nsd,node)=norm_node(1:nsd,node)+ph(1:nsd,inl)
  end do

end do ! end of quad loop

  end do ! end of element loop

  do i=1,nn
     temp=0.0
!*************************************************!
!  used to treat edge node                        !
!*************************************************!
!if(abs(xloc(1,i)).gt.1.49999) then
!     if((abs(norm_node(1,i)).gt.0) .and. (abs(norm_node(2,i)).gt.0)) norm_node(1,i)=0.0
!end if
!*************************************************!
     do isd=1,nsd
        temp=temp+norm_node(isd,i)**2
     end do
     temp=sqrt(temp)
     if(temp.lt.1.0e-4) then
        norm_node(:,i)=0.0
     else
        norm_node(:,i)=norm_node(:,i)/temp
     end if
  end do

open(120,file='spbcele.in',status='old')
   do i=1,ne_spbc
      read(120,110)spbcele(i)
   end do
close(120)
open(121,file='spbcnode.in',status='old')
   do i=1,nn_spbc
      read(121,110)spbcnode(i)
   end do
close(121)

110 format(I8)

open(776,file='pbnode.in',status='old')
do j=1,nn_pb
  read(776,'(I8,I8)')pbnode(1:2,j)
end do
close(776)

do i=1,ne_spbc
   ie=spbcele(i)
   icount=0
   do inl=1,nen
      node=ien(inl,ie)
      do j=1,nn_spbc
         if(node==spbcnode(j)) then
           icount=icount+1
           index_bcnode(icount,i)=node
         end if
      end do
   end do
end do

end subroutine normal_node
     










