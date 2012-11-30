

subroutine res_slipbc_3D(p,d,x,rngface,ien,spbcele,ne_local,ien_local,I_fluid)

  use fluid_variables,only:nn,ndf,nsd,ne,nen,neface,mapping,etype,lambda,vis_liq,ne_spbc,nnface
  use mpi_variables
  use interface_variables, only:vis_inter
  include 'mpif.h'
  
  real(8) p(ndf,nn),d(ndf,nn),x(nsd,nn)
  integer rngface(neface,ne),ien(nen,ne)
  integer spbcele(ne_spbc)      ! includes overlapping boundary elements
  integer ne_local,ien_local(ne_local)
  real(8) I_fluid(nn)

! right now only for hex element that the boundary face is a quadrilateral and has only 4 nodes.
  real(8) sq(0:2,1:4,1:4),sh(4)
  real(8) det

  integer iquad,nquad  !iquad=2, nquad=4 for quadrilateral element

  real(8) xq(2,4),wq(4)

  integer iq

  integer bcindex(6),num_slipbc
  real(8) xloc(3,4) !coordinates for the nodes on the face

  integer ieface,inface,inl,irng,ie,i,j,icount,jcount,iec,node(4),isd
  integer flag

  real(8) u(nsd,4),mu(4),I_node(4),res(nsd)!residual
  real(8) temp

  sq(0:2,1:4,1:4)=0.0
  sh(:)=0.0
  xq(:,:)=0.0
  wq(:)=0.0
  xloc(:,:)=0.0
  det=0.0
  u(:,:)=0.0
  mu(:)=0.0
  I_node(:)=0.0
  res(:)=0.0
  if(lambda.lt.0.0)goto 2000
  bcindex(:)=0
  bcindex(1)=6
  num_slipbc=1

  iquad=2
  nquad=4
  call quad2d4n(iquad,nquad,xq,wq,2,4)

  do iq=1,nquad
                  sq(0,1,iq) = (1 - xq(1,iq)) * (1 - xq(2,iq)) / 4
                  sq(0,2,iq) = (1 + xq(1,iq)) * (1 - xq(2,iq)) / 4
                  sq(0,3,iq) = (1 + xq(1,iq)) * (1 + xq(2,iq)) / 4
                  sq(0,4,iq) = (1 - xq(1,iq)) * (1 + xq(2,iq)) / 4
                
                  sq(1,1,iq) = - (1 - xq(2,iq)) / 4
                  sq(1,2,iq) = + (1 - xq(2,iq)) / 4
                  sq(1,3,iq) = + (1 + xq(2,iq)) / 4
                  sq(1,4,iq) = - (1 + xq(2,iq)) / 4
        
                  sq(2,1,iq) = - (1 - xq(1,iq)) / 4
                  sq(2,2,iq) = - (1 + xq(1,iq)) / 4
                  sq(2,3,iq) = + (1 + xq(1,iq)) / 4
                  sq(2,4,iq) = + (1 - xq(1,iq)) / 4
  end do

  do ie=1,ne_spbc
     iec=spbcele(ie)
     flag=0
     do j=1,ne_local
        if(iec==ien_local(j)) flag=1
     end do
     if (flag==1) then
        
	do ieface=1,neface  ! loop over all the faces for that elements
	   irng=rngface(ieface,iec)  ! irng indicates which bc the face is on. 0 for interior. 

	   if(irng==6) then  !slip bc id = 6
	      do inface=1,nnface  ! number of nodes per face
	         inl=mapping(ieface,inface,etype)
	         node(inface)=ien(inl,iec)
	         xloc(:,inface)=x(:,node(inface))
	         u(1:nsd,inface)=d(1:nsd,node(inface))
	         I_node(inface)=I_fluid(node(inface))
	         mu(inface)=(vis_liq+(vis_inter-vis_liq)*I_node(inface))/lambda
	      end do
	     do iq=1,nquad 
	        call shape_det_face(iq,xloc,sq,det,sh)
	        res(:)=0.0
	        do j=1,4
		   res(:)=res(:)+sh(j)*mu(j)*u(:,j)
	        end do
	        do j=1,4
		   p(1:nsd,node(j))=p(1:nsd,node(j))-sh(j)*abs(det)*wq(iq)*res(1:nsd)
	        end do
	        temp=0.0
!               do j=1,4
!                   temp=temp+sh(j)*mu(j)
!                end do
!		do j=1,4
!		   w(1:nsd,node(j))=w(1:nsd,node(j))+sh(j)*abs(det)*wq(iq)*temp
!		end do



	     end do ! end of iq over nquad

	     

	   endif  ! end of irng if


        end do ! end of ieface loop

     end if  ! end of flag

  end do ! end of loop over ie
call mpi_barrier(mpi_comm_world,ierror)
2000 continue


end subroutine res_slipbc_3D

