subroutine res_solid(x,x_curr,ien,p,xpre1,xpre2,kid,solid_acc,pre,solid_bcvel,solid_bcvel_old,solid_vel)
! subroutine to evaluate solid RHS residual after solid displacement is solved
! x  -- initial solid configuration
! x_curr --- current solid configuration
! disp --- solid displacement 
! ien --- solid connective matrix
! xpre1 --- last 1 time step solid position
! xpre2 --- last 2 time step solid position
use mpi_variables, only: myid
use solid_variables
implicit none

real(8) x(nsd_solid,nn_solid)
real(8) x_curr(nsd_solid,nn_solid)
real(8) disp(nsd_solid,nn_solid)
integer ien(ne_solid,nen_solid)
integer kid(nsd_solid,nn_solid)
real(8) xpre1(nsd_solid,nn_solid)
real(8) xpre2(nsd_solid,nn_solid)
real(8) disp1(nsd_solid,nn_solid)
real(8) disp2(nsd_solid,nn_solid)
real(8) solid_acc(nsd_solid,nn_solid)
real(8) pre(nn_solid)
real(8) solid_bcvel(nsd_solid,nn_solid)
real(8) solid_bcvel_old(nsd_solid,nn_solid)
real(8) solid_vel(nsd_solid,nn_solid)
!-------------------------------------------------
real(8) sq_solid(0:3,8,8)
integer iq
!real(8) :: xq_solid(nsdpad_solid,nquadpad_solid),wq(nquadpad_solid)
!-------------------------------------------------
integer inner
integer outer
parameter (inner = 50) ! solid equation inner 50 should be sufficient
parameter (outer = 10)  ! solid equation outer 5 should be sufficient
!-------------------------------------------------
real(8) w(nsd_solid,nn_solid)  ! pre-conditioner
real(8) p(nsd_solid,nn_solid)  ! residual vector
real(8) res
integer i
! Get sq for solid mesh to calculate sh in block_solid.f90
if (nsd_solid == 3) then
     do iq=1,nquad_solid
                if(nen_solid.eq.4) then
                  sq_solid(0,1,iq) = xq_solid(1,iq)
                  sq_solid(0,2,iq) = xq_solid(2,iq)
                  sq_solid(0,3,iq) = xq_solid(3,iq)
                  sq_solid(0,4,iq) = 1 - xq_solid(1,iq) - xq_solid(2,iq) - xq_solid(3,iq)
        else
                  sq_solid(0,1,iq) = (1 - xq_solid(1,iq))   &
                           * (1 - xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(0,2,iq) = (1 + xq_solid(1,iq))   &
                           * (1 - xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(0,3,iq) = (1 + xq_solid(1,iq))   &
                           * (1 + xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(0,4,iq) = (1 - xq_solid(1,iq))   &
                           * (1 + xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(0,5,iq) = (1 - xq_solid(1,iq))   &
                           * (1 - xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(0,6,iq) = (1 + xq_solid(1,iq))   &
                           * (1 - xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(0,7,iq) = (1 + xq_solid(1,iq))   &
                           * (1 + xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(0,8,iq) = (1 - xq_solid(1,iq))   &
                           * (1 + xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(1,1,iq) = - (1 - xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(1,2,iq) = + (1 - xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(1,3,iq) = + (1 + xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(1,4,iq) = - (1 + xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(1,5,iq) = - (1 - xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(1,6,iq) = + (1 - xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(1,7,iq) = + (1 + xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(1,8,iq) = - (1 + xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(2,1,iq) = - (1 - xq_solid(1,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(2,2,iq) = - (1 + xq_solid(1,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(2,3,iq) = + (1 + xq_solid(1,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(2,4,iq) = + (1 - xq_solid(1,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(2,5,iq) = - (1 - xq_solid(1,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(2,6,iq) = - (1 + xq_solid(1,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(2,7,iq) = + (1 + xq_solid(1,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(2,8,iq) = + (1 - xq_solid(1,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(3,1,iq) = - (1 - xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 8
                  sq_solid(3,2,iq) = - (1 + xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 8
                  sq_solid(3,3,iq) = - (1 + xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 8
                  sq_solid(3,4,iq) = - (1 - xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 8
                  sq_solid(3,5,iq) = + (1 - xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 8
                  sq_solid(3,6,iq) = + (1 + xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 8
                  sq_solid(3,7,iq) = + (1 + xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 8
                  sq_solid(3,8,iq) = + (1 - xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 8
        endif
          enddo
elseif (nsd_solid == 2) then
                  
      do iq=1,nquad_solid
                if(nen_solid==3) then
                  sq_solid(0,1,iq) = xq_solid(1,iq)
                  sq_solid(0,2,iq) = xq_solid(2,iq)
                  sq_solid(0,3,iq) = 1 - xq_solid(1,iq) - xq_solid(2,iq)
        elseif (nen_solid==4) then
                  sq_solid(0,1,iq) = (1 - xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 4
                  sq_solid(0,2,iq) = (1 + xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 4
                  sq_solid(0,3,iq) = (1 + xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 4
                  sq_solid(0,4,iq) = (1 - xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 4

                  sq_solid(1,1,iq) = - (1 - xq_solid(2,iq)) / 4
                  sq_solid(1,2,iq) = + (1 - xq_solid(2,iq)) / 4
                  sq_solid(1,3,iq) = + (1 + xq_solid(2,iq)) / 4
                  sq_solid(1,4,iq) = - (1 + xq_solid(2,iq)) / 4
                  
                  sq_solid(2,1,iq) = - (1 - xq_solid(1,iq)) / 4
                  sq_solid(2,2,iq) = - (1 + xq_solid(1,iq)) / 4
                  sq_solid(2,3,iq) = + (1 + xq_solid(1,iq)) / 4
                  sq_solid(2,4,iq) = + (1 - xq_solid(1,iq)) / 4
                  
        endif     
          enddo
end if
!------------------------------------------------------------------
disp(:,:)=x_curr(:,:)-x(:,:)
disp1(:,:)=xpre1(:,:) - x(:,:)
disp2(:,:)=xpre2(:,:) - x(:,:)


w(:,:)=0.0d0
!p(:,:)=0.0d0

call block_solid_res(x_curr,disp,w,p,ien,nsd_solid,nen_solid,ne_solid,&
		nn_solid,nquad_solid,wq_solid,sq_solid,disp1,disp2,solid_acc,pre,solid_bcvel,solid_bcvel_old,solid_vel)
!call setsolid_id(p,kid,nsd_solid)
disp1(:,:) = solid_vel(:,:) - solid_bcvel(:,:)

call getnorm(disp1,disp1,nsd_solid*nn_solid,res)
res=sqrt(res)
if (myid == 0) write(*,*) 'Norm of V^s - V^f', res

return
end 
