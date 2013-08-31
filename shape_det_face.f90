

subroutine shape_det_face(iq,xloc,sq,det,sh)

  integer iq
  real(8) xloc(3,4),sq(0:2,4,4)
  real(8) det,sh(4)

  real(8) x(2,4) ! transformed x on the same surface
  real(8) xr(2,2)
  real(8) nx(3),ny(3),nz(3)  ! x,y,z axis
  real(8) temp

  integer i,j,icount,jcount

  nx(:)=xloc(:,2)-xloc(:,1)

  temp=sqrt(nx(1)**2+nx(2)**2+nx(3)**2)
  nx(:)=nx(:)/temp

  ny(:)=xloc(:,3)-xloc(:,2)

  temp=sqrt(ny(1)**2+ny(2)**2+ny(3)**2)
  ny(:)=ny(:)/temp

  nz(1)=nx(2)*ny(3)-nx(3)*ny(2)
  nz(2)=nx(3)*ny(1)-nx(1)*ny(3)
  nz(3)=nx(1)*ny(2)-nx(2)*ny(1)

  temp=sqrt(nz(1)**2+nz(2)**2+nz(3)**2)
  nz(:)=nz(:)/temp

  ny(1)=nz(2)*nx(3)-nz(3)*nx(2)
  ny(2)=nz(3)*nx(1)-nz(1)*nx(3)
  ny(3)=nz(1)*nx(2)-nz(2)*nx(1)

  temp=sqrt(ny(1)**2+ny(2)**2+ny(3)**2)
  ny(:)=ny(:)/temp

  do i=1,4
     x(1,i)=(xloc(1,i)-xloc(1,1))*nx(1)+(xloc(2,i)-xloc(2,1))*nx(2)+(xloc(3,i)-xloc(3,1))*nx(3)
     x(2,i)=(xloc(1,i)-xloc(1,1))*ny(1)+(xloc(2,i)-xloc(2,1))*ny(2)+(xloc(3,i)-xloc(3,1))*ny(3)
  end do

        sh(1) = sq(0,1,iq)
        sh(2) = sq(0,2,iq)
        sh(3) = sq(0,3,iq)
        sh(4) = sq(0,4,iq)

        xr(1,1) = sq(1,1,iq) * x(1,1) + sq(1,2,iq) * x(1,2) &
                    + sq(1,3,iq) * x(1,3) + sq(1,4,iq) * x(1,4)
        xr(1,2) = sq(1,1,iq) * x(2,1) + sq(1,2,iq) * x(2,2) &
                    + sq(1,3,iq) * x(2,3) + sq(1,4,iq) * x(2,4)
        xr(2,1) = sq(2,1,iq) * x(1,1) + sq(2,2,iq) * x(1,2) &
                    + sq(2,3,iq) * x(1,3) + sq(2,4,iq) * x(1,4)
        xr(2,2) = sq(2,1,iq) * x(2,1) + sq(2,2,iq) * x(2,2) &
                    + sq(2,3,iq) * x(2,3) + sq(2,4,iq) * x(2,4)


    det = + (xr(1,1)*xr(2,2) - xr(2,1)*xr(1,2))
if(det.lt.0.0)write(*,*)'some thing wrong with boundary intergral'
end subroutine shape_det_face


