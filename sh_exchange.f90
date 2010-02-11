! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the shape function for both 2-D and 3-D case
! 3-D case has not been finished yet
! The shape function here also can be used as the cooeffiecints for the force distribution
! Except triangle, other shape function need be fixed manually,
! since the digital error can not be avoided
! Xingshi Wang 2008
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine sh_exchange(x,xe,nsd,nen,sh)
integer nen
integer nsd
real(8) x(nsd)
real(8) xe(nsd,nen) ! xyz coordinates of the element nodals
real(8) sh(nen)
real(16) sh1(nen)
real(16) me(nen,nen)
real(8) me1(nen,nen)
real(16) invme(nen,nen)
real(8) det
real(16) xp(nen)
integer indx(nen)
integer i ! loop variable
integer j
! Lapack variables
integer info
integer ipiv(nen)
real(8) work(nen)

        do i=1,nen
        sh1(i)=0
        do j=1,nen
        me(i,j)=0
        end do
        end do

        if (nsd .eq. 2) then  ! 2-D case
           if(nen .eq. 3) then ! 2-D triangle
           do i=1,nen
           me1(i,1)=1.0
           me1(i,2:(nsd+1))=xe(1:nsd,i)
           end do
           call determinant(me1,nen,nen,det)
          ! write(*,*) det
           sh(1)=(xe(1,2)*xe(2,3)-xe(1,3)*xe(2,2)+&
           (xe(2,2)-xe(2,3))*x(1)+&
           (xe(1,3)-xe(1,2))*x(2))/(det)

           sh(2)=(xe(1,3)*xe(2,1)-xe(1,1)*xe(2,3)+&
           (xe(2,3)-xe(2,1))*x(1)+&
           (xe(1,1)-xe(1,3))*x(2))/(det)
           sh(3)=(xe(1,1)*xe(2,2)-xe(1,2)*xe(2,1)+&
           (xe(2,1)-xe(2,2))*x(1)+&
           (xe(1,2)-xe(1,1))*x(2))/(det)

           sh(3)=1-sh(1)-sh(2)

           else if (nen .eq. 4) then ! 2-D quadrilateral


           do i=1,nen
           me(i,1)=1.0
           me(i,2)=xe(1,i)
           me(i,3)=xe(2,i)
           me(i,4)=xe(1,i)*xe(2,i)
           end do
           xp(1)=1
           xp(2)=x(1)
           xp(3)=x(2)
           xp(4)=x(1)*x(2)

           call MIGS(me,nen,invme,indx) ! get inverse

                call DGETRF(nen,nen,me,nen,ipiv,info)
                call DGETRI(nen,me,nen,ipiv,work,nen,info)
           write(*,*) 'difference', me(:,:)-invme(:,:)
           do j=1,nen
           do i=1,nen
 !          write(*,*) 'invme', invme(i,j)
           sh1(j)=sh1(j)+xp(i)*invme(i,j)
           sh(j)=sh1(j)
           end do
           !write(*,*) 'shape1', sh1(j)
           end do
           end if
        else
           if (nen .eq. 4) then ! 3-D tet case
                   call shx_tets(x,xe,nsd,nen,sh)
                   end if




        end if
        do i=1,nen
          if (sh(i)<0) then
            call fix_sh(sh,nen)
        write(*,*) '****errors in shape function, running fix munually***'
          do j=1,nen
!         write(*,*) 'xe=', xe(1,j), xe(2,j)
         write(*,*)'sh= ', sh(j)
!         write(*,*) 'xp= ', xp(j)
         write(*,*) 'sum of shape function should be 1', sum(sh)
          end do
          go to 1000
          end if

        end do

1000    continue
return
end

