subroutine bc_shape2d_node(sq,nsd,nen)
integer nsd
integer nen
real(8) sq(0:nsd,nen,nen)
real(8) xq(nsd,nen)
integer iq

if (nen == 3) then
        xq(1,1) = 1.0d0 ! integration point on  1
        xq(2,1) = 0.0d0
        
        xq(1,2)=0.0d0  ! integeration point on  2
        xq(2,2)=1.0d0
        
        xq(1,3)=0.0d0  ! integration point on   3
        xq(2,3)=0.0d0
end if

if (nen == 4) then
        xq(1,1)=-1.0d0 ! integration point on  1
        xq(2,1)=-1.0d0

        xq(1,2)=1.0d0  ! integration point on 2
        xq(2,2)=-1.0d0

        xq(1,3)=1.0d0 ! integration point on 3
        xq(2,3)=1.0d0

        xq(1,4)=-1.0d0 ! integration point on 4
        xq(2,4)=1.0d0
end if

do iq=1,nen
        if(nen==3) then
                  sq(0,1,iq) = xq(1,iq)
                  sq(0,2,iq) = xq(2,iq)
                  sq(0,3,iq) = 1 - xq(1,iq) - xq(2,iq)
        elseif (nen==4) then
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

        endif
end do
return
end
