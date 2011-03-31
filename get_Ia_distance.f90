!======================================
!reset Ia based on the distance between the center point and the interfacial point

!======================================

subroutine get_Ia_distance(ne_inter,inter_ele,x_inter,x_center,hg,I_fluid_center)

  use interface_variables
  use fluid_variables,only:nsd,ne

  integer ne_inter
  integer inter_ele(ne)
  real(8) x_inter(nsd,maxmatrix)
  real(8) x_center(nsd,ne)
  real(8) hg(ne)
  real(8) I_fluid_center(ne)

  integer i,j,icount,jcount
  real(8) mindis,hs,length
  real(8) xg(nsd)

  do i=1,ne_inter
     xg(1:nsd)=x_center(1:nsd,inter_ele(i))
     hs=hg(inter_ele(i))
     mindis=999.0
     do j=1,nn_inter
        length=(xg(1)-x_inter(1,j))**2+(xg(2)-x_inter(2,j))**2
        length=sqrt(length)
        if(length.lt.mindis) then
          mindis=length
        end if
     end do
     if(mindis.le.hs) then
        if(I_fluid_center(inter_ele(i)).lt.0.49999999)then
          I_fluid_center(inter_ele(i))=0.5*(1.0-mindis/hs)
        else if(I_fluid_center(inter_ele(i)).gt.0.50000001) then
          I_fluid_center(inter_ele(i))=0.5*(1.0+mindis/hs)
        end if
     end if
  end do
end subroutine get_Ia_distance

