!!calculate normal using B_Spline function!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_normal_Bspline(xp,xg,Ig_ini,infdomain_inter,corr_Ip,norm_inter,hg)

  use interface_variables
  use fluid_variables, only:nsd,ne,nn

  real(8) xp(nsd,maxmatrix),xg(nsd,ne)
  real(8) Ig_ini(ne)
  integer infdomain_inter(maxmatrix)
  real(8) corr_Ip(maxmatrix)
  real(8) norm_inter(nsd,maxmatrix)
  real(8) hg(ne)

  integer i,j,isd,jsd
  real(8) hs,Sp(nsd),temp,dx

  norm_inter(:,:) = 0.0
  do i=1,nn_inter
!     hs=hg(infdomain_inter(i))
     do j=1,ne ! loop over fluid nodes
	hs=hg(j)
	do isd=1,nsd ! loop over isd
	   do jsd=1,nsd ! loop over jsd,if jsd=isd, do the derivative
	      dx=xp(jsd,i)-xg(jsd,j)
	      if(jsd==isd) then
	         call B_Spline_1order(dx,hs,Sp(jsd))
	      else
	         call B_Spline_0order(dx,hs,Sp(jsd))
	      end if
	   end do ! get 1order derivative for isd
	   temp=1.0
	   do jsd=1,nsd
	      temp=temp*Sp(jsd)
	   end do 
	   norm_inter(isd,i)=norm_inter(isd,i)-Ig_ini(j)*temp
	end do ! end of loop over isd

     end do ! end of loop over fluid nodes
  end do ! end of loop over interfacial points
!write(*,*)'I_inter',norm_inter(1,1:nn_inter)
  do i=1,nn_inter
!     hs=hg(infdomain_inter(i))
     do j=1,nn_inter
	hs=hg(infdomain_inter(j))/1.0
	do isd=1,nsd
	   do jsd=1,nsd
	      dx=xp(jsd,i)-xp(jsd,j)
	      if(jsd==isd) then
	         call B_Spline_1order(dx,hs,Sp(jsd))
	      else
	         call B_Spline_0order(dx,hs,Sp(jsd))
	      end if
	   end do
	   temp=1.0
	   do jsd=1,nsd
	      temp=temp*Sp(jsd)
	   end do
	   norm_inter(isd,i)=norm_inter(isd,i)-corr_Ip(j)*temp
	end do

     end do
  end do
end subroutine get_normal_Bspline
















