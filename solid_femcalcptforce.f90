subroutine  solid_femcalcptforce
  use r_common
  use solid_variables
  implicit none

  !real*8 unitvector(nsd_solid)

  !real(kind=8) :: temp
  integer :: njc,ipt


  call r_timefun
  call r_load
!++++++++         
!     concentrated force
!++++++++         
  if (numfn .gt. 0) then
     call r_nodalf
  endif

  njc=0

!         do i=1,numgb
!ccccccccccccccccccccccccccccccccccccccccccc
!c     Fixed 1,2,3
!ccccccccccccccccccccccccccccccccccccccccccc
!            if (ndirgb(i) .eq. 111111) then
!               do j=1,numdir(i)
!                  nl=nodegb(i,j)
!
!     
!     (x-dir)
!     
!                  if (nxt(1,nl) .ne. 0) then
!                     njc=njc+1
!                     fix_con(nnd+njc)=-1
!	               temp = (dis(1,nl)-xtedis*nxt(1,nl))**2 +
!     $                    dis(2,nl)**2 + dis(3,nl)**2
!                     if (temp /= 0.0d0) then

!                        xsr=dsqrt((dis(1,nl)-xtedis*nxt(1,nl))**2+
!     $                    dis(2,nl)**2+dis(3,nl)**2)
!                  
!                        tension = xk*(xsr-xtedis+xstretch)

!                        unitvector(1)=(dis(1,nl)-
!     $                    xtedis*nxt(1,nl))/xsr
!                        unitvector(2)=dis(2,nl)/xsr          
!                        unitvector(3)=dis(3,nl)/xsr
                     
!                        predrf(nl)=predrf(nl)-
!     $                    tension*unitvector(1)-
!     $                    xvisc*vel_pt(1,nl)!
!	                  predrf(nl+nnd)=predrf(nl+nnd)-
!     $                    tension*unitvector(2)-
!     $                    xvisc*vel_pt(2,nl)
!                        predrf(nl+2*nnd)=predrf(nl+2*nnd)-
!     $                    tension*unitvector(3)-
!     $                    xvisc*vel_pt(3,nl)
!	               
!	               endif
 !                 endif
!
!
!     (y-dir)
!
!                  if (nxt(2,nl) .ne. 0) then
!                     njc=njc+1
!                     fix_con(nnd+njc)=-1
!				   
!				   temp = (dis(2,nl)-xtedis*nxt(2,nl))**2+
!     $                    dis(1,nl)**2+dis(3,nl)**2
!	               if (temp /= 0.0d0) then
                        
!					  xsr=dsqrt((dis(2,nl)-xtedis*nxt(2,nl))**2+
 !    $                    dis(1,nl)**2+dis(3,nl)**2)
                  
!                        tension = xk*(xsr-xtedis+xstretch)

!                        unitvector(1)=dis(1,nl)/xsr
!                        unitvector(2)=(dis(2,nl)-xtedis*nxt(2,nl))/xsr
!                        unitvector(3)=dis(3,nl)/xsr                     
                  
!                        predrf(nl)=predrf(nl)-
!     $                     xvisc*vel_pt(1,nl)-
!     $                     tension*unitvector(1)
!                        predrf(nl+nnd)=predrf(nl+nnd)-
!     $                     xvisc*vel_pt(2,nl)-
!     $                     tension*unitvector(2)
!                        predrf(nl+2*nnd)=predrf(nl+2*nnd)-
!     $                     xvisc*vel_pt(3,nl)-
!     $                     tension*unitvector(3)
!                     endif
!                  endif
!
!     z-dir
!
!                  if (nxt(3,nl) .ne. 0) then
!                     njc=njc+1
!                     fix_con(nnd+njc)=-1

! 	               temp = (dis(3,nl)-xtedis*nxt(3,nl))**2+
!     $                    dis(1,nl)**2+dis(2,nl)**2

!                     if (temp /= 0.0d0) then

!                        xsr=dsqrt((dis(3,nl)-xtedis*nxt(3,nl))**2+
 !    $                    dis(1,nl)**2+dis(2,nl)**2)
                  
!                        tension = xk*(xsr-xtedis+xstretch)
!
!                        unitvector(1)=dis(1,nl)/xsr
!                        unitvector(2)=dis(2,nl)/xsr
!                        unitvector(3)=(dis(3,nl)-xtedis*nxt(3,nl))/xsr                     
!                  
!                        predrf(nl)=predrf(nl)-
!     $                    xvisc*vel_pt(1,nl)-
!     $                    tension*unitvector(1)
!                        predrf(nl+nnd)=predrf(nl+nnd)-
!     $                    xvisc*vel_pt(2,nl)-
!     $                    tension*unitvector(2)
!                        predrf(nl+2*nnd)=predrf(nl+2*nnd)-
!     $                    xvisc*vel_pt(3,nl)-
!     $                    tension*unitvector(3)
!	               endif
!                  endif
!			 enddo
!            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Fixed 1
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!            if (ndirgb(i) .eq. 110111) then
!               do j=1,numdir(i)
!                  nl=nodegb(i,j)
!                  fix_con(nl)=-2
!                  predrf(nl)=0.0d0
!			 enddo
!            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Fixed 2
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!            if (ndirgb(i) .eq. 101111) then
!               do j=1,numdir(i)
!                  nl=nnd+nodegb(i,j)
!                  fix_con(nl)=-2
!                  predrf(nl)=0.0d0
!			 enddo
!            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Fixed 3
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!            if (ndirgb(i) .eq. 011111) then
!               do j=1,numdir(i)
!                  nl=2*nnd+nodegb(i,j)
!                  fix_con(nl)=-2
!                  predrf(nl)=0.0d0
!			 enddo
!            endif



!	    enddo


  do ipt = 1,nn_solid

     solid_force_FSI(1, ipt) = drf(ipt)            + predrf(ipt)
     solid_force_FSI(2, ipt) = drf(ipt+nn_solid)   + predrf(ipt+nn_solid)
     solid_force_FSI(3, ipt) = drf(ipt+2*nn_solid) + predrf(ipt+2*nn_solid)

  enddo

  return
end subroutine solid_femcalcptforce








