      subroutine  ibg_femcalcptforce(dnext_pt,
     $     dlptlocal_number,dlptlocal_head,dlptlocal_tail,
     $     acttype_con,fix_con,coord_pt,force_con,force_pt,
     $     accel_pt,vel_pt)
      implicit real*8 (a-h,o-z)
      include 'r_common' 
      include 'main_common'            

      integer dnext_pt( mn_point_alloc:mx_point_alloc )

      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail

      dimension force_con(ix:iz,mn_point_alloc:mx_point_alloc)
      dimension force_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      dimension acttype_con(mn_point_alloc:mx_point_alloc )
      dimension fix_con(mn_point_alloc:mx_point_alloc )
      dimension coord_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      dimension accel_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      dimension vel_pt(ix:iz, mn_point_alloc:mx_point_alloc )
      dimension unitvector(n_dim_space)

      if (dlptlocal_number .eq. 0) then
         return
      endif
      if (n_ibmfem .eq. 1) then
         call r_timefun
         call r_load
c++++++++         
c     concentrated force
c++++++++         
         if (numfn .gt. 0) then
            call r_nodalf
         endif

         do i=1,nnd
            fix_con(i)= 1.0
	   enddo

         njc=0

         do i=1,numgb
ccccccccccccccccccccccccccccccccccccccccccc
c     Fixed 1,2,3
ccccccccccccccccccccccccccccccccccccccccccc
            if (ndirgb(i) .eq. 111111) then
               do j=1,numdir(i)
                  nl=nodegb(i,j)
c     
                  if (nxt(1,nl) .ne. 0) then
                     njc=njc+1
                     fix_con(nnd+njc)=-1
c     
c     horizonal support (x-dir)
c
                     xsr=dsqrt((dis(1,nl)-xtedis*nxt(1,nl))**2+
     $                    dis(2,nl)**2+dis(3,nl)**2)
                  
                     tension = xk*(xsr-xtedis+xstretch)

                     unitvector(1)=(dis(1,nl)-
     $                    xtedis*nxt(1,nl))/xsr
                     unitvector(2)=dis(2,nl)/xsr          
                     unitvector(3)=dis(3,nl)/xsr
                     
                     predrf(nl)=predrf(nl)-
     $                    tension*unitvector(1)-
     $                    xvisc*vel_pt(ix,nl)
	               predrf(nl+nnd)=predrf(nl+nnd)-
     $                    tension*unitvector(2)-
     $                    xvisc*vel_pt(iy,nl)
                     predrf(nl+2*nnd)=predrf(nl+2*nnd)-
     $                    tension*unitvector(3)-
     $                    xvisc*vel_pt(iz,nl)
                  endif

c
c     vertical support (y-dir)
c
                  if (nxt(2,nl) .ne. 0) then
                     njc=njc+1
                     fix_con(nnd+njc)=-1

                     xsr=dsqrt((dis(2,nl)-xtedis*nxt(2,nl))**2+
     $                    dis(1,nl)**2+dis(3,nl)**2)
                  
                     tension = xk*(xsr-xtedis+xstretch)

                     unitvector(1)=dis(1,nl)/xsr
                     unitvector(2)=(dis(2,nl)-xtedis*nxt(2,nl))/xsr
                     unitvector(3)=dis(3,nl)/xsr                     
                  
                     predrf(nl)=predrf(nl)-
     $                    xvisc*vel_pt(ix,nl)-
     $                    tension*unitvector(1)
                     predrf(nl+nnd)=predrf(nl+nnd)-
     $                    xvisc*vel_pt(iy,nl)-
     $                    tension*unitvector(2)
                     predrf(nl+2*nnd)=predrf(nl+2*nnd)-
     $                    xvisc*vel_pt(iz,nl)-
     $                    tension*unitvector(3)
                  endif

c
c     z-dir
c
                  if (nxt(3,nl) .ne. 0) then
                     njc=njc+1
                     fix_con(nnd+njc)=-1

                     xsr=dsqrt((dis(3,nl)-xtedis*nxt(3,nl))**2+
     $                    dis(1,nl)**2+dis(2,nl)**2)
                  
                     tension = xk*(xsr-xtedis+xstretch)

                     unitvector(1)=dis(1,nl)/xsr
                     unitvector(2)=dis(2,nl)/xsr
                     unitvector(3)=(dis(3,nl)-xtedis*nxt(3,nl))/xsr                     
                  
                     predrf(nl)=predrf(nl)-
     $                    xvisc*vel_pt(ix,nl)-
     $                    tension*unitvector(1)
                     predrf(nl+nnd)=predrf(nl+nnd)-
     $                    xvisc*vel_pt(iy,nl)-
     $                    tension*unitvector(2)
                     predrf(nl+2*nnd)=predrf(nl+2*nnd)-
     $                    xvisc*vel_pt(iz,nl)-
     $                    tension*unitvector(3)
                  endif
			 enddo
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Fixed 1
ccccccccccccccccccccccccccccccccccccccccccccccccc
            if (ndirgb(i) .eq. 110111) then
               do j=1,numdir(i)
                  nl=nodegb(i,j)
                  fix_con(nl)=-2
                  predrf(nl)=0.0d0
			 enddo
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Fixed 2
ccccccccccccccccccccccccccccccccccccccccccccccccc
            if (ndirgb(i) .eq. 101111) then
               do j=1,numdir(i)
                  nl=nnd+nodegb(i,j)
                  fix_con(nl)=-2
                  predrf(nl)=0.0d0
			 enddo
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Fixed 3
ccccccccccccccccccccccccccccccccccccccccccccccccc
            if (ndirgb(i) .eq. 011111) then
               do j=1,numdir(i)
                  nl=2*nnd+nodegb(i,j)
                  fix_con(nl)=-2
                  predrf(nl)=0.0d0
			 enddo
            endif



	    enddo

      endif


      icon = dlptlocal_head

      ipt = icon

      do icon = dlptlocal_head, dlptlocal_tail
         ipt = icon   
         
         if ( (fix_con(icon) .eq. -1.0) ) then

            force_pt(ix, ipt) = 0.0d0
            force_pt(iy, ipt) = 0.0d0
            force_pt(iz, ipt) = 0.0d0

         elseif (fix_con(icon) .eq. -2.0)  then

            force_pt(ix, ipt) = 0.0d0
            force_pt(iy, ipt) = 0.0d0
            force_pt(iz, ipt) = 0.0d0

         elseif (n_ibmfem .eq. 1)  then !if FEM

            if(ipt. le. nptfem) then
            force_pt(ix, icon) = drf(icon)+predrf(icon)
            force_pt(iy, icon) = drf(icon+nnd)+predrf(icon+nnd)           
            force_pt(iz, icon) = drf(icon+2*nnd)+predrf(icon+2*nnd)
            else 
               force_pt(ix, icon) = force_con(ix,icon) -
     $              force_con(ix,icon-1)
               force_pt(iy, icon) = force_con(iy,icon) - 
     $              force_con(iy,icon-1)
               force_pt(iz, icon) = force_con(iz,icon) - 
     $              force_con(iz,icon-1)
            endif

c++++++
         elseif (n_ibmfem .eq. 0) then !if fiber

            force_pt(ix, ipt) = force_con(ix,icon) -
     $           force_con(ix,icon-1)
            force_pt(iy, ipt) = force_con(iy,icon) - 
     $           force_con(iy,icon-1)
            force_pt(iz, ipt) = force_con(iz,icon) - 
     $           force_con(iz,icon-1)
c Lucy comment this out, because cma is not defined in FEM
c            if (ipt .eq. nptfilea) then
c               force_pt(ix, ipt) = force_pt(ix,ipt)-
c     $              cma*accel_pt(ix,ipt)
c               force_pt(iy, ipt) = force_pt(iy,ipt)-
c     $              cma*accel_pt(iy,ipt)
c               force_pt(iz, ipt) = force_pt(iz,ipt)-
c     $              cma*accel_pt(iz,ipt)
c            endif
         endif

	enddo

      icon = dlptlocal_tail
      ipt = icon 

      return
      end








