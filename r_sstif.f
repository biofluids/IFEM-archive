c     element tangent stiffness matrix assemblage
c
      subroutine r_sstif(ocpp,ocuu,ocup,ne,w,toxj,body_force)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      include 'main_common'
      dimension ocuu(6,6),ocup(6),toxj(3,3),xac(3),xve(3)

c      if (iti .le. nmove) then
c         sdensit=0.0d0
c      else
         sdensit=sdensi
c      endif
ccccccccccc
c     i-u
ccccccccccc

      do ni=1,nis
         do i=1,3
            nu1=(i-1)*nnd+nea(ne,ni)
            mu1=(i-1)*nis+ni
ccccccccccc
c     fu
ccccccccccc
            call r_scalfu(fu,i,ni)
            predrf(nu1)=predrf(nu1)-fu*w
            do nk=1,nump
               call r_scalkup(fkup,ocup,i,nk,ni)
               xkup(mu1,nk,ne)=xkup(mu1,nk,ne)+fkup*w
		  enddo
		enddo
	enddo
c
      do i=1,nump
         call r_scalfp(fp,ocpp,i)
         xfp(i,ne)=xfp(i,ne)+fp*w
         do j=1,nump
            call r_scalkpp(fkpp,ocpp,i,j)
            xkpp(i,j,ne)=xkpp(i,j,ne)+fkpp*w
	   enddo
	enddo
cccccccccccccccccccccccccc
c     inertia forces
cccccccccccccccccccccccccc
      do i=1,3
         xac(i)=0.0d0
         xve(i)=0.0d0
         do k=1,nis
            ntem=nea(ne,k)
            xac(i)=xac(i)+h(k)*acm(i,ntem)
            xve(i)=xve(i)+h(k)*du(i,ntem)
	   enddo
	enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xac(1)=xac(1)+xmg1
      xac(2)=xac(2)+xmg2
	xac(3)=xac(3)+xmg3
	totalh=0
      do ni=1,nis
         nu1=nea(ne,ni)
         nv1=nnd+nea(ne,ni)
	   nw1=2*nnd+nea(ne,ni)

         predrf(nu1)=predrf(nu1)-w*sdensit*h(ni)*xac(1)
     $        -w*xviss*h(ni)*xve(1)
         predrf(nv1)=predrf(nv1)-w*sdensit*h(ni)*xac(2)
     $        -w*xviss*h(ni)*xve(2)
	   predrf(nw1)=predrf(nw1)-w*sdensit*h(ni)*xac(3)
     $        -w*xviss*h(ni)*xve(3)
	totalh=totalh+h(ni)
c	  body_force=body_force+predrf(nu1)
	enddo
		body_force=body_force+w*sdensit*totalh*xac(1)
c	write(*,*) 'body_force=',body_force, total_vol
      return
      end



