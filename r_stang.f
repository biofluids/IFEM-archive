      subroutine r_stang
      implicit real*8 (a-h,o-z)
      include 'r_common'
      real*8 x(3,9),toc(3,3),xto(3,3),xot(3,3),
     $     xj(3,3),xji(3,3),rs(3),toxj(3,3),toxji(3,3)
      real*8 xmj(3),xmi(3),dxmj(3,6),ddxmj(3,6,6),
     $     obc(6,6),ocuu(6,6),ocup(6)
      real*8 xfrtem(6,6),tem(6),ten(6),ttm(6)
c
      do i=1,3*nnd
         predrf2(i)=predrf(i)
         predrf(i)=0.0d0
         drf2(i)=drf(i)
	enddo
c     
	tot_vol=0.0
	body_force=0.0
      do ne=1,numel
c     position
         do j=1,nump
            xfp(j,ne)=0.0d0
            do i=1,nis
               xkup(i,j,ne)=0.0d0
               xkup(i+nis,j,ne)=0.0d0
			 xkup(i+2*nis,j,ne)=0.0d0
		  enddo
            do i=1,nump
               xkpp(i,j,ne)=0.0d0
		  enddo
	   enddo
c     
         do nos=1,nis
            ntem=nea(ne,nos)
            do noj=1,3
               x(noj,nos)=coor(ntem,noj)
               y(noj,nos)=coor(ntem,noj)+dis(noj,ntem)
		  enddo
	   enddo
c     gauss integration
         do lx=1,nint
           rs(1)=xg(lx,nint)
           do ly=1,nint
             rs(2)=xg(ly,nint)
		   do lz=1,nint
			 rs(3)=xg(lz,nint)
c     isoparametric interpolation
               call r_element(rs)
c     y-(r,s)
               call r_jacob(y,xj,xji,det)
c     x-(r,s)
               call r_jacob(x,toxj,toxji,todet)
c     derivative about ox
               call r_bdpd(toxji)
c     deformation gradient
               call r_stoxc(xto,xot,xj,xji,toxj,toxji,toc)
c     material j
               call r_smaterj(wto,toc,xmi,xmj,dxmj,ddxmj)
c     discretized pressure
               call r_spress(rs,ne)
c     continuous pressure
               call r_sbpress(dxmj,ddxmj,xmj)
c     material c
               call r_sboc(obc,ocpp,ocuu,ocup,xmj,dxmj,ddxmj)
c     strain
               call r_sstrain(toc,xto,lx,ly,lz,ne)
c     stress
               call r_spiola(ocpp,xmj,dxmj)
c     discretized pressure
               wp=wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)
               w=wp*todet
	tot_vol=tot_vol+w
               call r_sstif(ocpp,ocuu,ocup,ne,w,toxj,body_force)
               if (iflag .eq. 0) then
                  call r_scauchy(det,todet,xto,lx,ly,lz,ne)
               endif
			enddo
		enddo
	  enddo   
	enddo
	write(*,*) 'total volume=',tot_vol
c	write(*,*) 'body force=',body_force
	write(*,*) 'pre=',pre(1,1),'predrf=',predrf(1)

c
c     pressure condensation, inverse kpp
c
      do ne=1,numel
         do i=1,nump
            tem(i)=xfp(i,ne)
            do j=1,nump
               xfrtem(i,j)=xkpp(i,j,ne)
		  enddo
	   enddo

         call gaussj(xfrtem,nump,6,tem,1,1)

         do i=1,nump
            ttm(i)=0.0d0
            do k=1,nis
               do m=1,3
                  nu1=(m-1)*nnd+nea(ne,k)
                  mu1=(m-1)*nis+k
                  ttm(i)=ttm(i)+
     $              xkup(mu1,i,ne)*du(m,nea(ne,k))
                  predrf(nu1)=predrf(nu1)+
     $              xkup(mu1,i,ne)*tem(i)
			 enddo
	       enddo
		enddo
c
c     storage
c
         do i=1,nump
            ten(i)=0.0d0
            do j=1,nump
               ten(i)=ten(i)+xfrtem(i,j)*ttm(j)
		  enddo
            pre(i,ne)=-tem(i)-ten(i)
	   enddo
     
	enddo

	write(*,*) 'pre=',pre(1,1),'predrf=',predrf(1)
      return
      end
