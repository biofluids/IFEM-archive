c     element tangent stiffness matrix assemblage
c
      subroutine r_sstif(ocpp,ocuu,ocup,ne,w,toxj)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      include 'main_common'
      dimension ocuu(6,6),ocup(6),toxj(2,2),
     $     xac(2),xve(2)

      if (iti .le. nmove) then
         sdensit=0.0d0
      else
         sdensit=sdensi
      endif
ccccccccccc
c     i-u
ccccccccccc

      do 11 ni=1,nis
         do 12 i=1,2
            nu1=(i-1)*nnd+nea(ne,ni)
            mu1=(i-1)*nis+ni
ccccccccccc
c     fu
ccccccccccc
            call r_scalfu(fu,i,ni)
            predrf(nu1)=predrf(nu1)-fu*w
            do 10 nk=1,nump
               call r_scalkup(fkup,ocup,i,nk,ni)
               xkup(mu1,nk,ne)=xkup(mu1,nk,ne)+
     $              fkup*w
 10         continue
 12      continue
 11   continue
c
      do 13 i=1,nump
         call r_scalfp(fp,ocpp,i)
         xfp(i,ne)=xfp(i,ne)+fp*w
         do 14 j=1,nump
            call r_scalkpp(fkpp,ocpp,i,j)
            xkpp(i,j,ne)=xkpp(i,j,ne)+fkpp*w
 14      continue
 13   continue
cccccccccccccccccccccccccc
c     inertia forces
cccccccccccccccccccccccccc
      do 624 i=1,2
         xac(i)=0.0d0
         xve(i)=0.0d0
         do 625 k=1,nis
            ntem=nea(ne,k)
            xac(i)=xac(i)+h(k)*acm(i,ntem)
            xve(i)=xve(i)+h(k)*du(i,ntem)
 625     continue
 624  continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xac(1)=xac(1)+xmg1
      xac(2)=xac(2)+xmg2

      do 730 ni=1,nis
         nu1=nea(ne,ni)
         nv1=nnd+nea(ne,ni)
         predrf(nu1)=predrf(nu1)-w*sdensit*h(ni)*xac(1)
     $        -w*xviss*h(ni)*xve(1)
         predrf(nv1)=predrf(nv1)-w*sdensit*h(ni)*xac(2)
     $        -w*xviss*h(ni)*xve(2)
 730  continue
      return
      end



