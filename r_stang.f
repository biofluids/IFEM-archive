      Subroutine r_stang
      implicit real*8 (a-h,o-z)
      include 'r_common'
      dimension x(2,9),toc(3,3),xto(2,2),xot(2,2),
     $     xj(2,2),xji(2,2),rs(2),toxj(2,2),toxji(2,2)
      dimension xmj(3),xmi(3),dxmj(3,6),ddxmj(3,6,6),
     $     obc(6,6),ocuu(6,6),ocup(6)
      dimension xfrtem(6,6),tem(6),ten(6),ttm(6)
c
      do 19 i=1,2*nnd
         predrf2(i)=predrf(i)
         predrf(i)=0.0d0
         drf2(i)=drf(i)
c     drf(i)=0.0d0
 19   continue
c     
      do 99 ne=1,numel
c     position
         do 70 j=1,nump
            xfp(j,ne)=0.0d0
            do 71 i=1,nis
               xkup(i,j,ne)=0.0d0
               xkup(i+nis,j,ne)=0.0d0
 71         continue
            do 73 i=1,nump
               xkpp(i,j,ne)=0.0d0
 73         continue
 70      continue
c     
         do 27 nos=1,nis
            ntem=nea(ne,nos)
            do 26 noj=1,2
               x(noj,nos)=coor(ntem,noj)
               y(noj,nos)=coor(ntem,noj)+dis(noj,ntem)
 26         continue
 27      continue
c     gauss integration
         do 28 lx=1,nint
            rs(1)=xg(lx,nint)
            do 29 ly=1,nint
               rs(2)=xg(ly,nint)
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
               call r_sstrain(toc,xto,lx,ly,ne)
c     stress
               call r_spiola(ocpp,xmj,dxmj)
c     discretized pressure
               wp=wgt(lx,nint)*wgt(ly,nint)*thic 
               w=wp*todet
               call r_sstif(ocpp,ocuu,ocup,ne,w,toxj)
               if (iflag .eq. 0) then
                  call r_scauchy(det,todet,xto,lx,ly,ne)
               endif
   29       continue
   28    continue   
   99 continue

c
c     pressure condensation, inverse kpp
c
      do 89 ne=1,numel
         do 1 i=1,nump
            tem(i)=xfp(i,ne)
            do 2 j=1,nump
               xfrtem(i,j)=xkpp(i,j,ne)
 2          continue
 1       continue
c
         call gaussj(xfrtem,nump,6,tem,1,1)
c
         do 51 i=1,nump
            ttm(i)=0.0d0
            do 53 k=1,nis
               do 54 m=1,2
                  nu1=(m-1)*nnd+nea(ne,k)
                  mu1=(m-1)*nis+k
                  ttm(i)=ttm(i)+
     $              xkup(mu1,i,ne)*du(m,nea(ne,k))
                  predrf(nu1)=predrf(nu1)+
     $                 xkup(mu1,i,ne)*tem(i)
 54            continue
 53         continue
 51      continue
c
c     storage
c
         do 55 i=1,nump
            ten(i)=0.0d0
            do 56 j=1,nump
               ten(i)=ten(i)+
     $              xfrtem(i,j)*ttm(j)
 56         continue
            pre(i,ne)=-tem(i)-ten(i)

 55      continue
c     
 89   continue
      return
      end
