      subroutine calsts_hyp(xmk,disp,
     &                      sts_PK1,effstr,
     &                      mgk,ik)
c
c----------------------------------------------------c
c
c   Hyperelastic Material Subroutine
c
c** for rubber
c
c
c
c
c
c
c-----------------------------------------------------c
c
      implicit double precision (a-h,o-z)
      include 'parameter.h'
c
      real*8 xmk(2,maxGP)
      real*8 disp(2,maxNumnp)
      real*8 sts_PK1(3,3)
      real*8 defgra(3,3)
c
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      real*8 shpk(mnsch,maxGP),shpkdx(mnsch,maxGP),
     &       shpkdy(mnsch,maxGP),shpn(mnsch,maxNumnp)
c
      common /hyperelast/rho0,cc1,cc2,rambda
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /shapeK/shpk,shpkdx,shpkdy,shpn
c
c
        do ip = 1, 3
           do jp = 1, 3
              defgra(ip,jp) = 0.0d0
           enddo
        enddo
c
        mLoop = Lgp(ik)
        dgur = 0.0d0
c       
        do ir = 1, mLoop
          jp = Lmgp(ir,ik)
          defgra(1,1) = defgra(1,1) 
     &                + shpkdx(ir,ik)*disp(1,jp)
          defgra(1,2) = defgra(1,2) 
     &                + shpkdy(ir,ik)*disp(1,jp)
          defgra(2,1) = defgra(2,1) 
     &                + shpkdx(ir,ik)*disp(2,jp)
          defgra(2,2) = defgra(2,2) 
     &                + shpkdy(ir,ik)*disp(2,jp)
          dgur = dgur + shpk(ir,ik)*disp(1,jp)
        enddo
c        
c.....Plane Strain     
c        
        defgra(3,3) = 1.0
        defgra(1,1) = defgra(1,1) + 1.0
        defgra(2,2) = defgra(2,2) + 1.0

        dgjac = defgra(3,3)*defgra(1,1)*defgra(2,2)
     &         -defgra(3,3)*defgra(1,2)*defgra(2,1)

      gm11 = defgra(1,1)*defgra(1,1)+defgra(1,2)*defgra(1,2)
      gm12 = defgra(1,1)*defgra(2,1)+defgra(1,2)*defgra(2,2)
      gm21 = defgra(2,1)*defgra(1,1)+defgra(2,2)*defgra(1,2)
      gm22 = defgra(2,1)*defgra(2,1)+defgra(2,2)*defgra(2,2)
      gm33 = defgra(3,3)*defgra(3,3)
      gs11 = gm11*gm11+gm12*gm21
      gs12 = gm11*gm12+gm12*gm22
      gs21 = gm21*gm11+gm22*gm21
      gs22 = gm21*gm12+gm22*gm22
      gs33 = gm33*gm33
c
      t1 = gm11+gm22+gm33
      t3 = gm33*(gm11*gm22-gm12*gm21)
c
      cons1 = cc1+t1*cc2
      cons2 = cc1*t3**(1.0/3.0)+2*cc2*t3**(2.0/3.0)
     &        -rambda*dlog(t3)
c
      str11 = 2*(cons1*gm11-cc2*gs11-cons2)/dsqrt(t3)
      str12 = 2*(cons1*gm12-cc2*gs12      )/dsqrt(t3)
      str21 = 2*(cons1*gm21-cc2*gs21      )/dsqrt(t3)
      str22 = 2*(cons1*gm22-cc2*gs22-cons2)/dsqrt(t3)
      str33 = 2*(cons1*gm33-cc2*gs33-cons2)/dsqrt(t3)
c
      effstr = dsqrt(str11**2+str22**2+str33**2
     &        -str11*str22-str22*str33-str33*str11
     &        +3*(str12**2))
c
      dgjac = defgra(3,3)*defgra(1,1)*defgra(2,2)
     &       -defgra(3,3)*defgra(1,2)*defgra(2,1)
c
      dginv11 = defgra(2,2)*defgra(3,3)/dgjac
      dginv12 =-defgra(1,2)*defgra(3,3)/dgjac
      dginv21 =-defgra(2,1)*defgra(3,3)/dgjac
      dginv22 = defgra(1,1)*defgra(3,3)/dgjac
      dginv33 = 1/defgra(3,3)
c
      taus11 = dgjac*(dginv11*str11+dginv12*str21)
      taus12 = dgjac*(dginv11*str12+dginv12*str22)
      taus21 = dgjac*(dginv21*str11+dginv22*str21)
      taus22 = dgjac*(dginv21*str12+dginv22*str22)
      taus33 = dgjac*(dginv33*str33)

      sts_PK1(1,1) = taus11
      sts_PK1(2,1) = taus21
      sts_PK1(1,2) = taus12
      sts_PK1(2,2) = taus22
      sts_PK1(3,3) = taus33
      sts_PK1(1,3) = 0.0d0
      sts_PK1(2,3) = 0.0d0
      sts_PK1(3,1) = 0.0d0
      sts_PK1(3,2) = 0.0d0

      return
      end
