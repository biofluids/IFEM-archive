      subroutine getnqdf(rng)

      implicit none
      include "global.h"

      integer rng(neface,nec)
      integer i,ie,ieface

      nqdf = 0
      
      do ie = 1,nec
        do ieface = 1,neface
          if (fsurf(rng(ieface,ie)).eq.1) then
            nqdf = nqdf + 1
          end if
        end do
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine getfdrag(fdrag,dn,d,xloc,ien,rng,hn,hm,finv)
      
      implicit none
      include "global.h"

      real* 8 fdrag(nsdpad)
      real* 8 dn(ndf,nnc),d(ndf,nn_loc)
      real* 8 xloc(nsd,nn_loc)
      real* 8 sh(0:nsd,nen),det
c      real* 8 shrkf(0:nsd,maxconn,nqdf)
c      integer cnn(maxconn,nqdc),ncnn(nqdc)
      integer ien(nen,nec),rng(neface,nec)
      real* 8 hn(nnc),hm(nn_loc)

      integer qp,ie,ieface,inl,node,isd,jsd,ksd,ierr,iq
      real* 8 sigma(nsdpad,nsdpad),ui_j(nsdpad,nsdpad),press,mu
      real* 8 finv(nsd,nsd,nquad,nec)
      real* 8 dA,normal(nsdpad),x(nsdpad,nen)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

c      call grab_all(d,dn,ndf,hn,hm)
      call gather(d,dn,ndf,hn,hm)
      fdrag(:) = 0
      qp = 0

      do ie = 1,nec
        do ieface = 1,neface
          if (fsurf(rng(ieface,ie)).eq.1) then
c            qp = qp + 1
            sigma(:,:) = 0
            ui_j(:,:) = 0
            press     = 0
            mu = vis_liq
            do inl=1,nen
               node=ien(inl,ie)
               do isd=1,nsd
                  x(isd,inl) = xloc(isd,ien(inl,ie))
               enddo
            enddo

            do iq=1,nquad

              if (nen.eq.4) then
                 include "sh3d4n.h"
              else if (nen.eq.8) then
                 include "sh3d8n.h"
              end if
 
              do inl = 1,nen
                  node = ien(inl,ie)
                  do isd = 1,nsd
                     do jsd = 1,nsd
                        do ksd = 1,nsd
                           ui_j(isd,jsd) = 
c     &                       ui_j(isd,jsd) + sh(jsd,inl)*d(isd,node)
     +                ui_j(isd,jsd)+sh(ksd,inl)*d(isd,node)*finv(ksd,jsd,iq,ie)
                        enddo
                     end do
                  end do
                  press = press + sh(0,inl)*d(pdf,node)
               end do
            enddo

            do isd = 1,nsd
              do jsd = 1,nsd
                sigma(isd,jsd) = mu*(ui_j(isd,jsd)+ui_j(jsd,isd))
              end do
              sigma(isd,isd) = sigma(isd,isd) - press
            end do
            call getnormal(normal,dA,ie,ieface,xloc,ien)
            do isd = 1,nsd
              do jsd = 1,nsd
                fdrag(isd) = fdrag(isd)+sigma(isd,jsd)*normal(jsd)*dA
              end do
            end do
c            write(*,*) fdrag(1),fdrag(2),fdrag(3),sigma(1,1),sigma(2,2)
          end if
        end do
      end do

      call MPI_ALLREDUCE(fdrag, fdrag, nsd, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, MPI_COMM_WORLD,ierr)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getnormal(normal,dA,ie,ieface,xloc,ien)

      implicit none
      include "global.h"

      real* 8 normal(nsdpad), dA
      integer ie,ieface
      real* 8 xloc(nsd,nn_loc)
      integer ien(nen,nec)

      real* 8 xnode(nsdpad,4)
      integer inface,node
      real* 8 x12,x13,x14,y12,y13,y14,z12,z13,z14
      real* 8 nnorm

      do inface = 1,nnface
        node = ien(map(ieface,inface,etype),ie)
        xnode(1:nsd,inface) = xloc(1:nsd,node)
      end do
      
      x12 = xnode(1,2) - xnode(1,1)
      x13 = xnode(1,3) - xnode(1,1)
      y12 = xnode(2,2) - xnode(2,1)
      y13 = xnode(2,3) - xnode(2,1)
      z12 = xnode(3,2) - xnode(3,1)
      z13 = xnode(3,3) - xnode(3,1)

      dA = 0.0

      if (nnface.eq.4) then
        x14 = xnode(1,4)-xnode(1,1)
        y14 = xnode(2,4)-xnode(2,1)
        z14 = xnode(3,4)-xnode(3,1)
        normal(1) = z13*y14 - y13*z14
        normal(2) = x13*z14 - z13*x14
        normal(3) = y13*x14 - x13*y14
        dA = 0.5 * sqrt(normal(1)*normal(1)+
     &       normal(2)*normal(2)+normal(3)*normal(3))
      end if

      normal(1) = z12*y13 - y12*z13
      normal(2) = x12*z13 - z12*x13
      normal(3) = y12*x13 - x12*y13
      
      nnorm = sqrt(normal(1)*normal(1)+
     &     normal(2)*normal(2)+normal(3)*normal(3))
      normal(1:nsd) = normal(1:nsd)/nnorm
        
      dA = dA + 0.5 * nnorm

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fdragout(time,fdrag,dcomp,lcomp)

      implicit none
      include "global.h"

      real* 8 time, fdrag(nsdpad),deno
      integer dcomp,lcomp

      integer ifpd, ifpl
      character*8 FileStat

      ifpd = 12
      ifpl = 13

      if (time.eq.0.0) then
        FileStat = "NEW"
      else
        FileStat = "OLD"
      end if
      if (myid.eq.0) then

        open(ifpd, FILE = "drag.dat", STATUS = FileStat, POSITION = "APPEND")
        open(ifpl, FILE = "lift.dat", STATUS = FileStat, POSITION = "APPEND")

c... for cylinder only
        deno = 0.5*den_liq*ic(1)**2*4*1.5        

        write (ifpd, '(3F12.8)') time,fdrag(dcomp),fdrag(dcomp)/deno
        write (ifpl, '(3F12.8)') time,fdrag(lcomp),fdrag(lcomp)/deno
        
        close(ifpd)
        close(ifpl)
        
      end if

      return
      end


