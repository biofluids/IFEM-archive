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
      
      subroutine getfdrag(fdrag,d,xloc,ien,rng,hn,hm)
      
      implicit none
      include "global.h"

      real* 8 fdrag(nsdpad)
      real* 8 d(ndf,nn_loc)
      real* 8 xloc(nsd,nn_loc)
      integer ien(nen,nec),rng(neface,nec)
      real* 8 hn(nnc),hm(nn_loc)

      integer iq,ie,ieface,inl,node,isd,jsd,ierr
      real* 8 sigma(nsdpad,nsdpad),ui_j(nsdpad,nsdpad),press,mu
      real* 8 dA,normal(nsdpad),parentarea

      real* 8 eft,det
      real* 8 sh(0:nsdpad,nenpad)      
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)
      real* 8 x(nsdpad,nenpad)

      real* 8 surfacearea

      fdrag(:) = 0
      surfacearea = 0

      do ie = 1,nec
        do ieface = 1,neface
          if (fsurf(rng(ieface,ie)).eq.1) then
c           integrate over surface
            do inl=1,nen
              do isd=1,nsd
                x(isd,inl) = xloc(isd,ien(inl,ie))
              end do
            end do
            call getnormal(normal,dA,ie,ieface,xloc,ien)
            do iq = (nquad2d*(ieface-1)+1),(nquad2d*ieface)
              if (nen.eq.4) then
                include "sh3d4n.h"
                parentarea = 0.5
              else if (nen.eq.8) then
                include "sh3d8n.h"
                parentarea = 4.0
              end if

              eft = wq2d(iq) * dA/parentarea
              surfacearea = surfacearea + eft * 1.00
              sigma(:,:) = 0.0
              ui_j(:,:) = 0.0
              press = 0.0
              mu = vis_liq
              do inl = 1,nen
                do isd = 1,nsd
                  do jsd = 1,nsd
                    ui_j(isd,jsd) =
     &                   ui_j(isd,jsd) + sh(jsd,inl)*d(isd,ien(inl,ie))
                  end do
                end do
                press = press + sh(0,inl)*d(pdf,ien(inl,ie))
              end do
              do isd = 1,nsd
                do jsd = 1,nsd
                  sigma(isd,jsd) = mu*(ui_j(isd,jsd)+ui_j(jsd,isd))
                end do
                sigma(isd,isd) = sigma(isd,isd) - press
              end do
              do isd = 1,nsd
                do jsd = 1,nsd
                  fdrag(isd) = fdrag(isd)+sigma(isd,jsd)*normal(jsd)*eft
                end do
              end do
            end do
          end if
        end do
      end do

      call MPI_ALLREDUCE(fdrag, fdrag, nsd, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, MPI_COMM_WORLD,ierr)
      
      call MPI_ALLREDUCE(surfacearea, surfacearea, 1, 
     &     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
      
      if (myid.eq.0) then
        write(*,*) "forced surface area = ",surfacearea
      end if
      
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

      real* 8 time, fdrag(nsdpad)
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

        open(ifpd, FILE="drag.dat", STATUS=FileStat, POSITION="APPEND")
        open(ifpl, FILE="lift.dat", STATUS=FileStat, POSITION="APPEND")
        
        write (ifpd, '(2F12.8)') time,fdrag(dcomp)
        write (ifpl, '(2F12.8)') time,fdrag(lcomp)
        
        close(ifpd)
        close(ifpl)
        
      end if

      return
      end


