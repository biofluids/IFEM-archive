ccccccccccccccccccc
c  tecplotpost.f  c
c  G. Wagner      c
ccccccccccccccccccc

      subroutine tecplotpost
      include "global.h"
      include "malloc.h"

      integer ien(nen,nec)  
      real* 8 xn(nsd,nnc), dn(ndf,nnc), dfn(nnc), dd(ndf+4,nnc)
      integer plotmax, plotstep
      pointer (ienptr,ien),(xnptr, xn)
      pointer (dnptr,dn),(dfnptr,dfn),(ddptr,dd)
      character*24 title
      character*9 file
      character*16 et

      ienptr = malloc(nen*nec*isize)
      xnptr = malloc(nnc*nsd*fsize)
      dnptr = malloc(nnc*ndf*fsize)
      dfnptr = malloc(nnc*fsize)
      ddptr = malloc(nnc*(ndf+4)*fsize)

      ibase = ichar('0')        !! integer value for char '0'

      call error("read ien",-999,.false.)
      call readien(ien)

      call error("read x and rng",-999,.false.)
      call readx(xn)

      call error("write tecout.plt",-999,.false.)
      title = '3D Tet Mesh'
      call headerout(title,nsd,ndf)

      plotmax=nts
      plotstep=1
      file = "data.0000"
      write(*,*) plotmax
      do i = 0, plotmax, plotstep
c      do i=1,15
        i4 = idisk/1000
        i3 = (idisk-i4*1000)/100
        i2 = (idisk-i4*1000-i3*100)/10
        i1 = (idisk-i4*1000-i3*100-i2*10)/1
        i4 = i4 + ibase
        i3 = i3 + ibase
        i2 = i2 + ibase
        i1 = i1 + ibase
        file (6:6) = char(i4)
        file (7:7) = char(i3)
        file (8:8) = char(i2)
        file (9:9) = char(i1)

        write (*,*) file
        call postin(xn,dn,dfn,dd,file)
        call postout(idisk,xn,dd,nsd,ndf,nn,ne,nen)
        if (idisk == 0) call meshout(ien,nen,nec)
        idisk=idisk+1
      end do
      
      end
      


      subroutine tecplotpost2
      include "global.h"
      include "malloc.h"

      integer ien(nen,nec)  
      real* 8 xn(nsd,nnc), dn(ndf,nnc)
      integer plotmax, plotstep
      pointer (ienptr,ien),(xnptr, xn)
      pointer (dnptr,dn),(dfnptr,dfn),(ddptr,dd)
      character*24 title
      character*9 file
      character*16 et

      ienptr = malloc(nen*nec*isize)
      xnptr = malloc(nnc*nsd*fsize)
      dnptr = malloc(nnc*ndf*fsize)
c      dfnptr = malloc(nnc*fsize)
c      ddptr = malloc(nnc*(ndf+1)*fsize)

      ibase = ichar('0')        !! integer value for char '0'

      call error("read ien",-999,.false.)
      call readien(ien)

      call error("read x and rng",-999,.false.)
      call readx(xn)

      call error("write tecoutmesh.plt",-999,.false.)
      title = '3D Tet Mesh'
      call headerout2(title,nsd,ndf)
      write(*,*) 'idiskmsh',idiskmsh
        file = "mesh.0000"
c        idisk=idisk+1
        do i=1,nts+1
        i4 = idiskmsh/1000
        i3 = (idiskmsh-i4*1000)/100
        i2 = (idiskmsh-i4*1000-i3*100)/10
        i1 = (idiskmsh-i4*1000-i3*100-i2*10)/1
        i4 = i4 + ibase
        i3 = i3 + ibase
        i2 = i2 + ibase
        i1 = i1 + ibase
        file (6:6) = char(i4)
        file (7:7) = char(i3)
        file (8:8) = char(i2)
        file (9:9) = char(i1)

        write (*,*) file,'idiskmsh=',idiskmsh
        call postin2(xn,file)
        call postout2(idiskmsh,xn,nsd,ndf,nn,ne,nen)
        if (idiskmsh == 0) call meshout2(ien,nen,nec)
        idiskmsh=idiskmsh+1
        enddo
      end
      
