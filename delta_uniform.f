      subroutine delta_rkpm_uniform(deltaeachaxis,ipt,dlptlocal_number,
     +     coord_pt,deltaeachaxis)
      integer n,ipt,idim,ix,iz
      real* 8 arg(3),rad(3)
      integer dlptlocal_number,dlptlocal_head
      real* 8 coord_pt, deltaeachaxis(8,3)
      real* 8 xt3,xt0,xt1,xt2,xcont1,xarg
      
c
c     cubic spline
c dlptlocal_number = # of points in the influence domain
c ipt = current calculating point

c      ipt = dlptlocal_head

      do n = 1, dlptlocal_number ! loop over the influence domain
         do idim = ix, iz
            arg(idim) = realpart(coord_pt(idim,ipt))

            deltaeachaxis(3,idim)= (1.0d0/6.0d0*
     $           arg(idim)**3)
            
            xt3=deltaeachaxis(3,idim)
            
            xarg=arg(idim)+1.0d0
            xcont1=(27.0d0/17.0d0-30.0d0/17.0d0*xarg**2)/1.0d0
            
            deltaeachaxis(0,idim)= (1.0d0/6.0d0*
     $           (1.0d0-arg(idim))**3)
            xt0=deltaeachaxis(0,idim)
            
            xarg=arg(idim)-1.0d0
            xcont1=(27.0d0/17.0d0-30.0d0/17.0d0*xarg**2)/1.0d0
            
            deltaeachaxis(2,idim)= (2.0d0/3.0d0-
     $           (arg(idim)-1.0d0)**2*(0.5d0+arg(idim)/2.0d0)) 
            xt2=deltaeachaxis(2,idim)
            
            xarg=arg(idim)
            xcont1=(27.0d0/17.0d0-30.0d0/17.0d0*xarg**2)/1.0d0
            
            deltaeachaxis(1,idim)=(2.0d0/3.0d0-
     $           arg(idim)**2*(1.0d0-arg(idim)/2.0d0))
            xt1=deltaeachaxis(1,idim)

         enddo
      enddo

      return
      end


ccccccccccccccccccccccccccccccccccccccccc
      subroutine delta_original(deltaeachaxis,ipt,dlptlocal_number,coord_pt)
ccccccccccccccccccccccccccccccccccccccccc

      integer n,ipt,idim,ix,iz
      real* 8 arg(3),rad(3)
      integer dlptlocal_number,dlptlocal_head
      real* 8 coord_pt(3,nn), deltaeachaxis(8,3)
      real* 8 xt3,xt0,xt1,xt2,xcont1,xarg
c
c     original delta function
c

      do n = 1, dlptlocal_number
         do idim = ix, iz
            arg(idim)   = realpart(coord_pt(idim,ipt))
            rad(idim)   = sqrt(1.0d0+4.0d0*arg(idim)*(1.0d0-arg(idim)))
            deltaeachaxis(3,idim)= (1.0d0 + 2.0d0*arg(idim) 
     $           - rad(idim)) / 8.0d0
            deltaeachaxis(0,idim)= (3.0d0 - 2.0d0*arg(idim)
     $           - rad(idim)) / 8.0d0
            deltaeachaxis(2,idim)= 0.25d0 + (0.25d0 - 
     $           deltaeachaxis(0,idim))
            deltaeachaxis(1,idim)= 0.25d0 + (0.25d0 - 
     $           deltaeachaxis(3,idim))
         enddo
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccc
