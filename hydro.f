      subroutine hydrostatic(xloc,floc,p,rng,ien)
      use fluid_variables

      implicit none

      integer rng(neface,ne),ien(nen,ne)       
      real(8) xloc(nsd,nn),floc(nn),p(ndf,nn)

      real(8) x(nsdpad,10),f(10)
      real(8) xr(nsdpad,nsdpad),sh(10) 
      real(8) g,gr,ro,pp
      real(8) n1,n2,n3,c1,c2,c3,h
      integer inl,ie,isd,ieface,inface,node
      logical yes
c     integer ierr,io,status(MPI_STATUS_SIZE)

      yes = .false.
      if(interface(1)+999.gt.0.00001) then
         h = interface(1)
         gr=   gravity(1)
         yes = .true.
      endif
      if(interface(2)+999.gt.0.00001) then
         h = interface(2)
         gr=   gravity(2)
         yes = .true.
      endif
      if(interface(3)+999.gt.0.00001) then
         h = interface(3)
         gr=   gravity(3)
         yes = .true.
      endif
      if(.not.yes) return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ie=1,ne
         do ieface = 1,neface
            if(rng(ieface,ie).eq.hydro) then
           do inface=1,nnface
              inl = mapping(ieface,inface,etype)
              do isd = 1,nsd
                 x(isd,inface) = xloc(isd,ien(inl,ie))
              enddo
              f(inface) = floc(ien(inl,ie))
           enddo

           c1 = 0.0
           c2 = 0.0
           c3 = 0.0
           g  = 0.0
           do inface=1,nnface
              c1 = c1 + x(1,inface)
              c2 = c2 + x(2,inface)
              c3 = c3 + x(3,inface)
              g  = g  + f(  inface)
           enddo
           c1 = c1/nnface
           c2 = c2/nnface
           c3 = c3/nnface
           g  = g /nnface

           ro = ((1.0-g)*den_gas+g*den_liq)

cccccccc    since gravity is negative, then (-)(-) = +
           if(interface(1)+999.gt.0.00001) pp = ro*gr*(c1-h)
           if(interface(2)+999.gt.0.00001) pp = ro*gr*(c2-h)
           if(interface(3)+999.gt.0.00001) pp = ro*gr*(c3-h)

           if (nen.eq.4) then
              xr(1,1)=x(1,1)-x(1,3)
              xr(2,1)=x(2,1)-x(2,3)
              xr(3,1)=x(3,1)-x(3,3)
              xr(1,2)=x(1,2)-x(1,3)
              xr(2,2)=x(2,2)-x(2,3)
              xr(3,2)=x(3,2)-x(3,3)
              n1=0.5*(-xr(2,1)*xr(3,2)+xr(2,2)*xr(3,1))
              n2=0.5*(-xr(3,1)*xr(1,2)+xr(1,1)*xr(3,2))
              n3=0.5*(-xr(1,1)*xr(2,2)+xr(1,2)*xr(2,1))
              sh(1)= 0.3333333
              sh(2)= 0.3333333
              sh(3)= 0.3333333
           else
              xr(1,1)=0.5*(x(1,2)+x(1,3)-x(1,1)-x(1,4))
              xr(2,1)=0.5*(x(2,2)+x(2,3)-x(2,1)-x(2,4))
              xr(3,1)=0.5*(x(3,2)+x(3,3)-x(3,1)-x(3,4))
              xr(1,2)=0.5*(x(1,3)+x(1,4)-x(1,1)-x(1,2))
              xr(2,2)=0.5*(x(2,3)+x(2,4)-x(2,1)-x(2,2))
              xr(3,2)=0.5*(x(3,3)+x(3,4)-x(3,1)-x(3,2))
              n1=(-xr(2,1)*xr(3,2)+xr(2,2)*xr(3,1))
              n2=(-xr(3,1)*xr(1,2)+xr(1,1)*xr(3,2))
              n3=(-xr(1,1)*xr(2,2)+xr(1,2)*xr(2,1))
              sh(1)=0.25
              sh(2)=0.25
              sh(3)=0.25
              sh(4)=0.25
           endif

           n1 = abs(n1)
           n2 = abs(n2)
           n3 = abs(n3)

           do inface = 1,nnface
              inl = mapping(ieface,inface,etype)
              node = ien(inl,ie)
              p(1,node) = p(1,node) - sh(inface)*pp*n1
              p(2,node) = p(2,node) - sh(inface)*pp*n2
              p(3,node) = p(3,node) - sh(inface)*pp*n3
           enddo
            endif
         enddo
      enddo

      return
      end
