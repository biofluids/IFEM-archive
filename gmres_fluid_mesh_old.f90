module gmres_fluid_mesh
  implicit none


  private blockgmresm

contains

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gmresm(x,id,w,bg,dg,ien,z,v,zg,avg,sm,avloc,h,y,cc,ss)
  use fluid_variables, only: nsd,nn,ne,nen,kinner,kouter
  implicit none

  real(8) :: x(nsd,nn)
  !real(8) :: d(nsd,nn)
  real(8) :: bg(nsd*nn), dg(nsd*nn), w(nsd*nn)
  integer :: id(nsd,nn),ien(nen,ne)

  real(8) :: h(kinner+1,kinner)
  real(8) :: y(kinner+1)
  real(8) :: cc(kinner)
  real(8) :: ss(kinner)

  real(8) :: z(nsd*nn,kinner)
  real(8) :: v(nsd*nn,kinner+1)
  real(8) :: zg(nsd*nn), avg(nsd*nn),sm(nsd*nn)
  real(8) :: vloc(nsd,nn),avloc(nsd,nn)
  logical :: assemble
  real(8) :: eps, rnorm, rnorm0, order
  real(8) :: gam,hsave,ysave,tmpo
  integer :: i,j,k,ij,jj,i1,j1,k1,l,iqc,igmres

  !eps = 1.0e-12
  eps = 1.0e-24
  assemble = .true.

  !do iqc=1,nsd*nn
  sm(1:nsd*nn) = 1.0
  !enddo

  do iqc=1,nsd*nn
     if (w(iqc).lt.1e-6) w(iqc) = 1.0
     w(iqc) = 1.0 / sqrt(abs(w(iqc)))
  enddo

  !call fclear (v,nsd*nn*(kinner+1))
  v(1:nsd*nn,1:kinner+1) = 0.0d0

!  compute residual as r = W**(-1/2) * (b - A * d)

  do iqc=1,nsd*nn 
     v(iqc,1) = bg(iqc)
     v(iqc,1) = w(iqc) * v(iqc,1)
  enddo

  call getnorm(v(1,1),v(1,1),nsd*nn,rnorm0)
  rnorm0 = sqrt(rnorm0)

  if(rnorm0.le.eps) then
     write(6,102)
     write(7,102)
     return
  endif


!...outer GMRES loop (igmres)
  igmres = 0 
 10 igmres = igmres + 1

 !...convergence check
  call getnorm(v(1,1),v(1,1),nsd*nn,rnorm)
  rnorm = sqrt(rnorm)
  if (rnorm.le.eps.or.igmres.gt.kouter) goto 700

 !...first krylov vector
  do iqc=1,nsd*nn
     v(iqc,1) = v(iqc,1) / rnorm
  enddo

 !...construct krylov vectors 2:inner and hessenberg matrix
  do j=1,kinner

 !...insert contribution of preconditioner
     !write(*,*) v(:,j)
     do iqc=1,nsd*nn
        zg(iqc)   = sm(iqc) *  v(iqc,j)
        z(iqc,j) = zg(iqc)
        zg(iqc)   =  w(iqc) * zg(iqc)
     enddo
     call equal(zg,vloc,nsd*nn)
     !call fclear (avloc,nsd*nn)
     avloc(1:nsd,1:nn) = 0.0d0
     call blockgmresm(x,vloc,avloc,ien)
     call equal(avloc,avg,nsd*nn)
     call setid(avg,id,nsd)

     do iqc=1,nsd*nn
        avg(iqc) = w(iqc) * avg(iqc)
        v(iqc,j+1) = avg(iqc)
        !write(*,*) avg(iqc)
     enddo

     do i=1,j
        call getnorm(v(1,j+1),v(1,i),nsd*nn,tmpo)
        h(i,j) = tmpo
        do iqc=1,nsd*nn
           v(iqc,j+1) = v(iqc,j+1) - tmpo * v(iqc,i) 
        enddo
     enddo

     call getnorm(v(1,j+1),v(1,j+1),nsd*nn,tmpo)
     !write(*,*) "tmpo",tmpo,v(1,j+1),j
     tmpo = sqrt(tmpo)
     h(j+1,j) = tmpo
     do iqc=1,nsd*nn
        v(iqc,j+1) = v(iqc,j+1) / tmpo
     enddo
     !write(*,*) "axel",v(1,j+1),tmpo

  enddo         ! end inner loop

!...compute y(1:inner+1) from local hessenberg linear system
!...H_m * y = beta * e_1

 !...initialize reduced residual
  y(1         ) = rnorm
  y(2:kinner+1) = 0.0

 !...rotations
  do j=1,kinner
     j1 = j + 1

    !...previously computed rotations on colundf j
     do i=2,j
        i1 = i - 1
        hsave = h(i1,j)
        h(i1,j) = + cc(i1) * hsave + ss(i1) * h(i,j)
        h(i ,j) = - ss(i1) * hsave + cc(i1) * h(i,j)
     enddo
    !...new rotation on colundf j
       gam = sqrt(h(j,j)**2 + h(j1,j)**2)
       cc(j)=h(j,j)/gam
       ss(j)=h(j1,j)/gam
       h(j,j) = cc(j)*h(j,j) + ss(j)*h(j1,j)
    !...note: under-diagonal term h(j+1,j) becomes 0
            
       y(j1) = - ss(j) * y(j)
       y(j ) = + cc(j) * y(j)
    !...note: till now y(j+1) = 0
       rnorm = abs(y(j1))
    end do

!...back substitution
    j = kinner      !should reach here straight from rotation loop
    y(j) = y(j)/h(j,j)
    do jj=2,j
       k = j - jj + 1
       k1 = k + 1
       ysave = y(k)
       do l=k1,j
          ysave = ysave - h(k,l) * y(l)
       end do
       y(k) = ysave / h(k,k)
    end do

!...compute dg iterate dg = dg + Z_m * y
  j = kinner        !(PVM only fix)
  do jj=1,j
     tmpo = y(jj)
     do iqc=1,nsd*nn
        dg(iqc) = dg(iqc) + tmpo * z(iqc,jj)
     enddo
  enddo
!...if not done recover residual vector
  if (igmres.le.kouter) then
     do jj=1,j
        ij = j - jj + 2
        y(ij-1) = -ss(ij-1)*y(ij)
        y(ij) = cc(ij-1)*y(ij)
     end do
     do jj=1,j+1
        tmpo = y(jj)
        if (jj.eq.1) tmpo = tmpo - 1.0
        do iqc=1,nsd*nn
           v(iqc,1) = v(iqc,1) + tmpo*v(iqc,jj)
        enddo
     enddo
  endif

  goto 10

 700    continue

!...go back to unscaled system
  do iqc=1,nsd*nn
     dg(iqc) = w(iqc) * dg(iqc)
  enddo

  order = 0.43429448 * log(rnorm0/rnorm)
  write(6,101) order
  write(7,101) order

 101    format('Mesh     : convergence order = ', f5.2)
 102    format('Mesh     : zero residual')

  return
end subroutine gmresm


!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   S. Aliabadi                                                          c
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine blockgmresm(xloc,dloc,p,ien)
  use fluid_variables, only: nn,ne,nsd,nen,nsdpad,nenpad,landa_over_mu,nquad,sq,wq
  use global_constants
  implicit none

  integer,intent(in) :: ien(nen,ne)
  real(8),intent(in) :: xloc(nsd,nn)
  real(8),intent(in) :: dloc(nsd,nn)
  real(8) :: p(nsd,nn)

  real(8) :: x(nsdpad,nenpad)
  real(8) :: d(nsdpad,nenpad)

  real(8) :: eft0,det
  real(8) :: sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
  real(8) :: xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

  real(8) :: drx(nsdpad),dry(nsdpad),drz(nsdpad)
  real(8) :: ttt,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
  real(8) :: mu,la !,landa_over_mu
  integer :: inl, ie, isd, iq, node


  mu = 1.0
  la = mu * landa_over_mu
  do ie=1,ne
     do inl=1,nen
        do isd=1,nsd
           x(isd,inl) = xloc(isd,ien(inl,ie))
           d(isd,inl) = dloc(isd,ien(inl,ie))
        enddo
     enddo
     do iq=1,nquad
        if (nen.eq.4) then
           include "sh3d4n.h"
        else if (nen.eq.8) then
           include "sh3d8n.h"
        end if

        eft0 = abs(det) * wq(iq)
!         eft0 = wq(iq)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        CALCULATE k*delta(d)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.........  initialize variables
        do isd = 1,nsd
           drx(isd) = 0.0
           dry(isd) = 0.0
           drz(isd) = 0.0
        enddo

        do inl=1,nen
!c............... calculate the first derivative
           do isd=1,nsd
              drx(isd)=drx(isd)+sh(1,inl)*d(isd,inl)      
              dry(isd)=dry(isd)+sh(2,inl)*d(isd,inl)      
              drz(isd)=drz(isd)+sh(3,inl)*d(isd,inl)    
           enddo
        enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do inl=1,nen
           ph(1,inl) = sh(1,inl)*eft0
           ph(2,inl) = sh(2,inl)*eft0
           ph(3,inl) = sh(3,inl)*eft0
        enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.....Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ttt = la*(drx(xsd)+dry(ysd)+drz(zsd))
        txx = mu*(drx(xsd)+drx(xsd))
        tyx = mu*(dry(xsd)+drx(ysd))
        tzx = mu*(drz(xsd)+drx(zsd))
        txy = mu*(drx(ysd)+dry(xsd))
        tyy = mu*(dry(ysd)+dry(ysd))
        tzy = mu*(drz(ysd)+dry(zsd))
        txz = mu*(drx(zsd)+drz(xsd))
        tyz = mu*(dry(zsd)+drz(ysd))
        tzz = mu*(drz(zsd)+drz(zsd))
       
        do inl=1,nen
           node=ien(inl,ie)
!c.....Elastic Equation (calculate k*(delta(d))=p)
           p(xsd,node) = p(xsd,node) +          &
                           ph(xsd,inl) * ttt +  &
                           ph(xsd,inl) * txx +  &
                           ph(ysd,inl) * tyx +  &
                           ph(zsd,inl) * tzx 
           p(ysd,node) = p(ysd,node) +          &
                           ph(ysd,inl) * ttt +  &
                           ph(xsd,inl) * txy +  &
                           ph(ysd,inl) * tyy +  &
                           ph(zsd,inl) * tzy 
           p(zsd,node) = p(zsd,node) +          &
                           ph(zsd,inl) * ttt +  &
                           ph(xsd,inl) * txz +  &
                           ph(ysd,inl) * tyz +  &
                           ph(zsd,inl) * tzz 
        enddo
     enddo
  enddo
  return
end subroutine blockgmresm

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       S. K. ALIABADI
!       modified for linear elastic equation, Lucy Zhang 4/22/99
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine blockm(xloc, dloc, fluid_mesh_force, w, p, ien)
  use fluid_variables
  use global_constants
  implicit none

  integer,intent(in) :: ien(nen,ne)
  real(8) :: xloc(nsd,nn) !...initial mesh position
  real(8) :: dloc(nsd,nn) !...unknown vector, here 3 displacements
  !real(8) :: eloc(nsd,nn)
  real(8) :: fluid_mesh_force(nsd,nn)  !...influences the mesh movement
  real(8) :: x(nsdpad,nenpad), d(nsdpad,nenpad) !, e(nsdpad,nenpad)
  real(8) :: p(nsd,nn)
  real(8) :: fq(nsdpad)

  real(8) :: eft0,det
  real(8) :: sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
  real(8) :: xr(nsdpad,nsdpad), cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

  real(8) :: drx(nsdpad),dry(nsdpad),drz(nsdpad)
  real(8) :: ttt,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz

  real(8) :: mu,la
  real(8) :: erx,ery,erz

!....   stiffness matrix
  real(8) :: w(nsd,nn)  
  integer :: inl, ie, isd, iq, node


  mu = 1.0
  la = mu * landa_over_mu
  do ie=1,ne 
     do inl=1,nen
        node = ien(inl,ie)
        !e(isd,inl) = eloc(isd,node)
        x(1:nsd,inl)     = xloc(1:nsd,node)
        d(1:nsd,inl)     = dloc(1:nsd,node)
     enddo

     do iq=1,nquad
        if (nen.eq.4) then
           include "sh3d4n.h"
        else if (nen.eq.8) then
           include "sh3d8n.h"
        end if

        eft0 = abs(det) * wq(iq)
!.... no jacobian
!         eft0 = wq(iq)   
!............ Calculate local Stiffness Matrix k (kd=f)
        fq(1:nsd) = 0.0d0
        do inl=1,nen
           txx = sh(1,inl)**2
           tyy = sh(2,inl)**2
           tzz = sh(3,inl)**2
           ttt = txx + tyy + tzz 
           node = ien(inl,ie)

           fq(1:nsd) = fq(1:nsd) + sh(0,inl) * fluid_mesh_force(1:nsd,node)

           w(xsd,node)=w(xsd,node)+mu*(ttt+txx)*eft0
           w(ysd,node)=w(ysd,node)+mu*(ttt+tyy)*eft0
           w(zsd,node)=w(zsd,node)+mu*(ttt+tzz)*eft0  
         
           w(xsd,node)=w(xsd,node)+la*txx*eft0
           w(ysd,node)=w(ysd,node)+la*tyy*eft0
           w(zsd,node)=w(zsd,node)+la*tzz*eft0                  
        enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        CALCULATE kd-f
!c        ONLY CALCULATE  kd in local coordinates
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.........  initialize variables
        erx = 0.0
        ery = 0.0
        erz = 0.0

        do isd = 1,nsd
           drx(isd) = 0.0
           dry(isd) = 0.0
           drz(isd) = 0.0
        enddo

        do inl=1,nen
!          erx = erx + sh(0,inl)*e(1,inl)
!          ery = ery + sh(0,inl)*e(2,inl)
!          erz = erz + sh(0,inl)*e(3,inl)

!............... calculate the first derivative
           do isd=1,nsd
              drx(isd)=drx(isd)+sh(1,inl)*d(isd,inl)      
              dry(isd)=dry(isd)+sh(2,inl)*d(isd,inl)      
              drz(isd)=drz(isd)+sh(3,inl)*d(isd,inl)
           enddo
        enddo
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do inl=1,nen
           ph(1,inl) = sh(1,inl)*eft0
           ph(2,inl) = sh(2,inl)*eft0
           ph(3,inl) = sh(3,inl)*eft0
        enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.....Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ttt = la*(drx(xsd)+dry(ysd)+drz(zsd))
        txx = mu*(drx(xsd)+drx(xsd))
        tyx = mu*(dry(xsd)+drx(ysd))
        tzx = mu*(drz(xsd)+drx(zsd))
        txy = mu*(drx(ysd)+dry(xsd))
        tyy = mu*(dry(ysd)+dry(ysd))
        tzy = mu*(drz(ysd)+dry(zsd))
        txz = mu*(drx(zsd)+drz(xsd))
        tyz = mu*(dry(zsd)+drz(ysd))
        tzz = mu*(drz(zsd)+drz(zsd))

        do inl=1,nen
           node=ien(inl,ie)
!.....Elastic Equation (calculate residual: r=kd=p)
           p(xsd,node) = p(xsd,node) -  &
              ph(xsd,inl) * ttt -  &
              ph(xsd,inl) * txx -  &
              ph(ysd,inl) * tyx -  &
              ph(zsd,inl) * tzx    &
            + sh(0,inl) * fq(xsd)

           p(ysd,node) = p(ysd,node) -  &
              ph(ysd,inl) * ttt -  &
              ph(xsd,inl) * txy -  &
              ph(ysd,inl) * tyy -  &
              ph(zsd,inl) * tzy    &
            + sh(0,inl) * fq(ysd)

           p(zsd,node) = p(zsd,node) -  &
              ph(zsd,inl) * ttt -  &
              ph(xsd,inl) * txz -  &
              ph(ysd,inl) * tyz -  &
              ph(zsd,inl) * tzz    &
            + sh(0,inl) * fq(zsd)

        enddo
     enddo
  enddo
  return
end subroutine blockm


end module gmres_fluid_mesh