module gmres_fluid_ale
  implicit none

contains

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  modified by L. Zhang for Total ALE formulation, 7/21/99
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gmres_ale(x,d,do,id,w,bg,dg,hg,ien,z,v,zg,avg,sm,vloc,avloc,h,y,cc,ss,fext,finv,jac,jaco,convel,convelo)
  use fluid_variables
  implicit none

  real(8) x(nsd,nn)
  real(8) d(ndf,nn), do(ndf,nn),hg(ne)
  real(8) bg(ndf*nn), dg(ndf*nn), w(ndf*nn)
  integer id(ndf,nn),ien(nen,ne)
  !real(8) hn(nn),hm(nn)

  real(8) h(inner+1,inner)
  real(8) y(inner+1)
  real(8) cc(inner), ss(inner)
    
  real(8) finv(nsd,nsd,nquad,ne),jac(nquad,ne),jaco(nquad,ne)
  !real(8) refvel(nsd,nquad,ne),refvelo(nsd,nquad,ne)
  real(8) convel(nsd,nn)  
  real(8) convelo(nsd,nn) 
  real(8) z(ndf*nn,inner)
  real(8) v(ndf*nn,inner+1)
  real(8) zg(ndf*nn), avg(ndf*nn), sm(ndf*nn)
  real(8) vloc(ndf,nn),avloc(ndf,nn)
  logical assemble
  real(8) eps, rnorm, rnorm0, order
  real(8) gam,hsave,ysave,tmpo
  integer i,j,k,ij,jj,i1,j1,k1,l,iqc,igmres
  real(8) fext(nsd,nn)

  eps = 1.0e-12
  assemble = .true.

  do iqc=1,ndf*nn
     sm(iqc) = 1.0
  enddo

  if(iscaling.eq.0) then
     do iqc=1,ndf*nn
        w(iqc) = 1.0
     enddo
  else
     do iqc=1,ndf*nn
        if (w(iqc).lt.0.0) sm(iqc) = -1.0
        w(iqc) = 1.0 / sqrt(abs(w(iqc)))
     enddo
  endif
    
!...clear arrays
  call fclear (v,ndf*nn*(inner+1))
  call fclear (vloc,ndf*nn)
    
!...compute residual as r = W**(-1/2) * (b - A * d)
    
  do iqc=1,ndf*nn 
     v(iqc,1) = bg(iqc)
     v(iqc,1) = w(iqc) * v(iqc,1)
  enddo
  call getnorm(v(1,1),v(1,1),ndf*nn,rnorm0)
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
  call getnorm(v(1,1),v(1,1),ndf*nn,rnorm)
  rnorm = sqrt(rnorm)
  if (rnorm.le.eps.or.igmres.gt.outer) goto 700
    
!...first krylov vector
  do iqc=1,ndf*nn
     v(iqc,1) = v(iqc,1) / rnorm
  enddo
    
!...construct krylov vectors 2:inner and hessenberg matrix
  do j=1,inner
       
!...insert contribution of preconditioner
     do iqc=1,ndf*nn
        zg(iqc) = sm(iqc) * v(iqc,j)
        z(iqc,j) = zg(iqc)
     enddo
       
!...compute A * v_j
     do iqc=1,ndf*nn
        zg(iqc) = w(iqc) * zg(iqc)
     enddo
     call equal(zg,vloc,ndf*nn)
     call fclear (avloc,ndf*nn)
     !call equal(avloc,avg,ndf*nn)
     call blockgmres_ale(x,d,do,vloc,avloc,hg,ien,fext,finv,jac,jaco,convel,convelo)
     call equal(avloc,avg,ndf*nn)
     call setid(avg,id,ndf)
       
     do iqc=1,ndf*nn
        avg(iqc) = w(iqc) * avg(iqc)
        v(iqc,j+1) = avg(iqc)
     enddo
     do i=1,j
        call getnorm(v(1,j+1),v(1,i),ndf*nn,tmpo)
        h(i,j) = tmpo
        do iqc=1,ndf*nn
           v(iqc,j+1) = v(iqc,j+1) - tmpo * v(iqc,i) 
        enddo
     enddo
        
     call getnorm(v(1,j+1),v(1,j+1),ndf*nn,tmpo)
     tmpo = sqrt(tmpo)
     h(j+1,j) = tmpo
     do iqc=1,ndf*nn
        v(iqc,j+1) = v(iqc,j+1) / tmpo
     enddo
       
  enddo         ! end inner loop
    
!...compute y(1:inner+1) from local hessenberg linear system
!...H_m * y = beta * e_1
    
!...initialize reduced residual
  y(1) = rnorm
  do i=2,inner+1
     y(i) = 0.0
  enddo
    
!...rotations
  do j=1,inner
       
     j1 = j + 1
       
!...previously computed rotations on colundf j
     do i=2,j
        i1 = i - 1
        hsave = h(i1,j)
        h(i1,j) = + cc(i1) * hsave + ss(i1) * h(i,j)
        h(i ,j) = - ss(i1) * hsave + cc(i1) * h(i,j)
     end do
       
!...new rotation on colundf j
     gam = sqrt(h(j,j)**2 + h(j1,j)**2)
     cc(j) = h(j,j) / gam
     ss(j) = h(j1,j) / gam
     h(j,j) = cc(j) * h(j,j) + ss(j) * h(j1,j)
    !...note: under-diagonal term h(j+1,j) becomes 0
       
     y(j1) = - ss(j) * y(j)
     y(j ) = + cc(j) * y(j)
    !...note: till now y(j+1) = 0
     rnorm = abs(y(j1))
  end do
    
    
!...back substitution
  j = inner     !should reach here straight from rotation loop
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
  j = inner     !(PVM only fix)
  do jj=1,j
     tmpo = y(jj)
     do iqc=1,ndf*nn
        dg(iqc) = dg(iqc) + tmpo * z(iqc,jj)
     enddo
  end do
    
!...if not done recover residual vector
  if (igmres.le.outer) then
     do jj=1,j
        ij = j - jj + 2
        y(ij-1) = -ss(ij-1)*y(ij)
        y(ij) = cc(ij-1)*y(ij)
     end do
     do jj=1,j+1
        tmpo = y(jj)
        if (jj.eq.1) tmpo = tmpo - 1.0
        do iqc=1,ndf*nn
           v(iqc,1) = v(iqc,1) + tmpo*v(iqc,jj)
        enddo
     end do
  end if
    
  goto 10
    
 700    continue
    
!...go back to unscaled system
  do iqc=1,ndf*nn
     dg(iqc) = w(iqc) * dg(iqc)
  enddo
    
  order = 0.43429448 * log(rnorm0/rnorm)
  write(6,101) order
  write(7,101) order
    
 101    format('Flow     : convergence order = ', f5.2)
 102    format('Flow     : zero residual')
    
  return
end subroutine gmres_ale



!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   S. Aliabadi                                                          !
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine blockgmres_ale(xloc,dloc,doloc,qloc,p,hk,ien,fext,finv,jac,jaco,convel,convelo)
  use fluid_variables
  use run_variables
  use global_constants
  implicit none

  integer ien(nen,ne)
  real(8) xloc(nsd,nn) !,floc(nn)
  real(8) dloc(ndf,nn),doloc(ndf,nn)
  real(8) p(ndf,nn),qloc(ndf,nn),hk(ne)
    
  real(8) x(nsdpad,nenpad) !,f(nenpad)
  real(8) d(ndfpad,nenpad),do(ndfpad,nenpad),q(ndfpad,nenpad)

  real(8) eft0,det
  real(8) sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
  real(8) xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

  real(8) drs(ndfpad),qrt(ndfpad),qrs(ndfpad),drs_ref(nsdpad)
  real(8) drx(ndfpad),dry(ndfpad),drz(ndfpad)
  real(8) qrx(ndfpad),qry(ndfpad),qrz(ndfpad)
  real(8) qrx_hat,qry_hat,qrz_hat
  real(8) u,v,w,pp,wx,wy,wz
  real(8) txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
  real(8) pxx,pxy,pxz,pyx,pyy,pyz,pzx,pzy,pzz
    !real(8) temprefu,temprefv,temprefw
  real(8) finv(nsdpad,nsdpad,nquad,ne),jac(nquad,ne),jaco(nquad,ne)
  real(8) t(nsdpad,nsdpad),dtemp(nsdpad,nsdpad)
  real(8) hg,taum,tauc,vel,ree
  real(8) res_c,res_a(nsdpad),res_t(nsdpad)
  real(8) prs_c,prs_t(nsdpad)
  real(8) mu,nu,ro
  real(8) tempu,tempv,tempw,temp
  real(8) dtinv,oma,ama
  integer inl, ie, isd, idf, iq, node
  real(8) finv11,finv12,finv13,finv21,finv22,finv23,finv31
  real(8) finv32,finv33
  !real(8) refvel(nsd,nquad,ne),refvelo(nsd,nquad,ne)
  !real(8) refv(nsd),refvo(nsd)
  real(8) convel(nsd,nn)  
  real(8) convelo(nsd,nn) 
  real(8) convel_e(nsd,nen)
  real(8) conv(nsd)

  real(8) fext(nsd,nn)
  real(8) fnode(nsd,nen),fq(nsd)
    
  dtinv = 1.0/dt/alpha
  if(steady) dtinv = 0.0
  oma   = 1.0 -alpha
  ama   = 1.0 - oma

  do ie=1,ne
     do inl=1,nen
        do isd=1,nsd
           x(isd,inl)     = xloc(isd,ien(inl,ie))
           convel_e(isd,inl)  = ama*convel(isd,ien(inl,ie)) + oma*convelo(isd,ien(inl,ie))
           fnode(isd,inl) = fext(isd,ien(inl,ie))
        enddo
        do idf=1,ndf
           q(idf,inl) =  qloc(idf,ien(inl,ie))
           d(idf,inl) =  dloc(idf,ien(inl,ie))
           do(idf,inl) = doloc(idf,ien(inl,ie))
        enddo
        !f(inl) = floc(ien(inl,ie))
     enddo
     hg = hk(ie)
     do iq=1,nquad
        if (nen.eq.4) then
           include "sh3d4n.h"
        else if (nen.eq.8) then
           include "sh3d8n.h"
        end if
       
        eft0 = abs(det) * wq(iq) * alpha

        !do isd=1,nsd
        !   refv(isd)=refvel(isd,iq,ie)
        !   refvo(isd)=refvelo(isd,iq,ie)
        !   if(its.eq.1) refvo(isd)=refv(isd)
        !enddo

        do idf = 1,nsd
           drs(idf) = 0.0
           drx(idf) = 0.0
           dry(idf) = 0.0
           drz(idf) = 0.0
        enddo
        conv(1:3) = 0
        do inl=1,nen
           tempu= ama*d(udf,inl) + oma*do(udf,inl)
           tempv= ama*d(vdf,inl) + oma*do(vdf,inl)
           tempw= ama*d(wdf,inl) + oma*do(wdf,inl)
           drs(udf)=drs(udf)+sh(0,inl)*tempu            
           drx(udf)=drx(udf)+sh(1,inl)*tempu           
           dry(udf)=dry(udf)+sh(2,inl)*tempu          
           drz(udf)=drz(udf)+sh(3,inl)*tempu         
           drs(vdf)=drs(vdf)+sh(0,inl)*tempv            
           drx(vdf)=drx(vdf)+sh(1,inl)*tempv           
           dry(vdf)=dry(vdf)+sh(2,inl)*tempv          
           drz(vdf)=drz(vdf)+sh(3,inl)*tempv         
           drs(wdf)=drs(wdf)+sh(0,inl)*tempw            
           drx(wdf)=drx(wdf)+sh(1,inl)*tempw           
           dry(wdf)=dry(wdf)+sh(2,inl)*tempw          
           drz(wdf)=drz(wdf)+sh(3,inl)*tempw  
           fq(:)=fq(:)+sh(0,inl)*fnode(:,inl)
           do isd=1,nsd
              conv(isd)=conv(isd)+sh(0,inl)*convel_e(isd,inl)        
           enddo       
        end do
       
        do idf=1,ndf
           qrs(idf) = 0.0
           qrx(idf) = 0.0
           qry(idf) = 0.0
           qrz(idf) = 0.0
           qrt(idf) = 0.0
        enddo
        
        do inl=1,nen
           qrs(udf)=qrs(udf)+sh(0,inl)*q(udf,inl)      
           qrx(udf)=qrx(udf)+sh(1,inl)*q(udf,inl)      
           qry(udf)=qry(udf)+sh(2,inl)*q(udf,inl)      
           qrz(udf)=qrz(udf)+sh(3,inl)*q(udf,inl)      
           qrs(vdf)=qrs(vdf)+sh(0,inl)*q(vdf,inl)      
           qrx(vdf)=qrx(vdf)+sh(1,inl)*q(vdf,inl)      
           qry(vdf)=qry(vdf)+sh(2,inl)*q(vdf,inl)      
           qrz(vdf)=qrz(vdf)+sh(3,inl)*q(vdf,inl)      
           qrs(wdf)=qrs(wdf)+sh(0,inl)*q(wdf,inl)      
           qrx(wdf)=qrx(wdf)+sh(1,inl)*q(wdf,inl)      
           qry(wdf)=qry(wdf)+sh(2,inl)*q(wdf,inl)      
           qrz(wdf)=qrz(wdf)+sh(3,inl)*q(wdf,inl)      
        enddo

        !do isd = 1,nsd
        !   drs_ref(isd) = 0.0
        !enddo

        !drs_ref(udf) = ama*refv(udf) + oma*refvo(udf)
        !drs_ref(vdf) = ama*refv(vdf) + oma*refvo(vdf)
        !drs_ref(wdf) = ama*refv(wdf) + oma*refvo(wdf)

        do inl=1,nen
           qrt(udf)=qrt(udf)+sh(0,inl)*q(udf,inl)*dtinv
           qrt(vdf)=qrt(vdf)+sh(0,inl)*q(vdf,inl)*dtinv
           qrt(wdf)=qrt(wdf)+sh(0,inl)*q(wdf,inl)*dtinv
           qrs(pdf)=qrs(pdf)+sh(0,inl)*q(pdf,inl)      
           qrx(pdf)=qrx(pdf)+sh(1,inl)*q(pdf,inl)      
           qry(pdf)=qry(pdf)+sh(2,inl)*q(pdf,inl)      
           qrz(pdf)=qrz(pdf)+sh(3,inl)*q(pdf,inl)      
        end do
          
        u = drs(udf)
        v = drs(vdf)
        w = drs(wdf)
        pp= qrs(pdf)/alpha

        !wx = drs_ref(udf)
        !wy = drs_ref(vdf)
        !wz = drs_ref(wdf)

        wx = conv(udf)
        wy = conv(vdf)
        wz = conv(wdf)

        !fi = 0.0
        !do inl=1,nen
           !fi = fi + sh(0,inl)*f(inl)
        !enddo
          
!ccccccccccc
!         mu = (1.0-fi)*vis_gas+fi*vis_liq
!         ro = (1.0-fi)*den_gas+fi*den_liq

        mu = vis_liq
        ro = den_liq
        !ro = ro*jac(iq,ie)
!ccccccccccc
         
        nu = delta(4)*turb_kappa**2*hg**2   &
              * sqrt(2*drx(udf)**2+(dry(udf)+drx(vdf))**2   &
                    +2*dry(vdf)**2+(drz(udf)+drx(wdf))**2   &
                    +2*drz(wdf)**2+(drz(vdf)+dry(wdf))**2)  
        mu = mu + nu*ro                   

        res_c=0
        do inl=1,nen
           res_c = res_c+(sh(xsd,inl)*q(udf,inl)   &
                        + sh(ysd,inl)*q(vdf,inl)   &
                        + sh(zsd,inl)*q(wdf,inl))/alpha
           !do i=1,nsd
           !   do j=1,nsd
           !      res_c=res_c+sh(i,inl)*(ro*q(j,inl))*finv(i,j,iq,ie)/alpha
           !      !res_c = res_c + sh(i,inl)*(q(j,inl))*finv(i,j,iq,ie)/alpha
           !      res_c=res_c+2*sh(i,inl)*jac(iq,ie)*q(j,inl)*finv(i,j,iq,ie)/alpha
           !   enddo
           !enddo
        enddo
   
        if(stokes) then
           u = 0.0
           v = 0.0
           w = 0.0
        endif
          
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do isd = 1, nsd
           res_a(isd)=ro*(qrt(isd)+wx*qrx(isd)+wy*qry(isd)+wz*qrz(isd))
!     +  -fq(isd)
        enddo
!         if (ie.eq.10) then
!         write(*,*) res_a(1),res_a(2),res_a(3)
!         stop
!         endif

        !finv11 = finv(1,1,iq,ie)
        !finv12 = finv(1,2,iq,ie)
        !finv13 = finv(1,3,iq,ie)
        !finv21 = finv(2,1,iq,ie)
        !finv22 = finv(2,2,iq,ie)
        !finv23 = finv(2,3,iq,ie)
        !finv31 = finv(3,1,iq,ie)
        !finv32 = finv(3,2,iq,ie)
        !finv33 = finv(3,3,iq,ie)
        !qrx_hat=jac(iq,ie)*(finv11*qrx(pdf)+finv12*qry(pdf)+finv13*qrz(pdf))
        !qry_hat=jac(iq,ie)*(finv21*qrx(pdf)+finv22*qry(pdf)+finv23*qrz(pdf))
        !qrz_hat=jac(iq,ie)*(finv31*qrx(pdf)+finv32*qry(pdf)+finv33*qrz(pdf))

        res_t(xsd) = qrx(pdf)/alpha + res_a(xsd) 
        res_t(ysd) = qry(pdf)/alpha + res_a(ysd) 
        res_t(zsd) = qrz(pdf)/alpha + res_a(zsd)
        
        !res_t(xsd) = qrx_hat/alpha + res_a(xsd) 
        !res_t(ysd) = qry_hat/alpha + res_a(ysd) 
        !res_t(zsd) = qrz_hat/alpha + res_a(zsd)

!....    calculate the stabilization parameters, taum and tauc
        vel  = sqrt(u*u+v*v+w*w)
        ree  = vel*hg/mu/12.0
        if(steady.or.(.not.taudt)) then
           taum = 1.0/sqrt((2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
        else
           taum = 1.0/sqrt((2.0/dt)**2+(2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
        endif
        taum = delta(1)* taum
        tauc = delta(2)*hg*vel
        if(ree.lt.1.0) tauc = tauc*ree
       
!.....   Density optimization
        taum = taum/ro
        tauc = tauc*ro 
        do inl=1,nen
           ph(0,inl) = sh(0,inl)*eft0
           ph(1,inl) = sh(1,inl)*eft0
           ph(2,inl) = sh(2,inl)*eft0
           ph(3,inl) = sh(3,inl)*eft0
        enddo
          
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.....   Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !dtemp(1,1) = qrx(udf)
        !dtemp(2,1) = qrx(vdf)
        !dtemp(3,1) = qrx(wdf)
        !dtemp(1,2) = qry(udf)
        !dtemp(2,2) = qry(vdf)
        !dtemp(3,2) = qry(wdf)
        !dtemp(1,3) = qrz(udf)
        !dtemp(2,3) = qrz(vdf)
        !dtemp(3,3) = qrz(wdf)
        !do i = 1,nsd
        !   do j = 1,nsd
        !      t(i,j) = 0
        !      do k = 1,nsd
!       !       t(i,j) = t(i,j) + finv(i,k,iq,ie)*dtemp(k,j)
!   1   !                      + finv(j,k,iq,ie)*dtemp(k,i)
        !         t(i,j) = t(i,j) + dtemp(i,k)*finv(k,j,iq,ie)  &
        !                         + dtemp(j,k)*finv(k,i,iq,ie)
        !      enddo
        !   enddo
        !enddo

        txx = mu*(qrx(udf)+qrx(udf))
        tyx = mu*(qry(udf)+qrx(vdf))
        tzx = mu*(qrz(udf)+qrx(wdf))
        txy = mu*(qrx(vdf)+qry(udf))
        tyy = mu*(qry(vdf)+qry(vdf))
        tzy = mu*(qrz(vdf)+qry(wdf))
        txz = mu*(qrx(wdf)+qrz(udf))
        tyz = mu*(qry(wdf)+qrz(vdf))
        tzz = mu*(qrz(wdf)+qrz(wdf))

        !txx = mu * t(1,1)
        !tyx = mu * t(2,1)
        !tzx = mu * t(3,1)
        !txy = mu * t(1,2)
        !tyy = mu * t(2,2)
        !tzy = mu * t(3,2)
        !txz = mu * t(1,3)
        !tyz = mu * t(2,3)
        !tzz = mu * t(3,3)
          
        !pxx = jac(iq,ie)*(finv11 *txx+finv12 *tyx+finv13 *tzx)
        !pyx = jac(iq,ie)*(finv21 *txx+finv22 *tyx+finv23 *tzx)
        !pzx = jac(iq,ie)*(finv31 *txx+finv32 *tyx+finv33 *tzx)
        !pxy = jac(iq,ie)*(finv11 *txy+finv12 *tyy+finv13 *tzy)
        !pyy = jac(iq,ie)*(finv21 *txy+finv22 *tyy+finv23 *tzy)
        !pzy = jac(iq,ie)*(finv31 *txy+finv32 *tyy+finv33 *tzy)
        !pxz = jac(iq,ie)*(finv11 *txz+finv12 *tyz+finv13 *tzz)
        !pyz = jac(iq,ie)*(finv21 *txz+finv22 *tyz+finv23 *tzz)
        !pzz = jac(iq,ie)*(finv31 *txz+finv32 *tyz+finv33 *tzz)

        prs_t(udf) = res_t(udf) * taum
        prs_t(vdf) = res_t(vdf) * taum
        prs_t(wdf) = res_t(wdf) * taum
        prs_c     = res_c*tauc

        do inl=1,nen

           node = ien(inl,ie)
           temp = ro*(wx*ph(xsd,inl)+wy*ph(ysd,inl)+wz*ph(zsd,inl))
         
!.....      Continuty Equation
           p(pdf,node) = p(pdf,node)+ph(0,inl)*res_c
         
!.....      Momentum Equation (Euler Residual)
           p(udf,node) = p(udf,node)+ph(0,inl)*res_a(udf)
           p(vdf,node) = p(vdf,node)+ph(0,inl)*res_a(vdf)
           p(wdf,node) = p(wdf,node)+ph(0,inl)*res_a(wdf)

!.....      Momentum Equation (C: Diffusion terms)
!           Cauchy
           p(udf,node) = p(udf,node) - ph(xsd,inl)* pp + ph(xsd,inl)*txx + ph(ysd,inl)*tyx + ph(zsd,inl)*tzx 
           p(vdf,node) = p(vdf,node) - ph(ysd,inl)* pp + ph(xsd,inl)*txy + ph(ysd,inl)*tyy + ph(zsd,inl)*tzy 
           p(wdf,node) = p(wdf,node) - ph(zsd,inl)* pp + ph(xsd,inl)*txz + ph(ysd,inl)*tyz + ph(zsd,inl)*tzz 


!....       Piola Kirchoff Stress
           !p(udf,node) = p(udf,node) -             &
           !   ph(xsd,inl)*jac(iq,ie)*finv11* pp +  &
           !   ph(xsd,inl)*pxx +                    &
           !   ph(ysd,inl)*pyx +                    &
           !   ph(zsd,inl)*pzx 
           !p(vdf,node) = p(vdf,node) -             &
           !   ph(ysd,inl)*jac(iq,ie)*finv22* pp +  &
           !   ph(xsd,inl)*pxy +                    &
           !   ph(ysd,inl)*pyy +                    &
           !   ph(zsd,inl)*pzy 
           !p(wdf,node) = p(wdf,node) -             &
           !   ph(zsd,inl)*jac(iq,ie)*finv33* pp +  &
           !   ph(xsd,inl)*pxz +                    &
           !   ph(ysd,inl)*pyz +                    &
           !   ph(zsd,inl)*pzz 

!.....      Stablization with Tau_moment
           p(pdf,node) = p(pdf,node)         &
              + ph(xsd,inl)*prs_t(udf)       &
              + ph(ysd,inl)*prs_t(vdf)       &
              + ph(zsd,inl)*prs_t(wdf)
           p(udf,node) = p(udf,node) +prs_t(udf)*temp
           p(vdf,node) = p(vdf,node) +prs_t(vdf)*temp
           p(wdf,node) = p(wdf,node) +prs_t(wdf)*temp
         
!.....      Stablization with Tau_cont    
           p(udf,node) = p(udf,node) + ph(xsd,inl)*prs_c
           p(vdf,node) = p(vdf,node) + ph(ysd,inl)*prs_c
           p(wdf,node) = p(wdf,node) + ph(zsd,inl)*prs_c
         
        enddo

        if(stokes) goto 500
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.....   Non-linear Terms (from momentum eq.) 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do isd = 1, nsd
           res_a(isd)=ro*(qrs(1)*drx(isd)+qrs(2)*dry(isd)+qrs(3)*drz(isd))
        enddo
          
        do inl=1,nen
         
           node = ien(inl,ie)        
           temp=ro*(wx*ph(xsd,inl)+wy*ph(ysd,inl)+wz*ph(zsd,inl))
         
!.......    Additional nonlinear term form Galerkin
           p(udf,node) = p(udf,node)+ph(0,inl)*res_a(udf)
           p(vdf,node) = p(vdf,node)+ph(0,inl)*res_a(vdf)
           p(wdf,node) = p(wdf,node)+ph(0,inl)*res_a(wdf)
         
!.......    Additional nonlinear term form SUPG
           p(pdf,node) = p(pdf,node) + taum*(ph(xsd,inl)*res_a(udf)+ph(ysd,inl)*res_a(vdf)+ph(zsd,inl)*res_a(wdf)) 
           p(udf,node) = p(udf,node) + res_a(udf)*temp*taum
           p(vdf,node) = p(vdf,node) + res_a(vdf)*temp*taum
           p(wdf,node) = p(wdf,node) + res_a(wdf)*temp*taum

        enddo
          
 500          continue 
     enddo
  enddo
    
  return
end subroutine blockgmres_ale







!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       S. K. ALIABADI
!       modified by L. Zhang for total ALE formulation, 7/21/99
!       modified by A. Gerstenberger for IFEM: changed to updated ALE
!       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine block_ale(xloc,dloc,doloc,p,q,hk,ien,f_fluids,def,finv,jac,jaco,convel,convelo)
  use fluid_variables
  use run_variables
  use global_constants
  implicit none

  integer :: ien(nen,ne)
  real(8) :: xloc(nsd,nn) !,floc(nn)
  real(8) :: dloc(ndf,nn),doloc(ndf,nn)
  real(8) :: p(ndf,nn),q(ndf,nn),hk(ne)
    
  real(8) :: x(nsdpad,nenpad) !,f(nenpad)
  real(8) :: d(ndfpad,nenpad),do(ndfpad,nenpad)

  real(8) :: eft0,det,effd,effm,effc
  real(8) :: sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
  real(8) :: xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

  real(8) :: drt(ndfpad),drs(ndfpad)
  real(8) :: drx(ndfpad),dry(ndfpad),drz(ndfpad)
  real(8) :: drs_ref(nsdpad)
  real(8) :: u,v,w,pp,ug,wx,wy,wz !,fi

  real(8) :: dtemp(nsd,nsd),t(nsd,nsd)
  !real(8) drx_hat,dry_hat,drz_hat
  real(8) :: txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
  real(8) :: pxx,pxy,pxz,pyx,pyy,pyz,pzx,pzy,pzz
  real(8) :: hg,taum,tauc,vel,ree
  real(8) :: res_c,res_a(nsdpad),res_t(nsdpad)
  real(8) :: prs_c,prs_t(nsdpad)
  real(8) :: mu,nu,ro,g(nsdpad)
  real(8) :: tempu,tempv,tempw,temp
  !real(8) temprefu,temprefv,temprefw
  real(8) :: dtinv,oma,ama
  integer :: inl, ie, isd, idf, iq, node
  real(8) :: finv(nsd,nsd,nquad,ne),jac(nquad,ne),jaco(nquad,ne)
  real(8) :: def(nsd,nsd,nquad,ne)
  real(8) :: finv11,finv12,finv13,finv21,finv22,finv23,finv31
  real(8) :: finv32,finv33
  real(8) :: refvel(nsd,nquad,ne),refvelo(nsd,nquad,ne)
  real(8) :: convel(nsd,nn)  
  real(8) :: convelo(nsd,nn) 
  real(8) :: convel_e(nsd,nen)
  real(8) :: conv(nsd)
  !real(8) refv(nsd),refvo(nsd)
  real(8) :: dp1(nsd),dp3(nsd),dfinv(nsd,nsd,nsd)
  real(8) :: dp(nsd),dpp(nsd)

  real(8) :: f_fluids(nsd,nn)
  real(8) :: fnode(nsd,nen),fq(nsd)

  !integer i,j,k

  dtinv = 1.0/dt
  if(steady) dtinv = 0.0
  oma   = 1.0 -alpha
  ama   = 1.0 - oma

  e_loop: do ie=1,ne 
     do inl=1,nen
        node=ien(inl,ie)
        do isd=1,nsd
           x(isd,inl)         = xloc(isd,node)
           convel_e(isd,inl)  = ama*convel(isd,node) + oma*convelo(isd,node)
           fnode(isd,inl)     = f_fluids(isd,node)
        enddo
        do idf=1,ndf
           d(idf,inl)  =  dloc(idf,node)
           do(idf,inl) = doloc(idf,node)
        enddo
        !f(inl) = floc(ien(inl,ie))
     enddo
     hg = hk(ie)

     quad_loop: do iq=1,nquad
!...  calculate the shape function and the weight at quad point
        if (nen.eq.4) then
           include "sh3d4n.h"
        else if (nen.eq.8) then
           include "sh3d8n.h"
        end if

        eft0 = abs(det) * wq(iq)

        !do isd=1,nsd
        !   refv(isd) =refvel(isd,iq,ie)
        !   refvo(isd)=refvelo(isd,iq,ie)
        !   if (its.eq.1) refvo(isd)=refv(isd)
        !enddo

!...  initialize d, dd/dx, dd/dy, dd/dz, dd/dt
        do idf = 1,ndf
           drs(idf) = 0.0
           drx(idf) = 0.0
           dry(idf) = 0.0
           drz(idf) = 0.0
           drt(idf) = 0.0
        enddo

        fq(1:3) = 0
        conv(1:3) = 0
!... calculate vi, dvi/dxj
        do inl=1,nen
           tempu = ama*d(udf,inl) + oma*do(udf,inl)
           tempv = ama*d(vdf,inl) + oma*do(vdf,inl)
           tempw = ama*d(wdf,inl) + oma*do(wdf,inl)
           drs(udf)=drs(udf)+sh(0,inl)*tempu           
           drx(udf)=drx(udf)+sh(1,inl)*tempu           
           dry(udf)=dry(udf)+sh(2,inl)*tempu           
           drz(udf)=drz(udf)+sh(3,inl)*tempu           
           drs(vdf)=drs(vdf)+sh(0,inl)*tempv           
           drx(vdf)=drx(vdf)+sh(1,inl)*tempv           
           dry(vdf)=dry(vdf)+sh(2,inl)*tempv           
           drz(vdf)=drz(vdf)+sh(3,inl)*tempv           
           drs(wdf)=drs(wdf)+sh(0,inl)*tempw           
           drx(wdf)=drx(wdf)+sh(1,inl)*tempw           
           dry(wdf)=dry(wdf)+sh(2,inl)*tempw           
           drz(wdf)=drz(wdf)+sh(3,inl)*tempw 
           fq(1:3)=fq(1:3)+sh(0,inl)*fnode(1:3,inl)
           do isd=1,nsd
              conv(isd)=conv(isd)+sh(0,inl)*convel_e(isd,inl)        
           enddo
        enddo

        !do isd = 1,nsd
        !   drs_ref(isd) = 0.0
        !enddo

        !drs_ref(udf) = ama*refv(udf) + oma*refvo(udf)
        !drs_ref(vdf) = ama*refv(vdf) + oma*refvo(vdf)
        !drs_ref(wdf) = ama*refv(wdf) + oma*refvo(wdf)

!... calculate dvi/dt, p, dp/dxi
        do inl=1,nen
           drt(udf)=drt(udf)+sh(0,inl)*(d(udf,inl)-do(udf,inl))*dtinv
           drt(vdf)=drt(vdf)+sh(0,inl)*(d(vdf,inl)-do(vdf,inl))*dtinv
           drt(wdf)=drt(wdf)+sh(0,inl)*(d(wdf,inl)-do(wdf,inl))*dtinv
           drs(pdf)=drs(pdf)+sh(0,inl)*d(pdf,inl)      
           drx(pdf)=drx(pdf)+sh(1,inl)*d(pdf,inl)      
           dry(pdf)=dry(pdf)+sh(2,inl)*d(pdf,inl)      
           drz(pdf)=drz(pdf)+sh(3,inl)*d(pdf,inl)      
        end do
!... define u=v1, v=v2, w=v3, pp=p
        u = drs(udf)
        v = drs(vdf)
        w = drs(wdf)
        pp= drs(pdf)
          
        wx = conv(udf)
        wy = conv(vdf)
        wz = conv(wdf)

          !fi = 0.0
          !do inl=1,nen
          !   fi = fi + sh(0,inl)*f(inl)
          !enddo
!....  calculate liquid constant and gravity
!cccccccccccccc
!!        mu = (1.0-fi)*vis_gas+fi*vis_liq
!!        ro = (1.0-fi)*den_gas+fi*den_liq
        mu = vis_liq
        ro = den_liq
        !ro = ro*jac(iq,ie) !ro(hat)
        g(udf) = gravity(udf)
        g(vdf) = gravity(vdf)
        g(wdf) = gravity(wdf)
!cccccccccccccc

        nu = delta(4)*turb_kappa**2*hg**2                &
           * sqrt(2*drx(udf)**2+(dry(udf)+drx(vdf))**2   &
                 +2*dry(vdf)**2+(drz(udf)+drx(wdf))**2   &
                 +2*drz(wdf)**2+(drz(vdf)+dry(wdf))**2)
        mu = mu + nu*ro                   

        res_c=0

!....  calculate each term in the residual equation
        do inl=1,nen
         res_c = res_c+sh(xsd,inl)*d(udf,inl)   &
                      +sh(ysd,inl)*d(vdf,inl)   &
                      +sh(zsd,inl)*d(wdf,inl)
           !do i=1,nsd
           !   do j=1,nsd
           !      res_c = res_c +   sh(i,inl)*        ro*d(j,inl)*finv(i,j,iq,ie)
!                res_c = res_c +   sh(i,inl)*           d(j,inl)*finv(i,j,iq,ie)
           !      res_c = res_c + 2*sh(i,inl)*jac(iq,ie)*d(j,inl)*finv(i,j,iq,ie)
           !   enddo
           !enddo
        enddo

        !write(*,*) res_c


        if(stokes) then
           u = 0.0
           v = 0.0
           w = 0.0
        endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do isd = 1, nsd
           res_a(isd) = ro*(drt(isd)+wx*drx(isd)    &
                                    +wy*dry(isd)    &
                                    +wz*drz(isd)    &
                                    -g(isd))        &
                                    -fq(isd)
        enddo

        !finv11 = finv(1,1,iq,ie)
        !finv12 = finv(1,2,iq,ie)
        !finv13 = finv(1,3,iq,ie)
        !finv21 = finv(2,1,iq,ie)
        !finv22 = finv(2,2,iq,ie)
        !finv23 = finv(2,3,iq,ie)
        !finv31 = finv(3,1,iq,ie)
        !finv32 = finv(3,2,iq,ie)
        !finv33 = finv(3,3,iq,ie)

        !dpp(1)=drx(pdf)
        !dpp(2)=dry(pdf)
        !dpp(3)=drz(pdf)

        !dtemp(1,1) = drx(udf)
        !dtemp(2,1) = drx(vdf)
        !dtemp(3,1) = drx(wdf)
        !dtemp(1,2) = dry(udf)
        !dtemp(2,2) = dry(vdf)
        !dtemp(3,2) = dry(wdf)
        !dtemp(1,3) = drz(udf)
        !dtemp(2,3) = drz(vdf)
        !dtemp(3,3) = drz(wdf)

        !do isd=1,nsd
        !   do i=1,nsd
        !      do j=1,nsd
        !         dfinv(i,j,isd)=0
        !         do inl=1,nen
        !            dfinv(i,j,isd)=dfinv(i,j,isd)+sh(isd,inl)*finv(i,j,iq,ie)
        !         enddo
        !      enddo
        !   enddo
        !enddo

        !do i = 1,nsd
        !   dp1(i) = 0
        !   dp3(i) = 0
        !   do j = 1,nsd
        !      dp1(i) = dp1(i) + dpp(j)*jac(iq,ie)*finv(j,i,iq,ie)
!           do k = 1,nsd
!              dp2(k,i,j) = 0
!              do l = 1,nsd       
!             dp2(k,i,j) = dp2(k,i,j) + dtemp(k,l)*dfinv(l,i,j)
!   1               +dtemp(i,l)*dfinv(l,k,j)
!              enddo
!              dp3(i)=dp3(i)+mu*jac(iq,ie)*finv(j,k,iq,ie)*dp2(k,i,j)
!           enddo
        !   enddo
        !   dp(i) = dp1(i) - dp3(i)
        !enddo

        !res_t(xsd) = dp(1) + res_a(xsd) 
        !res_t(ysd) = dp(2) + res_a(ysd) 
        !res_t(zsd) = dp(3) + res_a(zsd)

        res_t(xsd) = drx(pdf) + res_a(xsd) 
        res_t(ysd) = dry(pdf) + res_a(ysd) 
        res_t(zsd) = drz(pdf) + res_a(zsd)

!.....  Stablization parameters, taum and tauc
        vel  = sqrt(u*u+v*v+w*w)
        ree  = vel*hg/mu/12.0
        if(steady.or.(.not.taudt)) then
           taum = 1.0/sqrt((2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
        else
           taum = 1.0/sqrt((2.0/dt)**2+(2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
        endif
        taum = delta(1)* taum
        tauc = delta(2)*hg*vel
        if(ree.lt.1.0) tauc = tauc*ree

         !...Density optimization
        taum = taum/ro
        tauc = tauc*ro 

        do inl=1,nen
           ph(0,inl) = sh(0,inl)*eft0
           ph(1,inl) = sh(1,inl)*eft0
           ph(2,inl) = sh(2,inl)*eft0
           ph(3,inl) = sh(3,inl)*eft0
        enddo
          
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!.....   Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !do i = 1,nsd
        !   do j = 1,nsd
        !      t(i,j) = 0
        !      do k = 1,nsd
        !         !t(i,j) = t(i,j) + finv(i,k,iq,ie)*dtemp(k,j)   &
        !        !                + finv(j,k,iq,ie)*dtemp(k,i)
        !         t(i,j) = t(i,j) + dtemp(i,k)*finv(k,j,iq,ie)    &
        !                     + dtemp(j,k)*finv(k,i,iq,ie)
        !      enddo
        !   enddo
        !enddo


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.....   Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c....   Calculate stress tau(ij) and pressure
        txx = mu*(drx(udf)+drx(udf))
        tyx = mu*(dry(udf)+drx(vdf))
        tzx = mu*(drz(udf)+drx(wdf))
        txy = mu*(drx(vdf)+dry(udf))
        tyy = mu*(dry(vdf)+dry(vdf))
        tzy = mu*(drz(vdf)+dry(wdf))
        txz = mu*(drx(wdf)+drz(udf))
        tyz = mu*(dry(wdf)+drz(vdf))
        tzz = mu*(drz(wdf)+drz(wdf))


        !txx = mu * t(1,1)
        !tyx = mu * t(2,1)
        !tzx = mu * t(3,1)
        !txy = mu * t(1,2)
        !tyy = mu * t(2,2)
        !tzy = mu * t(3,2)
        !txz = mu * t(1,3)
        !tyz = mu * t(2,3)
        !tzz = mu * t(3,3)
        !pxx = jac(iq,ie)*(finv11 *txx+finv12 *tyx+finv13 *tzx)
        !pyx = jac(iq,ie)*(finv21 *txx+finv22 *tyx+finv23 *tzx)
        !pzx = jac(iq,ie)*(finv31 *txx+finv32 *tyx+finv33 *tzx)
        !pxy = jac(iq,ie)*(finv11 *txy+finv12 *tyy+finv13 *tzy)
        !pyy = jac(iq,ie)*(finv21 *txy+finv22 *tyy+finv23 *tzy)
        !pzy = jac(iq,ie)*(finv31 *txy+finv32 *tyy+finv33 *tzy)
        !pxz = jac(iq,ie)*(finv11 *txz+finv12 *tyz+finv13 *tzz)
        !pyz = jac(iq,ie)*(finv21 *txz+finv22 *tyz+finv23 *tzz)
        !pzz = jac(iq,ie)*(finv31 *txz+finv32 *tyz+finv33 *tzz)


        prs_t(udf) = res_t(udf) * taum
        prs_t(vdf) = res_t(vdf) * taum
        prs_t(wdf) = res_t(wdf) * taum
        prs_c      = res_c*tauc

!.... calculate the residual at each degree of freedom
        do inl=1,nen
           node=ien(inl,ie)
           temp=ro*(wx*ph(xsd,inl)+wy*ph(ysd,inl)+wz*ph(zsd,inl))
!.....      Continuity Equation
           p(pdf,node) = p(pdf,node)-ph(0,inl)*res_c

!.....      Momentum Equation (Euler Residual)
           p(udf,node) = p(udf,node)-ph(0,inl)*res_a(udf)
           p(vdf,node) = p(vdf,node)-ph(0,inl)*res_a(vdf)
           p(wdf,node) = p(wdf,node)-ph(0,inl)*res_a(wdf)

!.....      Momentum Equation (C: Diffusion terms)
!.... Piola Kirchoff Stress
           !p(udf,node) = p(udf,node) +                   &
           ! ph(xsd,inl)* jac(iq,ie)*finv11 *pp -     &
           ! ph(xsd,inl)*pxx -                        &
           ! ph(ysd,inl)*pyx -                        &
           ! ph(zsd,inl)*pzx 
           !p(vdf,node) = p(vdf,node) +                   &
           ! ph(ysd,inl)* jac(iq,ie)*finv22 *pp -     &
           ! ph(xsd,inl)*pxy -                        &
           ! ph(ysd,inl)*pyy -                        &
           ! ph(zsd,inl)*pzy 
           !p(wdf,node) = p(wdf,node) +                   &
           ! ph(zsd,inl)* jac(iq,ie)*finv33 *pp -     &
           ! ph(xsd,inl)*pxz -                        &
           ! ph(ysd,inl)*pyz -                        &
           ! ph(zsd,inl)*pzz  
           !...Cauchy
           p(udf,node) = p(udf,node) + ph(xsd,inl)* pp -    &
                                       ph(xsd,inl)*txx -    &
                                       ph(ysd,inl)*tyx -    &
                                       ph(zsd,inl)*tzx 
           p(vdf,node) = p(vdf,node) + ph(ysd,inl)* pp -    &
                                       ph(xsd,inl)*txy -    &
                                       ph(ysd,inl)*tyy -    &
                                       ph(zsd,inl)*tzy 
           p(wdf,node) = p(wdf,node) + ph(zsd,inl)* pp -    &
                                       ph(xsd,inl)*txz -    &
                                       ph(ysd,inl)*tyz -    &
                                       ph(zsd,inl)*tzz 




!.....     Stablization with Tau_moment
           p(pdf,node) = p(pdf,node)           &
                - ph(xsd,inl)*prs_t(udf)       &
                - ph(ysd,inl)*prs_t(vdf)       &
                - ph(zsd,inl)*prs_t(wdf)
           p(udf,node) = p(udf,node) -prs_t(udf)*temp
           p(vdf,node) = p(vdf,node) -prs_t(vdf)*temp
           p(wdf,node) = p(wdf,node) -prs_t(wdf)*temp

!.....      Stablization with Tau_cont    
           p(udf,node) = p(udf,node) - ph(xsd,inl)*prs_c
           p(vdf,node) = p(vdf,node) - ph(ysd,inl)*prs_c
           p(wdf,node) = p(wdf,node) - ph(zsd,inl)*prs_c
        enddo
!   if (node.eq.10) write(*,*) p(1,node),p(2,node),p(3,node),p(4,node)  

!.....   Diagonal Preconditioner

        effd = mu*eft0*alpha
        effm = taum*eft0
        effc = tauc*eft0

        do inl=1,nen

           node=ien(inl,ie)
           ug = ro*(u*sh(1,inl)+v*sh(2,inl)+w*sh(3,inl))
           temp = alpha*ug + sh(0,inl)*dtinv*ro
         
           q(udf,node) = q(udf,node)+(alpha*ro*sh(0,inl)*drx(1)+temp)*(sh(0,inl)*eft0+ug*effm)
           q(vdf,node) = q(vdf,node)+(alpha*ro*sh(0,inl)*dry(2)+temp)*(sh(0,inl)*eft0+ug*effm) 
           q(wdf,node) = q(wdf,node)+(alpha*ro*sh(0,inl)*drz(3)+temp)*(sh(0,inl)*eft0+ug*effm)
           temp = sh(1,inl)**2+sh(2,inl)**2+sh(3,inl)**2
         
           q(udf,node) = q(udf,node)+(sh(1,inl)**2+temp)*effd
           q(vdf,node) = q(vdf,node)+(sh(2,inl)**2+temp)*effd
           q(wdf,node) = q(wdf,node)+(sh(3,inl)**2+temp)*effd
           q(pdf,node) = q(pdf,node)+temp*effm
           q(udf,node) = q(udf,node)+sh(1,inl)**2*effc
           q(vdf,node) = q(vdf,node)+sh(2,inl)**2*effc
           q(wdf,node) = q(wdf,node)+sh(3,inl)**2*effc
        enddo
     enddo quad_loop
  enddo e_loop
  return
end subroutine block_ale
    

end module gmres_fluid_ale