      subroutine hypo

      include "global.h"
      include "malloc.h"

      integer ien(nen,ne),rng(neface,ne)
      real* 8 xn(nsd,nn), x(nsd,nn), xna(nsd,nn),hg(ne)
      pointer (rngptr,rng),(ienptr,ien)
      pointer (xnptr,xn),(xptr,x),(hgptr,hg),(xnaptr,xna)

      integer id(ndf,nn) 
      real* 8 dn(ndf,nn),fn(ndf,nn)
      real* 8 d(ndf,nn),do(ndf,nn)
      real* 8 bg(ndf,nn),dg(ndf,nn),p(ndf,nn)
      real* 8 wg(ndf,nn), w(ndf,nn)
      pointer (idptr,id),(fnptr,fn),(dnptr,dn)
      pointer (dptr,d),(doptr,do),(uptr,u)
      pointer (bgptr,bg),(dgptr,dg),(pptr,p)
      pointer (wgptr,wg),(wptr,w)

      real* 8 disp(nsd,nn),dispn(nsd,nn), kdg(nsd,nn)
      real* 8 kp(nsd,nn), kw(nsd,nn), kid(nsd,nn)
      real* 8 kbg(nsd,nn), kzg(nsd*nn), kwg(nsd,nn)
      real* 8 e(nsd,nn), err(nsd,nn),xrefn(nsd,nn)
      real* 8 xref(nsd,nn),xrefold(nsd,nn)
      real* 8 xnold(nsd,nn),meshvel(nsd,nn)
      real* 8 meshveln(nsd,nn)
      real* 8 refvel(nsd,nquad,ne)
      real* 8 refvelo(nsd,nquad,ne)

      pointer (dispptr,disp),(dispnptr,dispn),(kdgptr,kdg)
      pointer (kpptr,kp), (kwptr,kw), (kidptr,kid)
      pointer (kbgptr,kbg), (kzgptr,kzg),(kwgptr,kwg)
      pointer (eptr,e),(errptr,err),(xrefnptr,xrefn)
      pointer (xrefptr,xref),(xnoldptr,xnold)
      pointer (xrefoldptr,xrefold)
      pointer (meshvelnptr,meshveln), (meshvelptr,meshvel)
      pointer (refvelptr,refvel)
      pointer (refveloptr,refvelo)

      real* 8 f(nsd,nsd,nquad,ne), finv(nsd,nsd,nquad,ne)
      real* 8 jac(nquad,ne),jacinv(nquad,ne),jaco(nquad,ne)
      pointer (fptr,f),(finvptr,finv)
      pointer (jacptr,jac),(jacinvptr,jacinv),(jacoptr,jaco)

      integer idf(nn) 
      real* 8 dnf(nn),fnf(nn)
      real* 8 df(nn),dof(nn)
      real* 8 bgf(nn),dgf(nn),pf(nn)
      real* 8 wgf(nn), wf(nn)
      pointer (idfptr,idf),(fnfptr,fnf),(dnfptr,dnf),(anfptr,anf)
      pointer (dfptr,df),(dofptr,dof),(afptr,af)
      pointer (bgfptr,bgf),(dgfptr,dgf),(pfptr,pf)
      pointer (wgfptr,wgf),(wfptr,wf)

      real* 8 dd(ndf+4,nn),hn(nn),hm(nn)
      pointer (ddptr,dd),(hnptr,hn),(hmptr,hm)

      real* 8 z(ndf*nn,inner),kz(ndf*nn,kinner)
      real* 8 v(ndf*nn,inner+1),kv(ndf*nn,kinner+1)
      real* 8 zg(ndf*nn), avg(ndf*nn), sm(ndf*nn)
      real* 8 vloc(ndf,nn),avloc(ndf,nn)
      pointer (zptr,z),(vptr,v),(kzptr,kz),(kvptr,kv)
      pointer (avgptr,avg),(zgptr,zg),(smptr,sm)
      pointer (vlocptr,vloc),(avlocptr,avloc)

      real* 8 zf(nn,inner)
      real* 8 vf(nn,inner+1)
      real* 8 zgf(nn), avgf(nn), smf(nn)
      real* 8 vlocf(nn),avlocf(nn)
      pointer (zfptr,zf),(vfptr,vf)
      pointer (avgfptr,avgf),(zgfptr,zgf),(smfptr,smf)
      pointer (vlocfptr,vlocf),(avlocfptr,avlocf)

      real* 8 h(inner+1,inner),kh(kinner,kinner)
      real* 8 y(inner+1),ky(kinner)
      real* 8 cc(inner), ss(inner),kcc(kinner),kss(kinner)
      pointer (hptr,h),(yptr,y),(ccptr,cc),(ssptr,ss)
      pointer (khptr,kh),(kyptr,ky),(kccptr,kcc),(kssptr,kss)

      logical assemble
      integer ierr,nqdf
      real* 8 diameter
      real* 8 fdrag(nsdpad)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("initialization",-999,.false.)
      tt = t_start
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("read ien",-999,.false.)
      ienptr = malloc(nen*ne*isize)
      call readien(ien)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  call nodesetup (cnn,mn) !!!!change me

      xnptr = malloc(nn*nsd*fsize)
      xptr = malloc(nn*nsd*fsize)
      rngptr = malloc(ne*neface*isize)
      xnaptr= malloc(nn*nsd*fsize)
      hgptr = malloc(ne*fsize)

      idptr = malloc(nn*ndf*isize)
      dnptr = malloc(nn*ndf*fsize)
      fnptr = malloc(nn*ndf*fsize)
      bgptr = malloc(nn*ndf*fsize)
      dgptr = malloc(nn*ndf*fsize)
      wgptr = malloc(nn*ndf*fsize)
      dptr = malloc(nn*ndf*fsize)
      uptr = malloc(nn*nsd*fsize)
      doptr = malloc(nn*ndf*fsize)
      pptr = malloc(nn*ndf*fsize)
      wptr = malloc(nn*ndf*fsize)

      idfptr = malloc(nn*isize)
      dnfptr = malloc(nn*fsize)
      fnfptr = malloc(nn*fsize)
      anfptr = malloc(nn*fsize)
      bgfptr = malloc(nn*fsize)
      dgfptr = malloc(nn*fsize)
      wgfptr = malloc(nn*fsize)
      dfptr = malloc(nn*fsize)
      afptr = malloc(nn*fsize)
      dofptr = malloc(nn*fsize)
      pfptr = malloc(nn*fsize)
      wfptr = malloc(nn*fsize)
      
      dispptr = malloc(nn*nsd*fsize)
      dispnptr = malloc(nn*nsd*fsize)
      kdgptr = malloc(nn*nsd*fsize)
      kpptr = malloc(nn*nsd*fsize)
      kwptr = malloc(nn*nsd*fsize)
      kidptr = malloc(nn*nsd*isize)
      kbgptr = malloc(nn*nsd*fsize)
      kwgptr = malloc(nn*nsd*fsize)
      kzgptr = malloc(nn*nsd*fsize)
      errptr = malloc(nn*nsd*fsize)
      eptr = malloc(nn*nsd*fsize)
      xrefnptr = malloc(nn*nsd*fsize)
      xrefptr = malloc(nn*nsd*fsize)
      xnoldptr=malloc(nn*nsd*fsize)
      xrefoldptr = malloc(nn*nsd*fsize)
      meshvelptr = malloc(nn*nsd*fsize)
      meshvelnptr = malloc(nn*nsd*fsize)
      refvelptr = malloc(nquad*ne*nsd*fsize)
      refveloptr = malloc(nquad*ne*nsd*fsize)
      fptr = malloc(nsd*nsd*nquad*ne*fsize)
      finvptr = malloc(nsd*nsd*nquad*ne*fsize)

      jacptr = malloc(nquad*ne*fsize)
      jacinvptr = malloc(nquad*ne*fsize)
      jacoptr = malloc(nquad*ne*fsize)

      ddptr = malloc((ndf+4)*nn*fsize)
      hnptr = malloc(nn*fsize)
      hmptr = malloc(nn*fsize)

      hptr = malloc((inner+1)*inner*fsize)
      yptr = malloc((inner+1)*fsize)
      ccptr = malloc(inner*fsize)
      ssptr = malloc(inner*fsize)

      khptr = malloc((kinner+1)*kinner*fsize)
      kyptr = malloc((kinner+1)*fsize)
      kccptr = malloc(kinner*fsize)
      kssptr = malloc(kinner*fsize)
      kzptr = malloc(ndf*nn*kinner*fsize)
      kvptr = malloc(ndf*nn*(kinner+1)*fsize)

      zptr = malloc(ndf*nn*inner*fsize)
      vptr = malloc(ndf*nn*(inner+1)*fsize)
      zgptr = malloc(ndf*nn*fsize)
      avgptr = malloc(ndf*nn*fsize)
      smptr = malloc(ndf*nn*fsize)
      vlocptr = malloc(ndf*nn*fsize)
      avlocptr = malloc(ndf*nn*fsize)

      zfptr = malloc(nn*inner*fsize)
      vfptr = malloc(nn*(inner+1)*fsize)
      zgfptr = malloc(nn*fsize)
      avgfptr = malloc(nn*fsize)
      smfptr = malloc(nn*fsize)
      vlocfptr = malloc(nn*fsize)
      avlocfptr = malloc(nn*fsize)

      diameter=1.5
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("read x and rng",-999,.false.)
      call readx(xref)
      call readrng(rng)
c      if (calcforce) then
c        call getnqdf(rng)
c      else
c        nqdf = 1
c      end if
c      if (nqdf.eq.0) nqdf=1
      call error("diskout",-999,.false.)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call fclear(disp,nsd*nn)
      call error("form d and id",-999,.false.)
      call formid(id,idf,rng,ien)
      call formd (xref,fn,fnf,rng,ien)
      call equal(fn ,d ,ndf*nn)
      call equal(fnf,df,nn)
      if(restart) then
         call error("restart",-999,.false.)
         call diskin(xref,d,dnf,dd)
      endif
      call setd(d ,fn ,id ,ndf)
      call setd(df,fnf,idf,1)
      call equal(df,anf,nn)
      call sharp(anf,nn,1)
      call error("diskout",-999,.false.)
      call diskout(disp,d,anf,dd)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if ((calcforce).and.(tt.eq.0.0)) then
c        fdrag(:) = 0
c        call fdragout(tt,fdrag,1,2)
c      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("local shape functions",-999,.false.)
      call shape

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("mesh information",-999,.false.)
      call vol(xref,ien,df)
         write(7,*)' '
         write(7,'(" Volume of gas............ = ",e15.8)') gas
         write(7,'(" Volume of liquid......... = ",e15.8)') liq
      call lenght(xref,ien,hg)
         write(7,'(" Minimum element lenght... = ",e15.8)') hmin
         write(7,'(" Maximum element lenght... = ",e15.8)') hmax
         write(7,'(" Minimum element volume... = ",e15.8)') vmin
         write(7,'(" Maximum element volume... = ",e15.8)') vmax

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (7,101) nq,nqf
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call equal(xref,xn,nsd*nn)
      call fclear(refvel,nsd*nquad*ne)  !not sure about this

      do ie=1,ne
         do iq=1,nquad
            jac(iq,ie)=1.0
         enddo
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Starts Time Loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do its=1,nts
         write (6,*) ' '
         write (6,*) ' TIME STEP = ', its
         write (7,*) ' '
         write (7,*) ' TIME STEP = ', its
         tt = tt + dt
         call equal(d,do,ndf*nn)
         call equal(df,dof,nn)
         call equal(xn,xnold,nsd*nn)     
         call equal(refvel,refvelo,nsd*nquad*ne)
         call equal(jac,jaco,nquad*ne)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Starts Mesh Update Iteration Loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         call error("form dd and kid",-999,.false.)
         call formidm(kid,rng,ien)
         call formdm (xref,disp,rng,ien,diameter)   
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do iit=1,40            !change if neessary
            call fclear(kp,nsd*nn)
            call fclear(kw,nsd*nn)
            call blockm(xref,e,disp,kw,kp,ien)

            call setid(kp,kid,nsd)
            call getnorm(kp,kp,nsd*nn,res_l)
            res_l= sqrt(res_l/nq)
            if (res_l.lt.epsilon) goto 200
            call fclear(kdg,nsd*nn)
            call gmresm(xref,kid,kw,kp,kdg,ien,hn,hm,kz,kv,kzg,
     &           avg,sm, avloc,kh,ky,kcc,kss)
            call getnorm(kdg,kdg,nsd*nn,del_l)
            del_l = sqrt(del_l/nq)
            call update(kp, disp, kdg, nsd,hn,hm)
            write(6,102)iit,del_l,del_g
            write(7,102)iit,del_l,del_g
         enddo
 200     continue

         call add(xn,xref,disp,nsd)
c.... calculate d(x)/d(xref), deformation gradient
         call defgrad(xn, xref,ien,f,jac,finv,jacinv)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
         do iit=1,nit
            call fclear(p,ndf*nn)
            call fclear(w,ndf*nn)
            call getfi(df,dof,af)
            call sharp(af,nn,1)
            if(hydro.gt.0) call hydrostatic(xref,af,p,rng,ien)
c.... Calculate mesh & referential velocities, vhat & w
            call velocity(xn,xnold,meshvel,refvel,
     +           finv,d,ien,hn,hm)
            call block(xref,d,do,af,p,w,hg,ien,f,finv,jac,jaco,
     +           refvel,refvelo)

            call setid(p,id,ndf)
            call getnorm(p,p,ndf*nn,res_l)
            res_l= sqrt(res_l/nq)
            call fclear(dg,ndf*nn)
            call gmres(xref,d,do,id,af,w,p,dg,hg,ien,hn,hm,z,v,zg,avg,
     &           sm,vloc,avloc,h,y,cc,ss,finv,jac,jaco,refvel,refvelo)
            call getnorm(dg,dg,ndf*nn,del_l)
            del_l = sqrt(del_l/nq)
            call update(p, d, dg, ndf,hn,hm)

            write(6,102)iit,res_l,del_l,res_g,del_g
            write(7,102)iit,res_l,del_l,res_g,del_g
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         enddo
c        if (calcforce) then
c          call getfdrag(fdrag,dn,d,x,ien,rng,hn,hm,finv)
c          call fdragout(tt,fdrag,1,2)
c        end if

      end do

      
      return
      
 101  format(/"Number of equations for Flow.........(nq) = ",i10,
     &     /"Number of equations for Interface...(nqf) = ",i10)
 102  format("Iteration",i3,':  ',4e14.7)

      
      end
      
      
