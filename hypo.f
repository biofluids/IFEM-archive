      subroutine hypo

      include "global.h"
      include "malloc.h"

      integer ien(nen,nec),rng(neface,nec)
      real* 8 xn(nsd,nnc), x(nsd,nn_loc), xna(nsd,nn),hg(nec)
      pointer (rngptr,rng),(ienptr,ien)
      pointer (xnptr,xn),(xptr,x),(hgptr,hg),(xnaptr,xna)

      integer id(ndf,nnc) 
      real* 8 dn(ndf,nnc),fn(ndf,nnc)
      real* 8 d(ndf,nn_loc),do(ndf,nn_loc)
      real* 8 bg(ndf,nnc),dg(ndf,nnc),p(ndf,nn_loc)
      real* 8 wg(ndf,nnc), w(ndf,nn_loc)
      pointer (idptr,id),(fnptr,fn),(dnptr,dn)
      pointer (dptr,d),(doptr,do),(uptr,u)
      pointer (bgptr,bg),(dgptr,dg),(pptr,p)
      pointer (wgptr,wg),(wptr,w)

      real* 8 disp(nsd,nn_loc),dispn(nsd,nnc), kdg(nsd,nnc)
      real* 8 kp(nsd,nn_loc), kw(nsd,nn_loc), kid(nsd,nnc)
      real* 8 kbg(nsd,nnc), kzg(nsd*nnc), kwg(nsd,nnc)
      real* 8 e(nsd,nn_loc), err(nsd,nnc),xrefn(nsd,nnc)
      real* 8 xref(nsd,nn_loc),xrefold(nsd,nn_loc)
      real* 8 xnold(nsd,nnc),meshvel(nsd,nn_loc)
      real* 8 meshveln(nsd,nnc)
      real* 8 refvel(nsd,nquad,nec)
      real* 8 refvelo(nsd,nquad,nec)

      pointer (dispptr,disp),(dispnptr,dispn),(kdgptr,kdg)
      pointer (kpptr,kp), (kwptr,kw), (kidptr,kid)
      pointer (kbgptr,kbg), (kzgptr,kzg),(kwgptr,kwg)
      pointer (eptr,e),(errptr,err),(xrefnptr,xrefn)
      pointer (xrefptr,xref),(xnoldptr,xnold)
      pointer (xrefoldptr,xrefold)
      pointer (meshvelnptr,meshveln), (meshvelptr,meshvel)
      pointer (refvelptr,refvel)
      pointer (refveloptr,refvelo)

      real* 8 f(nsd,nsd,nquad,nec), finv(nsd,nsd,nquad,nec)
      real* 8 jac(nquad,nec),jacinv(nquad,nec),jaco(nquad,nec)
      pointer (fptr,f),(finvptr,finv)
      pointer (jacptr,jac),(jacinvptr,jacinv),(jacoptr,jaco)

      integer idf(nnc) 
      real* 8 dnf(nnc),fnf(nnc)
      real* 8 df(nn_loc),dof(nn_loc)
      real* 8 bgf(nnc),dgf(nnc),pf(nn_loc)
      real* 8 wgf(nnc), wf(nn_loc)
      pointer (idfptr,idf),(fnfptr,fnf),(dnfptr,dnf),(anfptr,anf)
      pointer (dfptr,df),(dofptr,dof),(afptr,af)
      pointer (bgfptr,bgf),(dgfptr,dgf),(pfptr,pf)
      pointer (wgfptr,wgf),(wfptr,wf)

      real* 8 dd(ndf+4,nnc),hn(nnc),hm(nn_loc)
      pointer (ddptr,dd),(hnptr,hn),(hmptr,hm)

      real* 8 z(ndf*nnc,inner),kz(ndf*nnc,kinner)
      real* 8 v(ndf*nnc,inner+1),kv(ndf*nnc,kinner+1)
      real* 8 zg(ndf*nnc), avg(ndf*nnc), sm(ndf*nnc)
      real* 8 vloc(ndf,nn_loc),avloc(ndf,nn_loc)
      pointer (zptr,z),(vptr,v),(kzptr,kz),(kvptr,kv)
      pointer (avgptr,avg),(zgptr,zg),(smptr,sm)
      pointer (vlocptr,vloc),(avlocptr,avloc)

      real* 8 zf(nnc,inner)
      real* 8 vf(nnc,inner+1)
      real* 8 zgf(nnc), avgf(nnc), smf(nnc)
      real* 8 vlocf(nn_loc),avlocf(nn_loc)
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
      ienptr = malloc(nen*nec*isize)
      call readien(ien)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("communications",-999,.false.)
      call commsetup(ien)
c  call nodesetup (cnn,mn) !!!!change me

      xnptr = malloc(nnc*nsd*fsize)
      xptr = malloc(nn_loc*nsd*fsize)
      rngptr = malloc(nec*neface*isize)
      xnaptr= malloc(nn*nsd*fsize)
      hgptr = malloc(nec*fsize)

      idptr = malloc(nnc*ndf*isize)
      dnptr = malloc(nnc*ndf*fsize)
      fnptr = malloc(nnc*ndf*fsize)
      bgptr = malloc(nnc*ndf*fsize)
      dgptr = malloc(nnc*ndf*fsize)
      wgptr = malloc(nnc*ndf*fsize)
      dptr = malloc(nn_loc*ndf*fsize)
      uptr = malloc(nn_loc*nsd*fsize)
      doptr = malloc(nn_loc*ndf*fsize)
      pptr = malloc(nn_loc*ndf*fsize)
      wptr = malloc(nn_loc*ndf*fsize)

      idfptr = malloc(nnc*isize)
      dnfptr = malloc(nnc*fsize)
      fnfptr = malloc(nnc*fsize)
      anfptr = malloc(nnc*fsize)
      bgfptr = malloc(nnc*fsize)
      dgfptr = malloc(nnc*fsize)
      wgfptr = malloc(nnc*fsize)
      dfptr = malloc(nn_loc*fsize)
      afptr = malloc(nn_loc*fsize)
      dofptr = malloc(nn_loc*fsize)
      pfptr = malloc(nn_loc*fsize)
      wfptr = malloc(nn_loc*fsize)
      
      dispptr = malloc(nn_loc*nsd*fsize)
      dispnptr = malloc(nnc*nsd*fsize)
      kdgptr = malloc(nnc*nsd*fsize)
      kpptr = malloc(nn_loc*nsd*fsize)
      kwptr = malloc(nn_loc*nsd*fsize)
      kidptr = malloc(nnc*nsd*isize)
      kbgptr = malloc(nnc*nsd*fsize)
      kwgptr = malloc(nnc*nsd*fsize)
      kzgptr = malloc(nnc*nsd*fsize)
      errptr = malloc(nnc*nsd*fsize)
      eptr = malloc(nn_loc*nsd*fsize)
      xrefnptr = malloc(nnc*nsd*fsize)
      xrefptr = malloc(nn_loc*nsd*fsize)
      xnoldptr=malloc(nnc*nsd*fsize)
      xrefoldptr = malloc(nn_loc*nsd*fsize)
      meshvelptr = malloc(nn_loc*nsd*fsize)
      meshvelnptr = malloc(nnc*nsd*fsize)
      refvelptr = malloc(nquad*nec*nsd*fsize)
      refveloptr = malloc(nquad*nec*nsd*fsize)
      fptr = malloc(nsd*nsd*nquad*nec*fsize)
      finvptr = malloc(nsd*nsd*nquad*nec*fsize)

      jacptr = malloc(nquad*nec*fsize)
      jacinvptr = malloc(nquad*nec*fsize)
      jacoptr = malloc(nquad*nec*fsize)

      ddptr = malloc((ndf+4)*nnc*fsize)
      hnptr = malloc(nnc*fsize)
      hmptr = malloc(nn_loc*fsize)

      hptr = malloc((inner+1)*inner*fsize)
      yptr = malloc((inner+1)*fsize)
      ccptr = malloc(inner*fsize)
      ssptr = malloc(inner*fsize)

      khptr = malloc((kinner+1)*kinner*fsize)
      kyptr = malloc((kinner+1)*fsize)
      kccptr = malloc(kinner*fsize)
      kssptr = malloc(kinner*fsize)
      kzptr = malloc(ndf*nnc*kinner*fsize)
      kvptr = malloc(ndf*nnc*(kinner+1)*fsize)

      zptr = malloc(ndf*nnc*inner*fsize)
      vptr = malloc(ndf*nnc*(inner+1)*fsize)
      zgptr = malloc(ndf*nnc*fsize)
      avgptr = malloc(ndf*nnc*fsize)
      smptr = malloc(ndf*nnc*fsize)
      vlocptr = malloc(ndf*nn_loc*fsize)
      avlocptr = malloc(ndf*nn_loc*fsize)

      zfptr = malloc(nnc*inner*fsize)
      vfptr = malloc(nnc*(inner+1)*fsize)
      zgfptr = malloc(nnc*fsize)
      avgfptr = malloc(nnc*fsize)
      smfptr = malloc(nnc*fsize)
      vlocfptr = malloc(nn_loc*fsize)
      avlocfptr = malloc(nn_loc*fsize)

      diameter=1.5
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("read x and rng",-999,.false.)
      call readx(xrefn)
      call readrng(rng)
c      if (calcforce) then
c        call getnqdf(rng)
c      else
c        nqdf = 1
c      end if
c      if (nqdf.eq.0) nqdf=1
      call error("diskout",-999,.false.)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call readallx(xna)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call fclear(dispn,nsd*nnc)
      call error("form d and id",-999,.false.)
      call formid(id,idf,rng,ien,hn,hm)
      call formd (xrefn,fn,fnf,rng,ien,hn,hm)
      call equal(fn ,dn ,ndf*nnc)
      call equal(fnf,dnf,nnc)
      if(restart) then
         call error("restart",-999,.false.)
         call diskin(xrefn,dn,dnf,dd)
      endif
      call setd(dn ,fn ,id ,ndf)
      call setd(dnf,fnf,idf,1)
      call equal(dnf,anf,nnc)
      call sharp(anf,nnc,1)
      call error("diskout",-999,.false.)
      call diskout(dispn,dn,anf,dd)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if ((calcforce).and.(tt.eq.0.0)) then
c        fdrag(:) = 0
c        call fdragout(tt,fdrag,1,2)
c      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("local shape functions",-999,.false.)
      call shape
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("localization",-999,.false.)
      call gather (xref,xrefn,nsd,hn,hm)
      call gather (d,dn,ndf,hn,hm)
      call gather (df,dnf,1,hn,hm)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("mesh information",-999,.false.)
      call vol(xref,ien,df)
      if (myid.eq.0) then
         write(7,*)' '
         write(7,'(" Volume of gas............ = ",e15.8)') gas
         write(7,'(" Volume of liquid......... = ",e15.8)') liq
      endif
      call lenght(xref,ien,hg)
      if (myid.eq.0) then
         write(7,'(" Minimum element lenght... = ",e15.8)') hmin
         write(7,'(" Maximum element lenght... = ",e15.8)') hmax
         write(7,'(" Minimum element volume... = ",e15.8)') vmin
         write(7,'(" Maximum element volume... = ",e15.8)') vmax
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (myid.eq.0) write (7,101) nq,nqf
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call equal(xrefn,xn,nsd*nnc)
      call fclear(refvel,nsd*nquad*nec)  !not sure about this

      do ie=1,nec
         do iq=1,nquad
            jac(iq,ie)=1.0
         enddo
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Starts Time Loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do its=1,nts
         if (myid.eq.0) write (6,*) ' '
         if (myid.eq.0) write (6,*) ' TIME STEP = ', its
         if (myid.eq.0) write (7,*) ' '
         if (myid.eq.0) write (7,*) ' TIME STEP = ', its
         tt = tt + dt
         call equal(d,do,ndf*nn_loc)
         call equal(df,dof,nn_loc)
         call equal(xn,xnold,nsd*nnc)     
         call equal(refvel,refvelo,nsd*nquad*nec)
         call equal(jac,jaco,nquad*nec)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Starts Mesh Update Iteration Loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         call error("form dd and kid",-999,.false.)
         call formidm(kid,rng,ien,hn,hm)
         call formdm (xrefn,dispn,rng,ien,hn,hm,diameter)   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         call error("localization",-999,.false.)
         call gather (xref,xrefn,nsd,hn,hm)
         call gather (disp,dispn,nsd,hn,hm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do iit=1,40            !change if necessary
            call fclear(kp,nsd*nn_loc)
            call fclear(kw,nsd*nn_loc)
            call blockm(xref,e,disp,kw,kp,ien)
            assemble=.true.
            call scatter(kp,kbg,nsd,assemble,hn,hm)
            call scatter(kw,kwg,nsd,assemble,hn,hm)
            call setid(kbg,kid,nsd)
            call getnorm(kbg,kbg,nsd*nnc,res_l)
            res_l= sqrt(res_l/nq)
            if (res_l.lt.epsilon) goto 200
            call fclear(kdg,nsd*nnc)
            call gmresm(xref,kid,kwg,kbg,kdg,ien,hn,hm,kz,kv,kzg,
     &           avg,sm, avloc,kh,ky,kcc,kss)
            call getnorm(kdg,kdg,nsd*nnc,del_l)
            del_l = sqrt(del_l/nq)
            call update(kp, disp, kdg, nsd,hn,hm)
            if (myid.eq.0) write(6,102)iit,del_l,del_g
            if (myid.eq.0) write(7,102)iit,del_l,del_g
         enddo
 200     continue
         assemble=.false.
         call scatter(disp,dispn,nsd,assemble,hn,hm)
         call add(xn,xrefn,dispn,nsd)
         call gather(x,xn,nsd,hn,hm)
c.... calculate d(x)/d(xref), deformation gradient
         call defgrad(x, xref,ien,f,jac,finv,jacinv)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
         do iit=1,nit
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            call equal(dof,df,nn_loc)
c            call fclear(pf,nn_loc)
c            call fclear(wf,nn_loc)
c            call getu(d,do,u)
c            call blockf(xref,df,dof,u,pf,wf,hg,ien)
c            assemble=.true.
c            call scatter(pf,bgf,1,assemble,hn,hm)
c            call scatter(wf,wgf,1,assemble,hn,hm)
c            call setid(bgf,idf,1)
c            call getnorm(bgf,bgf,nnc,res_g)
c            res_g = sqrt(res_g/nqf)

c            call fclear(dgf,nnc)
c            call gmresf(xref,u,idf,wgf,bgf,dgf,hg,ien,hn,hm,
c     &           zf,vf,zgf,avgf,smf,vlocf,avlocf,h,y,cc,ss)
c            call getnorm(dgf,dgf,nnc,del_g)
c            del_g = sqrt(del_g/nqf)
c            call update(pf, df, dgf, 1,hn,hm)
c            call cut(df)
c            call locate(xref,ien,df)
c            call sharp(df,nn_loc,0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            call fclear(p,ndf*nn_loc)
            call fclear(w,ndf*nn_loc)
            call getfi(df,dof,af)
            call sharp(af,nn_loc,1)
            if(hydro.gt.0) call hydrostatic(xref,af,p,rng,ien)
c.... Calculate mesh & referential velocities, vhat & w
            call velocity(xn,xnold,x,meshveln,meshvel,refvel,
     +           finv,d,ien,hn,hm)
            call block(xref,d,do,af,p,w,hg,ien,f,finv,jac,jaco,
     +           refvel,refvelo)
            assemble=.true.
            call scatter(p,bg,ndf,assemble,hn,hm)
            call scatter(w,wg,ndf,assemble,hn,hm)
            call setid(bg,id,ndf)
            call getnorm(bg,bg,ndf*nnc,res_l)
            res_l= sqrt(res_l/nq)
            call fclear(dg,ndf*nnc)
            call gmres(xref,d,do,id,af,wg,bg,dg,hg,ien,hn,hm,z,v,zg,avg,
     &           sm,vloc,avloc,h,y,cc,ss,finv,jac,jaco,refvel,refvelo)
            call getnorm(dg,dg,ndf*nnc,del_l)
            del_l = sqrt(del_l/nq)
            call update(p, d, dg, ndf,hn,hm)

            if (myid.eq.0) write(6,102)iit,res_l,del_l,res_g,del_g
            if (myid.eq.0) write(7,102)iit,res_l,del_l,res_g,del_g
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         enddo
c        if (calcforce) then
c          call getfdrag(fdrag,dn,d,x,ien,rng,hn,hm,finv)
c          call fdragout(tt,fdrag,1,2)
c        end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..... Writing Output
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(mod(its,ntsbout).eq.0) then
            if(surf(0).gt.0) call dyno(xref,d,rng,ien)
            idisk = idisk+1
            assemble=.false.
            call scatter(d,dn,ndf,assemble,hn,hm)
            call scatter(df,dnf,1,assemble,hn,hm)
            call equal(dnf,anf,nnc)
            call sharp(anf,nnc,1)
            call diskout(dispn,dn,anf,dd)
         endif
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (myid.eq.0) then
         write(7,*)'Processors...............',numproc 
      endif
      
      return
      
 101  format(/"Number of equations for Flow.........(nq) = ",i10,
     &     /"Number of equations for Interface...(nqf) = ",i10)
 102  format("Iteration",i3,':  ',4e14.7)
      
      end
      
      
