c-------------------------------------*
c                                     *
c  ibm/fem version 1 poly/mech        *
c                                     *
c-------------------------------------*

      program ib1_main_iterate   
   



c==============================================================
c ***   begin 1 step   ***
c==============================================================
      do 1000 istep = nt0+1, nt0+n_step_run

         iti=istep

         write(6,*) istep
         
         klok = klok + 1
         
         if (n_iter_step .gt. 1) then
            if (dilocaldomain .eq. i_domain_soloist) then
               write(6,*) 'iteration # ',istep,
     $              ' -------------------------------'
               call flush(6)
            endif
         endif

c++++++++
c     fem
c++++++++ 
      if (n_ibmfem .eq. 1) then
         if (nrestart .ne. 1) then
            if (initdir .eq. 1) then
               xtime=tfun(ntfun)
               do 702 k=1,numgb
                  if (ndirgb(k) .eq. 121111) then
                     do 703 m=1,numdir(k) 
                        dis(1,nodegb(k,m))=
     $                       xtime*xindis(1,nodegb(k,m))              
 703                 continue
                  endif
c     
                  if (ndirgb(k) .eq. 112111) then
                     do 744 m=1,numdir(k) 
                        dis(2,nodegb(k,m))=
     $                       xtime*xindis(2,nodegb(k,m))              
 744                 continue
                  endif
 702           continue
            endif
         else
            if (iti .eq. nt0+1) then
               do 5545 i=1,nnd
                  read(25) dis(1,i),dis(2,i),
     $                 du(1,i),du(2,i),
     $                 acm(1,i),acm(2,i),
     $                 coord_pt(1,i),coord_pt(2,i),coord_pt(3,i),
     $                 vel_pt(1,i),vel_pt(2,i),vel_pt(3,i),
     $           accel_pt(1,i),accel_pt(2,i),accel_pt(3,i),
     $           prevel_pt(1,i),prevel_pt(2,i),prevel_pt(3,i)
 5545          continue

               do 3466 i=0,n_ss_alloc-1 
                  read(25) coord_ss(1,i),coord_ss(2,i),coord_ss(3,i)
                  read(25) flow_ss(i),pressure_ss(i)
 3466          continue

               do 6547 i=1,2*nnd
                  read(25) predrf(i),drf(i)
 6547          continue

               do 6533 ne=1,numel
                  do 6534 k=1,nump
                     nh1=(ne-1)*nump+k
                     read(25) prec(nh1)
 6534             continue
 6533          continue
            endif
         endif
         call r_stang
      endif

      if (n_ibmfem .ne. 1) then
         if (nrestart .eq. 1) then
            if (iti .eq. nt0+1) then
               do 5515 i=1,nptibm
                  read(25) coord_pt(1,i),coord_pt(2,i),coord_pt(3,i), 
     $                 vel_pt(1,i),vel_pt(2,i),vel_pt(3,i),
     $                 accel_pt(1,i),accel_pt(2,i),accel_pt(3,i),
     $                 prevel_pt(1,i),prevel_pt(2,i),prevel_pt(3,i)
 5515          continue
               do 3465 i=0,n_ss_alloc-1 
                  read(25) coord_ss(1,i),coord_ss(2,i),coord_ss(3,i)
                  read(25) flow_ss(i),pressure_ss(i)
 3465          continue
            endif
         endif
      endif

      if (nrestart .eq. 1) then
         if (n_ibmfem .eq. 1) then
            if (iti .eq. nt0+1) then
               read(25) vnorm,viter
            endif
         endif
      endif

c=============================
c ***   begin 1 iteration   ***
c=============================
      do 2000 iiterstep = 1, n_iter_step

         if (n_iter_step .gt. 1) then
            write(6,*) 'beginning iteration ', 
     $           iiterstep, ' of step ', istep
         endif
         
c++++++++ 
c     ** create connection lists, traversing dlist with key=idomain_ptcon.
c   asynchronous
c++++++++ 
         call create_conexchlists(
     $        dnext_pt, dlptlocal_number,
     $        dlptlocal_head, dlptlocal_tail,
     $        idomain_ptcon, dnext_con, dlconshare_number,
     $        dlconshare_head, dlconshare_tail)


c++++++++                   
c     ** assign masks for connection exchanges. mustfollowdacxlists.
c     asynchronous.
c++++++++                   


         call create_shareconmask(
     $        dlconshare_number,  dmasksharecon)
         
         call calctypeactivations(
     $        klok, bwriteon)

c++++++++
c     applying periodic boundary velocity
c     n_periodicvel equal to 1, apply periodic boundary velocity
c     equal to 0, do not apply
c++++++++

         if (n_periodicvel .eq. 1) then
            call periodicvel(u,v,w,dmnac1,dmnac2,
     $           dmnac3,dmxac1,dmxac2,dmxac3,dmnlc1,dmnlc2,
     $           dmnlc3,dmxlc1,dmxlc2,dmxlc3,info,istep,
     $           f1,f2,f3,p)

            call exchangefluiddata(dtid(dilocaldomain),dtidaxisneighbor,  
     $           ibd_fluid_vel,vel_fluid) 

         endif

c++++++++          
c     calculate point forces
c++++++++          
         call calcpointforces(
     $        dilocaldomain, dtid, dmasksharecon,
     $        dnext_pt,dlptlocal_number,
     $        dlptlocal_head,dlptlocal_tail,
     $        dnext_con, dlconshare_number, 
     $        dlconshare_head,dlconshare_tail,
     $        idomain_pt, idomain_ptcon,
     $        coord_pt, coord_ptfcu, 
     $        acttype_con,
     $        stiff0_con,
     $        rest0_con,
     $        fix_con,
     $        forces_fcu,
     $        stiffness_fcu, restlength_fcu, force_pt,
     $        accel_pt,vel_pt)

c++++++++
c     calculate the applied forces
c   n_dispforce equal to 1, use initial displacement.
c   n_dispforce equal to 0, use external excitation force
c++++++++
         if (n_dispforce .eq. 0) then
            
            call disturbance(u,v,w,istep,
     $           mn_pt_alloc, mx_pt_alloc,
     $           dnext_pt,dlptlocal_number, 
     $           dlptlocal_head, dlptlocal_tail, 
     $           dmnac1,dmnac2,dmnac3,dmxac1,dmxac2,dmxac3,
     $           dmnlc1,dmnlc2,dmnlc3,dmxlc1,dmxlc2,dmxlc3,
     $           vel_pt, coord_pt, force_pt, accel_pt)
         endif
 
         call interactviadelta(
     $        mn_pt_alloc, mx_pt_alloc,
     $        dnext_pt,
     $        dlptlocal_number, 
     $        dlptlocal_head, dlptlocal_tail,
     $        ib_spread_force,
     $        dmnac1, dmxac1,
     $        dmnac2, dmxac2,
     $        dmnac3, dmxac3,
     $        coord_pt, 
     $        force_pt, force_fluid)                       

c++++++++
c     ** exchange ghost cells for fluid force.   
c++++++++
         call exchangefluiddata(
     $        dtid(dilocaldomain), dtidaxisneighbor,
     $        ibd_fluid_force,
     $        force_fluid)

c==================================================
c ***   calculate sourcesink/divbalance flows ***
c==================================================

         call calcssbflow(pressure_ss, pressure_rsrvss, 
     $        resist_rsrvss, flow_ss, flow_bal)  

c==================================================
c ***   spread sourcesink/divbalance flows to the lattice ***
c==================================================

         call calcssbdelta(coord_ss,ice1_delss, ice2_delss, 
     $        ice3_delss,wt_delss)          

         call spreadssbflow(dilocaldomain,dmnac1, dmxac1, 
     $        dmnac2, dmxac2, dmnac3, dmxac3, dmnlc1, dmxlc1,
     $        dmnlc2, dmxlc2, dmnlc3, dmxlc3, flow_ss, wt_delss, 
     $        ice1_delss, ice2_delss, ice3_delss,flow_bal,         ! jsr !
     $        dbanybc,dmnbce1, dmxbce1, dmnbce2, dmxbce2, 
     $        dmnbce3, dmxbce3, div_fluid)

         
c++++++++
c     fluid solver
c++++++++     
         call  ibg_solvefluiddynamics(
     $           u,  v,  w,  p, 
     $           f1, f2, f3, q,   
     $           xc, yc, zc, qc,
     $           uc, vc, wc, pc,
     $           dmnac1,  dmnac2,  dmnac3,
     $           dmxac1,  dmxac2,  dmxac3,
     $           dmnlc1,  dmnlc2,  dmnlc3,
     $           dmxlc1,  dmxlc2,  dmxlc3,
     $           dmnam1,  dmnam2,  dmnam3,
     $           dmxam1,  dmxam2,  dmxam3,
     $           dmnmode1, dmnmode2, dmnmode3,
     $           dmxmode1, dmxmode2, dmxmode3,
     $           prdeno, qrfact, xcfact, ycfact, zcfact,
     $           n_table_fft,  table_fft, 
     $           n_work_fft,   workspace_fft,
     $           n_isys_fft,   isys_fft,
     $           info)


            call  ibg_checkcfdspeed(
     $           u, v, w,
     $           dmnac1,  dmnac2,  dmnac3,
     $           dmxac1,  dmxac2,  dmxac3,
     $           dmnlc1,   dmnlc2,   dmnlc3,
     $           dmxlc1,   dmxlc2,   dmxlc3,
     $           info)                                                   

            if (info .ne. 0) then
               write(6,*) 'main: ibg_checkcfdspeed'
               write(6,*) '  returned non-zero value.'
               write(6,*) '  cfl violation is a fatal error'

               call flush(6)
               call exit(1)
            endif
         
            call exchangefluiddata(
     $           dtid(dilocaldomain), dtidaxisneighbor,  
     $           ibd_fluid_vel,
     $           vel_fluid) 


c++++++++
c  **   local interpolation of fiber vels. asynchronous.
c++++++++
            call interactviadelta(
     $           mn_pt_alloc, mx_pt_alloc,
     $           dnext_pt,
     $           dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $           ib_interp_vel,
     $           dmnac1, dmxac1,
     $           dmnac2, dmxac2,
     $           dmnac3, dmxac3,
     $           coord_pt,
     $           vel_pt, vel_fluid)                  

            call interactviadelta(
     $           mn_marker_alloc, mx_marker_alloc,
     $           dnext_mk,
     $           dlmklocal_number, dlmklocal_head, dlmklocal_tail,
     $           ib_interp_vel,
     $           dmnac1, dmxac1,
     $           dmnac2, dmxac2,
     $           dmnac3, dmxac3,
     $           coord_mk,
     $           vel_mk, vel_fluid)

c++++++++
c     calculate the acceleration
c++++++++
            call calcaccel(klok,
     $           mn_pt_alloc,mx_pt_alloc,
     $           dnext_pt,
     $           dlptlocal_number,dlptlocal_head,dlptlocal_tail,
     $           vel_pt,prevel_pt,accel_pt)

c++++++++
c     update solid domain
c     ibm to rubber coordinate convertion
c++++++++
            if (n_ibmfem .eq. 1) then
               viter=0.0d0
               do 301 i=1,nnd
                  if (fix_con(i) .eq. -1) then
                     du(1,i)= 0.0d0
                     du(2,i)= 0.0d0
                  else
                     du(1,i)=vel_pt(1,i)
                     du(2,i)=vel_pt(3,i)
                     acm(1,i)=accel_pt(1,i)
                     acm(2,i)=accel_pt(3,i)
                  endif
                  dis(1,i)=dis(1,i)+du(1,i)
                  dis(2,i)=dis(2,i)+du(2,i)
 
                  if (klok .eq. 1) then
                     vnorm=vnorm+du(1,i)**2+
     $                    du(2,i)**2
                  else
                     viter=viter+du(1,i)**2+
     $                    du(2,i)**2
                  endif
 301           continue
               
               if (klok .eq. 1) then
                  vnorm=dsqrt(vnorm)
                  viter=1.0d0
               else
                  viter=dsqrt(viter)/vnorm
               endif

               write(*,*) viter

               if (niterf .eq. 0) then
                  if (viter .gt. threshold) then                     
                     do 45 i=1,nnd
                        dis(1,i)=dis(1,i)+du(1,i)
                        dis(2,i)=dis(2,i)+du(2,i)
 45                  continue
                     
                     do 33 ne=1,numel
                        do 34 k=1,nump
                           nh1=(ne-1)*nump+k
                           prec(nh1)=prec(nh1)+
     $                          pre(k,ne)
 34                     continue
 33                  continue
                  else
                     niterf=1
                  endif
               else
                  write(*,*) viter,viter
                  do 46 i=1,nnd
                     predrf(i)=predrf2(i)
                     predrf(i+nnd)=predrf2(i+nnd)
                     drf(i)=drf2(i)
                     drf(i+nnd)=drf2(i+nnd)
                     vel_pt(1,i)=0.0d0
                     vel_pt(3,i)=0.0d0
 46               continue
               endif
            endif

c++++++++
c ***   move the lagrangian things ***
c++++++++
            call movepoints(
     $           mn_pt_alloc, mx_pt_alloc,
     $           dnext_pt,
     $           dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $           vel_pt, coord_pt, acttypec_con,fix_con)

            call interpssbpressure(dilocaldomain,dmnac1, dmxac1,
     $           dmnac2, dmxac2, dmnac3, dmxac3,pressure_ss, 
     $           wt_delss,ice1_delss, ice2_delss, ice3_delss, 
     $           pressure_bal,dbanybc, dmnbce1, dmxbce1, dmnbce2, 
     $           dmxbce2,dmnbce3, dmxbce3, pressure_fluid)


c            call calccentroid(mn_cloud_sourcesink, 
c     $           mx_cloud_sourcesink, 
c     $           mn_marker_alloc, mx_marker_alloc, 
c     $           nmk_cloud(mn_cloud_sourcesink), 
c     $           mnmk_cloud(mn_cloud_sourcesink), 
c     $           mxmk_cloud(mn_cloud_sourcesink), 
c     $           coord_mk,coord_ss)

            call movemkpts(
     $           mn_marker_alloc, mx_marker_alloc,
     $           dnext_mk,
     $           dlmklocal_number,  dlmklocal_head, dlmklocal_tail,
     $           vel_mk,coord_mk,acttype_con)

 2000    continue

c++++++++++
c reapply velocity boundary condition
c++++++++++
         if (n_periodicvel .eq. 1) then
            call periodicvel(u,v,w,dmnac1,dmnac2,
     $           dmnac3,dmxac1,dmxac2,dmxac3,dmnlc1,dmnlc2,
     $           dmnlc3,dmxlc1,dmxlc2,dmxlc3,info,istep,
     $           f1,f2,f3,p)        
         endif
            
c++++++++++++++++
c     *** user file i/o - write out .xfdi file for current klok
c++++++++++++++++
c         if( (mod(klok, n_step_wr_ib_user_files) .eq. 0) .or.
c     $        (klok .eq. 1)) then
         
c     if(klok .eq. 1) then
c     currentstep = klok
c     else
c     currentstep = klok/n_step_wr_ib_user_files+1
c     endif    
         
c         if (mod(klok, n_step_wr_ib_user_files) .eq. 0) then

            
c            currentstep = klok/n_step_wr_ib_user_files
c            time_value(currentstep) = klok * td
c++++++++
c matlab output
c+++++++

            call exchangefluiddata(
     $           dtid(dilocaldomain), dtidaxisneighbor,  
     $           ibd_fluid_vel,
     $           vel_fluid) 
            if (n_matlab .eq. 1) then
               call  zmatlab(
     $              klok, td,
     $              ncloud_run,  mncloud_run, mxcloud_run, nmark_run,
     $              nmk_cloud,   mnmk_cloud,  mxmk_cloud,
     $              num_fiber, num_point,
     $              dnext_mk,
     $              dlfreemk_number, dlfreemk_head, dlfreemk_tail,
     $              dlmklocal_number, dlmklocal_head, dlmklocal_tail,
     $              mkpin,
     $              coord_mk,
     $              dnext_pt,
     $              dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $              pt_iptexp,
     $              coord_pt, attrib_fcu,
     $              force_con,force_pt,
     $              vel_pt,accel_pt,f1,f2,f3,tem,
     $              u,  v,  w, p,
     $              mnlatwr1,  mxlatwr1,
     $              mnlatwr2,  mxlatwr2,
     $              mnlatwr3,  mxlatwr3,
     $              dmnac1,    dmnac2,   dmnac3,
     $              dmxac1,    dmxac2,   dmxac3,
     $              dmnlc1,     dmnlc2,    dmnlc3,
     $              dmxlc1,     dmxlc2,    dmxlc3)
            endif
            
c         if (mod(klok, n_step_wr_ib_user_files) .eq. 0) then

         if ((mod(klok, n_step_wr_ib_user_files) .eq. 0) 
     $        .or. (klok .eq. 1)) then
            
            
            currentstep = klok/n_step_wr_ib_user_files
            time_value(currentstep) = klok * td

c++++++++
c     ibm :: tecplot and ensight output
c++++++++
            if (n_ibmfem .eq. 0) then
               
               if( n_tec_ens .eq. 1) then
                  call zibm_tec(
     $                 klok, td,
     $                 ncloud_run,  mncloud_run, mxcloud_run, nmark_run,
     $                 nmk_cloud,   mnmk_cloud,  mxmk_cloud,
     $                 num_fiber, num_point,
     $                 dnext_mk,
     $                 dlfreemk_number, dlfreemk_head, dlfreemk_tail,
     $                 dlmklocal_number, dlmklocal_head, dlmklocal_tail,
     $                 mkpin,
     $                 coord_mk,
     $                 dnext_pt,
     $                 dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $                 pt_iptexp,
     $                 coord_pt, attrib_fcu,
     $                 force_con,force_pt,
     $                 vel_pt,accel_pt,f1,f2,f3,tem,
     $                 u,  v,  w, p,
     $                 mnlatwr1,  mxlatwr1,
     $                 mnlatwr2,  mxlatwr2,
     $                 mnlatwr3,  mxlatwr3,
     $                 dmnac1,    dmnac2,   dmnac3,
     $                 dmxac1,    dmxac2,   dmxac3,
     $                 dmnlc1,     dmnlc2,    dmnlc3,
     $                 dmxlc1,     dmxlc2,    dmxlc3)
               else
                  call  zibm_ensmgeo(
     $                 klok, td,
     $                 ncloud_run,  mncloud_run, mxcloud_run, nmark_run,
     $                 nmk_cloud,   mnmk_cloud,  mxmk_cloud,
     $                 dnext_mk,
     $                 dlfreemk_number, dlfreemk_head, dlfreemk_tail,
     $                 dlmklocal_number, dlmklocal_head, dlmklocal_tail,
     $                 mkpin,
     $                 coord_mk)
                  call  zibm_ensgeo(
     $                 klok, td,
     $                 dnext_pt,
     $                 dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $                 pt_iptexp,
     $                 coord_pt, attrib_fcu)
                  call    zibm_ensfluid(
     $                 u,  v,  w, p,
     $                 klok,
     $                 mnlatwr1,  nskipwr,  mxlatwr1,
     $                 mnlatwr2,  nskipwr,  mxlatwr2,
     $                 mnlatwr3,  nskipwr,  mxlatwr3,
     $                 dmnac1,    dmnac2,   dmnac3,
     $                 dmxac1,    dmxac2,   dmxac3,
     $                 dmnlc1,     dmnlc2,    dmnlc3,
     $                 dmxlc1,     dmxlc2,    dmxlc3)
               endif
            else
c++++++++
c     fem :: tecplot and ensight output
c++++++++
               if( n_tec_ens .eq. 1) then
                  call zfem_tec(
     $                 klok, td,
     $                 ncloud_run,  mncloud_run, mxcloud_run, nmark_run,
     $                 nmk_cloud,   mnmk_cloud,  mxmk_cloud,
     $                 dnext_mk,
     $                 dlfreemk_number, dlfreemk_head, dlfreemk_tail,
     $                 dlmklocal_number, dlmklocal_head, dlmklocal_tail,
     $                 mkpin,
     $                 coord_mk,
     $                 dnext_pt,
     $                 dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $                 pt_iptexp,
     $                 coord_pt, attrib_fcu,
     $                 force_con,force_pt,
     $                 vel_pt,accel_pt,f1,f2,f3,tem,
     $                 u,  v,  w, p,
     $                 mnlatwr1,  nskipwr1,  mxlatwr1,
     $                 mnlatwr2,  nskipwr2,  mxlatwr2,
     $                 mnlatwr3,  nskipwr3,  mxlatwr3,
     $                 dmnac1,    dmnac2,   dmnac3,
     $                 dmxac1,    dmxac2,   dmxac3,
     $                 dmnlc1,     dmnlc2,    dmnlc3,
     $                 dmxlc1,     dmxlc2,    dmxlc3)
               else
c++++++++
c     marker file
c++++++++
                  call  zfem_ensmgeo(
     $                 klok, td,
     $                 ncloud_run,  
     $                 mncloud_run, mxcloud_run, nmark_run,
     $                 nmk_cloud,   mnmk_cloud,  mxmk_cloud,
     $                 dnext_mk,
     $                 dlfreemk_number, 
     $                 dlfreemk_head, 
     $                 dlfreemk_tail,
     $                 dlmklocal_number, 
     $                 dlmklocal_head, 
     $                 dlmklocal_tail,
     $                 mkpin,
     $                 coord_mk)
c++++++++
c     geometry file
c++++++++
                  call  zfem_ensgeo(
     $                 klok, td, dnext_pt,
     $                 dlptlocal_number, 
     $                 dlptlocal_head, 
     $                 dlptlocal_tail,
     $                 pt_iptexp, coord_pt, attrib_fcu)
c++++++++
c     v, p in fluid domain
c++++++++
                  call    zfem_ensfluid(
     $                 u,  v,  w, p, klok,
     $                 mnlatwr1,  nskipwr,  mxlatwr1,
     $                 mnlatwr2,  nskipwr,  mxlatwr2,
     $                 mnlatwr3,  nskipwr,  mxlatwr3,
     $                 dmnac1,    dmnac2,   dmnac3,
     $                 dmxac1,    dmxac2,   dmxac3,
     $                 dmnlc1,     dmnlc2,    dmnlc3,
     $                 dmxlc1,     dmxlc2,    dmxlc3)
                  
c++++++++
c     stress and strain in structure
c++++++++
                  call     zfem_ensstr(klok, 
     $                 dlptlocal_number,
     $                 dlptlocal_head, 
     $                 dlptlocal_tail)
                  
               endif               
            endif

            iitemp=klok/n_step_wr_ib_user_files
            nrc=223+mod(iitemp,2)

            rewind(nrc)
            write(nrc) istep
            do 6310 k = dmnac3 , dmxac3
               do 6320 j = dmnac2, dmxac2
                  do 6330 i = dmnac1, dmxac1
                     write(nrc) u( i, j, k),v( i, j, k),
     $                    w( i, j, k)
 6330             continue
 6320          continue
 6310       continue


            if (n_ibmfem .eq. 1) then
               do 6545 i=1,nnd
                  write(nrc) dis(1,i),dis(2,i),
     $                 du(1,i),du(2,i),
     $                 acm(1,i),acm(2,i),
     $                 coord_pt(1,i),coord_pt(2,i),coord_pt(3,i),
     $                 vel_pt(1,i),vel_pt(2,i),vel_pt(3,i),
     $                 accel_pt(1,i),accel_pt(2,i),accel_pt(3,i),
     $                 prevel_pt(1,i),prevel_pt(2,i),prevel_pt(3,i)
 6545          continue
               do 4465 i=0,n_ss_alloc-1 
                  write(nrc) coord_ss(1,i),coord_ss(2,i),coord_ss(3,i)
                  write(nrc) flow_ss(i),pressure_ss(i)
 4465          continue
            else
               do 6515 i=1,nptibm
                  write(nrc) coord_pt(1,i),coord_pt(2,i),coord_pt(3,i),
     $                 vel_pt(1,i),vel_pt(2,i),vel_pt(3,i),
     $                 accel_pt(1,i),accel_pt(2,i),accel_pt(3,i),
     $                 prevel_pt(1,i),prevel_pt(2,i),prevel_pt(3,i)
 6515          continue

               do 4466 i=0,n_ss_alloc-1 
                  write(nrc) coord_ss(1,i),coord_ss(2,i),coord_ss(3,i)
                  write(nrc) flow_ss(i),pressure_ss(i)
 4466          continue
            endif

            if (n_ibmfem .eq. 1) then
               do 6546 i=1,2*nnd
                  write(nrc) predrf(i),drf(i)
 6546          continue

               do 5533 ne=1,numel
                  do 5534 k=1,nump
                     nh1=(ne-1)*nump+k
                     write(nrc) prec(nh1)
 5534             continue
 5533          continue
               write(nrc) vnorm,viter
            endif
         endif
 1000 continue

      if ((n_ibmfem .eq. 0).and.(n_tec_ens .eq. 0)) then
         write(*, *) 'generate ensight case file'
         call zibm_enscase(time_value, currentstep)
      endif
      if((n_ibmfem .eq. 1).and. (n_tec_ens .eq. 0)) then
         write(*, *) 'generate ensight case file'
         call zfem_enscase(time_value, currentstep)
      endif

      naxx2=time()
      write(*,*) naxx1,naxx2,naxx2-naxx1

      stop
      end 
