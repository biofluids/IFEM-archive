      subroutine pospro(xm,dxm,dvm,xmk,dvmk,nxtime,nytime,
     &             LmTract,vTract,
     &             LmDispbcX,LmDispbcY,vDispbcX,vDispbcY,
     &             LmVelbcX,LmVelbcY,vVelbcX,vVelbcY,
     &             lnods,igp_elemind,igp_LocintInd,gp_loc2D, ! for FEM
     &             nelem,numnp,mgk,nndxdt,nndydt,nTract,
     &             nDispbcX,nDispbcY,nVelbcX,nVelbcY,
     &             iextf_LoadingType,ivel_loadingType,ifwPost)
c     
      implicit double precision (a-h,o-z)
c
      include 'parameter.h'
c
      real*8  xm(2,maxNumnp),dxm(2,maxNumnp),dvm(maxNumnp)
      real*8  xmk(2,maxGP),dvmk(maxGP)
      integer LmTract(2,maxTract)
      real*8  vTract(2,2,maxTract)
      integer LmDispbcX(maxDispbc),LmDispbcY(maxDispbc)
      real*8  vDispbcX(maxDispbc),vDispbcY(maxDispbc)
      integer LmVelbcX(maxVelbc),LmVelbcY(maxVelbc)
      real*8  vVelbcX(maxVelbc),vVelbcY(maxVelbc)
      integer nxtime(mxndt),nytime(mxndt)
      integer iextf_loadingType,ivel_loadingType,ifwEcho
c
      integer igp_elemind(maxGP),igp_locintInd(maxGP)
      integer lnods(maxNode,maxElem)
      real*8  gp_loc2D(2,maxIntElem)
c      
      	!** loc vars
      real*8  gp_weight2D(maxIntElem)
      real*8  shape(maxNode)	! the shape function values
c
c.......Definition of  common blocks
c
      real*8 el_E,el_V,el_G,el_K
      real*8 pl_k0,pl_H,pl_EP  
      real*8 r0,xc1,xc2,zimp
      integer iMaterType,nnode,iSplit
      integer imethRKPM,imethFEM,imeth,iLumping
      integer iIntMeth,iInter,iAdapt,iGL
      integer numadp,iad_stage,numnp
c
      character*120 linestr
c
c........Common Blocks............
c
      common /rkpm/ax,ay
      common /step/dt,nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp
      common /hyperelast/rho0,cc1,cc2,rambda
      common /elastic/ el_E, el_V, el_G, el_K
      common /plastic/ pl_k0, pl_H, pl_EP
      common /shear/r0,xc1,xc2,zimp
      common /algorithm/beta,gamma,theta,iLumping
c      
      common /ctrl/iMaterType,nnode
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iGL 
      common /adapt/iInter,iAdapt,numadp,iad_stage
      
      real*8 xyloc_of_elem4(2,4)
      data ((xyloc_of_elem4(idim,inode),idim=1,2),inode=1,4)
     ./
     .	-1.,-1.,
     .   1.,-1.,
     .   1., 1.,
     .  -1., 1.
     ./
c
c-------------------------------------------------------------
c
      real*8 xyloc_of_elem3(2,3)
      data ((xyloc_of_elem3(idim,inode),idim=1,2),inode=1,3)
     ./
     .	 1.,0.,
     .   0.,1.,
     .   0.,0.
     ./
c            
      integer itmp      
c            
      imethRKPM = 0
      imethFEM  = 1
c
c.....input data
c
      ifrInp = 10
c      
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) iMaterType
      write(ifwPost,'(1x,i5)') iMaterType
c      
c      itmp = iMaterType
c      iMaterType = itmp - itmp/10*10
c      
c      itmp  = itmp/10
c      imeth = itmp - itmp/10*10
c      write(ifwEcho,'(1x,a,i5)') iMaterType
c
c......input method: imeth  : (0: RKPM; 1: FEM )
c
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) imeth
      write(ifwEcho,'(1x,i5)') imeth
c      
c                    iAdapt : (0: without adaptivity) 
c                             (1: h-adaptivity)
c                             (2: p-adaptivity)
c
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) iAdapt
      write(ifwEcho,'(1x,i5)') iAdapt 
c
c  Number of point to refine
c
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) numadp
      write(ifwEcho,'(1x,i5)') numadp 
c
c.......input iInter: for interpolation function
c       for RKPM (1: linear; 11 : bilinear; 2: quadratuc )
c
      if (imeth .eq. imethRKPM) then
         call ReadALine(ifrInp,ifwEcho,linestr)
         read(linestr,*) iInter
         write(ifwEcho,'(1x,i5)') iInter
      endif
c
c......input 
c
		!** input nnode
      write(*,*) 'nnode=3, 4, 32(split2) or 34(split4)'
c      write(*,*) 'input nnode='
c      read(*,*) nnode
      call ReadALine(ifrinp,ifwEcho,linestr)
      read(linestr,*) nnode
      write(ifwEcho,'(1x,i5)') nnode
      write(*,*) 'Now nnode=',nnode
      if (nnode.eq.32) then
         iSplit=2
         nnode=3
      elseif (nnode.eq.34) then
         iSplit=4
         nnode=3    
      else
         iSplit=0
      endif
c
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) numnp
      write(ifwEcho,'(1x,i5)') numnp
      call checklim(numnp,maxNumnp,'numnp')
c
c##
c
      call ReadALine(ifrInp,ifwEcho,linestr)   
      read(linestr,*) nelem
      write(ifwEcho,'(1x,i5)') nelem
      call checklim(nelem,maxElem,'nelem')

      call ReadALine(ifrInp,ifwEcho,linestr)          
      read(linestr,*) nint      
      write(ifwEcho,'(1x,i5)') nint    
      if (nnode.eq.3) then
         nintElem=nint
      elseif (nnode.eq.4) then
         call checklim(nint,maxInt,'nint')
         nintElem=nint*nint
      endif    
      call checklim(nintElem,maxIntElem,'nintElem')  
      
      mgk =  nintElem * nelem
      if ( (nnode.eq.3) .and. (iSplit.ne.0) ) then
         mgk=mgk*iSplit
      endif
      
      call checklim(mgk,maxGP,'mgk(magGP)')

      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) dt,nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp
      write(ifwEcho,'(1x,e13.6,5i8)') dt,nstep,nOutputFrq,nPlotFrq,
     &                                iplot,iplot_gp

      call ReadALine(ifrInp,ifwEcho,linestr)          
      read(linestr,*) ax,ay
      write(ifwEcho,'(1x,2e13.6)') ax,ay


      call ReadALine(ifrInp,ifwEcho,linestr)             
      if ( iMaterType.eq.1 ) then 	
      			!** hyperelasticity
         read(linestr,*) rho0,cc1,cc2,rambda
         write(ifwEcho,'(1x,4e13.6)') rho0,cc1,cc2,rambda
         
      elseif ( iMaterType.eq.2 ) then	
      			!** elasticity-plasticity
         read(linestr,*) rho0,el_E,el_V,pl_k0,pl_EP
         write(ifwEcho,'(1x,5e13.6)') rho0,el_E,el_V,pl_k0,pl_EP
         
         	! shear modulus
         el_G=el_E / (2. * (1.+el_V) )
         	! volumn modulus
         el_K=el_E / (3. * (1.-2.*el_V) )
         
         	! hardening modulus
         pl_EP= el_E*pl_EP/(el_E-pl_EP)
         pl_H = 0.			! isotropic hardening
         
      elseif (iMaterType.eq.3) then	
      					! visco-plastic-plasticity
         read(linestr,*) rho0,el_E,el_V
         write(ifwEcho,'(1x,5e13.6)') rho0,el_E,el_V
         
         	! shear modulus
         el_G = el_E / (2. * (1. + el_V) )
         	! volumn modulus
         el_K = el_E / (3. * (1. - 2.*el_V) )
                                        
      else
         write(*,*) 'unknown iMaterType=',iMaterType
         stop	
      end if
c        
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) iIntMeth
      write(ifwEcho,'(1x,i8)') iIntMeth
c      
      if (iIntMeth .eq. 1) then	
      		!  predictor-corrector
         call ReadALine(ifrInp,ifwEcho,linestr)      
         read(linestr,*) beta,gamma
         write(ifwEcho,'(1x,2e13.6)') beta,gamma
      elseif (iIntMeth .eq. 2) then
     		!  central difference
                !  no time integration control parameter
      elseif (iIntMeth .eq. 3) then
      		! improved forward euler
         call ReadALine(ifrInp,ifwEcho,linestr)      
         read(linestr,*) theta
         write(ifwEcho,'(1x,e13.6)') theta
      else
         write(*,*) 'Unknown iIntMeth=',iIntMeth
         stop   
      endif                 

      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) nndxdt
      write(ifwEcho,'(1x,i8)') nndxdt
      call checklim(nndxdt,mxndt,'nndxdt(mxndt)')
      do ip = 1 , nndxdt
        call ReadALine(ifrInp,ifwEcho,linestr)
        read(linestr,*) nxtime(ip)
        write(ifwEcho,'(1x,i8)') nxtime(ip)
      enddo

      call ReadALine(ifrInp,ifwEcho,linestr)    
      read(linestr,*) nndydt
      write(ifwEcho,'(1x,i8)') nndydt
      call checklim(nndydt,mxndt,'nndydt(mxndt)')
      do ip = 1 , nndydt
        call ReadALine(ifrInp,ifwEcho,linestr)
        read(linestr,*) nytime(ip)
        write(ifwEcho,'(1x,i8)') nytime(ip)
      enddo

		!** read loading type
      call ReadALine(ifrInp,ifwEcho,linestr)                     
      read(linestr,*) iextf_LoadingType,ivel_LoadingType
      write(ifwEcho,'(1x,2i8)') iextf_LoadingType,ivel_LoadingType

		
		!** read traction BC              
      call ReadALine(ifrInp,ifwEcho,linestr)        
      read(linestr,*) nTract
      write(ifwEcho,'(1x,i8)') nTract
      call checklim(nTract,maxTract,'ntract')
      
      if (nTract.gt.0) then
      
            	!** read traction ratex, ratey
         call ReadALine(ifrInp,ifwEcho,linestr)                      
         read(linestr,*) tract_ratex,tract_ratey
         write(ifwEcho,'(1x,2e15.6)') tract_ratex, tract_ratey

         
         do iTract = 1 , nTract
      
     			!** read: ip1,ip2,fx1,fy1,fx2,fy2 
            call ReadALine(ifrInp,ifwEcho,linestr)     
            read(linestr,*) (LmTract(ip,iTract),ip=1,2),
     &             ( (vTract(idim,ip,iTract),idim=1,2), ip=1,2)	
            do ip=1,2
               vTract(1,ip,iTract)=vTract(1,ip,iTract)*tract_ratex
               vTract(2,ip,iTract)=vTract(2,ip,iTract)*tract_ratey
            enddo
            
            write(ifwEcho,'(1x,2i8,4e11.3)') 
     &             (LmTract(ip,iTract),ip=1,2),
     &             ( (vTract(idim,ip,iTract),idim=1,2), ip=1,2)	

            
         enddo
      endif
      
      		!** read Essen BC (disp & vel )
      call ReadALine(ifrInp,ifwEcho,linestr)               
      read(linestr,*) nDispbcX,nDispbcY, nVelbcX, nVelbcY
      write(ifwEcho,'(1x,4i8)') nDispbcX,nDispbcY,nVelbcX, nVelbcY
      call checklim(nDispbcX,maxDispbc,'nDispbcX')
      call checklim(nDispbcY,maxDispbc,'nDispbcY')
      call checklim(nVelbcX,maxVelbc,'nVelbcX')
      call checklim(nVelbcY,maxVelbc,'nVelbcY')
      
      if (nDispbcX.gt.0) then
                       
         do iDispbc=1,nDispbcX
         
            call ReadALine(ifrInp,ifwEcho,linestr)  
            read(linestr,*) LmDispbcX(iDispbc),vDispbcX(iDispbc)
            write(ifwEcho,'(1x,i8,e13.6)') 
     &           LmDispbcX(iDispbc),vDispbcX(iDispbc)
         enddo
      endif
      
      if (nDispbcY.gt.0) then
                        
         do iDispbc=1,nDispbcY
            call ReadALine(ifrInp,ifwEcho,linestr)  
            read(linestr,*) LmDispbcY(iDispbc),vDispbcY(iDispbc)
            write(ifwEcho,'(1x,i8,e13.6)') 
     &           LmDispbcY(iDispbc),vDispbcY(iDispbc)
         enddo
      endif
      
      		!** read velbc ratex, ratey
      if ( (nVelbcX.ne.0) .or. (nVelbcY.ne.0) ) then
         call ReadALine(ifrInp,ifwEcho,linestr)                        
         read(linestr,*) velbc_ratex,velbc_ratey
         write(ifwEcho,'(1x,2e15.6)') velbc_ratex, velbc_ratey
      endif
      
      if (nVelbcX.ne.0) then
           
         do iVelbc=1,nVelbcX    
            call ReadALine(ifrInp,ifwEcho,linestr)       
            read(linestr,*) LmVelbcX(iVelbc),vVelbcX(ivelbc)
            
            vVelbcX(ivelbc)=vVelbcX(ivelbc)*velbc_ratex

            write(ifwEcho,'(1x,i8,e15.6)') 
     &            LmVelbcX(iVelbc),vVelbcX(ivelbc)     
         enddo
      endif
      if (nVelbcY.ne.0) then
         
         do iVelbc=1,nVelbcY
            call ReadALine(ifrInp,ifwEcho,linestr)  
            read(linestr,*) LmVelbcY(iVelbc),vVelbcY(ivelbc)
     
            vVelbcY(ivelbc)=vVelbcY(ivelbc)*velbc_ratey 
    
            write(ifwEcho,'(1x,i8,e15.6)') 
     &            LmVelbcY(iVelbc),vVelbcY(ivelbc)

         enddo
      endif
      
      

		!** read nodal coords

      call ReadALine(ifrInp,ifwEcho,linestr)  
      read(linestr,*) coord_ratex,coord_ratey
      
      do ip=1,numnp
        call ReadALine(ifrInp,ifwEcho,linestr)  
        read(linestr,*) jp,xm(1,ip),xm(2,ip)
        
         xm(1,ip)=xm(1,ip)*coord_ratex
         xm(2,ip)=xm(2,ip)*coord_ratey

        
        write(ifwEcho,'(1x,i8,2e15.6)') jp,xm(1,ip),xm(2,ip)
        if (ip.ne.jp) then
           write(*,*) 'ip.ne.jp when read nodal coords'
           write(*,*) 'ip=',ip,' jp=',jp
           stop
        endif
      end do
      

           
           	!** read element connectivity info
      if (nnode.eq.4) then
         mnode=4
      elseif (nnode.eq.3) then
         if (iSplit.eq.0) then
            mnode=3
         else
            mnode=4
         endif
      endif
c
      do ie=1,nelem
        call ReadALine(ifrInp,ifwEcho,linestr)  
        read(linestr,*) je,
     &                  (lnods(inode,ie),inode=1,mnode)
        write(ifwEcho,'(1x,5i8)') je,(lnods(inode,ie),inode=1,mnode)
        if (ie.ne.je) then
           write(*,*) 'ie.ne.je when read element connectivity'
           write(*,*) 'ie=',ie,'je=',je
           stop
        endif
        
      end do
      
      if ( (iSplit.ne.0). and. (nnode.eq.3) ) then
         if ( iSplit.eq.2 ) then
            call split2tri(nelem,lnods,xm)
         elseif (iSplit.eq.4) then
            call split4tri(nelem,numnp,lnods,xm)               
         else
            write(*,*) 'Unknown number of iSplit=',iSplit
            stop
         endif
         
         	! output the new element conectivity
         do ie=1,nelem
            write(ifwEcho,'(1x,5i8)') ie,(lnods(inode,ie),inode=1,nnode)
         enddo   
         
      endif
c      
      do ip = 1, numnp
        dxm(1,ip) = ax*dabs(xm(1,nxtime(2))-xm(1,nxtime(1)))
        dxm(2,ip) = ay*dabs(xm(2,nytime(2))-xm(2,nytime(1)))
        dvm(ip) = 0.0d0
      end do
c
      kg = 0
      
      call gauss2D(gp_loc2D,gp_weight2D,nintElem,nnode)
      
      do ie=1,nelem

        xn1 = xm(1,lnods(1,ie))
        yn1 = xm(2,lnods(1,ie))
        xn2 = xm(1,lnods(2,ie))
        yn2 = xm(2,lnods(2,ie))
        xn3 = xm(1,lnods(3,ie))
        yn3 = xm(2,lnods(3,ie))
        
        if (nnode.eq.4) then
           xn4 = xm(1,lnods(4,ie))
           yn4 = xm(2,lnods(4,ie))
        endif

        iLocint=0
c        do ig=1,nint
c           do jg=1,nint
c              xsi=gp_loc1D(ig)
c              eta=gp_loc1D(jg)
        do iLocint=1,nintElem      
c              iLocint=iLocint+1
              kg=kg+1
              
              xsi=gp_loc2D(1,iLocint)
              eta=gp_loc2D(2,iLocint)
              
              if (imeth.eq.imethFEM) then
                 igp_elemind(kg)=ie
                 igp_LocintInd(kg)=iLocint
                 
c                 gp_locxy(1,kg)=xsi
c                 gp_locxy(2,kg)=eta
              endif
              
              if (nnode.eq.4) then
                 call evl_FEM_shape4(shape,xsi,eta)
                 xmk(1,kg)=shape(1)*xn1+shape(2)*xn2
     &                 +shape(3)*xn3+shape(4)*xn4
                 xmk(2,kg)=shape(1)*yn1+shape(2)*yn2
     &                 +shape(3)*yn3+shape(4)*yn4
                 call ajacob4(ajj,xn1,xn2,xn3,xn4,yn1,yn2,yn3,yn4,
     &                        xsi,eta)
                 dvmk(kg)=ajj*gp_weight2D(iLocint)
              elseif (nnode.eq.3) then
                 call evl_FEM_shape3(shape,xsi,eta)
                 xmk(1,kg)=shape(1)*xn1+shape(2)*xn2
     &                 +shape(3)*xn3
                 xmk(2,kg)=shape(1)*yn1+shape(2)*yn2
     &                 +shape(3)*yn3
                 call ajacob3(ajj,xn1,xn2,xn3,yn1,yn2,yn3,xsi,eta)
                 dvmk(kg)=ajj*gp_weight2D(iLocint)/2.
              endif
                  			
c              dvmk(kg)=ajj*gp_weight(ig)*gp_weight(jg) 

        enddo ! iLocint
c           enddo ! jg
c        enddo  ! ig
                

		! weight of nodes when by nodal integration
        if (nnode.eq.4) then
           do inode=1,nnode
              jnode=lnods(inode,ie)
              call ajacob4(ajj,xn1,xn2,xn3,xn4,yn1,yn2,yn3,yn4,
     &                    xyloc_of_elem4(1,inode),
     &                    xyloc_of_elem4(2,inode))
              dvm(jnode) = dvm(jnode)+ajj
           enddo   
        elseif(nnode.eq.3) then
           do inode=1,nnode
              jnode=lnods(inode,ie)
              call ajacob3(ajj,xn1,xn2,xn3,yn1,yn2,yn3,
     &                    xyloc_of_elem3(1,inode),
     &                    xyloc_of_elem3(2,inode))
              dvm(jnode)=dvm(jnode)+ajj/2.
           enddo
        else
        endif   
           
      end do	! endfor ie
      
      mgk=kg
      call checklim(mgk,maxGP,'mgk(magGP)')
c
 4040 format(5x,'#TITLE: DYNA2D RKPM-FEM ANALYSIS IN SOLIDS'//)
 4050 format(/5x,'#iProbType    1: hyper-elastic;',
     +       /5x,'#             2: elasto-plastic;',
     +       /5x,'#             3: elasto-visco-plastic;',
	    //5x, i5)
 4060 format(/5x,'#imeth      0: RKPM',
     +       /5x,'#           1: FEM',
     +       /5x,'#iLumping   0: normal lumping',
     +       /5x,'#           1: special lumping',
     +      //5x,i5, 5x, i5)
c
 4070 format(/5x,'#iadapt  0: without adaptivity;',
     +       /5x,'#        1: p-adaptivity;', 
     +       /5x,'#        2: h-adaptivity;', 
     +       /5x,'# numdap; ', 
     +       /5x,'# iGL    0: global', 
     +       /5x,'#        1: local', 
     +      //5x,3(i5,2x))
c
 4075 format(/5x,'#r0 xc1 xc2 zimp',
     +      //5x,3(e13.6,2x))
c
 4080 format(/5x,'#r0 , the magnitude of imperfection;',
     +      //5x, e14.6)
c
 4090 format(/5x,'#iInter    1: Linear;',
     +       /5x,'#         11: Bilinear;',
     +       /5x,'#          2: quadratic;',
     +      //5x, i5)
c
 4100 format(/5x,'#nnode     3: triangle;',
     +       /5x,'#          4: quadrilateral;',
     +       /5x,'#         32: quadrilateral split into 2 trangles;',
     +      //5x, i5)
c
 4110 format(/5x,'#NUMBER OF NODES',
     +      //5x,i5)
c
 4120 format(/5x,'#NUMBER OF ELEMENTS',
     +      //5x,i5)
c
 4130 format(/5x,'#Gauss pt along each direction of the cell',
     +      //5x,i5)
c
 4140 format(/5x,'#dt   nstep   nOutputFreq  nPlotFrq  iplot iplot_gp',
     +      //2x,e14.6,2x,5(i7,ix))
c
 4150 format(/5x,'# ax   ay ',
     +      //5x,2(e14.6,2x))
c
 4155 format(/5x,'# EP mater: rho,  E,  V,  k0, EP',
     +      //5x,5(e14.6,2x))
c
 4160 format(/5x,'# iIntMeth    1: predictor-corrector;',
     +       /5x,'#             2: central difference;',
     +       /5x,'#             3: improved Euler-thera;',
     +      //5x, i5)
c
 4170 format(/5x,'# beta, gamma, theta;',
     +      //5x,3(e14.6,2x) )
c
 4180 format(/5x,'# define the dilation parameters',
     +       /5x,'# x-dircetion')
 4190 format(/5x,'# y-direction')
c
 4110 format(/5x,'# Infomation about traction BC#',
     +       /5x,'# loading type:  1: static;',
     +       /5x,'#                2: heaviside;',
     +       /5x,'#                3: circular smoothed bilinear;',
     +       /5x,'#                4: impulse;',
     +      //5x,2(i5,1x))
c
 4120 format(/5x,'# Traction boundary: nTract:',
     +      //5x, i5)
c
 4140 format(/5x,'# Information about essential BC',
     +       /5x,'# (both displacement and velocity)',
     +       /5x,'# Essen BC: nDispbxX, nDispbcY nVelbcX, nVelbcY',
     +      //5x, 4(i5,2x))
c
 4150 format(/5x,'
	   )

      close(24)
      close(25)
      close(ifwPost)
c      
      end
c


