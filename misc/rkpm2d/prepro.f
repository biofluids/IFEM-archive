      subroutine prepro(xm,dxm,dvm,xmk,dvmk,nxtime,nytime,
     &             LmTract,vTract,
     &             LmDispbcX,LmDispbcY,vDispbcX,vDispbcY,
     &             LmVelbcX,LmVelbcY,vVelbcX,vVelbcY,
     &             lnods,igp_elemind,igp_LocintInd,gp_loc2D, ! for FEM
     &             nelem,numnp,mgk,nndxdt,nndydt,
     &             iextf_LoadingType,ivel_loadingType,ifwEcho)
c     
      implicit none
      include 'parameter.h'
c
      real*8  xm(2,maxNumnp),dxm(2,maxNumnp),dvm(maxNumnp),
     &        xmk(2,maxGP),dvmk(maxGP)
c
      integer LmTract(2,maxTract),
     &        LmDispbcX(maxDispbcX),LmDispbcY(maxDispbcY),
     &        LmVelbcX(maxVelbcX),LmVelbcY(maxVelbcY)
      real*8  vTract(2,2,maxTract),
     &        vDispbcX(maxDispbcX),vDispbcY(maxDispbcY),
     &        vVelbcX(maxVelbcX),vVelbcY(maxVelbcY)
c
      integer nxtime(mxndt),nytime(mxndt)
      integer iextf_loadingType,ivel_loadingType,ifwEcho
c
      integer igp_elemind(maxGP),igp_locintInd(maxGP),
     &        lgops(maxIntElem,maxElem),gpelem(maxGP)
      real*8  gp_loc2D(2,maxIntElem),
     &        gp_weight2D(maxIntElem),
     &        shape(maxNode)	     ! the shape function values
      integer lnods(maxNode,maxElem)
c      
      	!** loc vars
c
      integer numnp,nelem,mgk,nint,
     &        nndxdt,nndydt
c
      integer idim,inode,ip,ipt,
     &        itract,idispbc,ie,je,jp,kg,
     &        ivelbc,mnode,iLocint,jnode,ik,
     &        ifrInp,iSplit
c
      real*8 tract_ratex,tract_ratey,
     &       velbc_ratex,velbc_ratey,coord_ratex,coord_ratey,
     &       xn1,xn2,xn3,xn4,yn1,yn2,yn3,yn4,
     &       dvtt1,dvtt2,dvt,xsi,eta,ajj
c
c.......Definition of  common blocks
c
      real*8 ax,ay,afact1,afact2,
     &       dt,time,
     &       rho0,cc1,cc2,rambda,
     &       el_E,el_V,el_G,el_K,
     &       pl_k0,pl_H,pl_EP,
     &       r0,xc1,xc2,zimp,cr1,cr2, 
     &       beta,gamma,theta,
     &       temp,Temp_0,Cp
c
      integer nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp,
     &        iMaterType,nnode,nintElem,iLumping,
     &        nnc1,nnc2,
     &        imethRKPM,imethFEM,imeth,iIntMeth,iGL, 
     &        iInter,iAdapt,numadp,iad_stage,Istage,
     &        nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
c
       integer LmTractF(maxTF),nTractf,ip1,ip2
       real*8  bweight(maxTF),tempTF(maxTF)
c
       integer iCrack
       real*8  Cyplane,Cxmin,Cxmax
c
      character*120 linestr
c
c........Common Blocks............
c
      common /rkpm/ax,ay,afact1,afact2,nnc1,nnc2
      common /stepT/dt,time
      common /hyperelast/rho0,cc1,cc2,rambda
      common /elastic/el_E,el_V,el_G,el_K
      common /plastic/pl_k0,pl_H,pl_EP
      common /temperature/Temp_0,Cp
      common /shear/r0,xc1,xc2,zimp,cr1,cr2
      common /algorithm/beta,gamma,theta
      common /crack/Cyplane,Cxmin,Cxmax
c
      common /step/nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp
      common /ctrl/iMaterType,nnode,nintElem,iLumping
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter 
      common /adapt/iGL,iAdapt,numadp,iad_stage,Istage
      common /bound/nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
      common /force/nTractf,LmTractF
      common /connectK/lgops,gpelem
      common /Type1/iCrack
      common /Type2/bweight
c      
      real*8 xyloc_of_elem4(2,4)
      data ((xyloc_of_elem4(idim,inode),idim=1,2),inode=1,4)
     ./
     .	-1.,-1.,
     .   1.,-1.,
     .   1., 1.,
     .  -1., 1.
     ./
c
c------2x2 pt. Gaussian quadrature for 4-node elements
c
c			y	
c			^
c			|
c			|
c			|
c			|
c               4---------------3
c               |		|
c               |  4	    3	|
c               | 		|
c        	|		| -------------------> x
c		|		|
c		|  1	    2	|
c		|		|
c		1---------------2
c
c
c----------------------------------------------------
c
      real*8 xyloc_of_elem3(2,3)
      data ((xyloc_of_elem3(idim,inode),idim=1,2),inode=1,3)
     ./
     .	 1.,0.,
     .   0.,1.,
     .   0.,0.
     ./
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
      write(ifwEcho,'(1x,i5)') iMaterType
c      
c      iMaterType -----  1: hyper-elastic;
c                        2: elasto-plastic;
c                        3: visco-plastic;
c                        4: thermal-viscoelastic-plastic
c                        5: kaltoff
c      
c......input method: imeth  : (0: RKPM; 1: FEM )
c
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) imeth, iLumping
      write(ifwEcho,'(1x,2i5)') imeth, iLumping
c      
c                    iAdapt : (0: without adaptivity) 
c                             (1: h-adaptivity)
c                             (2: p-adaptivity)
c
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) iAdapt, numadp, iGL, iad_stage
      write(ifwEcho,'(1x,4i5)') iAdapt,numadp,iGL,iad_stage 
c
c   For shear-band simulation
c
      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) r0, xc1, xc2, zimp, cr1,cr2
      write(ifwEcho,'(1x,6e13.6)') r0,xc1,xc2,zimp,cr1,cr2
c
c.......input iInter: for interpolation function
c       for RKPM (1: linear; 11 : bilinear; 2: quadratuc )
c
         call ReadALine(ifrInp,ifwEcho,linestr)
         read(linestr,*) iInter
         write(ifwEcho,'(1x,i5)') iInter
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
         mgk = mgk*iSplit
      endif
      
      call checklim(mgk,maxGP,'mgk(magGP)')

      call ReadALine(ifrInp,ifwEcho,linestr)
      read(linestr,*) dt,nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp
      write(ifwEcho,'(1x,e13.6,5i8)') dt,nstep,nOutputFrq,nPlotFrq,
     &                                iplot,iplot_gp

      call ReadALine(ifrInp,ifwEcho,linestr)          
      read(linestr,*) ax,ay,afact1,afact2,nnc1,nnc2
      write(ifwEcho,'(1x,4e13.6,2i5)') ax,ay,afact1,afact2,nnc1,nnc2


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
         read(linestr,*) rho0,el_E,el_V,pl_k0
         write(ifwEcho,'(1x,4e13.6)') rho0,el_E,el_V,pl_k0
         
         	! shear modulus
         el_G = el_E / (2. * (1. + el_V) )
         	! volumn modulus
         el_K = el_E / (3. * (1. - 2.*el_V) )
                                        
      elseif (iMaterType.eq.4) then	
      			! thermo-elasto-visco-plasticity
         read(linestr,*) rho0,el_E,el_V,pl_k0,Temp_0,Cp
         write(ifwEcho,'(1x,6e13.6)') rho0,el_E,el_V,pl_k0,Temp_0,Cp
         
         	! shear modulus
         el_G = el_E / (2. * (1. + el_V) )
         	! volumn modulus
         el_K = el_E / (3. * (1. - 2.*el_V) )
      elseif (iMaterType.eq.5) then	
      			! thermo-elasto-visco-plasticity
         read(linestr,*) rho0,el_E,el_V,pl_k0,Temp_0,Cp
         write(ifwEcho,'(1x,6e13.6)') rho0,el_E,el_V,pl_k0,Temp_0,Cp
         
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
c......set the default
c
       theta = 0.50
       beta  = 0.00
       gamma = 0.50
c
      if (iIntMeth .eq. 1) then	
      		!  predictor-corrector
         call ReadALine(ifrInp,ifwEcho,linestr)      
         read(linestr,*) beta,gamma
         write(ifwEcho,'(1x,2e13.6)') beta,gamma
      elseif (iIntMeth .eq. 2) then
     		!  central difference
         call ReadALine(ifrInp,ifwEcho,linestr)      
         read(linestr,*) theta
         write(ifwEcho,'(1x,e13.6)') theta
      elseif (iIntMeth .eq. 3) then
      		! improved forward euler
         call ReadALine(ifrInp,ifwEcho,linestr)      
         read(linestr,*) theta
         write(ifwEcho,'(1x,e13.6)') theta
      else
         write(*,*) 'Unknown iIntMeth=',iIntMeth
         stop   
      endif                 
c
      call ReadAline(ifrInp,ifwEcho,linestr)
      read(linestr,*) iCrack
	    write(ifwEcho,'(1x,i8)') iCrack
c
      if (iCrack .eq. 1) then
          call ReadAline(ifrInp,ifwEcho,linestr)
	     read(linestr,*) Cyplane,Cxmin,Cxmax
	     write(ifwEcho,'(1x,3e14.8)') Cyplane,Cxmin,Cxmax
      endif
c
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
      read(linestr,*) iextf_loadingType,ivel_loadingType
      write(ifwEcho,'(1x,2i8)') iextf_loadingType,ivel_loadingType

		
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
c            
         enddo
      endif
c
c
c........This is for force calulation (assume all the arrangement
c        is sequential )
c
               !** read force calulation array
       call ReadALine(ifrInp,ifwEcho,linestr)
       read(linestr,*) nTractf
       write(ifwEcho,'(1x,i8)') nTractf
       call checklim(nTract,maxTract,'ntractf')
c
       if(nTractf .ne. 0) then
c
        do iTract = 1, nTractf
             call ReadALine(ifrInp,ifwEcho,linestr)
             read(linestr,*) LmTractF(iTract)
             write(ifwEcho,'(1x,i8)') LmTractF(iTract)
c
             bweight(iTract) = 0.0d0
c
        enddo
c
        endif
c
c      
      		!** read Essen BC (disp & vel )
      call ReadALine(ifrInp,ifwEcho,linestr)               
      read(linestr,*) nDispbcX,nDispbcY, nVelbcX, nVelbcY
      write(ifwEcho,'(1x,4i8)') nDispbcX,nDispbcY,nVelbcX, nVelbcY
      call checklim(nDispbcX,maxDispbcX,'nDispbcX')
      call checklim(nDispbcY,maxDispbcY,'nDispbcY')
      call checklim(nVelbcX,maxVelbcX,'nVelbcX')
      call checklim(nVelbcY,maxVelbcY,'nVelbcY')
      
      if (nDispbcX.gt.0) then
          
c
         do iDispbc=1,nDispbcX
         
            call ReadALine(ifrInp,ifwEcho,linestr)  
            read(linestr,*) LmDispbcX(iDispbc),vDispbcX(iDispbc)
            write(ifwEcho,'(1x,i8,e13.6)') 
     &           LmDispbcX(iDispbc),vDispbcX(iDispbc)
         enddo
      endif
c
c
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
c
      if (nVelbcY.ne.0) then
         
         do iVelbc=1,nVelbcY
            call ReadALine(ifrInp,ifwEcho,linestr)  
            read(linestr,*) LmVelbcY(iVelbc),vVelbcY(ivelbc)
     
            vVelbcY(ivelbc)=vVelbcY(ivelbc)*velbc_ratey 
    
            write(ifwEcho,'(1x,i8,e15.6)') 
     &            LmVelbcY(iVelbc),vVelbcY(ivelbc)

         enddo
      endif
c      
c
		!** read nodal coords
c
      call ReadALine(ifrInp,ifwEcho,linestr)  
      read(linestr,*) coord_ratex,coord_ratey
      write(ifwEcho,'(1x,2e15.6)') coord_ratex,coord_ratey
c
      do ip=1,numnp
        call ReadALine(ifrInp,ifwEcho,linestr)  
        read(linestr,*) jp,xm(1,ip),xm(2,ip),temp
        
         xm(1,ip)=xm(1,ip)*coord_ratex
         xm(2,ip)=xm(2,ip)*coord_ratey
c
        write(ifwEcho,'(1x,i8,2e15.6)') jp,xm(1,ip),xm(2,ip)
        if (ip.ne.jp) then
           write(*,*) 'ip.ne.jp when read nodal coords'
           write(*,*) 'ip=',ip,' jp=',jp
           stop
        endif
      end do
c      
c           
c..........!** read element connectivity info
c		
c
      if (nnode.eq.4) then
         mnode=4
      elseif (nnode.eq.3) then
         if (iSplit.eq.0) then
            mnode=3
         else
            mnode=4
         endif
      endif
      do ie=1,nelem
        call ReadALine(ifrInp,ifwEcho,linestr)  
        read(linestr,*) je,
     &                  (lnods(inode,je),inode=1,mnode)
        write(ifwEcho,'(1x,5i8)') je,(lnods(inode,je),inode=1,mnode)
      end do
c
c
cNote this change is made on July 23
c
c        if (ie.ne.je) then
c           write(*,*) 'ie.ne.je when read element connectivity'
c           write(*,*) 'ie=',ie,' je=',je
c           stop
c        endif
c        
c      
      if ((iSplit.ne.0). and. (nnode.eq.3)) then
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
        dvm(ip)   = 0.0d0
      end do
c
c      
      call gauss2D(gp_loc2D,gp_weight2D,nintElem,nnode)
c      
      kg = 0   ! total numerbing of Gauss quadrature pt.
      do ie= 1,nelem
c
        xn1 = xm(1,lnods(1,ie))
        yn1 = xm(2,lnods(1,ie))
        xn2 = xm(1,lnods(2,ie))
        yn2 = xm(2,lnods(2,ie))
        xn3 = xm(1,lnods(3,ie))
        yn3 = xm(2,lnods(3,ie))
c        
        if (nnode.eq.4) then
            xn4 = xm(1,lnods(4,ie))
            yn4 = xm(2,lnods(4,ie))
        endif

        do iLocint=1,nintElem      
              kg=kg+1
c
c....form the connectivity for Gauss pt.
c
	      lgops(iLocint,ie) = kg   !
              gpelem(kg)        = ie   !
c
c              
              xsi = gp_loc2D(1,iLocint)
              eta = gp_loc2D(2,iLocint)
c              
              if (imeth .eq. imethFEM) then
                 igp_elemind(kg)   = ie
                 igp_LocintInd(kg) = iLocint
              endif
c              
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
c      
c
        enddo ! iLocint
c                
c		! weight of nodes when by nodal integration
c
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
c      
       if(nTractf .ne. 0) then
c
	do ip = 1, nTractf
	   ipt        = LmTractF(ip)
	   tempTF(ip) = xm(1,ipt)
	enddo
c
	call permuteTF(tempTF,LmTractF,nTractf,maxTF)
c
        do iTract = 2, nTractf - 1
           ip1 = LmTractF(iTract-1)
           ip2 = LmTractF(iTract+1)
c
             bweight(iTract)  = 0.5*dabs(xm(1,ip2) - xm(1,ip1))
c
        enddo
c
        bweight(1) = 0.5 * abs(xm(1,LmTractF(2))
     &                    - xm(1,LmTractF(1)))
c
        bweight(nTractf) = 0.5 * abs(xm(1,LmTractF(nTractf))
     &                    - xm(1,LmTractF(nTractf-1)))
        endif
c
      mgk=kg
      call checklim(mgk,maxGP,'mgk(magGP)')
c
c---------------------------------------------
c    output original mesh
c---------------------------------------------
c
      if (1 .eq. 1) then
         dvtt1 = 0.0d0
         dvtt2 = 0.0d0
         dvt   = 0.0d0
c
c
	 write(25,'(1x,a)')  'VARIABLES="X","Y","V"'
         if (nnode.eq.4) then
	     write(25,'(1x,a,i6,a,i6,a)')
     &            'ZONE N =',numnp,', E =',nelem,
     &            'F = FEPOINT, ET = QUADRILATERAL'
	 elseif (nnode.eq.3) then
	     write(25,'(1x,a,i6,a,i6,a)')
     &            'ZONE N =',numnp,', E =',nelem,
     &            'F = FEPOINT, ET = TRIANGLE'
	 endif
c
         do ip = 1, numnp
            dvtt1 = dvtt1 + dvm(ip)
            write(25,9201)xm(1,ip),xm(2,ip),dvm(ip)
         enddo
         do ie=1,nelem
	    write(25,'(1x,4i8)') (lnods(inode,ie),inode=1,nnode)
	 enddo
c
c     
         do ik = 1, mgk
            dvtt2 = dvtt2 + dvmk(ik)
c            write(24,9201)xmk(1,ik),xmk(2,ik),dvmk(ik)
            write(24,9201)xmk(1,ik),xmk(2,ik)
         enddo
         print *,'total mass-zero',dvtt1,dvtt2
      endif
c
 9200    format(i6,2(1x,e12.5))
 9201    format(3(1x,e12.5))
 9202    format(2(1x,e12.5))
c
      close(24)
      close(25)
c
      close(ifwEcho)
      return      
      end
c
