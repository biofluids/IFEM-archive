      subroutine write2d(xm,xs,disp,vel,acc,
     &           disp_node,vel_node,acc_node,           
     &           lnods,nelem,numnp,mgk,  
     &           istep,iOutput,fhead)
c
c------------------------------------------------------------
c
c    The output subroutine for 2-D code
c
c
c    do weighted average outside the subroutine
c
c   September, 1998
c
c-----------------------------------------------------------
c
       implicit none
       include 'parameter.h'
c
       real*8  xm(2,maxNumnp),xs(2,maxNumnp)
       real*8  disp(2,maxNumnp),vel(2,maxNumnp),acc(2,maxNumnp)
       real*8  disp_node(2,maxNumnp),
     &         vel_node(2,maxNumnp),
     &         acc_node(2,maxNumnp)
       real*8  stsnpt(4,maxNumnp),
     &         effEnpt(maxNumnp),effSnpt(maxNumnp),
     &         effTemp(maxNumnp),effRate(maxNumnp),
     &         effMax(maxNumnp)
       integer nelem,numnp,mgk
       integer lnods(maxNode,maxElem)
c
       character*60 fhead,fwname,tmpstr1,tmpstr2,tmpstr3
       character*26 suffix
       character*1  single
c
       integer istep,iOutput,
     &         ifwout,ist,
     &         ie,ip,ipt,inode,ists,j
c
       real*8  dt,time
       integer nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp,
     &         iMaterType,nnode,nintElem,iLumping,
     &         imethRKPM,imethFEM,imeth,iIntMeth,iInter,
     &         iGL,iAdapt,numadp,iad_stage,Istage
c
       common /stepT/dt,time
       common /step/nstep,nOutputFrq,nPlotFrq,iplot,iplot_gp
       common /ctrl/iMaterType,nnode,nintElem,iLumping
       common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter
       common /adapt/iGL,iAdapt,numadp,iad_stage,Istage
       common /output/stsnpt,effEnpt,effSnpt,effTemp,effRate,effMax
c
c
       suffix = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
c
c
c           Output the intermid result
c
       ifwout = 51
c
c.......! output to files every nOutputFrq steps or at the end step
c
	  iOutput = iOutput + 1
	  write(*,*) 'iOutput=',iOutput
	  write(*,3384)istep,dt,time
cc##
	  ist    = Istage + 1 
	  single = suffix(ist:ist)
cc#
	 if ( 1 .eq. 1) then
	  call Int2Str(iOutput,tmpstr1)
          call appstr(single,tmpstr1,tmpstr2)
          call appstr('.gpsts',tmpstr2,tmpstr3)
          call appstr(fhead,tmpstr3,fwname)
	  write(*,*) 'fwname=',fwname
	  open(ifwout,file=fwname)
	    if (iMaterType .eq. 1) then
                write(ifwout,'(1x,a)') ' VARIABLES =
     &          "X","Y","S1"'
            elseif(iMaterType .eq. 2) then
                write(ifwout,'(1x,a)') ' VARIABLES =
     &          "X","Y","EE","SS","S1","S2","S4"'
            elseif(iMaterType .eq. 3) then
                write(ifwout,'(1x,a)') ' VARIABLES =
     &          "X","Y","EE","SS","S1","S2","S4"'
            elseif(iMaterType .eq. 4) then
                write(ifwout,'(1x,a)') ' VARIABLES =
     &          "X","Y","EE","SS","TT","S1"'
            elseif(iMaterType .eq. 5) then
                write(ifwout,'(1x,a)') ' VARIABLES =
     &          "X","Y","SS","EE","TT","S1","S12","RR"'
            else
                write(ifwout,'(1x,a)') ' VARIABLES =
     &          "X","Y","E1","S1"'
	    endif
c
            if (nnode.eq.3) then
	        write(ifwout,'(1x,a,i6,a,i6,a)')
     &         'ZONE N =',numnp,', E =',nelem,
     &         ', F = FEPOINT, ET = TRIANGLE'
            elseif (nnode.eq.4) then
	        write(ifwout,'(1x,a,i6,a,i6,a)')
     &         'ZONE N =',numnp,', E =',nelem,
     &         ', F = FEPOINT, ET = QUADRILATERAL'
            endif
c
	  do ipt = 1, numnp
	     if (iMaterType .eq. 1) then
                 write(ifwout,'(1x,3e15.7)') 
     &           xs(1,ipt),xs(2,ipt),
     &           effSnpt(ipt)
             elseif(iMaterType .eq. 2) then
                 write(ifwout,'(1x,7e15.7)') 
     &           xs(1,ipt),xs(2,ipt),
     &           effEnpt(ipt),effSnpt(ipt),
     &           stsnpt(1,ipt),stsnpt(2,ipt),
     &           stsnpt(4,ipt)
             elseif(iMaterType .eq. 3) then
                 write(ifwout,'(1x,7e15.7)') 
     &           xs(1,ipt),xs(2,ipt),
     &           effEnpt(ipt),effSnpt(ipt),
     &           stsnpt(1,ipt),stsnpt(2,ipt),
     &           stsnpt(4,ipt)
             elseif(iMaterType .eq. 4) then
                 write(ifwout,'(1x,6e15.7)') 
     &           xs(1,ipt),xs(2,ipt),
     &           effEnpt(ipt),effSnpt(ipt),
     &           effTemp(ipt),stsnpt(1,ipt)
             elseif(iMaterType .eq. 5) then
                 write(ifwout,'(1x,8e15.7)') 
     &           xs(1,ipt),xs(2,ipt),
     &           effSnpt(ipt),effEnpt(ipt),
     &           effTemp(ipt),stsnpt(1,ipt),
     &           stsnpt(4,ipt),effRate(ipt)
	     else
                 write(ifwout,'(1x,4e15.7)') 
     &           xs(1,ipt),xs(2,ipt),
     &           effEnpt(ipt),effSnpt(ipt)
	     endif
c
c...................................................
c
          enddo
c
	    if (nnode .eq. 4) then
               do ie=1,nelem
	           write(ifwout,'(1x,4i8)')
     &             (lnods(inode,ie),inode=1,nnode)
               enddo
	    elseif(nnode .eq. 8) then
               do ie=1,nelem
	           write(ifwout,'(1x,8i8)')
     &             (lnods(inode,ie),inode=1,nnode)
               enddo
	    endif
cc##
	  close(ifwout)
          endif
cc##
	  if (1 .eq. 0) then   ! output deformed shape
	     call Int2Str(iOutput,tmpstr1)
	     call appstr('.newxy',tmpstr1,tmpstr2)
	     call appstr(fhead,tmpstr2,fwname)
	     write(*,*) 'fwname=',fwname
	     open(ifwout,file=fwname)
	     do ip = 1, numnp
		write(ifwout,'(1x,2e15.7)') xs(1,ip),xs(2,ip)
             enddo
	     close(ifwout)
           endif
c
c...........Output the deformed shape in TECPLOT format.........
c
           if ( 1 .eq. 0) then
            call Int2Str(iOutput,tmpstr1)
            call appstr(single,tmpstr1,tmpstr2)
            call appstr('.tecplt',tmpstr2,tmpstr3)
            call appstr(fhead,tmpstr3,fwname)
            write(*,*) 'fwname=',fwname
            open(ifwout,file=fwname)
            write(ifwout,'(1x,a)')  'VARIABLES="X","Y"'
c
            if (nnode.eq.3) then
	        write(ifwout,'(1x,a,i6,a,i6,a)')
     &         'ZONE N =',numnp,', E =',nelem,
     &         ', F = FEPOINT, ET = TRIANGLE'
            elseif (nnode.eq.4) then
	        write(ifwout,'(1x,a,i6,a,i6,a)')
     &         'ZONE N =',numnp,', E =',nelem,
     &         ', F = FEPOINT, ET = QUADRILATERAL'
            endif
c
            do ip = 1,numnp
               write(ifwout,'(1x,2e15.6)')xs(1,ip),xs(2,ip)
            enddo
c
	    if (nnode .eq. 3) then
               do ie=1,nelem
	           write(ifwout,'(1x,3i8)')
     &             (lnods(inode,ie),inode=1,nnode)
               enddo
	    elseif(nnode .eq. 4) then
               do ie=1,nelem
	           write(ifwout,'(1x,4i8)')
     &             (lnods(inode,ie),inode=1,nnode)
               enddo
	    endif
	    close(ifwout)
        endif
cc##
cc##   
	  !Output the displacement of each nodes
        if (1.eq.0) then
	   call Int2Str(iOutput,tmpstr1)
	   call appstr('.disp',tmpstr1,tmpstr2)
	   call appstr(fhead,tmpstr2,fwname)
	   write(*,*) 'fwname=',fwname
	   open(ifwout,file=fwname)
c
	   do ip = 1, numnp
	      write(ifwout,'(1x,2e12.5)')
     &             (xs(j,ip)-xm(j,ip), j=1,2)
           enddo
	   close(ifwout)
        endif
c
c         Output the velocity of each nodes
c
	if (1 .eq. 0) then
	   call Int2Str(iOutput,tmpstr1)
	   call appstr('.vel',tmpstr1,tmpstr2)
	   call appstr(fhead,tmpstr2,fwname)
	   write(*,*) 'fwname=',fwname
	   open(ifwout,file=fwname)
c
	   do ip = 1, numnp
	      write(ifwout,'(1x,2e12.5)') (vel_node(j,ip), j=1,2)
           enddo
	   close(ifwout)
        endif
c
c      Output the acceleration of each nodes
c
        if (1.eq.0 ) then
	   call Int2Str(iOutput,tmpstr1)
	   call appstr('.acc',tmpstr1,tmpstr2)
	   call appstr(fhead,tmpstr2,fwname)
	   write(*,*) 'fwname=',fwname
	   open(ifwout,file=fwname)
c
	   do ip = 1, numnp
	      write(ifwout,'(1x,2e12.5)') (acc_node(j,ip), j=1,2)
           enddo
	   close(ifwout)
        endif
c
 3384 format(i7,3(1x,12e12.5))
c
      return
      end
c
c
