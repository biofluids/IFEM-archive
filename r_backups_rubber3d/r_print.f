      subroutine r_print
      implicit real*8 (a-h,o-z)
      include 'r_common'
      include 'main_common'
c
      write(3,105) (iti+1)/n_step_wr_ib_user_files,
     $     (iti)*dt_solid
ccccccccccccccccccccccccccccccc
c     print all nodes
ccccccccccccccccccccccccccccccc
      do ip=1,npr
c     displacement
         if (npdis .eq. 1) then
            xtt=dis(ndprint(ip),nprint(ip))+
     $           coor(nprint(ip),ndprint(ip))
         else           
            xtt=dis(ndprint(ip),nprint(ip))
         endif
         vtt=du(ndprint(ip),nprint(ip))
         write(3,104) ip,(iti+1)/n_step_wr_ib_user_files,xtt
 104     format(1x, 'disy',i3,'(',i4,')=',e23.7,';'/)
         write(3,115) ip,(iti+1)/n_step_wr_ib_user_files,vtt
 115     format(1x, 'vely',i3,'(',i4,')=',e23.7,';'/)
	enddo
 105  format(1x, 'time(',i4,')=',e23.7,';'/)
c     output displacement
      na=iti/nina-nai
      if (na .ne. 0) then
         do i=1,nnd
            write(22,*) i,dis(1,i),dis(2,i),dis(3,i)
	   enddo
ccccccccccccccccccccccccccccc
c     Print
ccccccccccccccccccccccccccccc
         do ni=1,ncop
            write(24,*) ni,prec(ni)
	   enddo
cccccccccccccccccccccccccccccc
c     print reaction forces
cccccccccccccccccccccccccccccc
         if (nreact .eq. 1) then
            do k=1,nrtp
               write(3,607) ndraf(k),nraf(k),iti+1,raf(nraf(k),ndraf(k))
 607           format(1x, 'rforce',i1,i3,'(',i4,')=',e23.7,';'/)
		  enddo
         endif
         write(22,*)
         nai=nai+1     
      endif
      return
      end
