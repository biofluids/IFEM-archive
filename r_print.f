      subroutine r_print
      implicit real*8 (a-h,o-z)
      include 'r_common'
      include 'main_common'
c
      write(3,105) (iti+1)/n_step_wr_ib_user_files,
     $     (iti)*dt_solid*unit_time
ccccccccccccccccccccccccccccccc
c     print all nodes
ccccccccccccccccccccccccccccccc
      do 1 ip=1,npr
c     displacement
         if (npdis .eq. 1) then
            xtt=unit_length*(dis(ndprint(ip),nprint(ip))+
     $           coor(nprint(ip),ndprint(ip)))
         else           
            xtt=dis(ndprint(ip),nprint(ip))*unit_length
         endif
         vtt=du(ndprint(ip),nprint(ip))*unit_velocity
         write(3,104) ip,(iti+1)/n_step_wr_ib_user_files,xtt
 104     format(1x, 'disy',i3,'(',i4,')=',e23.7,';'/)
         write(3,115) ip,(iti+1)/n_step_wr_ib_user_files,vtt
 115     format(1x, 'vely',i3,'(',i4,')=',e23.7,';'/)
 1    continue
 105  format(1x, 'time(',i4,')=',e23.7,';'/)
c     output displacement
      na=iti/nina-nai
      if (na .ne. 0) then
         do 4 i=1,nnd
            write(22,*) i,dis(1,i),dis(2,i)
 4       continue
ccccccccccccccccccccccccccccc
c     Print
ccccccccccccccccccccccccccccc
         do 33 ni=1,ncop
            write(24,*) ni,prec(ni)
 33      continue
cccccccccccccccccccccccccccccc
c     print reaction forces
cccccccccccccccccccccccccccccc
         if (nreact .eq. 1) then
            do 602 k=1,nrtp
               write(3,607) ndraf(k),nraf(k),iti+1,raf(nraf(k),ndraf(k))
 607           format(1x, 'rforce',i1,i3,'(',i4,')=',e23.7,';'/)
 602        continue
         endif
         write(22,*)
         nai=nai+1     
      endif
      return
      end
