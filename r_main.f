      subroutine r_main
      implicit real*8 (a-h,o-z)
      include 'r_common'
      include 'main_common'      

      include 'ibd0_exchange_pars.fh'
      include 'ibg_variables_cfd.fh'
      include 'ibg_variables_cloud.fh'
      include 'ibg_variables_domain.fh'
      include 'ibg_variables_fluid.fh'
      include 'ibg_variables_marker.fh'
      include 'ibg_variables_point.fh'
      include 'ibg_variable_equivalences.fh'
      include 'iba_application_parameters.fh'
      include 'ibg_change_me_ptcon_var_decl.fh'
      include 'ibg_change_me_ptcon_var_common_equiv.fh'
      include 'ibg_parameters_run.fh' 
c
      x13=1.0d0/3.0d0
      x23=2.0d0/3.0d0
      x43=4.0d0/3.0d0
      x53=5.0d0/3.0d0
      x73=7.0d0/3.0d0
      x83=8.0d0/3.0d0
      x49=4.0d0/9.0d0
      x109=10.0d0/9.0d0
c
      do 1 i=1,nnda
         write(20,*) i,coor(i,1),coor(i,2)
 1    continue
c
      do 2 i=1,numela
         write(21,*) i,(nea(i,j),j=1,nis)
 2    continue
c
      close(20)
      close(21)
c
      call r_sinit
c
      do 10 i=1,nnd
         drf(i)=0.0d0
         drf(i+nnd)=0.0d0
         do 11 j=1,2
            dis(j,i)=0.0d0
            diso(j,i)=0.0d0
            prevel_pt(j,i)=0.0d0
            prevel2_pt(j,i)=0.0d0
            velt(j,i)=0.0d0
            acm(j,i)=0.0d0
            acm1(j,i)=0.0d0
 11      continue
 10   continue
      write(3,105) 
 105  format(1x,'time(1)=0.0d0;')
cccccccccccccccccccccccccccccccccccccc
      do 41 ip=1,npr
         xtt=0.0d0
         write(3,107) ip,xtt
 107     format(1x,'disy',i3,'(1)=',e23.7,';'/)
 41   continue
cccccccccccccccccccccccccccccccccccccc
      nai=0
      na=0
      ncop=nump*numel
c
      do 55 i=1,ncop
         prec(i)=0.0d0
         preco(i)=0.0d0
 55   continue
c
      return
      end

