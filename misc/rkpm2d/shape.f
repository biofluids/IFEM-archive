      subroutine RKPMshape(shp,shpd,shpdd,b,bd,bdd,
     &           cpt,cjt,dcjt,wjt,iInter)
c
       implicit none
c
      real*8  cpt(2),cjt(2),dcjt(2)
      real*8  b(*),bd(2,*),bdd(3,*)
      real*8  shpd(2),shpdd(3)
      real*8  shp,wjt
      integer iInter
c
c
      if (iInter .eq. 1 ) then
	 call RKPMshape_La(shp,shpd,b,bd,cpt,cjt,dcjt,wjt)
      elseif(iInter .eq. 11) then
	 call RKPMshape_bl(shp,shpd,shpdd,b,bd,bdd,cpt,cjt,dcjt,wjt)
      elseif(iInter .eq. 2)  then
	 call RKPMshape_qd(shp,shpd,shpdd,b,bd,bdd,cpt,cjt,dcjt,wjt)
      else
	 call RKPMshape_Lb(shp,shpd,shpdd,b,bd,bdd,cpt,cjt,dcjt,wjt)
      endif
c
      return
      end	!ends RKPMshape
c
