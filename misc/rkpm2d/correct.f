      subroutine correct(b,bd,bdd,cpt,cjp,dcjp,
     &                   dwjp,nep,iInter,
     &                   n_support,Lmap)
c
      implicit none
      include 'parameter.h'
c
      integer nep,iInter,n_support,Lmap(mnsch)
      real*8 b(*),bd(2,*),bdd(3,*),cpt(2)
      real*8 cjp(2,nep),dcjp(2,nep),dwjp(nep)
c
      if (iInter .eq. 1) then
	 call correct_La(b,bd,cpt,cjp,dcjp,dwjp,nep,
     &                   n_support,Lmap)
      elseif(iInter .eq. 11) then
	 call correct_bl(b,bd,bdd,cpt,cjp,dcjp,dwjp,nep,
     &                   n_support,Lmap)
      elseif(iInter .eq. 2)  then
	 call correct_qd(b,bd,bdd,cpt,cjp,dcjp,dwjp,nep,
     &                   n_support,Lmap)
      else
	 call correct_Lb(b,bd,bdd,cpt,cjp,dcjp,dwjp,nep,
     &                   n_support,Lmap)
      endif
c
c
      return
      end


