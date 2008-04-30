!      subroutine shx_tets(x,xe,nsd,nen,sh,shg)
      subroutine shx_tets(x,xe,nsd,nen,sh)
	integer nen
      integer nsd
      real(8) x(nsd)
	real(8) xe(nsd,nen) ! xyz coordinates of the element nodals
	real(8) sh(nen)
	real(8) shg(nsd,nen)

	real(8) x1, x2, x3, x4
	real(8) y1, y2, y3, y4
	real(8) z1, z2, z3, z4
	real(8) xx, yy, zz
	real(8) norm

	x1=xe(1,1)
	x2=xe(1,2)
	x3=xe(1,3)
	x4=xe(1,4)

	
	y1=xe(2,1)
	y2=xe(2,2)
	y3=xe(2,3)
	y4=xe(2,4)

	
	z1=xe(3,1)
	z2=xe(3,2)
	z3=xe(3,3)
	z4=xe(3,4)

	xx=x(1)
	yy=x(2)
	zz=x(3)

	norm=-x2*y3*z4+x2*z3*y4+x3*y2*z4-x3*z2*y4-x4*y2*z3+x4*z2*y3 &
          +x1*y3*z4-x1*z3*y4-x3*y1*z4+x3*z1*y4+x4*y1*z3-x4*z1*y3 &
     	 -x1*y2*z4+x1*z2*y4+x2*y1*z4-x2*z1*y4-x4*y1*z2+x4*z1*y2 &
     	 +x1*y2*z3-x1*z2*y3-x2*y1*z3+x2*z1*y3+x3*y1*z2-x3*z1*y2

      norm=1.0/norm
!==========================================================================
      sh(1)=-(x2*y3*z4-x2*z3*y4-x3*y2*z4+x3*z2*y4+x4*y2*z3-x4*z2*y3) &
           *norm 
	sh(1)=sh(1)+(y3*z4-z3*y4-y2*z4+z2*y4+y2*z3-z2*y3)*xx*norm
	sh(1)=sh(1)-(x3*z4-z3*x4-x2*z4+z2*x4+x2*z3-z2*x3)*yy*norm
	sh(1)=sh(1)+(x3*y4-y3*x4-x2*y4+y2*x4+x2*y3-y2*x3)*zz*norm
!===========================================================================
	sh(2)=(x1*y3*z4-x1*z3*y4-x3*y1*z4+x3*z1*y4+x4*y1*z3-x4*z1*y3)*norm
	sh(2)=sh(2)-xx*(y3*z4-z3*y4-y1*z4+z1*y4+y1*z3-z1*y3)*norm
	sh(2)=sh(2)+yy*(x3*z4-z3*x4-x1*z4+z1*x4+x1*z3-z1*x3)*norm
	sh(2)=sh(2)-zz*(x3*y4-y3*x4-x1*y4+y1*x4+x1*y3-y1*x3)*norm
!===========================================================================
      sh(3)=-(x1*y2*z4-x1*z2*y4-x2*y1*z4+x2*z1*y4+x4*y1*z2-x4*z1*y2) &
     	  *norm
      sh(3)=sh(3)+xx*(y2*z4-z2*y4-y1*z4+z1*y4+y1*z2-z1*y2)*norm
	sh(3)=sh(3)-yy*(x2*z4-z2*x4-x1*z4+z1*x4+x1*z2-z1*x2)*norm
	sh(3)=sh(3)+zz*(x2*y4-y2*x4-x1*y4+y1*x4+x1*y2-y1*x2)*norm
!============================================================================
      sh(4)=(x1*y2*z3-x1*z2*y3-x2*y1*z3+x2*z1*y3+x3*y1*z2-x3*z1*y2)*norm
	sh(4)=sh(4)-xx*(y2*z3-z2*y3-y1*z3+z1*y3+y1*z2-z1*y2)*norm
	sh(4)=sh(4)+yy*(x2*z3-z2*x3-x1*z3+z1*x3+x1*z2-z1*x2)*norm
	sh(4)=sh(4)-zz*(x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2)*norm
!============================================================================
!N_{,x_1}      
	shg(1,1)=(y3*z4-z3*y4-y2*z4+z2*y4+y2*z3-z2*y3)*norm
	shg(1,2)=-(y3*z4-z3*y4-y1*z4+z1*y4+y1*z3-z1*y3)*norm
	shg(1,3)=(y2*z4-z2*y4-y1*z4+z1*y4+y1*z2-z1*y2)*norm
	shg(1,4)=-(y2*z3-z2*y3-y1*z3+z1*y3+y1*z2-z1*y2)*norm
!N_{,x_2}
      shg(2,1)=-(x3*z4-z3*x4-x2*z4+z2*x4+x2*z3-z2*x3)*norm
      shg(2,2)=(x3*z4-z3*x4-x1*z4+z1*x4+x1*z3-z1*x3)*norm
	shg(2,3)=-(x2*z4-z2*x4-x1*z4+z1*x4+x1*z2-z1*x2)*norm      
      shg(2,4)=(x2*z3-z2*x3-x1*z3+z1*x3+x1*z2-z1*x2)*norm
!N_{,x_3}
      shg(3,1)=(x3*y4-y3*x4-x2*y4+y2*x4+x2*y3-y2*x3)*norm
      shg(3,2)=-(x3*y4-y3*x4-x1*y4+y1*x4+x1*y3-y1*x3)*norm
	shg(3,3)=(x2*y4-y2*x4-x1*y4+y1*x4+x1*y2-y1*x2)*norm
	shg(3,4)=-(x2*y3-y2*x3-x1*y3+y1*x3+x1*y2-y1*x2)*norm	
	
	return
	end

