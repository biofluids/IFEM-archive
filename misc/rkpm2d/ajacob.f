      subroutine ajacob4(ajj,x1,x2,x3,x4,y1,y2,y3,y4,
     &                  xai,eta)

c*** cal the determinant of J matrix for 4-node element
c                          
c                       dx         dx
c                     ------      -----
c                      dxsi        deta
c            det (                         )                  
c                       dy         dy
c                     ------      -----
c                      dxsi        deta


      implicit double precision (a-h,o-z)

      aj11=(1.-eta)*(x2-x1)/4.+(1.+eta)*(x3-x4)/4.
      aj12=(1.-xai)*(x4-x1)/4.+(1.+xai)*(x3-x2)/4.
      aj21=(1.-eta)*(y2-y1)/4.+(1.+eta)*(y3-y4)/4.
      aj22=(1.-xai)*(y4-y1)/4.+(1.+xai)*(y3-y2)/4.
      ajj = aj11*aj22-aj12*aj21

      end	!ends ajacob4
      
      
      subroutine ajacob3(ajj,x1,x2,x3,y1,y2,y3,
     &                  xai,eta)

c*** cal the determinant of J matrix for 3-node element
c                          
c                       dx         dx
c                     ------      -----
c                      dxsi        deta
c            det (                         )                  
c                       dy         dy
c                     ------      -----
c                      dxsi        deta


      implicit double precision (a-h,o-z)

      aj11=x1-x3
      aj12=x2-x3
      aj21=y1-y3
      aj22=y2-y3
      ajj = aj11*aj22-aj12*aj21

      end	!ends ajacob3
