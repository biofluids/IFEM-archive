subroutine r_spiola_elastic(det,xot,ge,iq,ne,cstr_element)
  use solid_variables, only: nsd_solid,nen_solid,ne_solid,nquadpad_solid
  use r_common
  implicit none 

  real(8),intent(in) :: det
  real(8) xot(nsd_solid,nsd_solid)
  real(8) :: cstr_element(nsd_solid*2),ge(nsd_solid*2,ne_solid,nquadpad_solid)  !...Cauchy stress
  integer :: isd,jsd,ksd,iq,ne
  real(8) :: cstr(nsd_solid,nsd_solid)
  real(8) :: tempC, tempD, tempE
! Xingshi May 2008 
!======================================
  real(8) detf ! determninat of deformation gradient matrix
  real(8) xot_temp(nsd_solid,nsd_solid)
  detf=0.0d0
  xot_temp(:,:)=0.0d0
  xot_temp(:,:)=xot(:,:)
  call determinant(xot_temp,nsd_solid,nsd_solid,detf)
!=====================================
  cstr(:,:) = 0.0d0
  
  threedim: if (nsd_solid .eq. 3) then 

  tempC= young_mod*(1-Poisson)/((1+Poisson)*(1-2*Poisson))
  tempD= Poisson/(1-Poisson)
  tempE= (1-2*Poisson)/(2*(1-Poisson))
  
  cstr_element(1)= tempC*(ge(1,ne,iq)+tempD*ge(2,ne,iq)+tempD*ge(3,ne,iq))
  cstr_element(2)= tempC*(tempD*ge(1,ne,iq)+ge(2,ne,iq)+tempD*ge(3,ne,iq))
  cstr_element(3)= tempC*(tempD*ge(1,ne,iq)+tempD*ge(2,ne,iq)+ge(3,ne,iq))
  cstr_element(4)= tempC*tempE*ge(4,ne,iq)
  cstr_element(5)= tempC*tempE*ge(5,ne,iq)
  cstr_element(6)= tempC*tempE*ge(6,ne,iq)

  do isd=1,3
    cstr(isd,isd)=cstr_element(isd)
  enddo
  cstr(2,3)=cstr_element(4)
  cstr(3,2)=cstr_element(4)
  cstr(1,3)=cstr_element(5)
  cstr(3,1)=cstr_element(5)
  cstr(1,2)=cstr_element(6)
  cstr(2,1)=cstr_element(6)

  do isd = 1,3
     do jsd = 1,3
        PK1str_tens(isd,jsd) = 0.0d0
        do ksd = 1,3
           PK1str_tens(isd,jsd) = PK1str_tens(isd,jsd) + det*xot(isd,ksd)*cstr(ksd,jsd)
        enddo
     enddo
  enddo

  endif  threedim

  twodim: if (nsd_solid .eq. 2) then 
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Mickael Plane stress
	cstr_element(1)= (young_mod/(1-Poisson**2))*(ge(1,ne,iq)+Poisson*ge(2,ne,iq))
	cstr_element(2)= (young_mod/(1-Poisson**2))*(ge(2,ne,iq)+Poisson*ge(1,ne,iq))
    cstr_element(3)= (young_mod/(1+Poisson))*ge(3,ne,iq) 
    cstr_element(4)= 0.0 
! write(*,*) 'xot', xot(:,:)
  do isd=1,2
    cstr(isd,isd)=cstr_element(isd)
  enddo
  cstr(1,2)= cstr_element(3)
  cstr(2,1)= cstr_element(3)
 
! write(*,*)'stress', cstr(:,:)
	do isd = 1,nsd_solid
     do jsd = 1,nsd_solid
        PK1str_tens(isd,jsd) = 0.0d0
        do ksd = 1,nsd_solid
!           PK1str_tens(isd,jsd) = PK1str_tens(isd,jsd) + det*xot(isd,ksd)*cstr(ksd,jsd)
! Use detf instead of det which is the detminant of jacobi to inintial
           PK1str_tens(isd,jsd) = PK1str_tens(isd,jsd) + xot(isd,ksd)*cstr(ksd,jsd)/detf
	    enddo
     enddo
  enddo
! write(*,*) 'PK1str' , PK1str_tens(:,:)
! write(*,*) 'DET F', det
  endif twodim

end subroutine
