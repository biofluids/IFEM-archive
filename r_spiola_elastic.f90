subroutine r_spiola_elastic(det,xot,ge,iq,ne,cstr_element)
  use solid_variables, only: nen_solid,ne_solid,nquadpad_solid
  use r_common
  implicit none 

  real(8),intent(in) :: det
  real(8) xot(3,3)
  real(8) :: cstr_element(6),ge(6,ne_solid,nquadpad_solid)  !...Cauchy stress
  
  integer :: isd,jsd,ksd,iq,ne
  real(8) :: cstr(3,3)
 
  cstr(:,:) = 0.0d0

  cstr_element(1:6)=young_mod*ge(1:6,ne,iq)
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

end subroutine