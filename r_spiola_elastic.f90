subroutine r_spiola_elastic(det,xot,ge,iq,ne,cstr_element)
  use solid_variables, only: nsd_solid,nen_solid,ne_solid,nquadpad_solid
  use r_common
  implicit none 

  real(8),intent(in) :: det
  real(8) xot(nsd_solid,nsd_solid)
  real(8) :: cstr_element(nsd_solid*2),ge(nsd_solid*2,ne_solid,nquadpad_solid)  !...Cauchy stress
  
  integer :: isd,jsd,ksd,iq,ne
  real(8) :: cstr(nsd_solid,nsd_solid)
 
  cstr(:,:) = 0.0d0

  cstr_element(1:nsd_solid*2)=young_mod*ge(1:nsd_solid*2,ne,iq)

  threedim: if (nsd_solid .eq. 3) then 

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

  do isd=1,2
    cstr(isd,isd)=cstr_element(isd)
  enddo
  cstr(1,2)=cstr_element(3)
  cstr(2,1)=cstr_element(3)

  do isd = 1,nsd_solid
     do jsd = 1,nsd_solid
        PK1str_tens(isd,jsd) = 0.0d0
        do ksd = 1,nsd_solid
           PK1str_tens(isd,jsd) = PK1str_tens(isd,jsd) + det*xot(isd,ksd)*cstr(ksd,jsd)
        enddo
     enddo
  enddo

  endif twodim

end subroutine