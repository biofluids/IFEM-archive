      subroutine  interpssbpressure(! accumulate lattice press.->sourcesink/bal! sub !
     .     dilocaldomain,       ! sub !
     .     dmnac1, dmxac1,      ! sub !
     .     dmnac2, dmxac2,      ! sub !
     .     dmnac3, dmxac3,      ! sub !
      
     .     pressure_ss,         ! sub !
     .     wt_delss,            ! sub !
     .     ice1_delss, ice2_delss, ice3_delss, ! sub !
     .     pressure_bal,        ! sub !
     .     dnbc,                ! sub !
     .     dmnbce1, dmxbce1,    ! sub !
     .     dmnbce2, dmxbce2,    ! sub !
     .     dmnbce3, dmxbce3,    ! sub !
     .     pressure_fluid)      ! sub !

      implicit real*8 (a-h,o-z)

c----------------------------------------------------------------------
c     this       subroutine  is part of a distributed memory implementation
c
c  this       subroutine  computes the pressure at points in space by
c  taking weighted averages of the neighboring grid pressures.
c----------------------------------------------------------------------

c---------------------------------------------------------------------- ! doc !
c           subroutine  interpssbpressure
c
c     this       subroutine  is for a distributed memory implementation.
c
c  this routine accumulates the lattice weights associated with
c  a sourcesink of pressure *pressure_ss* at *coord_ss* from the
c  nearby planes of *cfdpressure*, from cells within the local domain.
c  the total weight accumulated is *pressure_ss*. the portion of
c  this *pressure_ss* accumulated from a particular plane k depends
c  on the z-distance between the sourcesink and the plane. if *coord_ss(iz)*
c  is sufficiently far from k, then the contribution from *pressure_fluid* is 0.0,
c  otherwise, the portion is determined by del(z-distance).
c  accumulation from each lattice on the plane is similarly determined
c  by the x-distance from the lattice to *coord_ss(ix)*
c  and the y-distance from the lattice to *coord_ss(iy)*.
c---------------------------------------------------------------------- ! doc !

c data ---------------

      include 'main_common'      
c      include 'ibd0_parameters.fh'
      include 'ibd0_nonzerocfddiv_pars.fh'

c input

      integer    dilocaldomain

      dimension  wt_delss(
     .          mn_wt_delss1:mx_wt_delss1,
     .          mn_wt_delss2:mx_wt_delss2,
     .          mn_wt_delss3:mx_wt_delss3,
     .          mn_ss_alloc:mx_ss_alloc )

      dimension  ice1_delss( mn_wt_delss1:mx_wt_delss1,
     .          mn_ss_alloc:mx_ss_alloc )
      dimension  ice2_delss( mn_wt_delss2:mx_wt_delss2,
     .          mn_ss_alloc:mx_ss_alloc )
      dimension  ice3_delss( mn_wt_delss3:mx_wt_delss3,
     .          mn_ss_alloc:mx_ss_alloc )


      integer    dnbc
      integer    dmnbce1, dmxbce1
      integer    dmnbce2, dmxbce2
      integer    dmnbce3, dmxbce3


c     fluid array - dimension using global coordinates
      integer    dmnac1, dmxac1
      integer    dmnac2, dmxac2
      integer    dmnac3, dmxac3
      dimension  pressure_fluid( dmnac1:dmxac1, dmnac2:dmxac2,
     .                                            dmnac3:dmxac3)

      dimension  pressure_ss( mn_ss_alloc:mx_ss_alloc )
c
c initialize pressure_ss = 0.0 for all sourcesinks
c
      do 1 iss = mn_ss_exp, mx_ss_exp
         pressure_ss(iss) = 0.0
    1 continue
 

c
c add the weighted lattice pressures to the existing values for source/sink
c
c     outside this neighborhood the values of *wt_delss* are zero, and
c     no contribution to pressure_ss is made by any values of pressure_fluid
c     outside this neighborhood.
c

      do 10 iss = mn_ss_exp, mx_ss_exp

      do 2 iwt3 = mn_wt_delss3, mx_wt_delss3                            ! 3-d !
                   ilce3 = ice3_delss(iwt3, iss)                        ! 3-d !
      do 2 iwt2 = mn_wt_delss2, mx_wt_delss2
                   ilce2 = ice2_delss(iwt2, iss)
      do 2 iwt1 = mn_wt_delss1, mx_wt_delss1
                   ilce1 = ice1_delss(iwt1, iss)


                pressure_ss(iss)
     .        = pressure_ss(iss)
     .      + (   wt_delss(      iwt1, iwt2, iwt3, iss)
     .          * pressure_fluid(ilce1,ilce2,ilce3)    )
    2 continue

   10 continue


c ------------------------------------------------
c initialize local sum of pressure_bal to 0.
c for all local lattice points in the balance planes (if any),
c   set pressure_bal = pressure_bal + pressure_fluid
c these planes balance the fluid created/destroyed by the sources and sinks
c ------------------------------------------------

      pressure_bal = 0.0



        do 3 ilce3 = dmnbce3, dmxbce3
        do 3 ilce2 = dmnbce2, dmxbce2
        do 3 ilce1 = dmnbce1, dmxbce1

           pressure_bal = pressure_bal +
     .                   wt_lc_bal * pressure_fluid(ilce1,ilce2,ilce3)

    3   continue







c
c 3. assign and subtract the reference pressure from the source/sink pressures
c

      refpressure = pressure_bal

      do 5 iss = mn_ss_exp, mx_ss_exp
      pressure_ss(iss) = pressure_ss(iss) - refpressure
    5 continue

      pressure_bal = pressure_bal - refpressure ! == 0. ! ????????????? ! ??? !


      return
      end                                                               ! end !
c======================================================================

