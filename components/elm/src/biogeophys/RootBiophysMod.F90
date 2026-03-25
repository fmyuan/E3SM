module RootBiophysMod

#include "shr_assert.h"
  !-------------------------------------------------------------------------------------- 
  ! DESCRIPTION:
  ! module contains subroutine for root biophysics
  !
  ! HISTORY
  ! created by Jinyun Tang, Mar 1st, 2014
  ! added variable DTB option for Zeng-Decker, Michael A. Brunke, Aug. 25, 2016
  implicit none
  private
  public :: init_vegrootfr
  public :: init_rootprof
  integer, parameter :: zeng_2001_root = 0 !the zeng 2001 root profile function
  integer, parameter :: spruce_root = 1    !the spruce linear root profile 

  integer :: root_prof_method              !select the type of root profile parameterization   
  !-------------------------------------------------------------------------------------- 

contains

  !-------------------------------------------------------------------------------------- 
  subroutine init_rootprof()
    !
    !DESCRIPTION
    ! initialize methods for root profile calculation
    implicit none

    root_prof_method = zeng_2001_root
    ! Keep the Zeng/Spruce mixed path active. init_vegrootfr applies the
    ! SPRUCE profile only to HUMHOL peatland topounit classes.

  end subroutine init_rootprof

  !-------------------------------------------------------------------------------------- 
  subroutine init_vegrootfr(bounds, nlevsoi, nlevgrnd, nlev2bed, rootfr)
    !
    !DESCRIPTION
    !initialize plant root profiles
    !
    ! USES
    use shr_kind_mod   , only : r8 => shr_kind_r8   
    use shr_assert_mod , only : shr_assert
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use decompMod      , only : bounds_type
    use abortutils     , only : endrun         
    use TopounitType   , only : top_pp
    use VegetationType , only : veg_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                     ! bounds
    integer,           intent(in) :: nlevsoi                    ! number of hydactive layers
    integer,           intent(in) :: nlevgrnd                   ! number of soil layers
    integer,           intent(in) :: nlev2bed(bounds%begc: )    ! number of layers to bedrock
    real(r8),          intent(out):: rootfr(bounds%begp: , 1: ) !
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'init_vegrootfr'  ! subroutine name
    real(r8) :: rootfr_zeng(bounds%begp:bounds%endp, 1:nlevsoi)
    real(r8) :: rootfr_spruce(bounds%begp:bounds%endp, 1:nlevsoi)
    integer  :: p, t
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(rootfr) == (/bounds%endp, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    rootfr(bounds%begp:bounds%endp, 1:nlevgrnd) = 0._r8

    select case (root_prof_method)
    case (zeng_2001_root)
       rootfr_zeng(bounds%begp:bounds%endp, 1:nlevsoi) = zeng2001_rootfr(bounds, nlevsoi, nlev2bed)
       rootfr_spruce(bounds%begp:bounds%endp, 1:nlevsoi) = spruce_rootfr(bounds, nlevsoi, nlev2bed)
       do p = bounds%begp,bounds%endp
          t = veg_pp%topounit(p)
          if (veg_pp%wtcol(p) <= 0._r8) cycle
          if (is_peatland_topounit(t)) then
             rootfr(p,1:nlevsoi) = rootfr_spruce(p,1:nlevsoi)
          else
             rootfr(p,1:nlevsoi) = rootfr_zeng(p,1:nlevsoi)
          endif
       enddo
    case (spruce_root)
       rootfr(bounds%begp:bounds%endp, 1 : nlevsoi) = spruce_rootfr(bounds,nlevsoi, nlev2bed)

       !case (jackson_1996_root)
       !jackson root, 1996, to be defined later
       !rootfr(bounds%begp:bounds%endp, 1 : ubj) = jackson1996_rootfr(bounds,
       !ubj, pcolumn, ivt, zi)
       !case (schenk_jackson_2002_root)
       !schenk and Jackson root, 2002, to be defined later
       !rootfr(bounds%begp:bounds%endp, 1 : ubj) = schenk2002_rootfr(bounds,
       !ubj, pcolumn, ivt, zi)        
    case default
       call endrun(subname // ':: a root fraction function must be specified!')
    end select

    call check_rootfr_nonnegative(bounds, nlevsoi, rootfr, subname)
  end subroutine init_vegrootfr

  !--------------------------------------------------------------------------------------
  subroutine check_rootfr_nonnegative(bounds, nlevsoi, rootfr, context)
    !
    ! DESCRIPTION
    ! Abort early if the initialized root profile contains negative layer
    ! fractions. Negative root fractions are not physically meaningful and
    ! can later appear as hydrology or carbon balance failures.
    !
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use decompMod      , only : bounds_type
    use abortutils     , only : endrun
    use elm_varctl     , only : iulog
    use elm_varcon     , only : namep
    use VegetationType , only : veg_pp
    use TopounitType   , only : top_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer,          intent(in) :: nlevsoi
    real(r8),         intent(in) :: rootfr(bounds%begp: , 1: )
    character(len=*), intent(in) :: context
    !
    ! !LOCAL VARIABLES:
    integer :: p, lev, t
    !------------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       if (veg_pp%wtcol(p) <= 0._r8) cycle
       do lev = 1, nlevsoi
          if (rootfr(p,lev) < 0._r8) then
             t = veg_pp%topounit(p)
             write(iulog,*) 'negative root fraction diagnostic'
             write(iulog,*) 'patch index               = ', p
             write(iulog,*) 'column index              = ', veg_pp%column(p)
             write(iulog,*) 'topounit index            = ', t
             write(iulog,*) 'topounit topo_grc_ind     = ', top_pp%topo_grc_ind(t)
             write(iulog,*) 'patch weight in column    = ', veg_pp%wtcol(p)
             write(iulog,*) 'pft type                  = ', veg_pp%itype(p)
             write(iulog,*) 'layer                     = ', lev
             write(iulog,*) 'root fraction             = ', rootfr(p,lev)
             write(iulog,*) 'root fraction profile     = ', rootfr(p,1:nlevsoi)
             write(iulog,*) 'topounit peat_depth       = ', top_pp%peat_depth(t)
             call endrun(decomp_index=p, elmlevel=namep, &
                  msg=trim(context)//':: negative root fraction')
          end if
       end do
    end do
  end subroutine check_rootfr_nonnegative

  !--------------------------------------------------------------------------------------
  logical function is_peatland_topounit(t)
    !
    ! DESCRIPTION
    ! The multi-topounit peatlands surface uses topo_grc_ind 1-3 for fen,
    ! bog hollow, and bog hummock, and 4 for upland. Root profile selection
    ! should follow topounit class, not peat depth, because peat depth can be
    ! zero or missing for an otherwise peatland topounit.
    !
    use elm_varctl  , only : use_humhol
    use TopounitType, only : top_pp
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: t
    !------------------------------------------------------------------------

    is_peatland_topounit = use_humhol .and. &
         top_pp%topo_grc_ind(t) >= 1 .and. &
         top_pp%topo_grc_ind(t) <= 3
  end function is_peatland_topounit

  !--------------------------------------------------------------------------------------   
  function zeng2001_rootfr(bounds, ubj, njbed) result(rootfr)
    !
    ! DESCRIPTION
    ! compute root profile for soil water uptake
    ! using equation from Zeng 2001, J. Hydrometeorology
    !
    ! USES
    use shr_kind_mod   , only : r8 => shr_kind_r8   
    use shr_assert_mod , only : shr_assert
    use shr_log_mod    , only : errMsg => shr_log_errMsg   
    use decompMod      , only : bounds_type
    use pftvarcon      , only : noveg, roota_par, rootb_par  !these pars shall be moved to here and set as private in the future
    use elm_varctl     , only : use_var_soil_thick
    use VegetationType , only : veg_pp
    use ColumnType     , only : col_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                  ! bounds
    integer           , intent(in)    :: ubj                     ! ubnd
    integer           , intent(in)    :: njbed(bounds%begc: )    ! nlev2bed
    !
    ! !RESULT
    real(r8) :: rootfr(bounds%begp:bounds%endp , 1:ubj ) !
    !
    ! !LOCAL VARIABLES:
    integer :: p, lev, c, nlevbed
    real    :: totrootfr
    !------------------------------------------------------------------------

    !(computing from surface, d is depth in meter):
    ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
    ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
    ! beta & d_obs given in Zeng et al. (1998).   

    rootfr(bounds%begp:bounds%endp, 1:ubj) = 0._r8

    do p = bounds%begp,bounds%endp   

       if (veg_pp%itype(p) /= noveg .and. .not.veg_pp%is_fates(p)) then
          c = veg_pp%column(p)
	  nlevbed = njbed(c)
	  totrootfr = 0._r8
          do lev = 1, ubj-1
             rootfr(p,lev) = .5_r8*( exp(-roota_par(veg_pp%itype(p)) * col_pp%zi(c,lev-1))  &
                  + exp(-rootb_par(veg_pp%itype(p)) * col_pp%zi(c,lev-1))  &
                  - exp(-roota_par(veg_pp%itype(p)) * col_pp%zi(c,lev  ))  &
                  - exp(-rootb_par(veg_pp%itype(p)) * col_pp%zi(c,lev  )) )
	     if(lev <= nlevbed) then
                totrootfr = totrootfr + rootfr(p,lev)
	     end if
          end do
          rootfr(p,ubj) = .5_r8*( exp(-roota_par(veg_pp%itype(p)) * col_pp%zi(c,ubj-1))  &
               + exp(-rootb_par(veg_pp%itype(p)) * col_pp%zi(c,ubj-1)) )
          totrootfr = totrootfr + rootfr(p, ubj)

          ! Adjust layer root fractions if nlev2bed < nlevsoi
          if (use_var_soil_thick .and. nlevbed < ubj) then
             do lev = 1, nlevbed
                rootfr(p,lev) = rootfr(p,lev) / totrootfr
             end do
             rootfr(p,nlevbed+1:ubj) = 0.0_r8
          endif
       else
          rootfr(p,1:ubj) = 0._r8
       endif

    enddo
    return

  end function zeng2001_rootfr

  function spruce_rootfr(bounds, ubj, njbed) result(rootfr)
    !
    ! DESCRIPTION
    ! compute root profile for soil water uptake
    ! using equation from spruce field observation 
    !
    ! USES
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_assert_mod , only : shr_assert
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use decompMod      , only : bounds_type
    use pftvarcon      , only : noveg, roota_par, rootb_par, graminoid
    use elm_varctl     , only : use_var_soil_thick
    use VegetationType , only : veg_pp
    use ColumnType     , only : col_pp
    use TopounitType   , only : top_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                  ! bounds
    integer           , intent(in)    :: ubj                     ! ubnd
    integer           , intent(in)    :: njbed(bounds%begc: )    ! nlev2bed
    !
    ! !RESULT
    real(r8) :: rootfr(bounds%begp:bounds%endp , 1:ubj ) !
    !
    ! !LOCAL VARIABLES:
    integer, parameter :: peatland_evergreen_needleleaf_boreal_pft = 3
    integer, parameter :: broadleaf_deciduous_boreal_tree_pft = 10
    integer :: p, lev, c, nlevbed, t
    real    :: totrootfr
    real(r8) :: cumdist,cumdist_last,roota_slope,rootb_intercept
    logical :: is_peat_topounit
    !------------------------------------------------------------------------

    !(computing from surface, d is peat depth  in centimeter):
    ! Y = root_slope*d+root_intercept  
    ! 

    rootfr(bounds%begp:bounds%endp, 1:ubj) = 0._r8

    do p = bounds%begp,bounds%endp
       if (veg_pp%itype(p) /= noveg .and. .not.veg_pp%is_fates(p)) then
          c = veg_pp%column(p)
          t = veg_pp%topounit(p)
          is_peat_topounit = is_peatland_topounit(t)
          ! For this profile, roota_par/rootb_par store the observed
          ! SPRUCE linear slope/intercept for peatland-tuned PFTs.
          rootb_intercept = rootb_par(veg_pp%itype(p))
          roota_slope     = roota_par(veg_pp%itype(p))
          if (is_peat_topounit) then
             if (veg_pp%itype(p) == broadleaf_deciduous_boreal_tree_pft) then
                rootb_intercept = rootb_par(peatland_evergreen_needleleaf_boreal_pft)
                roota_slope     = roota_par(peatland_evergreen_needleleaf_boreal_pft)
             else if (graminoid(veg_pp%itype(p)) == 1._r8) then
                ! Graminoid default roota/rootb values are Zeng-profile
                ! parameters. Under the SPRUCE linear profile, use a 100 cm
                ! zero-intercept cumulative rooting depth instead.
                rootb_intercept = 0._r8
                roota_slope     = 1._r8
             end if
          end if
          nlevbed = njbed(c)
          totrootfr = 0._r8
          cumdist_last = rootb_intercept
          do lev = 1, ubj
             cumdist = min(cumdist_last + col_pp%dz(c,lev) * roota_slope, 1.0_r8)
             if (cumdist > 0._r8) then
               rootfr(p,lev) = max(cumdist, 0._r8) - max(cumdist_last, 0._r8)
             endif
             if (lev .eq. 1 .and. rootb_intercept .gt. 0) then
               rootfr(p,lev) = rootfr(p,lev) + rootb_intercept
             end if
             cumdist_last = cumdist
             if(lev <= nlevbed) then
                totrootfr = totrootfr + rootfr(p,lev)
             end if
          end do

          ! Adjust layer root fractions if nlev2bed < nlevsoi
          if (use_var_soil_thick .and. nlevbed < ubj) then
             do lev = 1, nlevbed
                rootfr(p,lev) = rootfr(p,lev) / totrootfr
             end do
             rootfr(p,nlevbed+1:ubj) = 0.0_r8
          endif
       else
          rootfr(p,1:ubj) = 0._r8
       endif
    enddo
    return

  end function spruce_rootfr

end module RootBiophysMod
