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
    integer  :: p, t
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(rootfr) == (/bounds%endp, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    rootfr(bounds%begp:bounds%endp, 1:nlevgrnd) = 0._r8

    select case (root_prof_method)
    case (zeng_2001_root)
       rootfr_zeng(bounds%begp:bounds%endp, 1:nlevsoi) = zeng2001_rootfr(bounds, nlevsoi, nlev2bed)
       do p = bounds%begp,bounds%endp
          t = veg_pp%topounit(p)
          if (veg_pp%wtcol(p) <= 0._r8) cycle
          rootfr(p,1:nlevsoi) = rootfr_zeng(p,1:nlevsoi)
       enddo
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
             call endrun(decomp_index=p, elmlevel=namep, &
                  msg=trim(context)//':: negative root fraction')
          end if
       end do
    end do
  end subroutine check_rootfr_nonnegative

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

end module RootBiophysMod
