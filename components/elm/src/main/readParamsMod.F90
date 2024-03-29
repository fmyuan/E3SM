module readParamsMod

  !-----------------------------------------------------------------------
  !
  ! Read parameters
  ! module used to read parameters for individual modules
  !
  use elm_varctl   , only: use_cn, use_century_decomp
  use elm_varctl   , only: use_lch4, use_fates
  implicit none
  save
  private
  !
  public :: readSharedParameters
  public :: readPrivateParameters
  
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readSharedParameters ()
    !
    implicit none
    !-----------------------------------------------------------------------

    call CNParamsSharedReadFile()

  end subroutine readSharedParameters

  !-----------------------------------------------------------------------

  subroutine CNParamsSharedReadFile ()
    !
    ! read CN and BGC shared parameters
    !

    use SharedParamsMod       , only : ParamsReadShared

    use elm_varctl              , only : paramfile, iulog
    use spmdMod                 , only : masterproc
    use fileutils               , only : getfil
    use ncdio_pio               , only : ncd_pio_closefile, ncd_pio_openfile, &
                                         file_desc_t, ncd_inqdid, ncd_inqdlen                                       
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !OTHER LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNParamsSharedReadFile'
    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'readParamsMod.F90::'//trim(subname)//' :: reading CN '//&
          'and BGC parameter file'
    end if

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid) 
    call ncd_inqdlen(ncid,dimid,npft) 

    !
    ! some parameters (eg. organic_max) are used in non-CN, non-BGC cases
    !
    call ParamsReadShared(ncid)


  end subroutine CNParamsSharedReadFile

  !-----------------------------------------------------------------------
  subroutine readPrivateParameters
    ! read CN and BGC shared parameters
    !
    use AllocationMod          , only : readCNAllocParams    
    use SoilLittDecompMod              , only : readSoilLittDecompParams
    use DecompCascadeBGCMod    , only : readDecompBGCParams
    use DecompCascadeCNMod     , only : readDecompCNParams
    use PhenologyMod             , only : readPhenolParams
    use MaintenanceRespMod               , only : readMaintenanceRespParams
    use NitrogenDynamicsMod           , only : readNitrogenDynamicsParams
    use GapMortalityMod          , only : readGapMortParams 
    use NitrifDenitrifMod      , only : readNitrifDenitrifParams
    use SoilLittVertTranspMod    , only : readSoilLittVertTranspParams
    use CH4Mod                   , only : readCH4Params
    use elm_varctl               , only : paramfile, iulog, use_hydrstress
    use spmdMod                  , only : masterproc
    use fileutils                , only : getfil
    use ncdio_pio                , only : ncd_pio_closefile, ncd_pio_openfile, &
                                          file_desc_t, ncd_inqdid, ncd_inqdlen
    use PhotosynthesisMod        , only : params_inst
    
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !OTHER LOCAL VARIABLES:
    character(len=32)  :: subname = 'readPrivateParameters'
    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'readParamsMod.F90::'//trim(subname)//' :: reading CN '//&
          'and BGC parameter file'
    end if

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid) 
    call ncd_inqdlen(ncid,dimid,npft)
    
    !
    ! populate each module with private parameters
    !       

    if ( (use_cn .or. use_fates) ) then

       call readCNAllocParams(ncid)
       !
         call readSoilLittDecompParams(ncid)
         if (use_century_decomp) then
            call readDecompBGCParams(ncid)
         else
            call readDecompCNParams(ncid)
         end if
       
         call readNitrifDenitrifParams(ncid)

         call readSoilLittVertTranspParams(ncid)
       
         if (use_lch4) then
            call readCH4Params (ncid)
         end if
      !
      call readNitrogenDynamicsParams (ncid)
    end if

    if (use_cn) then
         call readPhenolParams(ncid)
         call readMaintenanceRespParams (ncid)
         call readGapMortParams (ncid)
    end if
    !
    ! Biogeophysics
    !
    if(use_hydrstress)then
      call params_inst%readParams( ncid )
    endif

    !
    ! close CN params file
    !
    call ncd_pio_closefile(ncid)


 end subroutine readPrivateParameters
end module readParamsMod
