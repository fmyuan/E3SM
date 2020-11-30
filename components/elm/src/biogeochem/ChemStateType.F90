module ChemStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  
  implicit none
  save
  private
  !----------------------------------------------------
  ! column chemical state variables structure
  !----------------------------------------------------
  type, public :: chemstate_type

     real(r8), pointer :: soil_pH(:,:)    ! soil pH (-nlevsno+1:nlevgrnd)

     ! Data that must be saved for chemistry model (via alquimia)
     ! Sizes are set by alquimia
     ! State variables [col x layer]: water_density, porosity, temperature, aqueous_pressure
     ! [col x layer x num_primary]: total_mobile, total_immobile
     ! [col x layer x num_minerals]: mineral_volume_fraction, mineral_specific_surface_area
     ! [col x layer x num_surface_sites]: surface_site_density
     ! [col x layer x num_ion_exchange_sites]: cation_exchange_capacity
     ! [col x layer x num_aux_ints]: aux_ints
     ! [col x layer x num_aux_doubles]: aux_doubles
     ! Question: Is there a problem if these are not c doubles?
     real(r8), pointer :: water_density(:,:)
     !  real(r8), pointer :: porosity(:,:)    ! Redundant with soilstate_type%watsat_col
     !  real(r8), pointer :: temperature(:,:) ! Redundant with columnenergystate%t_soisno
     real(r8), pointer :: aqueous_pressure(:,:)

     real(r8), pointer :: total_mobile(:,:,:) 
     real(r8), pointer :: total_immobile(:,:,:) 
     real(r8), pointer :: mineral_volume_fraction(:,:,:) 
     real(r8), pointer :: mineral_specific_surface_area(:,:,:) 
     real(r8), pointer :: surface_site_density(:,:,:) 
     real(r8), pointer :: cation_exchange_capacity(:,:,:) 
     real(r8), pointer :: aux_doubles(:,:,:) 
     integer,  pointer :: aux_ints(:,:,:)
    
  contains
    procedure, public  :: Init         
    procedure, private :: InitAllocate   
    procedure, public  :: Restart
  end type chemstate_type
  
  contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    ! use ExternalModelInterfaceMod, only : EMI_Init_EM
    ! use ExternalModelConstants   , only : EM_ID_ALQUIMIA
    use clm_varctl               , only : use_alquimia

    implicit none

    class(chemstate_type)         :: this
    type(bounds_type), intent(in) :: bounds  

    ! Maybe it's better to initialize alquimia here?
    ! That way we guarantee that sizes are set right before we allocate the data
    ! But this happens before temperature and moisture are initialized so we should not transfer any of that data to alquimia at this point
    ! Also currently initialized as part of biogeophys which happens before decomp cascade is initialized. Could be a problem for pool mapping
    ! Problem: Calling EMI here introduces a circular dependency because EMI uses CLM_instmod and this Init is called from within CLM_instmod
    ! if (use_alquimia) then
    !   call EMI_Init_EM(EM_ID_ALQUIMIA)
    ! endif

    call this%InitAllocate(bounds)

  end subroutine Init
  
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use elm_varpar            , only : nlevdecomp_full
    ! Sizes were set by alquimia as part of initialization
    use elm_varctl            , only : use_alquimia
    use elm_varpar            , only : alquimia_num_primary, alquimia_num_minerals,&
                                       alquimia_num_surface_sites, alquimia_num_ion_exchange_sites, &
                                       alquimia_num_aux_doubles, alquimia_num_aux_ints
    !
    ! !ARGUMENTS:
    class(chemstate_type)            :: this
    type(bounds_type), intent(in)    :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------
    begc = bounds%begc;
    endc = bounds%endc
    lbj  = 1;
    ubj  = nlevdecomp_full
    
    allocate(this%soil_pH(begc:endc, lbj:ubj))

    ! Data for chemistry model (via alquimia)
    ! State variables [col x layer]: water_density, porosity*, temperature*, aqueous_pressure
    ! [col x layer x num_primary]: total_mobile, total_immobile
    ! [col x layer x num_minerals]: mineral_volume_fraction, mineral_specific_surface_area
    ! [col x layer x num_surface_sites]: surface_site_density
    ! [col x layer x num_ion_exchange_sites]: cation_exchange_capacity
    ! [col x layer x num_aux_ints]: aux_ints
    ! [col x layer x num_aux_doubles]: aux_doubles
    if(use_alquimia) then
      allocate(this%water_density(begc:endc,lbj:ubj))
      allocate(this%aqueous_pressure(begc:endc,lbj:ubj))

      allocate(this%total_mobile(begc:endc,lbj:ubj,1:alquimia_num_primary))
      allocate(this%total_immobile(begc:endc,lbj:ubj,1:alquimia_num_primary))
      allocate(this%mineral_volume_fraction(begc:endc,lbj:ubj,1:alquimia_num_minerals))
      allocate(this%mineral_specific_surface_area(begc:endc,lbj:ubj,1:alquimia_num_minerals))
      allocate(this%surface_site_density(begc:endc,lbj:ubj,1:alquimia_num_surface_sites))
      allocate(this%cation_exchange_capacity(begc:endc,lbj:ubj,1:alquimia_num_ion_exchange_sites))
      allocate(this%aux_ints(begc:endc,lbj:ubj,1:alquimia_num_aux_ints))
      allocate(this%aux_doubles(begc:endc,lbj:ubj,1:alquimia_num_aux_doubles))
    endif

  end subroutine InitAllocate


  subroutine Restart (this,  bounds, ncid, flag )

    use restUtilMod     , only : restartvar
    use ncdio_pio       , only : file_desc_t,ncd_double, ncd_int
    use clm_varpar      , only : alquimia_num_primary, alquimia_num_minerals,&
                                 alquimia_num_surface_sites, alquimia_num_ion_exchange_sites, &
                                 alquimia_num_aux_doubles, alquimia_num_aux_ints
    use clm_varctl               , only : use_alquimia
    use clm_varpar            , only : nlevdecomp_full

    implicit none
    !
    ! !ARGUMENTS:
    class (chemstate_type)     :: this
    type(bounds_type) , intent(in)     :: bounds 
    type(file_desc_t) , intent(inout)  :: ncid   ! netcdf id
    character(len=*)  , intent(in)     :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    integer :: ii
    character(len=256) :: nc_varname, var_longname, alq_poolname
    real(r8), pointer :: real2d(:,:)
    real(r8) , pointer :: int2d(:,:) ! Restart system doesn't actually support 2D integer arrays for some reason. Workaround is cast to real and back

    ! TODO: Check on read that number and order of variables is correct.
    !       - Make model fail if it fails to read an expected variable (i.e., not enough values stored)
    !       - At end of expected list, check if there is another in the netCDF file (too many variables stored)
    !       - See if long_name can be read from file and compared with expected long_name
    !       In either of these cases, restart does not match current reaction network spec and model should fail
    if(use_alquimia) then
      alq_poolname=''
      do ii=1,alquimia_num_primary
        !!! TOTAL_MOBILE !!!
        ! call c_f_string_ptr(name_list(ii),alq_poolname) ! Need to get metadata from alquimia somehow... EMI will not pass character data
        ! Generate field name as ALQUIMIA_MOBILE_01, ALQUIMIA_MOBILE_02, ...
        write(nc_varname,'(a,i2.2)') 'ALQUIMIA_MOBILE_',ii
        var_longname = 'Alquimia total mobile '//trim(alq_poolname)
        real2d => this%total_mobile(:,:,ii)

        call restartvar(ncid=ncid, flag=flag, varname=nc_varname, xtype=ncd_double,   &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name=var_longname, units='M', &
            interpinic_flag='interp', readvar=readvar, data=real2d)

        write(nc_varname,'(a,i2.2)') 'ALQUIMIA_IMMOBILE_',ii
        var_longname = 'Alquimia total immobile '//trim(alq_poolname)
        real2d => this%total_immobile(:,:,ii)

        call restartvar(ncid=ncid, flag=flag, varname=nc_varname, xtype=ncd_double,   &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name=var_longname, units='mol/m^3', &
            interpinic_flag='interp', readvar=readvar, data=real2d)

      enddo ! End of primary species loop

        
      do ii=1,alquimia_num_minerals
        !!! mineral_volume_fraction !!!
        ! Generate field name as ALQUIMIA_MINERAL_01, ALQUIMIA_MINERAL_02, ...
        write(nc_varname,'(a,i2.2)') 'ALQUIMIA_MINERAL_VF_',ii
        var_longname = 'Alquimia mineral volume fraction '//trim(alq_poolname)
        real2d => this%mineral_volume_fraction(:,:,ii)

        call restartvar(ncid=ncid, flag=flag, varname=nc_varname, xtype=ncd_double,   &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name=var_longname, units='[m^3 mineral/m^3 bulk]', &
            interpinic_flag='interp', readvar=readvar, data=real2d)

        !!! Mineral specific surface areas !!!
        ! Generate field name as ALQUIMIA_MINERAL_SSA_01, ALQUIMIA_MINERAL_SSA_02, ...
        write(nc_varname,'(a,i2.2)') 'ALQUIMIA_MINERAL_SSA_',ii
        var_longname = 'Alquimia mineral specific surface area '//trim(alq_poolname)
        real2d => this%mineral_specific_surface_area(:,:,ii)

        call restartvar(ncid=ncid, flag=flag, varname=nc_varname, xtype=ncd_double,   &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name=var_longname, units='[m^2 mineral/m^3 bulk]', &
            interpinic_flag='interp', readvar=readvar, data=real2d)

      enddo ! End of mineral species loop
      
      
      do ii=1,alquimia_num_surface_sites
        !!! surface site density !!!
        ! call c_f_string_ptr(name_list(ii),alq_poolname)
        ! Generate field name as ALQUIMIA_SURFACE_SITE_DENS_01, ALQUIMIA_SURFACE_SITE_DENS_02, ...
        write(nc_varname,'(a,i2.2)') 'ALQUIMIA_SURFACE_SITE_DENS_',ii
        var_longname = 'Alquimia surface site density '//trim(alq_poolname)
        real2d => this%surface_site_density(:,:,ii)

        call restartvar(ncid=ncid, flag=flag, varname=nc_varname, xtype=ncd_double,   &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name=var_longname, units='moles/m^3 bulk', &
            interpinic_flag='interp', readvar=readvar, data=real2d)

      enddo ! End of surface site densities


      do ii=1,alquimia_num_ion_exchange_sites
        !!! surface site density !!!
        ! call c_f_string_ptr(name_list(ii),alq_poolname)
        ! Generate field name as ALQUIMIA_CEC_01, ALQUIMIA_CEC_02, ...
        write(nc_varname,'(a,i2.2)') 'ALQUIMIA_CEC_',ii
        var_longname = 'Alquimia cation exchange capacity '//trim(alq_poolname)
        real2d => this%cation_exchange_capacity(:,:,ii)

        call restartvar(ncid=ncid, flag=flag, varname=nc_varname, xtype=ncd_double,   &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name=var_longname, units='moles/m^3 bulk', &
            interpinic_flag='interp', readvar=readvar, data=real2d)

      enddo ! End of ion exchange sites


      ! Aux doubles. These don't have metadata
      do ii=1,alquimia_num_aux_doubles
        write(nc_varname,'(a,i2.2)') 'ALQUIMIA_AUX_DOUBLE_',ii
        var_longname = ''
        real2d => this%aux_doubles(:,:,ii)

        call restartvar(ncid=ncid, flag=flag, varname=nc_varname, xtype=ncd_double,   &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name=var_longname, units='-', &
            interpinic_flag='interp', readvar=readvar, data=real2d)

      enddo ! End of aux doubles


      ! Aux integers. These don't have metadata
      ! Restart system only supports 1D ints so I am casting this to real
      allocate(int2d(bounds%begc:bounds%endc,1:nlevdecomp_full))
      do ii=1,alquimia_num_aux_ints
        write(nc_varname,'(a,i2.2)') 'ALQUIMIA_AUX_INT_',ii
        var_longname = ''
        if(flag == 'write') int2d = real(this%aux_ints(:,:,ii))

        call restartvar(ncid=ncid, flag=flag, varname=nc_varname, xtype=ncd_int,   &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name=var_longname, units='-', &
            interpinic_flag='interp', readvar=readvar, data=int2d)

        if(flag == 'read') this%aux_ints(:,:,ii) = int(int2d)

      enddo ! End of aux ints
      deallocate(int2d)

      call restartvar(ncid=ncid, flag=flag, varname='ALQUIMIA_WATER_DENSITY', xtype=ncd_double,   &
        dim1name='column', dim2name='levgrnd', switchdim=.true., &
        long_name='alquimia water density', units='kg/m^3', &
        interpinic_flag='interp', readvar=readvar, data=this%water_density)

      call restartvar(ncid=ncid, flag=flag, varname='ALQUIMIA_AQUEOUS_PRESSURE', xtype=ncd_double,   &
        dim1name='column', dim2name='levgrnd', switchdim=.true., &
        long_name='alquimia aqueous pressure', units='Pa', &
        interpinic_flag='interp', readvar=readvar, data=this%aqueous_pressure)
    endif

  end subroutine Restart
end module ChemStateType
