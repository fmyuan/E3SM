module ncdio_nf90Mod

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: ncdio_nf90Mod
  !
  ! !DESCRIPTION:
  ! Generic interfaces to read/write fields to netcdf files for ELM
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8, i4=>shr_kind_i4
  use shr_infnan_mod , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_sys_mod    , only : shr_sys_abort
  use elm_varcon     , only : spval,ispval
  use elm_varcon     , only : grlnd, nameg, namet, namel, namec, namep, nameCohort
  use elm_varctl     , only : iulog
  use spmdMod        , only : masterproc, iam
  use netcdf
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  public :: check_var          ! determine if variable is on netcdf file
  public :: check_att          ! check if attribute is on file
  public :: check_dim          ! validity check on dimension
  public :: ncd_pio_openfile   ! open a file (dummy module name for convienence)
  public :: ncd_pio_createfile ! create a new file (dummy module name for convienence)
  public :: ncd_pio_closefile  ! close a file (dummy module name for convienence)
  !public :: ncd_nf90_init     ! called from elm_comp (TODO: needed?)
  public :: ncd_putatt         ! put attribute
  public :: ncd_getatt         ! get attribute
  public :: ncd_defdim         ! define dimension
  public :: ncd_defvar         ! define variables
  public :: ncd_enddef         ! end define mode
  public :: ncd_inqdid         ! inquire dimension id
  public :: ncd_inqdname       ! inquire dimension name
  public :: ncd_inqdlen        ! inquire dimension length
  public :: ncd_inqfdims       ! inquire file dimnesions 
  public :: ncd_inqvid         ! inquire variable id
  public :: ncd_inqvname       ! inquire variable name
  public :: ncd_inqvdims       ! inquire variable ndims
  public :: ncd_inqvdids       ! inquire variable dimids
  public :: ncd_inqvdlen       ! inquire variable dimension size
  public :: ncd_inqvdname      ! inquire variable dimension name
  public :: ncd_io             ! write local data

  integer,parameter,public :: ncd_int       = nf90_int
  integer,parameter,public :: ncd_log       =-nf90_int
  integer,parameter,public :: ncd_float     = nf90_real
  integer,parameter,public :: ncd_double    = nf90_double
  integer,parameter,public :: ncd_char      = nf90_char
  integer,parameter,public :: ncd_global    = nf90_global
  integer,parameter,public :: ncd_write     = nf90_write
  integer,parameter,public :: ncd_nowrite   = nf90_nowrite
  integer,parameter,public :: ncd_clobber   = nf90_clobber
  integer,parameter,public :: ncd_noclobber = nf90_noclobber
  integer,parameter,public :: ncd_nofill    = nf90_nofill
  integer,parameter,public :: ncd_unlimited = nf90_unlimited

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  !
  interface ncd_pio_openfile
     module procedure ncd_nf90_openfile
  end interface

  interface ncd_pio_closefile
     module procedure ncd_nf90_closefile
  end interface

  interface ncd_pio_createfile
     module procedure ncd_nf90_createfile
  end interface

  interface ncd_defvar
     module procedure ncd_defvar_bynf
     module procedure ncd_defvar_bygrid
  end interface

  interface ncd_putatt
     module procedure ncd_putatt_int
     module procedure ncd_putatt_real
     module procedure ncd_putatt_char
  end interface

  interface ncd_getatt
     module procedure ncd_getatt_char
  end interface ncd_getatt

  interface ncd_io 
     module procedure ncd_io_char_var0_start

     module procedure ncd_io_0d_log
     module procedure ncd_io_1d_log
     
     module procedure ncd_io_0d_text
     module procedure ncd_io_1d_text
     module procedure ncd_io_2d_text

     module procedure ncd_io_0d_int
     module procedure ncd_io_1d_int
     module procedure ncd_io_2d_int
     module procedure ncd_io_3d_int
     module procedure ncd_io_4d_int

     module procedure ncd_io_0d_double
     module procedure ncd_io_1d_double
     module procedure ncd_io_2d_double
     module procedure ncd_io_3d_double
     module procedure ncd_io_4d_double

  end interface

  interface ncd_inqvdlen
     module procedure ncd_inqvdlen_byName
  end interface

  interface ncd_inqvdname
     module procedure ncd_inqvdname_byName
  end interface ncd_inqvdname

  integer , parameter  , public  :: max_string_len = 256     ! length of strings
  real(r8), parameter  , public  :: fillvalue = 1.e36_r8     ! fill value for netcdf fields

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ncd_nf90_openfile(fileid, fname, mode)
    !
    ! !DESCRIPTION:
    ! Open a NetCDF fortran file
    !
    ! !ARGUMENTS:
    integer            , intent(inout) :: fileid ! Output file handle
    character(len=*)   , intent(in)    :: fname  ! Input filename to open
    integer            , intent(in)    :: mode   ! file mode: nf90_nowrite, nf90_write
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    status = nf90_open(trim(fname), mode, fileid)

    if(status /= nf90_NOERR) then
       call shr_sys_abort('ncd_nf90_openfile ERROR: Failed to open file -' &
        // trim(fname))
    else
       if (masterproc) &
          write(iulog,*) 'Opened existing file ', trim(fname)
    end if

  end subroutine ncd_nf90_openfile

  !-----------------------------------------------------------------------
  subroutine ncd_nf90_closefile(fileid)
    !
    ! !DESCRIPTION:
    ! Close a NetCDF fortran file
    !
    ! !ARGUMENTS:
    integer, intent(inout) :: fileid   ! file handle to close

    ! ! Local variables
    integer :: status

    !-----------------------------------------------------------------------

    status = nf90_close(fileid)
    if (status /= nf90_NOERR) then
       call shr_sys_abort('ncd_nf90_closefile Failed!' &
          // trim(nf90_strerror(status)) )
    end if

  end subroutine ncd_nf90_closefile

  !-----------------------------------------------------------------------
  subroutine ncd_nf90_createfile(fileid, fname)
    !
    ! !DESCRIPTION:
    ! Create a new NetCDF file with PIO
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer          , intent(inout)  :: fileid  ! file descriptor
    character(len=*) ,  intent(in)    :: fname   ! File name to create
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    status = nf90_create(trim(fname), nf90_CLOBBER, fileid)

    if(status /= nf90_NOERR) then
       call shr_sys_abort( ' ncd_nf90_createfile Failed to open file to write: ' &
         //trim(fname) )
    else
       if (masterproc) &
          write(iulog,*) 'Opened file ', trim(fname),  ' to write', fileid
    end if

  end subroutine ncd_nf90_createfile
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine check_var(ncid, varname, vardesc, readvar, print_err)
    !
    ! !DESCRIPTION:
    ! Check if variable is on netcdf file
    !
    ! !ARGUMENTS:
    integer            , intent(inout)   :: ncid      ! file descriptor
    character(len=*)   , intent(in)      :: varname   ! Varible name to check
    character(len=*)   , intent(out)     :: vardesc   ! Varible description
    logical            , intent(out)     :: readvar   ! If variable exists or not
    logical, optional  , intent(in)      :: print_err ! If should print about error
    !
    ! !LOCAL VARIABLES:
    integer :: varid
    integer :: status     ! return value
    logical :: log_err    ! if should log error
    character(len=*),parameter :: subname='check_var' ! subroutine name
    !-----------------------------------------------------------------------


    if ( present(print_err) )then
       log_err = print_err
    else
       log_err = .true.
    end if
    readvar = .true.

    status = nf90_inq_varid (ncid, varname, varid)
    if (status /= nf90_NOERR) then
       readvar = .false.
       if (masterproc .and. log_err) &
            write(iulog,*) subname//': variable ',trim(varname),' is not on dataset'
    end if
    vardesc = 'not-checked'  !

  end subroutine check_var

  !-----------------------------------------------------------------------
  subroutine check_att(ncid, varid, attrib, att_found)
    !
    ! !DESCRIPTION:
    ! Check if attribute for a variable (by varid) or global att by (nf90_global) is on file
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer           ,intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id, or, nf90_global
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    logical           ,intent(out)   :: att_found ! true if the attribute was found
    !
    ! !LOCAL VARIABLES:
    integer :: status

    character(len=*), parameter :: subname = 'check_att'
    !-----------------------------------------------------------------------

    att_found = .true.
    status = nf90_inquire_attribute(ncid, varid, trim(attrib))
    if (status /= nf90_NOERR) then
       att_found = .false.
    end if

  end subroutine check_att

  !-----------------------------------------------------------------------
  subroutine check_dim(ncid, dimname, value)
    !
    ! !DESCRIPTION:
    ! Validity check on dimension
    !
    ! !ARGUMENTS:
    integer          , intent(in) :: ncid      ! PIO file handle
    character(len=*) , intent(in) :: dimname   ! Dimension name
    integer          , intent(in) :: value     ! Expected dimension size
    !
    ! !LOCAL VARIABLES:
    integer :: dimid, dimlen    ! temporaries
    integer :: status           ! error code      
    character(len=*),parameter :: subname='check_dim' ! subroutine name
    !-----------------------------------------------------------------------

    status = nf90_inq_dimid (ncid, trim(dimname), dimid)
    if (status /= nf90_NOERR) then
       write(iulog,*) subname//' ERROR: to inquire dimension ', trim(dimname)
       call shr_sys_abort(errMsg(__FILE__, __LINE__))
    end if

    status = nf90_inquire_dimension (ncid, dimid, len=dimlen)
    if (dimlen /= value) then
       write(iulog,*) subname//' ERROR: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for dimension ',trim(dimname)
       call shr_sys_abort(errMsg(__FILE__, __LINE__))
    end if

  end subroutine check_dim

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ncd_inqdid(ncid,name,dimid,dimexist)
    !
    ! !DESCRIPTION:
    ! inquire on a dimension id
    !
    ! !ARGUMENTS:
    integer          ,intent(inout):: ncid      ! netcdf file id
    character(len=*) , intent(in)  :: name      ! dimension name
    integer          , intent(out) :: dimid     ! dimension id
    logical,optional , intent(out) :: dimexist  ! if this dimension exists or not
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    status = nf90_inq_dimid(ncid,name,dimid)
    if ( present(dimexist) )then
       if ( status == nf90_NOERR)then
          dimexist = .true.
       else
          dimexist = .false.
       end if
    end if

  end subroutine ncd_inqdid

  !-----------------------------------------------------------------------
  subroutine ncd_inqdlen(ncid,dimid,dlen,name)
    !
    ! !DESCRIPTION:
    ! dimension length inquire, by id or name
    !
    ! !ARGUMENTS:
    integer           , intent(inout) :: ncid       ! netcdf file id
    integer           , intent(inout) :: dimid      ! dimension id
    integer           , intent(out)   :: dlen       ! dimension len
    character(len=*), optional, intent(in) :: name  ! dimension name
    !
    ! !LOCAL VARIABLES:
    integer :: status
    logical :: dimexist
    !-----------------------------------------------------------------------

    dlen = -1
    dimexist = .true.
    if ( present(name) )then
       call ncd_inqdid(ncid,name,dimid, dimexist)
    end if
    if (dimexist) then
       status = nf90_inquire_dimension(ncid,dimid,len=dlen)
        if (status /= nf90_NOERR) then
           call shr_sys_abort('ncd_nf90_inquire_dimension ERROR')
        end if
    else
        call shr_sys_abort('ncd_nf90_inquire_dimension ERROR')
    end if

  end subroutine ncd_inqdlen

  !-----------------------------------------------------------------------
  subroutine ncd_inqdname(ncid,dimid,dname)
    !
    ! !DESCRIPTION:
    ! inquire dim name by id
    !
    ! !ARGUMENTS:
    integer           , intent(in) :: ncid      ! netcdf file id
    integer           , intent(in) :: dimid     ! dimension id
    character(len=*)  , intent(out):: dname     ! dimension name
    !
    ! !LOCAL VARIABLES:
    integer :: status
    !-----------------------------------------------------------------------

    status = nf90_inquire_dimension(ncid,dimid,name=dname)
    if (status /= nf90_NOERR) then
       call shr_sys_abort('ncd_nf90_inquire_dimension ERROR' &
       // trim(nf90_strerror(status)))
    end if

  end subroutine ncd_inqdname

  !-----------------------------------------------------------------------
  subroutine ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    ! this subroutine is ELM specific, especially the dimension name(s)
    ! for earth grids: ni/nj, lon/lat, lsmlon/lsmlat, gridcell (unstructured)
    !
    ! !ARGUMENTS:
    integer           , intent(inout):: ncid
    logical           , intent(out)  :: isgrid2d
    integer           , intent(out)  :: ni
    integer           , intent(out)  :: nj
    integer           , intent(out)  :: ns
    !
    ! !LOCAL VARIABLES:
    integer  :: dimid                             ! netCDF id
    integer  :: status                            ! error status
    character(len=32) :: subname = 'ncd_inqfdims' ! subroutine name
    !-----------------------------------------------------------------------

    ni = 0
    nj = 0
    isgrid2d = .true.

    status = nf90_inq_dimid (ncid, 'lon', dimid)
    if (status == nf90_NOERR) then
       status = nf90_inquire_dimension(ncid, dimid, len=ni)
       status = nf90_inq_dimid (ncid, 'lat', dimid)
       if (status == nf90_NOERR) then
          status = nf90_inquire_dimension(ncid, dimid, len=nj)
       else
          nj = 1
          isgrid2d = .false.
       endif
    !
    else
       status = nf90_inq_dimid (ncid, 'lsmlon', dimid)
       if (status == nf90_NOERR) then
          status = nf90_inquire_dimension(ncid, dimid, len=ni)
          status = nf90_inq_dimid (ncid, 'lsmlat', dimid)
          if (status == nf90_NOERR) then
            status = nf90_inquire_dimension(ncid, dimid, len=nj)
          else
            nj = 1
            isgrid2d = .false.
          endif
       !
       else
          status = nf90_inq_dimid (ncid, 'ni', dimid)
          if (status == nf90_NOERR) then
             status = nf90_inquire_dimension(ncid, dimid, len=ni)
             status = nf90_inq_dimid (ncid, 'nj', dimid)
             if (status == nf90_NOERR) then
               status = nf90_inquire_dimension(ncid, dimid, len=nj)
             else
               nj = 1
               isgrid2d = .false.
             endif
          !
          else
             status = nf90_inq_dimid (ncid, 'gridcell', dimid)
             if (status == nf90_NOERR) then
                status = nf90_inquire_dimension(ncid, dimid, len=ni)
                if (status == nf90_NOERR) nj = 1
                isgrid2d = .false.

             end if
          end if

       endif
    end if
   
    if (ni == 0 .or. nj == 0) then
       write(iulog,*) trim(subname),' ERROR: ni,nj = ',ni,nj,' cannot be zero '
       call shr_sys_abort(errMsg(__FILE__, __LINE__))
    end if


    ns = ni*nj

  end subroutine ncd_inqfdims

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ncd_inqvid(ncid,name,varid,readvar)
    !
    ! !DESCRIPTION:
    ! Inquire on a variable ID by name
    !
    ! !ARGUMENTS:
    integer           , intent(inout) :: ncid      ! netcdf file id
    character(len=*)  , intent(in)    :: name      ! variable name
    integer           , intent(out)   :: varid     ! variable id
    logical, optional , intent(out)   :: readvar   ! does variable exist
    !
    ! !LOCAL VARIABLES:
    integer :: status
    character(len=*),parameter :: subname='ncd_inqvid' ! subroutine name
    !-----------------------------------------------------------------------

    varid = -1
    if (present(readvar)) then
       readvar = .false.
       status = nf90_inq_varid(ncid,name, varid)
       if (status /= nf90_NOERR) then
          if (masterproc) write(iulog,*) subname//': variable ',trim(name),' is not on dataset'
          readvar = .false.
       else
          readvar = .true.
       end if

    else
       status = nf90_inq_varid(ncid,name,varid)

       if (status /= nf90_NOERR) then
         call shr_sys_abort('ncd_inqvdims ERROR ' &
          // trim(nf90_strerror(status)))
       end if
    endif

  end subroutine ncd_inqvid

  !-----------------------------------------------------------------------
  subroutine ncd_inqvdims(ncid, ndims, varname)
    !
    ! !DESCRIPTION:
    ! inquire variable dimension numbers by variable name
    !
    ! !ARGUMENTS:
    integer           , intent(in)   :: ncid      ! netcdf file id
    integer           , intent(out)  :: ndims     ! variable ndims
    character(len=*)  , intent(inout):: varname   ! variable name
    !
    ! !LOCAL VARIABLES:
    integer :: status
    integer :: varid
    character(len=*),parameter :: subname='ncd_inqvdims' ! subroutine name
    !-----------------------------------------------------------------------

    ndims = -1
    status = nf90_inq_varid(ncid,trim(varname),varid)
    if (status /= nf90_NOERR) then
       call shr_sys_abort('ncd_inqvdims ERROR - varid for ' // trim(varname) &
       // trim(nf90_strerror(status)))
    end if

    status = nf90_inquire_variable(ncid, varid, ndims = ndims)
    if (status /= nf90_NOERR) then
       call shr_sys_abort('ncd_inqvdims ERROR - ndims for ' // trim(varname) &
       // trim(nf90_strerror(status)))
    end if

  end subroutine ncd_inqvdims

  !-----------------------------------------------------------------------
  subroutine ncd_inqvname(ncid,varid,vname)
    !
    ! !DESCRIPTION:
    ! inquire variable name by its id
    !
    ! !ARGUMENTS:
    integer           , intent(in)   :: ncid      ! netcdf file id
    integer           , intent(in)   :: varid     ! variable id
    character(len=*)  , intent(out)  :: vname     ! variable vname
    !
    ! !LOCAL VARIABLES:
    integer :: status
    character(len=*),parameter :: subname='ncd_inqvname' ! subroutine name
    !-----------------------------------------------------------------------

    vname = ''
    status = nf90_inquire_variable(ncid, varid, name = vname)
    if (status /= nf90_NOERR) then
       call shr_sys_abort('ncd_inqvname ERROR' //trim(vname) &
       // trim(nf90_strerror(status)))
    end if

  end subroutine ncd_inqvname

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ncd_inqvdids(ncid, varname, dids, error_code)
    !
    ! !DESCRIPTION:
    ! inquire variable' all dimension ids by name
    !
    ! !ARGUMENTS:
    integer          ,intent(in)   :: ncid      ! netcdf file id
    character(len=*) ,intent(in)   :: varname   ! variable descriptor
    integer          ,intent(out)  :: dids(:)   ! variable dids
    integer          ,intent(out)  :: error_code
    !
    ! !LOCAL VARIABLES:
    integer     :: status
    integer     :: varid
    character(len=*),parameter :: subname='ncd_inqvdids' ! subroutine name
    !-----------------------------------------------------------------------

    dids(:) = -1
    status = nf90_inq_varid(ncid, trim(varname), varid)
    if (status == nf90_NOERR) then
       status = nf90_inquire_variable(ncid, varid, dimids = dids)
        if (status /= nf90_NOERR) then
          if (masterproc) write(iulog,*) subname// &
           ' ncd_inqvdids ERROR' //trim(varname) &
           // trim(nf90_strerror(status))
        end if
    end if
    error_code = status

  end subroutine ncd_inqvdids

  !-------------------------------------------------------------------------
  subroutine ncd_inqvdlen_byName(ncid,varname,dimnum,dlen,err_code)
    !
    ! !DESCRIPTION:
    ! inquire size of one of a variable's dimensions, given a variable name
    !
    ! If the variable has n dimensions, then dimnum should be between 1 and n; this routine
    ! returns the size of the dimnum'th dimension.
    !
    ! If there is an error condition, dlen will be -1, and err_code will hold the error
    ! code; possible error codes are:
    ! 0: no error
    ! 1: dimnum out of range
    ! 11: variable not found
    !
    ! !ARGUMENTS:
    integer,intent(inout) :: ncid      ! netcdf file id
    character(len=*)  ,intent(in)    :: varname   ! variable name
    integer           ,intent(in)    :: dimnum    ! dimension number to query
    integer           ,intent(out)   :: dlen      ! length of the dimension
    integer           ,intent(out)   :: err_code  ! error code (0 means no error)
    !
    ! !LOCAL VARIABLES:
    integer            :: status
    integer            :: varid, vdnums
    integer, pointer   :: vdids(:)

    integer, parameter :: dlen_invalid = -1
    integer, parameter :: error_variable_not_found = 11
    character(len=*),parameter :: subname='ncd_inqvdlen_byName' ! subroutine name
    !-----------------------------------------------------------------------

    dlen = dlen_invalid
    err_code = error_variable_not_found

    status = nf90_inq_varid (ncid, varname, varid)
    if (status == nf90_NOERR) then
       status = nf90_inquire_variable(ncid, varid, ndims = vdnums)
       if(status == nf90_NOERR) then
         allocate(vdids(vdnums))
         status = nf90_inquire_variable(ncid, varid, dimids = vdids)
         if (status == nf90_NOERR) then
           err_code = nf90_inquire_dimension(ncid,vdids(dimnum),len=dlen)
         endif
         deallocate(vdids)
       end if

    end if

  end subroutine ncd_inqvdlen_byName

  !-----------------------------------------------------------------------
  subroutine ncd_inqvdname_byName(ncid,varname,dimnum,dname,err_code)
    !
    ! !DESCRIPTION:
    ! Inquire name of one of a variable's dimensions, given a variable name
    !
    ! If the variable has n dimensions, then dimnum should be between 1 and n; this
    ! routine returns the name of the dimnum'th dimension.
    !
    ! If there is an error condition, dname will be ' ', and err_code will hold the error
    ! code; possible error codes are:
    ! 0: no error
    ! 1: dimnum out of range
    ! 11: variable not found
    !
    ! !ARGUMENTS:
    integer           ,intent(inout) :: ncid      ! netcdf file id
    character(len=*)  ,intent(in)    :: varname   ! variable name
    integer           ,intent(in)    :: dimnum    ! dimension number to query
    character(len=*)  ,intent(out)   :: dname     ! name of the dimension
    integer           ,intent(out)   :: err_code  ! error code (0 means no error)
    !
    ! !LOCAL VARIABLES:
    integer            :: status
    integer            :: varid, vdnums
    integer, pointer   :: vdids(:)
    character(len=*), parameter :: dname_invalid = ' '
    integer, parameter :: error_variable_not_found = 11

    character(len=*), parameter :: subname = 'ncd_inqvdname_byName'
    !-----------------------------------------------------------------------

    dname = dname_invalid
    err_code = error_variable_not_found

    status = nf90_inq_varid (ncid, varname, varid)
    if (status == nf90_NOERR) then
       status = nf90_inquire_variable(ncid, varid, ndims = vdnums)
       if(status == nf90_NOERR) then
         allocate(vdids(vdnums))
         status = nf90_inquire_variable(ncid, varid, dimids = vdids)
         if (status == nf90_NOERR) then
           err_code = nf90_inquire_dimension(ncid, vdids(dimnum), name=dname)
         endif
         deallocate(vdids)
       end if
    end if

  end subroutine ncd_inqvdname_byName
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------

  subroutine ncd_putatt_int(ncid,varid,attrib,value)
    !
    ! !DESCRIPTION:
    ! put integer attributes
    !
    ! !ARGUMENTS:
    integer           ,intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    integer           ,intent(in)    :: value     ! netcdf attrib value
    !
    ! !LOCAL VARIABLES:
    integer :: status
    character(len=*),parameter :: subname='ncd_putatt_int' ! subroutine name
    !-----------------------------------------------------------------------

    status = nf90_put_att(ncid,varid,trim(attrib),value)

  end subroutine ncd_putatt_int

  !-----------------------------------------------------------------------
  subroutine ncd_putatt_char(ncid,varid,attrib,value)
    !
    ! !DESCRIPTION:
    ! put character attributes
    !
    ! !ARGUMENTS:
    integer,intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    character(len=*)  ,intent(in)    :: value     ! netcdf attrib value
    !
    ! !LOCAL VARIABLES:
    integer :: status
    character(len=*),parameter :: subname='ncd_putatt_char' ! subroutine name
    !-----------------------------------------------------------------------

    status = nf90_put_att(ncid,varid,trim(attrib),value)

  end subroutine ncd_putatt_char

  !-----------------------------------------------------------------------
  subroutine ncd_putatt_real(ncid,varid,attrib,value,xtype)
    !
    ! !DESCRIPTION:
    ! put real attributes
    !
    ! !ARGUMENTS:
    integer,intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    real(r8)          ,intent(in)    :: value     ! netcdf attrib value
    integer           ,intent(in)    :: xtype     ! netcdf data type
    !
    ! !LOCAL VARIABLES:
    integer :: status
    real*4  :: value4
    character(len=*),parameter :: subname='ncd_putatt_real' ! subroutine name
    !-----------------------------------------------------------------------

    value4 = value

    if (xtype == nf90_double) then
       status = nf90_put_att(ncid,varid,trim(attrib),value)
    else
       status = nf90_put_att(ncid,varid,trim(attrib),value4)
    endif

  end subroutine ncd_putatt_real

  !-----------------------------------------------------------------------
  subroutine ncd_getatt_char(ncid,varid,attrib,value)
    !
    ! !DESCRIPTION:
    ! get a character attribute
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer,intent(inout) :: ncid      ! netcdf file id
    integer           ,intent(in)    :: varid     ! netcdf var id
    character(len=*)  ,intent(in)    :: attrib    ! netcdf attrib
    character(len=*)  ,intent(out)   :: value     ! netcdf attrib value
    !
    ! !LOCAL VARIABLES:
    integer :: status
    
    character(len=*), parameter :: subname = 'ncd_getatt_char'
    !-----------------------------------------------------------------------
    
    status = nf90_get_att(ncid,varid,trim(attrib),value)

  end subroutine ncd_getatt_char

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ncd_defdim(ncid,dname,dlen,dimid)
    !
    ! !DESCRIPTION:
    ! define dimension
    !
    ! !ARGUMENTS:
    integer, intent(in) :: ncid      ! netcdf file id
    character(len=*)  , intent(in) :: dname     ! netcdf dimension name
    integer           , intent(in) :: dlen      ! netcdf dimension length
    integer           , intent(out):: dimid     ! netcdf dimension id
    !
    ! !LOCAL VARIABLES:
    integer :: status
    character(len=*),parameter :: subname='ncd_defdim' ! subroutine name
    !-----------------------------------------------------------------------

    status = nf90_def_dim(ncid,trim(dname),dlen,dimid)

  end subroutine ncd_defdim

  !-----------------------------------------------------------------------
  subroutine ncd_defvar_bynf(ncid, varname, xtype, ndims, dimid, varid, &
       long_name, standard_name, units, cell_method, missing_value,     &
       fill_value, imissing_value, ifill_value, comment, flag_meanings, &
       flag_values, nvalid_range )
    !
    ! !DESCRIPTION:
    !  Define a netcdf variable, by known dimension length and ids
    !
    ! !ARGUMENTS:
    integer            , intent(inout) :: ncid                  ! netcdf file id
    character(len=*)   , intent(in)  :: varname                 ! variable name
    integer            , intent(in)  :: xtype                   ! external type
    integer            , intent(in)  :: ndims                   ! number of dims
    integer            , intent(inout) :: varid                 ! returned var id
    integer            , intent(in), optional :: dimid(:)       ! dimids
    character(len=*)   , intent(in), optional :: long_name      ! attribute
    character(len=*)   , intent(in), optional :: standard_name  ! attribute
    character(len=*)   , intent(in), optional :: units          ! attribute
    character(len=*)   , intent(in), optional :: cell_method    ! attribute
    character(len=*)   , intent(in), optional :: comment        ! attribute
    character(len=*)   , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)           , intent(in), optional :: missing_value  ! attribute for real
    real(r8)           , intent(in), optional :: fill_value     ! attribute for real
    integer            , intent(in), optional :: imissing_value ! attribute for int
    integer            , intent(in), optional :: ifill_value    ! attribute for int
    integer            , intent(in), optional :: flag_values(:) ! attribute for int
    integer            , intent(in), optional :: nvalid_range(2)! attribute for int
    !
    ! !LOCAL VARIABLES:
    integer :: n                   ! indices
    integer :: ldimid(4)           ! local dimid 
    integer :: dimid0(1)           ! local dimid
    integer :: status              ! error status 
    integer :: lxtype              ! local external type (in case logical variable)
    character(len=128) :: dimname  ! temporary
    character(len=256) :: str      ! temporary
    character(len=*),parameter :: subname='ncd_defvar_bynf' ! subroutine name
    !-----------------------------------------------------------------------

    varid = -1

    dimid0 = 0
    ldimid(:) = 0
    if (present(dimid)) then
       ldimid(1:ndims) = dimid(1:ndims)
    else   ! ndims must be zero if dimid not present
       if (ndims /= 0) then
          write(iulog,*) subname//' ERROR: dimid not supplied and ndims ne 0 ',trim(varname),ndims
          call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif
    endif

    if ( xtype == ncd_log )then
       lxtype = ncd_int
    else
       lxtype = xtype
    end if

    if (ndims >  0) then 
       status = nf90_inquire_dimension(ncid,ldimid(ndims),name=dimname)
    end if

    ! Define variable
    if (present(dimid)) then
       status = nf90_def_var(ncid,trim(varname),lxtype,dimid(1:ndims),varid)
    else
       status = nf90_def_var(ncid,trim(varname),lxtype,dimid0        ,varid)
    endif
    !
    ! Add attributes
    !
    if (present(long_name)) then
       call ncd_putatt(ncid, varid, 'long_name', trim(long_name))
    end if
    if (present(standard_name)) then
       if (standard_name /= ' ') then
          call ncd_putatt(ncid, varid, 'standard_name', trim(standard_name))
       end if
    end if
    if (present(flag_values)) then
       status = nf90_put_att(ncid,varid,'flag_values',flag_values)
       if ( .not. present(flag_meanings)) then
          write(iulog,*) 'Error in defining variable = ', trim(varname)
          call shr_sys_abort(" ERROR:: flag_values set -- but not flag_meanings"//errMsg(__FILE__, __LINE__))
       end if
    end if
    if (present(flag_meanings)) then
       if ( .not. present(flag_values)) then
          write(iulog,*) 'Error in defining variable = ', trim(varname)
          call shr_sys_abort(" ERROR:: flag_meanings set -- but not flag_values"//errMsg(__FILE__, __LINE__) )
       end if
       if ( size(flag_values) /= size(flag_meanings) ) then
          write(iulog,*) 'Error in defining variable = ', trim(varname)
          call shr_sys_abort(" ERROR:: flag_meanings and flag_values dimension different"//errMsg(__FILE__, __LINE__))
       end if
       str = flag_meanings(1)
       do n = 1, size(flag_meanings)
          if ( index(flag_meanings(n), ' ') /= 0 )then
             write(iulog,*) 'Error in defining variable = ', trim(varname)
             call shr_sys_abort(" ERROR:: flag_meanings has an invalid space in it"//errMsg(__FILE__, __LINE__) )
          end if
          if ( n > 1 ) str = trim(str)//" "//flag_meanings(n)
       end do
       status = nf90_put_att(ncid,varid,'flag_meanings', trim(str) )
    end if
    if (present(comment)) then
       call ncd_putatt(ncid, varid, 'comment', trim(comment))
    end if
    if (present(units)) then
       call ncd_putatt(ncid, varid, 'units', trim(units))
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call ncd_putatt(ncid, varid, 'cell_methods', trim(str))
    end if
    if (present(fill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', fill_value, lxtype)
    end if
    if (present(missing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', missing_value, lxtype)
    end if
    if (present(ifill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', ifill_value)
    end if
    if (present(imissing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', imissing_value)
    end if
    if (present(nvalid_range)) then
       status = nf90_put_att(ncid,varid,'valid_range', nvalid_range )
    end if
    if ( xtype == ncd_log )then
       status = nf90_put_att(ncid,varid,'flag_values',     (/0, 1/) )
       status = nf90_put_att(ncid,varid,'flag_meanings',  "FALSE TRUE" )
       status = nf90_put_att(ncid,varid,'valid_range',    (/0, 1/) )
    end if

  end subroutine ncd_defvar_bynf

  !-----------------------------------------------------------------------
  subroutine ncd_defvar_bygrid(ncid, varname, xtype,           &
       dim1name, dim2name, dim3name, dim4name, dim5name,       &
       long_name, standard_name, units, cell_method,           &
       missing_value, fill_value, imissing_value, ifill_value, &
       switchdim, comment, flag_meanings, flag_values, nvalid_range )
    !
    ! !DESCRIPTION:
    !  Define a netcdf variable, by flexible dimension names (up to 5)
    !
    ! !ARGUMENTS:
    integer            , intent(inout) :: ncid                 ! netcdf file id
    character(len=*)   , intent(in)  :: varname                 ! variable name
    integer            , intent(in)  :: xtype                   ! external type
    character(len=*)   , intent(in), optional :: dim1name       ! dimension name
    character(len=*)   , intent(in), optional :: dim2name       ! dimension name
    character(len=*)   , intent(in), optional :: dim3name       ! dimension name
    character(len=*)   , intent(in), optional :: dim4name       ! dimension name
    character(len=*)   , intent(in), optional :: dim5name       ! dimension name
    character(len=*)   , intent(in), optional :: long_name      ! attribute
    character(len=*)   , intent(in), optional :: standard_name  ! attribute
    character(len=*)   , intent(in), optional :: units          ! attribute
    character(len=*)   , intent(in), optional :: cell_method    ! attribute
    character(len=*)   , intent(in), optional :: comment        ! attribute
    character(len=*)   , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)           , intent(in), optional :: missing_value  ! attribute for real
    real(r8)           , intent(in), optional :: fill_value     ! attribute for real
    integer            , intent(in), optional :: imissing_value ! attribute for int
    integer            , intent(in), optional :: ifill_value    ! attribute for int
    logical            , intent(in), optional :: switchdim      ! true=> permute dim1 and dim2 for output
    integer            , intent(in), optional :: flag_values(:)  ! attribute for int
    integer            , intent(in), optional :: nvalid_range(2)  ! attribute for int
    !
    ! !LOCAL VARIABLES:
    integer :: n              ! indices
    integer :: ndims          ! dimension counter
    integer :: dimid(5)       ! dimension ids
    integer :: varid          ! variable id
    integer :: itmp           ! temporary
    character(len=256) :: str ! temporary
    character(len=*),parameter :: subname='ncd_defvar_bygrid' ! subroutine name
    !-----------------------------------------------------------------------

    dimid(:) = 0

    ! Determine dimension ids for variable

    if (present(dim1name)) call ncd_inqdid(ncid, dim1name, dimid(1))
    if (present(dim2name)) call ncd_inqdid(ncid, dim2name, dimid(2))
    if (present(dim3name)) call ncd_inqdid(ncid, dim3name, dimid(3))
    if (present(dim4name)) call ncd_inqdid(ncid, dim4name, dimid(4))
    if (present(dim5name)) call ncd_inqdid(ncid, dim5name, dimid(5))

    ! Permute dim1 and dim2 if necessary

    if (present(switchdim)) then
       itmp = dimid(2)
       dimid(2) = dimid(1)
       dimid(1) = itmp
    end if

    ! Define variable

    ndims = 0
    if (present(dim1name)) then
       do n = 1, size(dimid)
          if (dimid(n) /= 0) ndims = ndims + 1
       end do
    end if

    call ncd_defvar_bynf(ncid,varname,xtype,ndims,dimid,varid,   &
         long_name=long_name, standard_name=standard_name,       &
         units=units, cell_method=cell_method,                   &
         missing_value=missing_value, fill_value=fill_value,     &
         imissing_value=imissing_value, ifill_value=ifill_value, &
         comment=comment, flag_meanings=flag_meanings,           &
         flag_values=flag_values, nvalid_range=nvalid_range )

  end subroutine ncd_defvar_bygrid
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ncd_enddef(ncid)
    !
    ! !DESCRIPTION:
    !
    ! !ARGUMENTS:
    integer,intent(inout) :: ncid      ! netcdf file id
    !
    ! !LOCAL VARIABLES:
    integer :: status   ! error status
    character(len=*),parameter :: subname='ncd_enddef' ! subroutine name
    !-----------------------------------------------------------------------

    status = nf90_enddef(ncid)
    if (status /= nf90_NOERR) then
       call shr_sys_abort('ncd_enddef ERROR' &
       // trim(nf90_strerror(status)))
    end if

  end subroutine ncd_enddef
  !-----------------------------------------------------------------------

  !------------------------------------------------------------------------
  subroutine ncd_io_char_var0_start(varname, data, flag, ncid, start )
    !
    ! !DESCRIPTION:
    ! netcdf I/O of character array with start indices input
    !
    ! !ARGUMENTS:
    integer          , intent(inout) :: ncid             ! netcdf file id
    character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname          ! local variable name
    character(len=*) , intent(inout) :: data             ! raw data for this index
    integer          , intent(in)    :: start(:)         ! output bounds
    !
    ! !LOCAL VARIABLES:
    integer :: status               ! error code
    integer :: varid
    logical :: readvar
    character(len=*),parameter :: subname='ncd_io_char_var0_start'
    !-----------------------------------------------------------------------

    call ncd_inqvid(ncid, trim(varname), varid, readvar)
    if (readvar) then
       if (flag == 'read') then
          status = nf90_get_var(ncid, varid, data, start=start)
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED ' //trim(varname) &
               // trim(nf90_strerror(status)))
          end if
       elseif (flag == 'write') then
          status = nf90_put_var(ncid, varid, data, start=start)
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io writeFAILED  ' //trim(varname) &
               // trim(nf90_strerror(status)))
          end if
       end if
    else
       if (masterproc) write(iulog,*) subname//': variable ',trim(varname),' is not on dataset'
    end if

  end subroutine ncd_io_char_var0_start

  !------------------------------------------------------------------------
  subroutine ncd_io_0d_log(varname, data, flag, ncid, readvar)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of logical variable
    !
    ! !ARGUMENTS:
    integer            , intent(inout) :: ncid      ! netcdf file id
    character(len=*)   , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*)   , intent(in)    :: varname   ! variable name
    logical            , intent(inout) :: data      ! raw data
    logical, optional  , intent(out)   :: readvar   ! was var read?
    !
    ! !LOCAL VARIABLES:
    integer           :: varid              ! netCDF variable id
    integer           :: start(1), count(1) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    integer           :: idata
    integer, pointer  :: idata1d(:)         ! Temporary integer data to send to file
    character(len=*),parameter :: subname='ncd_io_0d_log'
    !-----------------------------------------------------------------------

    start(:) = 0
    count(:) = 0

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
       if (varpresent) then
          status = nf90_get_var(ncid, varid, idata)
          if ( idata == 0 )then
             data = .false.
          else if ( idata == 1 )then
             data = .true.
          else
             call shr_sys_abort(' ERROR: bad integer value for logical data'//errMsg(__FILE__, __LINE__))
          end if
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       start(1) = 1 ; count(1) = 1
       call ncd_inqvid(ncid, varname, varid)
       allocate(idata1d(1))
       if ( data )then
          idata1d(1) = 1
       else
          idata1d(1) = 0
       end if
       status = nf90_put_var(ncid, varid, idata1d, start=start, count=count)
       deallocate(idata1d)

    endif   ! flag

  end subroutine ncd_io_0d_log

  !------------------------------------------------------------------------

  subroutine ncd_io_1d_log(varname, data, flag, ncid, readvar, nt)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of logical variable
    !
    ! !ARGUMENTS:
    integer            , intent(inout) :: ncid      ! netcdf file id
    character(len=*)   , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*)   , intent(in)    :: varname   ! variable name
    logical            , intent(inout) :: data(:)   ! raw data
    logical, optional  , intent(out)   :: readvar   ! was var read?
    integer, optional  , intent(in)    :: nt        ! time sample index
    !
    ! !LOCAL VARIABLES:
    integer           :: varid              ! netCDF variable id
    integer           :: start(2), count(2) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    integer           :: idata
    integer, pointer  :: idata1d(:)         ! Temporary integer data to send to file
    character(len=*),parameter :: subname='ncd_io_1d_log'
    !-----------------------------------------------------------------------

    start(:) = 0
    count(:) = 0

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
       if (varpresent) then
          status = nf90_get_var(ncid, varid, idata)
          data(:)=.False.
          data   = (idata == 1)
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       start(1) = 1 ; count(1) = size(data)
       start(2) = 1 ; count(2) = 1       ! need checking here - why 2D index?
       if (present(nt)) start(2) = nt    ! why time index in dim2?
       allocate(idata1d(size(data))) 
       where( data )
          idata1d = 1
       elsewhere
          idata1d = 0
       end where
       call ncd_inqvid(ncid, varname, varid)
       status = nf90_put_var(ncid, varid, idata1d, start=start, count=count)
       deallocate( idata1d )

    endif   ! flag

  end subroutine ncd_io_1d_log
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  subroutine ncd_io_0d_text(varname, data, flag, ncid, readvar)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of text variable of 1 char
    !
    ! !ARGUMENTS:
    integer         ,           intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    character(len=*),           intent(inout) :: data         ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    !
    ! !LOCAL VARIABLES:
    integer           :: varid              ! netCDF variable id
    integer           :: start(1), count(1) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_0d_text'
    !-----------------------------------------------------------------------

    start(1) = 1
    count(1) = len(data)

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
       if (varpresent) then
          data   = ' '
          status = nf90_get_var(ncid, varid, data, start=start, count=count)
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid(ncid, varname, varid)

       status = nf90_put_var(ncid, varid, data, start=start, count=count)

    endif
    !
    if (status /= nf90_NOERR .and. .not.present(readvar)) then
       call shr_sys_abort('ncd_io read/write ERROR' //trim(varname) &
           // trim(nf90_strerror(status)))
    end if

  end subroutine ncd_io_0d_text

  !-----------------------------------------------------------------------

  subroutine ncd_io_1d_text(varname, data, flag, ncid, readvar, nt, starts, counts)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of text (char *)  variable of 1d
    !
    ! !ARGUMENTS:
    integer         ,         intent(inout)   :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    character(len=*),           intent(inout) :: data(:)      ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    integer         , optional, intent(in)    :: starts(1)    ! starts index
    integer         , optional, intent(in)    :: counts(1)    ! counts

    ! !LOCAL VARIABLES:
    character(len=32) :: dimname            ! temporary
    integer           :: varid              ! varid
    integer           :: ndims              ! ndims for var
    integer           :: dims(1)            ! dim sizes
    integer           :: dids(1)            ! dim ids
    integer           :: start(2), count(2) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_1d_text'
    !-----------------------------------------------------------------------

    start(1) = 1
    count(1) = len(data(1))                     ! dim1 is for string length (char*)
    start(2) = 1
    count(2) = size(data,dim=1)
    if (present(starts)) start(2) = starts(1)   ! dim2 is for rest dims
    if (present(counts)) count(2) = counts(1)

    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    status = nf90_inquire_variable(ncid, varid, ndims = ndims)
    call ncd_inqvdids(ncid, varname, dids, status)
    call ncd_inqdname(ncid,dids(ndims),dimname)
    if ('time' == trim(dimname) .and. present(nt)) then
       !this will override starts/counts if input as well
       start(ndims) = nt; count(ndims) = 1
    end if

    if (flag == 'read') then
       if (varpresent) then
          status = nf90_get_var(ncid, varid, data, start=start, count=count)
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       status = nf90_put_var(ncid, varid, data, start=start, count=count)
       if (status /= nf90_NOERR) then
          call shr_sys_abort('ncd_io read/write ERROR' //trim(varname) &
            // trim(nf90_strerror(status)))
       end if

    else if (masterproc) then
       write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
       call shr_sys_abort(errMsg(__FILE__, __LINE__))
    endif

    !
  end subroutine ncd_io_1d_text

  !-----------------------------------------------------------------------

  subroutine ncd_io_2d_text(varname, data, flag, ncid, readvar, nt, starts, counts)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of text variable
    !
    ! !ARGUMENTS:
    integer         ,         intent(inout)   :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    character(len=*),           intent(inout) :: data(:,:)    ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    integer         , optional, intent(in)    :: starts(2)    ! starts index
    integer         , optional, intent(in)    :: counts(2)    ! counts
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: dimname    ! temporary
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: dims(2)    ! dim sizes
    integer           :: dids(2)    ! dim ids
    integer           :: start(3), count(3) ! output bounds
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_2d_text'
    !-----------------------------------------------------------------------

    start(1) = 1
    count(1) = len(data(1,1))                         ! dim1 is for string length (char*)
    start(2) = 1
    count(2) = size(data,dim=1)
    start(3) = 1
    count(3) = size(data,dim=2)
    if (present(starts)) start(2:3) = starts(1:2)   ! this is for rest dims
    if (present(counts)) count(2:3) = counts(1:2)

    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    status = nf90_inquire_variable(ncid, varid, ndims = ndims)
    call ncd_inqvdids(ncid, varname, dids, status)
    call ncd_inqdname(ncid,dids(ndims),dimname)
    if ('time' == trim(dimname) .and. present(nt)) then
       !this will override starts/counts if input as well
       start(ndims) = nt; count(ndims) = 1
    end if

    if (flag == 'read') then
       if (varpresent) then
          status = nf90_get_var(ncid, varid, data, start=start, count=count)
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then
       status = nf90_put_var(ncid, varid, data, start=start, count=count)
       if (status /= nf90_NOERR) then
          call shr_sys_abort('ncd_io read/write ERROR' //trim(varname) &
           // trim(nf90_strerror(status)))
       end if

    else if (masterproc) then
       write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
       call shr_sys_abort(errMsg(__FILE__, __LINE__))
    endif

    !
  end subroutine ncd_io_2d_text
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ncd_io_0d_int(varname, data, flag, ncid, readvar, nt)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    integer         ,           intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    integer(i4)     ,           intent(inout) :: data         ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    !
    ! !LOCAL VARIABLES:
    integer           :: varid              ! netCDF variable id
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_0d_int'
    !-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, trim(varname), varid, readvar=varpresent)
       if (varpresent) then
          status = nf90_get_var(ncid, varid, data)
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid  (ncid, trim(varname), varid)

       status = nf90_put_var(ncid, varid, data)

    endif

  end subroutine ncd_io_0d_int

  !-----------------------------------------------------------------------
  subroutine ncd_io_0d_double(varname, data, flag, ncid, readvar, nt)
    !
    ! !DESCRIPTION:
    ! netcdf I/O of global variable
    !
    ! !ARGUMENTS:
    integer         ,           intent(inout) :: ncid         ! netcdf file id
    character(len=*),           intent(in)    :: flag         ! 'read' or 'write'
    character(len=*),           intent(in)    :: varname      ! variable name
    real(r8)        ,           intent(inout) :: data         ! raw data
    logical         , optional, intent(out)   :: readvar      ! was var read?
    integer         , optional, intent(in)    :: nt           ! time sample index
    !
    ! !LOCAL VARIABLES:
    integer           :: varid              ! netCDF variable id
    integer           :: status             ! error code
    logical           :: varpresent         ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_0d_double'
    !-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, trim(varname), varid, readvar=varpresent)
       if (varpresent) then
          status = nf90_get_var(ncid, varid, data)
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid  (ncid, trim(varname), varid)

       status = nf90_put_var(ncid, varid, data)

    endif

end subroutine ncd_io_0d_double

  !-----------------------------------------------------------------------

  subroutine ncd_io_1d_int(varname, data, dim1name, flag, ncid, nt, &
             readvar, cnvrtnan2fill, starts, counts)
    !
    ! !DESCRIPTION:
    ! netcdf I/O for 1d 
    !
    ! !ARGUMENTS:
    integer          , intent(inout)         :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    integer(i4)      , pointer               :: data(:)       ! local decomposition data
    character(len=*) , optional, intent(in)  :: dim1name      ! dimension 1 name
    integer          , optional, intent(in)  :: nt            ! time sample index
    logical          , optional, intent(out) :: readvar       ! true => variable is on initial dataset (read only)
    logical          , optional, intent(in)  :: cnvrtnan2fill ! true => convert any NaN's to _FillValue (spval)
    integer          , optional, intent(in)  :: starts(1)     ! starts index
    integer          , optional, intent(in)  :: counts(1)     ! counts
    !
    ! Local Variables
    character(len=32)                :: dimname    ! temporary
    integer                          :: varid      ! varid
    integer                          :: ndims      ! ndims for var
    integer                          :: i, dlen(4) ! (up to 4) dim index/sizes
    integer                          :: dids(4)    ! dim ids
    integer                          :: start(4)   ! netcdf start index
    integer                          :: count(4)   ! netcdf count index
    integer, pointer                 :: itemp(:,:,:,:)
    integer, pointer                 :: itemp1d(:)
    integer                          :: status     ! error code  
    logical                          :: varpresent ! if true, variable is on tape
    character(len=*),parameter       :: subname='ncd_io_1d_int' ! subroutine name
    !-----------------------------------------------------------------------

    start(:) = 0
    count(:) = 0
    if (present(starts)) then
      start(1)=starts(1)
    else
      start(1)=1
    endif
    if (present(counts)) then
      count(1)=counts(1)
    else
      count(1)=1
    endif

    !
    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    if (varpresent) then
       status = nf90_inquire_variable(ncid, varid, ndims = ndims)
       call ncd_inqvdids(ncid, varname, dids(1:ndims), status)

       dlen(:)=1
       do i = 1, ndims
          call ncd_inqdname(ncid,dids(i),dimname)
          call ncd_inqvdlen_byName(ncid,trim(varname),i,dlen(i),status)

          if ('time' == trim(dimname) .and. present(nt)) then
             !this will override starts/counts if input as well
             ! i.e. only read/write 1 timestep-sliced data
             start(i) = nt
             count(i) = 1
          elseif (.not. present(counts)) then
             count(i) = dlen(i)
          end if

       end do

       if (flag == 'read') then
          if (ndims == 1) then
             status= nf90_get_var(ncid, varid, data, start=start, count=count)
          else if (ndims == 2) then
             ! it implies read 2D data and then flatten (column-first)
             ! so be aware of start/count if as input
             allocate( itemp(dlen(1),dlen(2),1,1) )
             allocate( itemp1d(dlen(1)*dlen(2)) )
             status= nf90_get_var(ncid, varid, itemp)
             itemp1d = reshape(itemp, shape(itemp1d))
             data = itemp1d(start(1):start(1)+count(1)-1)  !
             deallocate(itemp,itemp1d)
          else if (ndims == 3) then
             ! it implies read 2D data and then flatten (column-first)
             ! so be aware of start/count if as input
             allocate( itemp(dlen(1),dlen(2),dlen(3),1) )
             allocate( itemp1d(dlen(1)*dlen(2)*dlen(3)) )
             status= nf90_get_var(ncid, varid, itemp)
             itemp1d = reshape(itemp, shape(itemp1d))
             data = itemp1d(start(1):start(1)+count(1)-1)  !
             deallocate(itemp,itemp1d)
          else if (ndims == 4) then
             ! it implies read 4D data and then flatten (column-first)
             ! so be aware of start/count if as input
             allocate( itemp(dlen(1),dlen(2),dlen(3),dlen(4)) )
             allocate( itemp1d(dlen(1)*dlen(2)*dlen(3)*dlen(4)) )
             status= nf90_get_var(ncid, varid, itemp)
             itemp1d = reshape(itemp, shape(itemp1d))
             data = itemp1d(start(1):start(1)+count(1)-1)  !
             deallocate(itemp,itemp1d)
          else
             write(iulog,*) trim(subname),' ERROR: ndims over 4 NOT supported yet'
             call shr_sys_abort(errMsg(__FILE__, __LINE__))
          end if

          if (present(readvar)) then
            readvar = varpresent
          elseif (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       elseif (flag == 'write') then
           status= nf90_put_var(ncid, varid, data, start=start, count=count)

           if (status /= nf90_NOERR) then
              call shr_sys_abort('ncd_io read/write ERROR' //trim(varname) &
                 // trim(nf90_strerror(status)))
          end if

       else if (masterproc) then
           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
           call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif

       !
    !
    elseif (masterproc .and. .not.present(readvar)) then
       write(iulog,*) subname//' ERROR: non-existed variable in file id',trim(varname), ncid
       call shr_sys_abort(errMsg(__FILE__, __LINE__))

    endif

  end subroutine ncd_io_1d_int
  !-----------------------------------------------------------------------

  subroutine ncd_io_1d_double(varname, data, dim1name, flag, ncid, nt, &
             readvar, cnvrtnan2fill, starts, counts)
    !
    ! !DESCRIPTION:
    ! netcdf I/O for 1d real
    !
    ! !ARGUMENTS:
    integer, intent(inout)                   :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    real(r8)         , pointer               :: data(:)       ! local decomposition data
    character(len=*) , optional, intent(in)  :: dim1name      ! dimension 1 name
    integer          , optional, intent(in)  :: nt            ! time sample index
    logical          , optional, intent(out) :: readvar       ! true => variable is on initial dataset (read only)
    logical          , optional, intent(in)  :: cnvrtnan2fill ! true => convert any NaN's to _FillValue (spval)
    integer          , optional, intent(in)  :: starts(1)     ! starts index
    integer          , optional, intent(in)  :: counts(1)     ! counts
    !
    ! Local Variables
    character(len=32)                :: dimname    ! temporary
    integer                          :: varid      ! varid
    integer                          :: ndims      ! ndims for var
    integer                          :: i, dlen(4) ! (up to 4) dim index/size
    integer                          :: dids(4)    ! dim ids
    integer                          :: start(4)   ! netcdf start index
    integer                          :: count(4)   ! netcdf count index
    real(r8)         , pointer       :: rtemp(:,:,:,:)
    real(r8)         , pointer       :: rtemp1d(:)
    integer                          :: status     ! error code
    logical                          :: varpresent ! if true, variable is on tape
    character(len=*),parameter       :: subname='ncd_io_1d_double' ! subroutine name
    !-----------------------------------------------------------------------

    start(:) = 0
    count(:) = 0

    if ( present(cnvrtnan2fill) )then
       if (.not. cnvrtnan2fill) then
          call shr_sys_abort(' ERROR: cnvrtnan2fill present but NOT set to true -- MUST set it to TRUE if used'//&
               errMsg(__FILE__, __LINE__))
       endif
    end if

    if (present(starts)) then
      start(1)=starts(1)
    else
      start(1)=1
    endif
    if (present(counts)) then
      count(1)=counts(1)
    else
      count(1)=1   ! it will be over-ride below
    endif

    !
    dlen(:)=1
    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    if (varpresent) then
       status = nf90_inquire_variable(ncid, varid, ndims = ndims)
       call ncd_inqvdids(ncid, varname, dids(1:ndims), status)
       do i = 1, ndims
          call ncd_inqdname(ncid,dids(i),dimname)
          call ncd_inqvdlen_byName(ncid,trim(varname),i,dlen(i),status)

          if ('time' == trim(dimname) .and. present(nt)) then
             !this will override starts/counts if input as well
             ! i.e. only read/write 1 timestep-sliced data
             start(i) = nt
             count(i) = 1
          elseif (.not. present(counts)) then
             count(i) = dlen(i)
          end if
       end do

       if (flag == 'read') then
          if (ndims == 1) then
             status= nf90_get_var(ncid, varid, data, start=start, count=count)
          else if (ndims == 2) then
             ! it implies read 2D data and then flatten (column-first)
             ! so be aware of start/count if as input
             allocate(rtemp(dlen(1),dlen(2),1,1) )
             allocate(rtemp1d(dlen(1)*dlen(2)) )
             status= nf90_get_var(ncid, varid, rtemp)
             rtemp1d = reshape(rtemp, shape(rtemp1d))
             data = rtemp1d(start(1):start(1)+count(1)-1)
             deallocate(rtemp,rtemp1d)
          else if (ndims == 3) then
             ! it implies read 3D data and then flatten
             ! so be aware of start/count if as input
             allocate(rtemp(dlen(1),dlen(2),dlen(3),1) )
             allocate(rtemp1d(dlen(1)*dlen(2)*dlen(3)) )
             status= nf90_get_var(ncid, varid, rtemp)
             rtemp1d = reshape(rtemp, shape(rtemp1d))
             data = rtemp1d(start(1):start(1)+count(1)-1)
             deallocate(rtemp,rtemp1d)
          else if (ndims == 4) then
             ! it implies read 4D data and then flatten
             ! so be aware of start/count if as input
             allocate(rtemp(dlen(1),dlen(2),dlen(3),dlen(4)) )
             allocate(rtemp1d(dlen(1)*dlen(2)*dlen(3)*dlen(4)) )
             status= nf90_get_var(ncid, varid, rtemp)
             rtemp1d = reshape(rtemp, shape(rtemp1d))
             data = rtemp1d(start(1):start(1)+count(1)-1)
             deallocate(rtemp,rtemp1d)
          else
             write(iulog,*) trim(subname),' ERROR: ndims over 4 NOT supported yet'
             call shr_sys_abort(errMsg(__FILE__, __LINE__))
          end if

          if ( present(cnvrtnan2fill) )then
             do i = 1, size(data)
                if ( data(i) == spval ) data(i) = nan
             end do
          end if

          if (present(readvar)) then
            readvar = varpresent
          elseif (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       elseif (flag == 'write') then
          if ( present(cnvrtnan2fill) )then
             do i = 1, size(data)
                if ( isnan(data(i)) ) data(i) = spval
             end do
          end if
          status= nf90_put_var(ncid, varid, data, start=start, count=count)
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io writeFAILED ' //trim(varname) &
                // trim(nf90_strerror(status)))
          end if

       else if (masterproc) then
          write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
          call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif

       !
    !
    elseif (masterproc .and. .not.present(readvar)) then
       write(iulog,*) subname//' ERROR: non-existed variable in file id',trim(varname), ncid
       call shr_sys_abort(errMsg(__FILE__, __LINE__))

    endif

  end subroutine ncd_io_1d_double

  !-----------------------------------------------------------------------

  subroutine ncd_io_2d_int(varname, data, dim1name, lowerb2, upperb2, &
       flag, ncid, nt, readvar, switchdim, cnvrtnan2fill, starts, counts)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 2d integer (single)
    !
    ! !ARGUMENTS:
    integer          , intent(inout)         :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    integer(i4)      , pointer               :: data(:,:)     ! local decomposition input data
    character(len=*) , optional, intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)            :: nt              ! time sample index
    integer, optional, intent(in)            :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical, optional, intent(out)           :: readvar         ! true => variable is on initial dataset (read only)
    logical, optional, intent(in)            :: switchdim       ! true=> permute dim1 and dim2 for output
    logical, optional, intent(in)            :: cnvrtnan2fill   ! true => convert any NaN's to _FillValue (spval)
    integer, optional, intent(in)            :: starts(2)       ! starts index
    integer, optional, intent(in)            :: counts(2)       ! counts
    !
    ! !LOCAL VARIABLES:
    integer , pointer :: temp(:,:)
    integer           :: ndim1,ndim2
    character(len=32) :: dimname    ! temporary      
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: varid      ! varid
    integer           :: n,i,j      ! indices
    integer           :: dlen       ! dim sizes
    integer           :: dids(2)    ! dim ids
    integer           :: start(2)   ! netcdf start index
    integer           :: count(2)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    integer           :: lb1,lb2
    integer           :: ub1,ub2
    character(len=*),parameter :: subname='ncd_io_2d_int' ! subroutine name
    !-----------------------------------------------------------------------

    start(:)=0
    count(:)=0
    if (present(starts)) then
      start(1:2)=starts(1:2)
    else
      start(1:2)=1
    endif
    if (present(counts)) then
      count(1:2)=counts(1:2)
    else
      count(1:2)=1
    endif

    lb1 = lbound(data, dim=1)
    ub1 = ubound(data, dim=1)
    lb2 = lbound(data, dim=2)
    ub2 = ubound(data, dim=2)

    if (present(switchdim)) then
       if (present(lowerb2)) lb2 = lowerb2
       if (present(upperb2)) ub2 = upperb2
       allocate(temp(lb2:ub2,lb1:ub1))
    end if

    !
    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    if (varpresent) then
       status = nf90_inquire_variable(ncid, varid, ndims = ndims)
       if (ndims /= 2) then
          write(iulog,*) trim(subname),' ERROR: ndims must be 2'
          call shr_sys_abort(errMsg(__FILE__, __LINE__))
       end if
       call ncd_inqvdids(ncid, varname, dids, status)
       do i = 1, ndims
          call ncd_inqdname(ncid,dids(i),dimname)
          call ncd_inqvdlen_byName(ncid,trim(varname),i,dlen,status)

          if ('time' == trim(dimname) .and. present(nt)) then
             !this will override starts/counts if input as well
             ! i.e. only read/write 1 timestep-sliced data
             start(i) = nt
             count(i) = 1
          elseif (.not. present(counts)) then
             count(i) = dlen
          end if
       end do

       if (flag == 'read') then
          if (present(switchdim)) then
             status= nf90_get_var(ncid, varid, temp, start=start, count=count)
             do j = lb2,ub2
                do i = lb1,ub1
                   data(i,j) = temp(j,i)
                end do
             end do
          else
             status= nf90_get_var(ncid, varid, data, start=start, count=count)
          end if
          !
          if (present(readvar)) then
            readvar = varpresent
          elseif (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if
       !
       elseif (flag == 'write') then
          if (present(switchdim)) then
             do j = lb2,ub2
                do i = lb1,ub1
                   temp(j,i) = data(i,j)
                end do
             end do
             status= nf90_put_var(ncid, varid, temp, start=start, count=count)
          else
             status= nf90_put_var(ncid, varid, data, start=start, count=count)
          end if
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io writeFAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       else if (masterproc) then
           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
           call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif

    !
    elseif (masterproc .and. .not.present(readvar)) then
       write(iulog,*) subname//' ERROR: non-existed variable in file id',trim(varname), ncid
       call shr_sys_abort(errMsg(__FILE__, __LINE__))

    endif

    if (present(switchdim)) then
       deallocate(temp)
    end if

  end subroutine ncd_io_2d_int

  !-----------------------------------------------------------------------

  subroutine ncd_io_2d_double(varname, data, dim1name, lowerb2, upperb2, &
       flag, ncid, nt, readvar, switchdim, cnvrtnan2fill, starts, counts)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 2d real
    !
    ! !ARGUMENTS:
    integer, intent(inout)                   :: ncid            ! netcdf file id
    character(len=*) , intent(in)            :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname         ! variable name
    real(r8)         , pointer               :: data(:,:)       ! local decomposition input data
    character(len=*) , optional, intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)            :: nt              ! time sample index
    integer, optional, intent(in)            :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical, optional, intent(out)           :: readvar         ! true => variable is on initial dataset (read only)
    logical, optional, intent(in)            :: switchdim       ! true=> permute dim1 and dim2 for output
    logical, optional, intent(in)            :: cnvrtnan2fill   ! true => convert any NaN's to _FillValue (spval)
    integer, optional, intent(in)            :: starts(2)       ! starts index
    integer, optional, intent(in)            :: counts(2)       ! counts
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: temp(:,:)
    integer           :: ndim1,ndim2       
    character(len=32) :: dimname    ! temporary      
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: varid      ! varid
    integer           :: n,i,j      ! indices
    integer           :: dlen       ! dim size
    integer           :: dids(2)    ! dim ids
    integer           :: start(2)   ! netcdf start index
    integer           :: count(2)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    integer           :: lb1,lb2
    integer           :: ub1,ub2
    character(len=*),parameter :: subname='ncd_io_2d_double' ! subroutine name
    !-----------------------------------------------------------------------

    start(:)=0
    count(:)=0
    if (present(starts)) then
      start(1:2)=starts(1:2)
    else
      start(1:2)=1
    endif
    if (present(counts)) then
      count(1:2)=counts(1:2)
    else
      count(1:2)=1
    endif

    if ( present(cnvrtnan2fill) )then
       if (.not. cnvrtnan2fill) then
          call shr_sys_abort( ' ERROR: cnvrtnan2fill present but NOT set to true -- MUST set it to TRUE if used'//&
               errMsg(__FILE__, __LINE__))
       endif
    end if

    lb1 = lbound(data, dim=1)
    ub1 = ubound(data, dim=1)
    lb2 = lbound(data, dim=2)
    ub2 = ubound(data, dim=2)

    if (present(switchdim)) then
       if (present(lowerb2)) lb2 = lowerb2
       if (present(upperb2)) ub2 = upperb2
       allocate(temp(lb2:ub2,lb1:ub1))
    end if

    !
    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    if (varpresent) then
       status = nf90_inquire_variable(ncid, varid, ndims = ndims)
       if (ndims /= 2) then
          write(iulog,*) trim(subname),' ERROR: ndims must be 2'
          call shr_sys_abort(errMsg(__FILE__, __LINE__))
       end if
       call ncd_inqvdids(ncid, varname, dids, status)
       do i = 1, ndims
          call ncd_inqdname(ncid,dids(i),dimname)
          call ncd_inqvdlen_byName(ncid,trim(varname),i,dlen,status)

          if ('time' == trim(dimname) .and. present(nt)) then
             !this will override starts/counts if input as well
             ! i.e. only read/write 1 timestep-sliced data
             start(i) = nt
             count(i) = 1
          elseif (.not. present(counts)) then
             count(i) = dlen
          end if
       end do

       if (flag == 'read') then
          if (present(switchdim)) then
             status= nf90_get_var(ncid, varid, temp, start=start, count=count)
             do j = lb2,ub2
                do i = lb1,ub1
                   data(i,j) = temp(j,i)
                end do
             end do
          else
             status= nf90_get_var(ncid, varid, data, start=start, count=count)
          end if
          !
          if ( present(cnvrtnan2fill) )then
             do j = lb2,ub2
                do i = lb1,ub1
                   if ( data(i,j) == spval )then
                      data(i,j) = nan
                   end if
                end do
             end do
          end if
          !
          if (present(readvar)) then
            readvar = varpresent
          elseif (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if
       !
       elseif (flag == 'write') then
          if ( present(cnvrtnan2fill) )then
             do j = lb2,ub2
                do i = lb1,ub1
                   if ( isnan(data(i,j)) )then
                      data(i,j) = spval
                   end if
                end do
             end do
          end if
          !
          if (present(switchdim)) then
             do j = lb2,ub2
                do i = lb1,ub1
                   temp(j,i) = data(i,j)
                end do
             end do
             status= nf90_put_var(ncid, varid, temp, start=start, count=count)
          else
             status= nf90_put_var(ncid, varid, data, start=start, count=count)
          end if
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io writeFAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       else if (masterproc) then
           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
           call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif

       !
    !
    elseif (masterproc .and. .not.present(readvar)) then
       write(iulog,*) subname//' ERROR: non-existed variable in file id',trim(varname), ncid
       call shr_sys_abort(errMsg(__FILE__, __LINE__))

    endif

    if (present(switchdim)) then
       deallocate(temp)
    end if

  end subroutine ncd_io_2d_double

  !-----------------------------------------------------------------------

  subroutine ncd_io_3d_int(varname, data, dim1name, flag, ncid, nt, &
       readvar, starts, counts)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 3d integer (single)
    !
    ! !ARGUMENTS:
    integer          , intent(inout)         :: ncid            ! netcdf file id
    character(len=*) , intent(in)            :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname         ! variable name
    integer(i4)      , pointer               :: data(:,:,:)     ! local decomposition input data
    character(len=*) , optional, intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)            :: nt              ! time sample index
    logical, optional, intent(out)           :: readvar         ! true => variable is on initial dataset (read only)
    integer, optional, intent(in)            :: starts(3)       ! starts index
    integer, optional, intent(in)            :: counts(3)       ! counts
    !
    ! !LOCAL VARIABLES:
    integer           :: ndim1,ndim2
    character(len=32) :: dimname    ! temporary
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: varid      ! varid
    integer           :: i, dlen    ! dim size
    integer           :: dids(3)    ! dim ids
    integer           :: start(3)   ! netcdf start index
    integer           :: count(3)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_3d_int' ! subroutine name
    !-----------------------------------------------------------------------

    start(:)=0
    count(:)=0
    if (present(starts)) then
      start(1:3)=starts(1:3)
    else
      start(1:3)=1
    endif
    if (present(counts)) then
      count(1:3)=counts(1:3)
    else
      count(1:3)=1
    endif

    !
    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    if (varpresent) then
       status = nf90_inquire_variable(ncid, varid, ndims = ndims)
       if (ndims /= 3) then
          write(iulog,*) trim(subname),' ERROR: ndims must be 3'
          call shr_sys_abort(errMsg(__FILE__, __LINE__))
       end if
       call ncd_inqvdids(ncid, varname, dids, status)
       do i = 1, ndims
          call ncd_inqdname(ncid,dids(i),dimname)
          call ncd_inqvdlen_byName(ncid,trim(varname),i,dlen,status)

          if ('time' == trim(dimname) .and. present(nt)) then
             !this will override starts/counts if input as well
             ! i.e. only read/write 1 timestep-sliced data
             start(i) = nt
             count(i) = 1
          elseif (.not. present(counts)) then
             count(i) = dlen
          end if
       end do

       if (flag == 'read') then
          status= nf90_get_var(ncid, varid, data, start=start, count=count)
          !
          if (present(readvar)) then
            readvar = varpresent
          elseif (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if
       !
       elseif (flag == 'write') then
          status= nf90_put_var(ncid, varid, data, start=start, count=count)
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io writeFAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       else if (masterproc) then
           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
           call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif

       !
    !
    elseif (masterproc .and. .not.present(readvar)) then
       write(iulog,*) subname//' ERROR: non-existed variable in file id',trim(varname), ncid
       call shr_sys_abort(errMsg(__FILE__, __LINE__))

    endif

  end subroutine ncd_io_3d_int

  !-----------------------------------------------------------------------

  subroutine ncd_io_3d_double(varname, data, dim1name, flag, ncid, nt, &
       readvar, starts, counts)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 3d real data
    !
    ! !ARGUMENTS:
    integer          , intent(inout)         :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    real(r8)         , pointer               :: data(:,:,:)   ! local decomposition input data
    character(len=*) , optional, intent(in)  :: dim1name      ! dimension 1 name
    integer, optional, intent(in)            :: nt            ! time sample index
    logical, optional, intent(out)           :: readvar       ! true => variable is on initial dataset (read only)
    integer, optional, intent(in)            :: starts(3)     ! starts index
    integer, optional, intent(in)            :: counts(3)     ! counts
    !
    ! !LOCAL VARIABLES:
    integer           :: ndim1,ndim2
    character(len=32) :: dimname    ! temporary
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: varid      ! varid
    integer           :: n,i,j      ! indices
    integer           :: dlen       ! dim size
    integer           :: dids(3)    ! dim ids
    integer           :: start(3)   ! netcdf start index
    integer           :: count(3)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_3d_double' ! subroutine name
    !-----------------------------------------------------------------------

    start(:)=0
    count(:)=0
    if (present(starts)) then
      start(1:3)=starts(1:3)
    else
      start(1:3)=1
    endif
    if (present(counts)) then
      count(1:3)=counts(1:3)
    else
      count(1:3)=1
    endif

    !
    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    if (varpresent) then
       status = nf90_inquire_variable(ncid, varid, ndims = ndims)
       if (ndims /= 3) then
          write(iulog,*) trim(subname),' ERROR: ndims must be 3'
          call shr_sys_abort(errMsg(__FILE__, __LINE__))
       end if
       call ncd_inqvdids(ncid, varname, dids, status)
       do i = 1, ndims
          call ncd_inqdname(ncid,dids(i),dimname)
          call ncd_inqvdlen_byName(ncid,trim(varname),i,dlen,status)

          if ('time' == trim(dimname) .and. present(nt)) then
             !this will override starts/counts if input as well
             ! i.e. only read/write 1 timestep-sliced data
             start(i) = nt
             count(i) = 1
          elseif (.not. present(counts)) then
             count(i) = dlen
          end if
       end do

       if (flag == 'read') then
          status= nf90_get_var(ncid, varid, data, start=start, count=count)
          !
          if (present(readvar)) then
            readvar = varpresent
          elseif (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if
       !
       elseif (flag == 'write') then
          status= nf90_put_var(ncid, varid, data, start=start, count=count)
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io writeFAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       else if (masterproc) then
           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
           call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif

       !
    !
    elseif (masterproc .and. .not.present(readvar)) then
       write(iulog,*) subname//' ERROR: non-existed variable in file id',trim(varname), ncid
       call shr_sys_abort(errMsg(__FILE__, __LINE__))

    endif

  end subroutine ncd_io_3d_double

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine ncd_io_4d_int(varname, data, dim1name, flag, ncid, nt, &
       readvar, starts, counts)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 4d integer (single)
    !
    ! !ARGUMENTS:
    integer          , intent(inout)         :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    integer(i4)      , pointer               :: data(:,:,:,:) ! local decomposition input data
    character(len=*) , optional, intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)            :: nt              ! time sample index
    logical, optional, intent(out)           :: readvar         ! true => variable is on initial dataset (read only)
    integer, optional, intent(in)            :: starts(4)       ! starts index
    integer, optional, intent(in)            :: counts(4)       ! counts
    !
    ! !LOCAL VARIABLES:
    integer           :: ndim1,ndim2
    character(len=32) :: dimname    ! temporary
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: varid      ! varid
    integer           :: i, dlen    ! dim size
    integer           :: dids(4)    ! dim ids
    integer           :: start(4)   ! netcdf start index
    integer           :: count(4)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_4d_int' ! subroutine name
    !-----------------------------------------------------------------------

    start(:)=0
    count(:)=0
    if (present(starts)) then
      start(1:4)=starts(1:4)
    else
      start(1:4)=1
    endif
    if (present(counts)) then
      count(1:4)=counts(1:4)
    else
      count(1:4)=1
    endif

    !
    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    if (varpresent) then
       status = nf90_inquire_variable(ncid, varid, ndims = ndims)
       if (ndims /= 4) then
          write(iulog,*) trim(subname),' ERROR: ndims must be 4'
          call shr_sys_abort(errMsg(__FILE__, __LINE__))
       end if
       call ncd_inqvdids(ncid, varname, dids, status)
       do i = 1, ndims
          call ncd_inqdname(ncid,dids(i),dimname)
          call ncd_inqvdlen_byName(ncid,trim(varname),i,dlen,status)

          if ('time' == trim(dimname) .and. present(nt)) then
             !this will override starts/counts if input as well
             ! i.e. only read/write 1 timestep-sliced data
             start(i) = nt
             count(i) = 1
          elseif (.not. present(counts)) then
             count(i) = dlen
          end if
       end do

       if (flag == 'read') then
          status= nf90_get_var(ncid, varid, data, start=start, count=count)
          !
          if (present(readvar)) then
            readvar = varpresent
          elseif (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       !
       elseif (flag == 'write') then
          status= nf90_put_var(ncid, varid, data, start=start, count=count)
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io writeFAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       else if (masterproc) then
           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
           call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif

       !
    !
    elseif (masterproc .and. .not.present(readvar)) then
       write(iulog,*) subname//' ERROR: non-existed variable in file id',trim(varname), ncid
       call shr_sys_abort(errMsg(__FILE__, __LINE__))

    endif

  end subroutine ncd_io_4d_int

  !-----------------------------------------------------------------------

  subroutine ncd_io_4d_double(varname, data, dim1name, flag, ncid, nt, &
       readvar, starts, counts)
    !
    ! !DESCRIPTION:
    ! Netcdf i/o of 4d real data
    !
    ! !ARGUMENTS:
    integer          , intent(inout)         :: ncid          ! netcdf file id
    character(len=*) , intent(in)            :: flag          ! 'read' or 'write'
    character(len=*) , intent(in)            :: varname       ! variable name
    real(r8)         , pointer               :: data(:,:,:,:) ! local decomposition input data
    character(len=*) , optional, intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)            :: nt              ! time sample index
    logical, optional, intent(out)           :: readvar         ! true => variable is on initial dataset (read only)
    integer, optional, intent(in)            :: starts(4)       ! starts index
    integer, optional, intent(in)            :: counts(4)       ! counts
    !
    ! !LOCAL VARIABLES:
    integer           :: ndim1,ndim2
    character(len=32) :: dimname    ! temporary
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: varid      ! varid
    integer           :: n,i,j      ! indices
    integer           :: dlen       ! dim size
    integer           :: dids(4)    ! dim ids
    integer           :: start(4)   ! netcdf start index
    integer           :: count(4)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    character(len=*),parameter :: subname='ncd_io_4d_double' ! subroutine name
    !-----------------------------------------------------------------------

    start(:)=0
    count(:)=0
    if (present(starts)) then
      start(1:4)=starts(1:4)
    else
      start(1:4)=1
    endif
    if (present(counts)) then
      count(1:4)=counts(1:4)
    else
      count(1:4)=1
    endif

    !
    call ncd_inqvid(ncid, varname, varid, readvar=varpresent)
    if (varpresent) then
       status = nf90_inquire_variable(ncid, varid, ndims = ndims)
       if (ndims /= 4) then
          write(iulog,*) trim(subname),' ERROR: ndims must be 4'
          call shr_sys_abort(errMsg(__FILE__, __LINE__))
       end if
       call ncd_inqvdids(ncid, varname, dids, status)
       do i = 1, ndims
          call ncd_inqdname(ncid,dids(i),dimname)
          call ncd_inqvdlen_byName(ncid,trim(varname),i,dlen,status)

          if ('time' == trim(dimname) .and. present(nt)) then
             !this will override starts/counts if input as well
             ! i.e. only read/write 1 timestep-sliced data
             start(i) = nt
             count(i) = 1
          elseif (.not. present(counts)) then
             count(i) = dlen
          end if
       end do

       if (flag == 'read') then
          status= nf90_get_var(ncid, varid, data, start=start, count=count)
          !
          if (present(readvar)) then
            readvar = varpresent
          elseif (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io read FAILED ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if
       !
       elseif (flag == 'write') then
          status= nf90_put_var(ncid, varid, data, start=start, count=count)
          if (status /= nf90_NOERR) then
             call shr_sys_abort('ncd_io writeFAILED  ' //trim(varname) &
                // ' with error - ' //trim(nf90_strerror(status)))
          end if

       else if (masterproc) then
           write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
           call shr_sys_abort(errMsg(__FILE__, __LINE__))
       endif

       !
    !
    elseif (masterproc .and. .not.present(readvar)) then
       write(iulog,*) subname//' ERROR: non-existed variable in file id',trim(varname), ncid
       call shr_sys_abort(errMsg(__FILE__, __LINE__))

    endif

  end subroutine ncd_io_4d_double

  !------------------------------------------------------------------------

end module ncdio_nf90Mod
