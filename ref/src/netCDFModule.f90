!------------------------------------------------------------------
! Wrappers over NetCDF API to simplify kernel I/O code.
!------------------------------------------------------------------

module netCDFModule

  use netcdf

  implicit none

  character(len=15) :: routineName

contains

  !------------------------------------------------------------------
  ! open_file
  !
  ! Opens a netCDF file
  !------------------------------------------------------------------
  subroutine open_file(filename, mode, ncFileID)

    character(len=*), intent( in) :: filename
    character(len=1), intent( in) :: mode
    integer,          intent(out) :: ncFileID

    routineName="open_file"

    if (mode =="r") then
      call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncFileID))
    elseif (mode=="w") then
      call nc_check(nf90_create(trim(filename), NF90_CLOBBER, ncFileID))
    else
      write(*,*) "ERROR: Unsupported mode for open_file"
      stop 1
    end if

  end subroutine open_file

  !------------------------------------------------------------------
  ! close_file
  !
  ! Closes a netCDF file
  !------------------------------------------------------------------
  subroutine close_file(ncFileID)

    integer, intent(in) :: ncFileID

    routineName="close_file"

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

  end subroutine close_file

  !------------------------------------------------------------------
  ! define_on
  !
  ! Turns define mode on
  !------------------------------------------------------------------
  subroutine define_on(ncFileID)

    integer, intent(in) :: ncFileID

    routineName="define_on"

    call nc_check(nf90_redef(ncfileID))

  end subroutine define_on

  !------------------------------------------------------------------
  ! define_off
  !
  ! Turns define mode off
  !------------------------------------------------------------------
  subroutine define_off(ncFileID)

    integer, intent(in) :: ncFileID

    routineName="define_off"

    call nc_check(nf90_enddef(ncfileID))

  end subroutine define_off

  !------------------------------------------------------------------
  ! define_dim
  !
  ! Defines a dimensions
  !------------------------------------------------------------------
  subroutine define_dim(ncFileID, dimName, dimSize, dimID)

    integer,          intent( in) :: ncFileID
    character(len=*), intent( in) :: dimName
    integer,          intent( in) :: dimSize
    integer,          intent(out) :: dimID

    routineName="define_dim"

    call nc_check(nf90_def_dim(ncid=ncFileID, name=TRIM(dimName), len=dimSize, dimid = dimID))

  end subroutine define_dim

  !------------------------------------------------------------------
  ! define_var_2d_real
  !
  ! Defines a real 2D variable
  !------------------------------------------------------------------
  subroutine define_var_2d_real(ncFileID, varName, dim1ID, dim2ID, varID)

    integer,          intent( in) :: ncFileID
    character(len=*), intent( in) :: varName
    integer,          intent( in) :: dim1ID, dim2ID
    integer,          intent(out) :: varID

    routineName="define_var_2d_real"

    call nc_check(nf90_def_var(ncid=ncFileID, name=TRIM(varName), &
                  xtype=nf90_double, dimids=(/dim1ID, dim2ID/),   &
                  varid=varID))

  end subroutine define_var_2d_real

  !------------------------------------------------------------------
  ! define_var_3d_real
  !
  ! Defines a real 3D variable
  !------------------------------------------------------------------
  subroutine define_var_3d_real(ncFileID, varName, dim1ID, dim2ID, dim3ID, varID)

    integer,          intent( in) :: ncFileID
    character(len=*), intent( in) :: varName
    integer,          intent( in) :: dim1ID, dim2ID, dim3ID
    integer,          intent(out) :: varID

    routineName="define_var_3d_real"

    call nc_check(nf90_def_var(ncid=ncFileID, name=TRIM(varName),       &
                  xtype=nf90_double, dimids=(/dim1ID, dim2ID, dim3ID/), &
                  varid=varID))

  end subroutine define_var_3d_real

  !------------------------------------------------------------------
  ! write_var_2d_real
  !
  ! Writes a real 2D variable
  !------------------------------------------------------------------
  subroutine write_var_2d_real(ncFileID, varID, data)

    integer, intent(in) :: ncFileID
    integer, intent(in) :: varID
    real,    intent(in) :: data(:,:)

    routineName="write_var_2d_real"

    call nc_check(nf90_put_var(ncFileID, varID, data))

  end subroutine write_var_2d_real

  !------------------------------------------------------------------
  ! write_var_3d_real
  !
  ! Writes a real 3D variable
  !------------------------------------------------------------------
  subroutine write_var_3d_real(ncFileID, varID, data)

    integer, intent(in) :: ncFileID
    integer, intent(in) :: varID
    real,    intent(in) :: data(:,:,:)

    routineName="write_var_3d_real"

    call nc_check(nf90_put_var(ncFileID, varID, data))

  end subroutine write_var_3d_real

  !------------------------------------------------------------------
  ! read_global_logical
  !
  ! Reads a global logical
  !------------------------------------------------------------------
  subroutine read_global_logical(ncFileID, globalName, data)

    integer,          intent(   in) :: ncFileID
    character(len=*), intent(   in) :: globalName
    logical,          intent(inout) :: data

    integer :: val

    routineName="read_global_logical"

    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, TRIM(globalName), val))
    if (val == 1) then
      data = .true.
    else
      data = .false.
    end if

  end subroutine read_global_logical

  !------------------------------------------------------------------
  ! read_global_int
  !
  ! Reads a global integer
  !------------------------------------------------------------------
  subroutine read_global_int(ncFileID, globalName, data)

    integer,          intent(   in) :: ncFileID
    character(len=*), intent(   in) :: globalName
    integer,          intent(inout) :: data

    routineName="read_global_int"

    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, TRIM(globalName), data))

  end subroutine read_global_int

  !------------------------------------------------------------------
  ! read_global_real
  !
  ! Reads a global real
  !------------------------------------------------------------------
  subroutine read_global_real(ncFileID, globalName, data)

    integer,          intent(   in) :: ncFileID
    character(len=*), intent(   in) :: globalName
    real,             intent(inout) :: data

    routineName="read_global_real"

    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, TRIM(globalName), data))

  end subroutine read_global_real

  !------------------------------------------------------------------
  ! read_dimension
  !
  ! Reads a dimension from NetCDF file.
  !------------------------------------------------------------------
  subroutine read_dimension(ncFileID, dimName, data)

    integer,          intent(   in) :: ncFileID
    character(len=*), intent(   in) :: dimName
    integer,          intent(inout) :: data

    integer :: dimID

    routineName="read_dimension"

    call nc_check(nf90_inq_dimid(ncFileID, TRIM(dimName), dimID))
    call nc_check(nf90_inquire_dimension(ncFileID, dimID, len=data))

  end subroutine read_dimension

  !------------------------------------------------------------------
  ! read_3d_real
  !
  ! Reads a 3D real array from NetCDF file.
  !------------------------------------------------------------------
  subroutine read_3d_real(ncFileID, varName, data)

    integer,          intent(   in) :: ncFileID
    character(len=*), intent(   in) :: varName
    real,             intent(inout) :: data(:,:,:)

    integer :: varID

    routineName="read_3d_real"

    call nc_check(nf90_inq_varid(ncFileID, TRIM(varName), varID))
    call nc_check(nf90_get_var(ncFileID, varID, data))

  end subroutine read_3d_real

  !------------------------------------------------------------------
  ! read_2d_real
  !
  ! Reads a 2D real array from NetCDF file.
  !------------------------------------------------------------------
  subroutine read_2d_real(ncFileID, varName, data)

    integer,          intent(   in) :: ncFileID
    character(len=*), intent(   in) :: varName
    real,             intent(inout) :: data(:,:)

    integer :: varID

    routineName="read_2d_real"

    call nc_check(nf90_inq_varid(ncFileID, TRIM(varName), varID))
    call nc_check(nf90_get_var(ncFileID, varID, data))

  end subroutine read_2d_real

  !------------------------------------------------------------------
  ! write_global_character
  !
  ! Reads a global character
  !------------------------------------------------------------------
  subroutine write_global_character(ncFileID, globalName, data)

    integer,          intent(in) :: ncFileID
    character(len=*), intent(in) :: globalName
    character(len=*), intent(in) :: data

    routineName="write_global_character"

    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, TRIM(globalName), TRIM(data)))

  end subroutine write_global_character

  !------------------------------------------------------------------
  ! write_global_logical
  !
  ! Reads a global logical
  !------------------------------------------------------------------
  subroutine write_global_logical(ncFileID, globalName, data)

    integer,          intent(in) :: ncFileID
    character(len=*), intent(in) :: globalName
    logical,          intent(in) :: data

    routineName="write_global_logical"

    if (data) then
      call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, TRIM(globalName), 1))
    else
      call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, TRIM(globalName), 0))
    end if

  end subroutine write_global_logical

  !------------------------------------------------------------------
  ! write_global_int
  !
  ! Reads a global integer
  !------------------------------------------------------------------
  subroutine write_global_int(ncFileID, globalName, data)

    integer,          intent(in) :: ncFileID
    character(len=*), intent(in) :: globalName
    integer,          intent(in) :: data

    routineName="write_global_int"

    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, TRIM(globalName), data))

  end subroutine write_global_int

  !------------------------------------------------------------------
  ! write_global_real
  !
  ! Reads a global real
  !------------------------------------------------------------------
  subroutine write_global_real(ncFileID, globalName, data)

    integer,          intent(in) :: ncFileID
    character(len=*), intent(in) :: globalName
    real,             intent(in) :: data

    routineName="write_global_real"

    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, TRIM(globalName), data))

  end subroutine write_global_real

  !------------------------------------------------------------------
  ! nc_check
  !
  ! Checks return status from a NetCDF API call.  If an error was
  ! returned, print the message and abort the program.
  !------------------------------------------------------------------
  subroutine nc_check(istatus)

    integer, intent (in) :: istatus

    character(len=512) :: error_msg

    ! if no error, nothing to do here.  we are done.
    if( istatus == nf90_noerr) return

    error_msg = nf90_strerror(istatus)

    write(*, '(A9,A,A2,A)') "Error in ", TRIM(routineName), ": ", error_msg
    stop 1

  end subroutine nc_check

end module netCDFModule
