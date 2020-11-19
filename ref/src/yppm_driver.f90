!------------------------------------------------------------------
! yppm_driver
!
! Driver program to run the yppm kernel
!------------------------------------------------------------------
program yppm_driver

  use OMP_LIB
  use yppm_core_mod, ONLY: yppm
  use netCDFModule
  use ieee_arithmetic

  implicit none

  ! State variables
  integer :: js, je, isd, ied, jsd, jed
  integer :: npx
  integer :: npy
  integer :: grid_type, ord_in
  real    :: lim_fac
  logical :: nested, regional
  real, allocatable :: cry(:,:)
  real, allocatable :: q(:,:)
  real, allocatable :: fy2(:,:)
  real, allocatable :: dya(:,:)

  ! Driver variables
  integer   :: nthreads                              ! # of OpenMP threads
  integer   :: narg                                  ! # of command line arguments
  integer*8 :: count_start, count_end, count_rate    ! Timer start/stop
  integer   :: ret                                   ! Return status
  logical   :: underflow_support, gradual, underflow ! Underflow control vars
  real      :: fptest                                ! Underflow control var
  integer, external :: print_affinity                ! External subroutine

  ! Input configuration variables
  character(len=64) :: namelist_file = "test_input/yppm_0.0.1.nl"
  character(len=64) :: input_file, output_file
  integer           :: nl_unit

  ! Input namelists
  namelist /io/       input_file, output_file

  ! Get the number of arguments
  narg = command_argument_count()
  if (narg /= 1) then
    write(*,*) "Usage: yppm <namelist_file>"
    stop 1
  end if

  ! Get the namelist file name from the argument list
  call get_command_argument(1, namelist_file)

  ! Open the namelist file
  open(newunit=nl_unit, file=TRIM(namelist_file), form='formatted', status='old')

  ! Read the data IO settings from the namelist
  read(nl_unit, nml=io)

  ! Get OMP_NUM_THREADS value
  nthreads = omp_get_max_threads()

  ! Force thie program to run in abrupt underflow mode
  ! Intel's -ftz option flushes denormalized values to zero for non-SSE instructions anyway
  underflow_support = ieee_support_underflow_control(fptest)
  if (underflow_support) then
    write(*,*) "Underflow control supported for the default real kind"
  else
    write(*,*) "No underflow control support"
    stop 1
  end if
  call ieee_set_underflow_mode(.false.)
  call ieee_get_underflow_mode(gradual)
  if (.not. gradual) then 
    write(*,*) "Able to set abrupt underflow mode"
  else    
    write(*,*) "Error setting abrupt underflow mode"
    stop 1
  end if

  ! Print out configuration settings
  write(*, '(A,A)') 'Input file = ', TRIM(input_file)
  write(*, '(A,I0)') 'nthreads = ', nthreads
  ret = print_affinity(0)

  ! Read the input state from the NetCDF input file
  call read_state(TRIM(input_file))

  ! Print the input state
  call print_state("Input State")

  ! Get the start time
  call system_clock(count_start, count_rate)

  ! Run the kernel
  call yppm(fy2, q, cry, ord_in, isd, ied, isd, ied, js, je, jsd, jed, npx, npy, dya, &
            nested, grid_type, lim_fac, regional)

  ! Get the stop time
  call system_clock(count_end, count_rate)  

  ! Print the output state
  call print_state("Output State")

  ! Write the output state to the NetCDF output file
  call write_state(TRIM(output_file))

  ! Write timing information
  write(*,*)
  write(*,'(A,F12.9)') "Total time (in seconds)=" ,  ((count_end - count_start) * 1.0) / count_rate

  ! Deallocate the state variables
  call deallocate_state()

contains

  !------------------------------------------------------------------
  ! print_state
  !
  ! prints statistics for the kernel state variables
  !------------------------------------------------------------------
  subroutine print_state(msg)

    character(len=*) :: msg

    write(*,*)
    write(*,'(A5,A115)') "TEST ", repeat("=",115)
    write(*,'(A5,A20)') "TEST ", msg
    write(*,'(A5,A115)') "TEST ", repeat("=",115)
    write(*,'(A5,A15,5A20)') "TEST ", "Variable", "Min", "Max", "First", "Last", "RMS"
    write(*,'(A5,A115)') "TEST ", repeat("-",115)

    call print_2d_variable("fy2", fy2)
    call print_2d_variable("q",   q)
    call print_2d_variable("cry", cry)
    call print_2d_variable("dya", dya)

    write(*,'(A5,A115)') "TEST ", repeat("-",115)
    write(*,*)

  end subroutine print_state

  !------------------------------------------------------------------
  ! print_2d_variable
  !
  ! prints statistics for a 2d state variable
  !------------------------------------------------------------------
  subroutine print_2d_variable(name, data)

    character(len=*) :: name
    real             :: data(:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we don't have to know start/end indices here
    write(*,'(A5, A15,5E20.6)') "TEST ", name, minval(DBLE(data)), maxval(DBLE(data)), DBLE(data(1,1)), &
                            DBLE(data(size(data,1), size(data,2))),            &
                            sqrt(sum(DBLE(data)**2) / size(data))

  end subroutine print_2d_variable

  !------------------------------------------------------------------
  ! deallocate_state
  !
  ! Deallocate state variables
  !------------------------------------------------------------------
  subroutine deallocate_state

    if (allocated(cry)) then
      deallocate(cry)
    end if
    if (allocated(q)) then
      deallocate(q)
    end if
    if (allocated(dya)) then
      deallocate(dya)
    end if
    if (allocated(fy2)) then
      deallocate(fy2)
    end if

  end subroutine deallocate_state

  !------------------------------------------------------------------
  ! allocate_state
  !
  ! Aallocate state variables
  !------------------------------------------------------------------
  subroutine allocate_state

    call deallocate_state()

    allocate(cry(isd:ied, js:je + 1))
    allocate(q(isd:ied, jsd:jed))
    allocate(dya(isd:ied, jsd:jed))
    allocate(fy2(isd:ied, js:je + 1))

  end subroutine allocate_state

  !------------------------------------------------------------------
  ! read_state
  !
  ! Read state from NetCDF file
  !------------------------------------------------------------------
  subroutine read_state(filename)

    character(len=*), intent(in) :: filename

    ! netCDF variables
    integer :: ncFileID

    ! Local variables
    integer :: nid, njd, njp1

    ! Open file for read only
    call open_file(filename, "r", ncFileID)

    ! Read global attributes
    call read_global_int(ncFileID, "isd", isd)
    call read_global_int(ncFileID, "ied", ied)
    call read_global_int(ncFileID, "jsd", jsd)
    call read_global_int(ncFileID, "jed", jed)
    call read_global_int(ncFileID, "js", js)
    call read_global_int(ncFileID, "je", je)
    call read_global_int(ncFileID, "ord_in", ord_in)
    call read_global_int(ncFileID, "npx", npx)
    call read_global_int(ncFileID, "npy", npy)
    call read_global_int(ncFileID, "grid_type", grid_type)
    call read_global_real(ncFileID, "lim_fac", lim_fac)
    call read_global_logical(ncFileID, "regional", regional)
    call read_global_logical(ncFileID, "nested", nested)

    ! Read the model dimensions
    call read_dimension(ncFileID, "nid", nid)
    call read_dimension(ncFileID, "njd", njd)
    call read_dimension(ncFileID, "njp1", njp1)

    ! Check to make sure state dimensions matches state indices
    if ((nid   /= (ied-isd+1)) .OR. (njd   /= (jed-jsd+1)) .OR. &
        (njp1 /= (je-js+2))) then
      write(*,*) "Dimensions and indices of input data are inconsistent"
      stop 1
    end if

    ! Allocate and initialize the state
    call allocate_state()

    ! Read the state variables
    call read_2d_real(ncFileID, "fy2", fy2)
    call read_2d_real(ncFileID, "q", q)
    call read_2d_real(ncFileID, "cry", cry)
    call read_2d_real(ncFileID, "dya", dya)

    ! Close the NetCDF file
    call close_file(ncFileID)

  end subroutine read_state

  !------------------------------------------------------------------
  ! write_state
  !
  ! Write state to NetCDF file
  !------------------------------------------------------------------
  subroutine write_state(filename)

    character(len=*), intent(in) :: filename

    ! General netCDF variables
    integer :: ncFileID
    integer :: nidDimID, njdDimID, njp1DimID
    integer :: fy2VarID, qVarID, cryVarID, dyaVarID

    ! Local variables
    character(len=8)      :: crdate  ! Needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime  ! Needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone  ! Needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values  ! Needed by F90 DATE_AND_TIME intrinsic
    character(len=19)     :: timestr ! String representation of clock

    ! Open new file, overwriting previous contents
    call open_file(filename, "w", ncFileID)

    ! Write Global Attributes
    call DATE_AND_TIME(crdate,crtime,crzone,values)
    write(timestr,'(i4,2(a,i2.2),1x,i2.2,2(a,i2.2))') &
          values(1), '/', values(2), '/', values(3), values(5), ':', &
          values(6), ':', values(7)
    call write_global_character(ncFileID, "creation_date", timestr)
    call write_global_character(ncFileID, "kernel_name", "yppm")
    call write_global_int(ncFileID, "isd", isd)
    call write_global_int(ncFileID, "ied", ied)
    call write_global_int(ncFileID, "jsd", jsd)
    call write_global_int(ncFileID, "jed", jed)
    call write_global_int(ncFileID, "js", js)
    call write_global_int(ncFileID, "je", je)
    call write_global_int(ncFileID, "ord_in", ord_in)
    call write_global_int(ncFileID, "npx", npx)  
    call write_global_int(ncFileID, "npy", npy)
    call write_global_int(ncFileID, "grid_type", grid_type)
    call write_global_real(ncFileID, "lim_fac", lim_fac)
    call write_global_logical(ncFileID, "regional", regional)
    call write_global_logical(ncFileID, "nested", nested)

    ! Define the i, j dimensions
    call define_dim(ncFileID, "nid", ied-isd+1, nidDimID)
    call define_dim(ncFileID, "njd", jed-jsd+1, njdDimID)

    ! Define the jp1 dimension
    call define_dim(ncFileID, "njp1", je-js+2, njp1DimID)

    ! Define the fields
    call define_var_2d_real(ncFileID, "fy2",  nidDimID, njp1DimID, fy2VarID)
    call define_var_2d_real(ncFileID, "q",    nidDimID, njdDimID,  qVarID)
    call define_var_2d_real(ncFileID, "cry",  nidDimID, njp1DimID, cryVarID)
    call define_var_2d_real(ncFileID, "dya",  nidDimID, njdDimID,  dyaVarID)

    ! Leave define mode so we can fill
    call define_off(ncFileID)

    ! Fill the variables
    call write_var_2d_real(ncFileID, fy2VarID, fy2)
    call write_var_2d_real(ncFileID, qVarID, q)
    call write_var_2d_real(ncFileID, cryVarID, cry)
    call write_var_2d_real(ncFileID, dyaVarID, dya)

    ! Close the NetCDF file
    call close_file(ncFileID)

  end subroutine write_state

end program yppm_driver
