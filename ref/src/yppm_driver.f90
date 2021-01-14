!------------------------------------------------------------------
! yppm_driver
!
! Driver program to run the yppm kernel
!------------------------------------------------------------------
program yppm_driver

  use OMP_LIB
  use yppm_core_mod 
  use ieee_arithmetic

  implicit none

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

end program yppm_driver
