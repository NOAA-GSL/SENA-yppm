!------------------------------------------------------------------
! yppm_driver
!
! Driver program to run the yppm kernel
!------------------------------------------------------------------
program yppm_driver

    use OMP_LIB
    use yppm_core_mod
#ifdef ENABLE_GPTL
  use gptl
#endif
  
    implicit none
  
    ! Driver variables
    integer            :: nthreads                              ! # of OpenMP threads
    integer            :: narg                                  ! # of command line arguments
    integer*8          :: count_start, count_end, count_rate    ! Timer start/stop
    integer            :: ret                                   ! Return status
    integer, external  :: print_affinity                ! External subroutine
  
    ! Input configuration variables
    character(len=128) :: namelist_file = "test_input/yppm_0.0.1.nl"
    character(len=128) :: input_file, output_file
    integer            :: nl_unit, interpFactor
  
    ! Input namelists
    namelist /io/        input_file, output_file, interpFactor
    namelist /debug/     do_profile                ! Defined in yppm_core_mod
  
    ! Set defaults.  Null strings should cause a quick error.
    input_file = ""
    output_file = ""
    interpFactor = 0
  
    ! Get the number of arguments
    narg = command_argument_count()
    if (narg /= 1) then
      write (*, *) "Usage: yppm <namelist_file>"
      stop 1
    end if
  
    ! Get the namelist file name from the argument list
    call get_command_argument(1, namelist_file)
  
    ! Open the namelist file
    open (newunit=nl_unit, file=TRIM(namelist_file), form='formatted', status='old')
  
    ! Read the data IO settings from the namelist
    read (nl_unit, nml=io)
  
    ! Read the debug settings from the namelist
    read(nl_unit, nml=debug)
  
    ! Get OMP_NUM_THREADS value
    nthreads = omp_get_max_threads()
  
    ! Initialize GPTL if enabled
#ifdef ENABLE_GPTL
  if (do_profile == 1) then
    ret = GPTLinitialize()
  end if
#endif

    ! Print out configuration settings
    write (*, '(A,A)') 'Input file = ', TRIM(input_file)
    write (*, '(A,I0)') 'nthreads = ', nthreads
    ret = print_affinity(0)

    ! Read the input state from the NetCDF input file
    call read_state(TRIM(input_file))
    
    ! Print the input state
    call print_state("Input State - Original")
  
    ! Interpolate the data according to the interpolation factor.
    if (interpFactor > 0) then
      call interpolate_state(interpFactor)
      ! Print the input state
      call print_state("Input State - Interpolated")
    elseif (interpFactor == 0) then
      ! Do nothing.
    else
      print *,"Error, InterpFactor less than zero."
      stop 1
    endif
  
    ! Get the start time
    call system_clock(count_start, count_rate)
  
#ifdef ENABLE_GPTL
  if (do_profile == 1) then
     ret = gptlstart('kernel')
  end if
#endif
  
    ! Run the kernel
    call yppm(fy2, q, cry, ord_in, isd, ied, isd, ied, js, je, jsd, jed, npx, npy, dya, &
              nested, grid_type, lim_fac, regional)
  
#ifdef ENABLE_GPTL
  if (do_profile == 1) then
    ret = gptlstop('kernel')
  end if
#endif
  
    ! Get the stop time
    call system_clock(count_end, count_rate)
  
    ! Print the output state
    call print_state("Output State")
  
    ! Write the output state to the NetCDF output file
    call write_state(TRIM(output_file))
  
    ! Write timing information
    write (*, *)
    write (*, '(A, F12.9)') "Total time (in seconds)=", ((count_end - count_start) * 1.0) / count_rate
  
    ! Deallocate the state variables
    call deallocate_state()
  
    ! Turn off GPTL if enabled
#ifdef ENABLE_GPTL
  if (do_profile == 1) then
    ret = gptlpr(0)
    ret = gptlfinalize()
  end if
#endif
  
end program yppm_driver