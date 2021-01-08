  !------------------------------------------------------------------
  ! interpolateTest
  !
  ! Unit test to test interpolation software.  The origArray is
  ! expanded by interpFactor to interpolateArray.  The origArray is
  ! filled with 1.0(s), interpolation should create a larger
  ! array, also filled with 1.0(s) if successful.
  !------------------------------------------------------------------

program interpolateTest

  use interpolate

  implicit none

  integer :: i, j, k ! do loop indexes
  integer :: interpFactor ! also a do loop index
  integer :: odims(4), idims(4) ! lower/upper bounds for 2D arrays
  real, allocatable :: origArray(:, :) ! original array
  real, allocatable :: interpArray(:, :) ! interpolated array
  real :: epsilon
  parameter(epsilon=1e-5)

  do i = 3, 7 ! create these 2D arrays
    do interpFactor = 1, 5 ! create these interpFactor(s)

      ! Create the [i,i+i] 2D array.
      if (allocated(origArray)) then
        deallocate (origArray)
      endif
      allocate (origArray(i, i + i))
      odims(1) = 1
      odims(2) = i
      odims(3) = 1
      odims(4) = i + i ! non-sqare matrix

      ! Populate the originalArray with 1.0
      origArray = 1.0

      ! Create the interpolateArray.
      call interpolate_allocate(odims, interpFactor, idims)
      if (allocated(interpArray)) then
        deallocate (interpArray)
      endif
      allocate (interpArray(idims(1):idims(2), idims(3):idims(4)))
      interpArray = 2.0 ! populate the interpArray with 2.0 (2.0 != 1.0)

      ! Interpolate the interpolateArray.
      call interpolate_array(origArray, odims, interpArray, idims, interpFactor)

      ! The inperpolateArray should be all 1.0, +/- epsilon.
      ! Print success or failure.  Exit with non-zero error code if falure.
      do j = idims(3), idims(4)
        do k = idims(1), idims(2)
          if (abs(interpArray(k, j) - 1.0) .gt. epsilon) then
            print *, "Failure"
            print *, k, j, interpArray(k, j), interpArray(k, j) - 1.0
            stop -10
          endif
        enddo
      enddo

    enddo
  enddo

  print *, "Success"
  stop 0

end program interpolateTest
