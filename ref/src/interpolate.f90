module interpolate
  implicit none

contains

subroutine interpolate_allocate(odims,interpFactor,idims)

  ! Subroutine for calculating the lower and upper bounds of an
  ! interpolated array with regard to interpFactor.  The user is
  ! expoected to call this routine, allocate the interpolated array,
  ! then call subroutine: interpolate_arrray to fill the 
  ! interpolated array.
  !
  ! Subroutine parameters:
  ! originalArray, 2D array to be interpolated
  ! odims, 4 elements integer array, lower/upper subscripts of original array
  ! interpFactor, interpolation factor
  ! odims, 4 elements integer array, lower/upper subscripts of interpolated array
  !
  ! For example:
  ! dimension interpArray(:,:)
  ! call interpolate_allocate(odims,interpFactor,idims)
  ! l1=idims(1)
  ! u1=idims(2)
  ! l2=idims(3)
  ! l2=idims(4)
  ! allocate (interpArray(l1:u1,l2:u2))
  ! call interpolate_array(originalArray, odims, interpolatedArray, idims)

  integer, intent(in) :: odims(4)
  integer, intent(in) :: interpFactor
  integer, intent(out):: idims(4)

  idims(1) = odims(1)
  idims(3) = odims(3)
  idims(2) = (odims(1)-1) + (odims(2)-odims(1)+1) + (odims(2)-odims(1)) * interpFactor
  idims(4) = (odims(3)-1) + (odims(4)-odims(3)+1) + (odims(4)-odims(3)) * interpFactor

end subroutine interpolate_allocate

subroutine interpolate_array(originalArray, odims, interpolatedArray, idims, interpFactor)

  ! Subroutine to interpolate interpolatedArray with regard to
  ! interpolatedArray and interpFactor.
  !
  ! We'll assume that the interpolation factor is greater or equal to one.
  ! An interpolation factor of 1 means 1 new interpolated element.
  ! So for example, a 3x3 matrix with interpFactor=1 becomes a 5x5 matrix
  ! [ x x x ]     ==> [x o x o x]
  ! [ x x x ]         [o o o o o]
  ! [ x x x ]         [x o x o x]
  !                   [o o o o o]
  !                   [x o x o x]
  !
  ! And an interpFactor=2, means a 3x3 matrix becomes a 7x7 matrix.
  ! [ x x x ]     ==> [x o o x o o x]
  ! [ x x x ]         [o o o o o o o]
  ! [ x x x ]         [o o o o o o o]
  !                   [x o o x o o o]
  !                   [o o o o o o o]
  !                   [o o o o o o o]
  !                   [x o o x o o x]
  ! We'll be interpolating on the horizonal, then inperpolating on 
  ! the vertical.  That way, all the "o"s are interpolated.
  !
  ! Subroutine parameters:
  ! originalArray, 2D array to interpolate
  ! interpolatedArray, resultant 2D array
  ! interpFactor, integer, number of points to add between original points
 
  integer, intent(in) :: odims(4)
  real, intent(in)    :: originalArray(odims(1):odims(2),odims(3):odims(4))
  integer, intent(in) :: idims(4)
  real, intent(inout) :: interpolatedArray(idims(1):idims(2),idims(3):idims(4))
  integer, intent(in) :: interpFactor

  ! Locals:
  integer :: icount
  integer :: i,j,k ! loop indexes
  integer :: ol1,ol2,ou1,ou2
  integer :: il1,iu1,il2,iu2
  real    :: w1, w2
 
  ! Unlock the subscript bounds of the originalArray.
  ol1 = odims(1)
  ou1 = odims(2) 
  ol2 = odims(3)
  ou2 = odims(4)

  ! Unlock the subscript bounds of the interpolatedArray.
  il1 = idims(1)
  iu1 = idims(2) 
  il2 = idims(3)
  iu2 = idims(4)

  ! Intersperse the originalArray points into the interpolatedArray.
  icount = 0
  do j=il2,iu2,interpFactor+1
     do i= il1,iu1,interpFactor+1
        interpolatedArray(i,j) = originalArray(ol1+modulo(icount,ou1-ol1+1), &
                                              (ol2+    icount/(ou1-ol1+1)))
        icount = icount+1
     enddo
  enddo

  ! Horizontal interpolation.  This interpolates the "o"s in the rows with "x"s.
  do j=il2,iu2,interpFactor+1
     do i=il1,iu1-1,interpFactor+1 ! -1 to skip the last iteration
        do k=1,interpFactor
           w1 = (interpFactor-k+1)/dble(interpFactor+1)
           w2 = k/dble(interpFactor + 1)
           interpolatedArray(i+k,j) = interpolatedArray(i,j) * w1 + &
                                      interpolatedArray(i+interpFactor+1,j) * w2
        enddo
     enddo
  enddo

  ! Vertical interpolation.  This interpolates all the "o"s in the columns.
  do j=il2,iu2-1,interpFactor+1 ! -1 to skip the last iteration
     do i=il1,iu1
        do k=1,interpFactor
           w1 = (interpFactor-k+1)/dble(interpFactor+1)
           w2 = k/dble(interpFactor + 1)
           interpolatedArray(i,j+k) = interpolatedArray(i,j) * w1  + &
                                      interpolatedArray(i,j+interpFactor+1) * w2
        enddo
     enddo
  enddo

  end subroutine interpolate_array

end module interpolate
