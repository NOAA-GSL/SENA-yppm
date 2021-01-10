module interpolate
  implicit none

contains

  subroutine interpolate_allocate(odims, interpFactor, fdims)
  implicit none

    ! Subroutine for calculating the lower and upper bounds of an
    ! interpolated array with regard to interpFactor.  The user is
    ! expected to call this routine, allocate the interpolated array,
    ! then call subroutine: interpolate_arrray to fill the
    ! interpolated array.
    !
    ! Subroutine parameters:
    ! originalArray, 2D array to be interpolated
    ! odims, 4 elements integer array, lower/upper subscripts of original array
    ! interpFactor, interpolation factor
    ! odims, 4 elements integer array, lower/upper subscripts of interpolated
    ! array
    !
    ! For example:
    ! dimension interpArray(:,:)
    ! call interpolate_allocate(odims,interpFactor,fdims)
    ! l1=fdims(1)
    ! u1=fdims(2)
    ! l2=fdims(3)
    ! l2=fdims(4)
    ! allocate (interpArray(l1:u1,l2:u2))
    ! call interpolate_array(originalArray, odims, interpolatedArray, fdims,
    ! interpFactor)

    integer, intent(in) :: odims(4)
    integer, intent(in) :: interpFactor
    integer, intent(out):: fdims(4)

    fdims(1) = odims(1)
    fdims(3) = odims(3)
    fdims(2) = (odims(1)-1) + (odims(2)-odims(1)+1) + (odims(2)-odims(1))*interpFactor
    fdims(4) = (odims(3)-1) + (odims(4)-odims(3)+1) + (odims(4)-odims(3))*interpFactor

  end subroutine interpolate_allocate

  subroutine interpolate_array(oa, odims, f, fdims, interpFactor)

    ! Bilinear interpolation of f with respect to oa.
    ! See Wikipedia's description of bilinear interpolation.
    !    https://en.wikipedia.org/wiki/Bilinear_interpolation
    !
    ! This subroutine follows the "Alternative algorithm" in the
    ! article refererenced above.
    ! Referring to the square in the article:
    !   Q11 = (x1,y1) ! bottom left
    !   Q12 = (x1,y2) ! top left
    !   Q21 = (x2,y1) ! botton right
    !   Q22 = (x2,y2) ! top right
    !
    ! f's size has been calculated in subroutine interpolate_allocate
    ! f has been allocated in the call tree above this subroutine
    !
    ! Subroutine Parameters:
    ! oa          - original array, 2D
    ! odims       - low/high subscripts for oa
    ! f           - interpolated array, 2D
    ! fdims       - low/high subscripts for f 
    ! interpFactor- interpolation factor

    integer, intent(in)    :: odims(4) ! 
    real,    intent(in)    :: oa(odims(1):odims(2), odims(3):odims(4))
    integer, intent(in)    :: fdims(4)
    real,    intent(inout) :: f(fdims(1):fdims(2), fdims(3):fdims(4))
    integer, intent(in)    :: interpFactor

    ! Here are some illustrations of interpFactor.
    ! We are assuming that the grid points are equally spaced.
    !
    ! An interpolation factor of 1 means 1 new interpolated element.
    ! So for example, a 3x3 matrix with interpFactor=1 becomes a 5x5 matrix
    !
    ! x's are data from oa, the o's are points to be interpolated 
    !
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
    !                   [x o o x o o x]
    !                   [o o o o o o o]
    !                   [o o o o o o o]
    !                   [x o o x o o x]

    ! Locals
    integer :: x1,x2,y1,y2         ! indexes of the interpolated point
    integer :: ix,iy               ! indexes of the interpolated point
    real    :: a0,a1,a2,a3         ! coefficients
    integer :: icount              ! 1D counter
    real    :: fQ11,fQ12,fQ21,fQ22 ! function values at the corners

    integer :: ol1,ou1,ol2,ou2     ! lower and upper bounds of 2D array oa
    integer :: fl1,fu1,fl2,fu2     ! lower and upper bounds of 2d array f
    integer :: i,j                 ! do loop indexes to loop over the squares

    ! Unlock the subscript bounds of the originalArray.
    ol1 = odims(1)
    ou1 = odims(2)
    ol2 = odims(3)
    ou2 = odims(4)

    ! Unlock the subscript bounds of the interpolatedArray.
    fl1 = fdims(1)
    fu1 = fdims(2)
    fl2 = fdims(3)
    fu2 = fdims(4)

    ! Intersperse the oa points into f.
    icount = 0
    do j = fl2, fu2, interpFactor + 1
      do i = fl1, fu1, interpFactor + 1
        f(i, j) = oa(ol1 + modulo(icount, ou1-ol1+1), &
                                 (ol2 + icount/(ou1-ol1+1)))
        icount = icount + 1
      enddo
    enddo

    ! Loop over all the squares.
    do j = fl2, fu2-1, interpFactor + 1
      do i = fl1, fu1-1, interpFactor + 1
        ! Loop over all the interpolated points
        do iy=j,j+interpFactor +1
           do ix=i,i+interpFactor +1
             ! Find the indices of the corner points
             x1 = i
             x2 = i + interpFactor +1
             y1 = j
             y2 = j + interpFactor +1
             ! Skip the corner points
             if (((ix==x1) .and. (iy==y1)) .or. &
                 ((ix==x1) .and. (iy==y2)) .or. &
                 ((ix==x2) .and. (iy==y1)) .or. &
                 ((ix==x2) .and. (iy==y2))) then
               cycle 
             endif
             ! Find the value of the corner points
             fQ11 = f(x1,y1)
             fQ12 = f(x1,y2)
             fQ21 = f(x2,y1)
             fQ22 = f(x2,y2)
             ! Calulate the values of a0,a1,a2,a3.
             a0 = (fQ11*x2*y2)/((x1-x2)*(y1-y2)) + &
                  (fQ12*x2*y1)/((x1-x2)*(y2-y1)) + &
                  (fQ21*x1*y2)/((x1-x2)*(y2-y1)) + &
                  (fQ22*x1*y1)/((x1-x2)*(y1-y2)) 
             a1 = (fQ11*y2)/((x1-x2)*(y2-y1)) + &
                  (fQ12*y1)/((x1-x2)*(y1-y2)) + &
                  (fQ21*y2)/((x1-x2)*(y1-y2)) + &
                  (fQ22*y1)/((x1-x2)*(y2-y1)) 
             a2 = (fQ11*x2)/((x1-x2)*(y2-y1)) + &
                  (fQ12*x2)/((x1-x2)*(y1-y2)) + &
                  (fQ21*x1)/((x1-x2)*(y1-y2)) + &
                  (fQ22*x1)/((x1-x2)*(y2-y1))
             a3 = fQ11/((x1-x2)*(y1-y2)) + &
                  fQ12/((x1-x2)*(y2-y1)) + &
                  fQ21/((x1-x2)*(y2-y1)) + &
                  fQ22/((x1-x2)*(y1-y2)) 
             ! Store the interpolated point.
             f(ix,iy) = a0 + a1*ix + a2*iy + a3*ix*iy
           enddo ! x loop
         enddo   ! y loop
       enddo     ! i loop
     enddo       ! j loop
 
  end subroutine interpolate_array

end module interpolate
 
