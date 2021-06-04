!***********************************************************************
!*
!* This file contains derivative work from the Geophysical Fluid Dynamics
!* Laboratory GitHub project for the FV3 dynamical core
!* (https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere). It is a stand alone
!* subset of the original specifically licensed for use under the
!* Apache License (https://www.apache.org/licenses/LICENSE-2.0).
!* The application of the Apache License to this standalone work does not imply
!* nor shall it be construed as any change to the original source license.
!*
!***********************************************************************

module yppm_core_mod

  use netCDFModule
  use interpolate
#ifdef ENABLE_GPTL
 use gptl
#endif

  implicit none

  private 
  public yppm, print_state, deallocate_state, read_state, write_state, interpolate_state
  public js, je, isd, ied, jsd, jed, npx, npy, grid_type, ord_in, lim_fac, nested, regional, cry, q, fy2, dya, do_profile

  integer :: do_profile = 0  ! Flag for enabling profiling at runtime

  real, parameter :: ppm_fac = 1.5   !< nonlinear scheme limiter: between 1 and 2
  real, parameter :: r3 = 1. / 3.
  real, parameter :: near_zero = 1.E-25
  ! scheme 2.1: perturbation form

  real, parameter :: s11 = 11. / 14.
  real, parameter :: s14 =  4. /  7.
  real, parameter :: s15 =  3. / 14.

  !----------------------------------------------------
  ! volume-conserving cubic with 2nd drv=0 at end point:
  !----------------------------------------------------
  ! Non-monotonic
  real, parameter :: c1 = -2. / 14.
  real, parameter :: c2 = 11. / 14.
  real, parameter :: c3 =  5. / 14.

  !----------------------
  ! PPM volume mean form:
  !----------------------
  real, parameter :: p1 = 7. / 12.     ! 0.58333333
  real, parameter :: p2 = -1. / 12.

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

contains

  !------------------------------------------------------------------
  ! yppm
  !
  ! Top-level routine for yppm kernel
  !------------------------------------------------------------------
  subroutine yppm(flux, q, c, jord, ifirst, ilast, isd, ied, js, je, jsd, jed, npx, npy, dya, nested, grid_type, lim_fac, regional)

    real,    intent(out) :: flux(ifirst:ilast, js:je + 1) !<  Flux
    real,    intent( in) :: q(ifirst:ilast, jsd:jed)
    real,    intent( in) :: c(isd:ied, js:je + 1)         !< Courant number
    integer, intent( in) :: jord
    integer, intent( in) :: ifirst, ilast                 !< Compute domain
    integer, intent( in) :: isd, ied, js, je, jsd, jed
    integer, intent( in) :: npx, npy
    real,    intent( in) :: dya(isd:ied, jsd:jed)
    logical, intent( in) :: nested
    integer, intent( in) :: grid_type
    real,    intent( in) :: lim_fac
    logical, intent( in) :: regional

    ! Local:
    real    :: dm(ifirst:ilast, js - 2:je + 2)
    real    :: al(ifirst:ilast, js - 1:je + 2)
    real    :: bl(ifirst:ilast, js - 1:je + 1)
    real    :: br(ifirst:ilast, js - 1:je + 1)
    real    :: b0(ifirst:ilast, js - 1:je + 1)
    real    :: dq(ifirst:ilast, js - 3:je + 2)
    real    :: fx0(ifirst:ilast)
    real    :: fx1(ifirst:ilast)
    real    :: xt1(ifirst:ilast)
    logical :: smt5(ifirst:ilast, js - 1:je + 1)
    logical :: smt6(ifirst:ilast, js - 1:je + 1)
    logical :: hi5(ifirst:ilast)
    logical :: hi6(ifirst:ilast)
    real    :: x0, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2, r1
    integer :: i, j, js1, je3, je1, mord
    integer :: ret

#ifdef ENABLE_GPTL
   if (do_profile == 1) then
     ret = gptlstart('yppm')
   end if
#endif

    if (.not. (nested .or. regional) .and. grid_type < 3) then
      ! Cubed-sphere:
      js1 = max(3, js - 1)
      je3 = min(npy - 2, je + 2)
      je1 = min(npy - 3, je + 1)
    else
      ! Nested grid OR Doubly periodic domain:
      js1 = js - 1
      je3 = je + 2
      je1 = je + 1
    endif

    mord = abs(jord)

    if (jord < 8) then

      do j = js1, je3
        do i = ifirst, ilast
          al(i, j) = p1 * (q(i, j - 1) + q(i, j)) + p2 * (q(i, j - 2) + q(i, j + 1))
        enddo
      enddo

      if (.not. (nested .or. regional) .and. grid_type < 3) then
        if (js == 1) then
          do i = ifirst, ilast
            al(i, 0) = c1 * q(i, -2) + c2 * q(i, -1) + c3 * q(i, 0)
            al(i, 1) = 0.5 * (((2. * dya(i, 0) + dya(i, -1)) * q(i, 0) - dya(i, 0) * q(i, -1)) / (dya(i, -1) + dya(i, 0)) &
                            + ((2. * dya(i, 1) + dya(i, 2))  * q(i, 1) - dya(i, 1) * q(i,  2)) / (dya(i,  1) + dya(i, 2)))
            al(i, 2) = c3 * q(i, 1) + c2 * q(i, 2) + c1 * q(i, 3)
          enddo
        endif
        if ((je + 1) == npy) then
          do i = ifirst, ilast
            al(i, npy - 1) = c1 * q(i, npy - 3) + c2 * q(i, npy - 2) + c3 * q(i, npy - 1)
            al(i, npy) = 0.5 * (((2. * dya(i, npy - 1) + dya(i, npy - 2)) * q(i, npy - 1) - dya(i, npy - 1) * &
                               q(i, npy - 2)) / (dya(i, npy - 2) + dya(i, npy - 1)) + &
                                ((2. * dya(i, npy) + dya(i, npy + 1)) * q(i, npy) -   &
                               dya(i, npy) * q(i, npy + 1)) / (dya(i, npy) + dya(i, npy + 1)))
            al(i, npy + 1) = c3 * q(i, npy) + c2 * q(i, npy + 1) + c1 * q(i, npy + 2)
          enddo
        endif
      endif

      if (jord < 0) then
        do j = js - 1, je + 2
          do i = ifirst, ilast
            al(i, j) = max(0., al(i, j))
          enddo
        enddo
      endif

      if (mord == 1) then
        do j = js - 1, je + 1
          do i = ifirst, ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j + 1) - q(i, j)
            b0(i, j) = bl(i, j) + br(i, j)
            smt5(i, j) = abs(lim_fac * b0(i, j)) < abs(bl(i, j) - br(i, j))
          enddo
        enddo
        do j = js, je + 1
!DEC$ VECTOR ALWAYS
          do i = ifirst, ilast
            if (c(i, j) > 0.) then
              fx1(i) = (1. - c(i, j)) * (br(i, j - 1) - c(i, j) * b0(i, j - 1))
              flux(i, j) = q(i, j - 1)
            else
              fx1(i) = (1. + c(i, j)) * (bl(i, j) + c(i, j) * b0(i, j))
              flux(i, j) = q(i, j)
            endif
            if (smt5(i, j - 1) .or. smt5(i, j)) flux(i, j) = flux(i, j) + fx1(i)
          enddo
        enddo

      elseif (mord == 2) then   ! Perfectly linear scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7

        do j = js, je + 1
!DEC$ VECTOR ALWAYS
          do i = ifirst, ilast
            xt = c(i, j)
            if (xt > 0.) then
              qtmp = q(i, j - 1)
              flux(i, j) = qtmp + (1. - xt) * (al(i, j) - qtmp - xt * (al(i, j - 1) + al(i, j) - (qtmp + qtmp)))
            else
              qtmp = q(i, j)
              flux(i, j) = qtmp + (1. + xt) * (al(i, j) - qtmp + xt * (al(i, j) + al(i, j + 1) - (qtmp + qtmp)))
            endif
          enddo
        enddo

      elseif (mord == 3) then

        do j = js - 1, je + 1
          do i = ifirst, ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j + 1) - q(i, j)
            b0(i, j) = bl(i, j) + br(i, j)
            x0 = abs(b0(i, j))
            xt = abs(bl(i, j) - br(i, j))
            smt5(i, j) = x0 < xt
            smt6(i, j) = 3. * x0 < xt
          enddo
        enddo
        do j = js, je + 1
          do i = ifirst, ilast
            fx1(i) = 0.
            xt1(i) = c(i, j)
            hi5(i) = smt5(i, j - 1) .and. smt5(i, j)
            hi6(i) = smt6(i, j - 1) .or. smt6(i, j)
          enddo
          do i = ifirst, ilast
            if (xt1(i) > 0.) then
              if (hi6(i)) then
                fx1(i) = br(i, j - 1) - xt1(i) * b0(i, j - 1)
              elseif (hi5(i)) then ! both up-downwind sides are noisy; 2nd order, piece-wise linear
                fx1(i) = sign(min(abs(bl(i, j - 1)), abs(br(i, j - 1))), br(i, j - 1))
              endif
              flux(i, j) = q(i, j - 1) + (1. - xt1(i)) * fx1(i)
            else
              if (hi6(i)) then
                fx1(i) = bl(i, j) + xt1(i) * b0(i, j)
              elseif (hi5(i)) then ! both up-downwind sides are noisy; 2nd order, piece-wise linear
                fx1(i) = sign(min(abs(bl(i, j)), abs(br(i, j))), bl(i, j))
              endif
              flux(i, j) = q(i, j) + (1. + xt1(i)) * fx1(i)
            endif
          enddo
        enddo

      elseif (mord == 4) then

        do j = js - 1, je + 1
          do i = ifirst, ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j + 1) - q(i, j)
            b0(i, j) = bl(i, j) + br(i, j)
            x0 = abs(b0(i, j))
            xt = abs(bl(i, j) - br(i, j))
            smt5(i, j) = x0 < xt
            smt6(i, j) = 3. * x0 < xt
          enddo
        enddo
        do j = js, je + 1
          do i = ifirst, ilast
            xt1(i) = c(i, j)
            hi5(i) = smt5(i, j - 1) .and. smt5(i, j)
            hi6(i) = smt6(i, j - 1) .or. smt6(i, j)
            hi5(i) = hi5(i) .or. hi6(i)
          enddo
!DEC$ VECTOR ALWAYS
          do i = ifirst, ilast
            if (xt1(i) > 0.) then
              fx1(i) = (1. - xt1(i)) * (br(i, j - 1) - xt1(i) * b0(i, j - 1))
              flux(i, j) = q(i, j - 1)
            else
              fx1(i) = (1. + xt1(i)) * (bl(i, j) + xt1(i) * b0(i, j))
              flux(i, j) = q(i, j)
            endif
            if (hi5(i)) flux(i, j) = flux(i, j) + fx1(i)
          enddo
        enddo

      else  ! mord=5,6,7
        if (mord == 5) then
          do j = js - 1, je + 1
            do i = ifirst, ilast
              bl(i, j) = al(i, j) - q(i, j)
              br(i, j) = al(i, j + 1) - q(i, j)
              b0(i, j) = bl(i, j) + br(i, j)
              smt5(i, j) = bl(i, j) * br(i, j) < 0.
            enddo
          enddo
        else
          do j = js - 1, je + 1
            do i = ifirst, ilast
              bl(i, j) = al(i, j) - q(i, j)
              br(i, j) = al(i, j + 1) - q(i, j)
              b0(i, j) = bl(i, j) + br(i, j)
              smt5(i, j) = 3. * abs(b0(i, j)) < abs(bl(i, j) - br(i, j))
            enddo
          enddo
        endif

        do j = js, je + 1
!DEC$ VECTOR ALWAYS
          do i = ifirst, ilast
            if (c(i, j) > 0.) then
              fx1(i) = (1. - c(i, j)) * (br(i, j - 1) - c(i, j) * b0(i, j - 1))
              flux(i, j) = q(i, j - 1)
            else
              fx1(i) = (1. + c(i, j)) * (bl(i, j) + c(i, j) * b0(i, j))
              flux(i, j) = q(i, j)
            endif
            if (smt5(i, j - 1) .or. smt5(i, j)) flux(i, j) = flux(i, j) + fx1(i)
          enddo
        enddo

      endif

      return

    else
! Monotonic constraints:
! ord = 8: PPM with Lin's PPM fast monotone constraint
! ord > 8: PPM with Lin's modification of Huynh 2nd constraint

      do j = js - 2, je + 2
        do i = ifirst, ilast
          xt = 0.25 * (q(i, j + 1) - q(i, j - 1))
          dm(i, j) = sign(min(abs(xt), max(q(i, j - 1), q(i, j), q(i, j + 1)) - q(i, j), &
                              q(i, j) - min(q(i, j - 1), q(i, j), q(i, j + 1))), xt)
        enddo
      enddo
      do j = js1, je1 + 1
        do i = ifirst, ilast
          al(i, j) = 0.5 * (q(i, j - 1) + q(i, j)) + r3 * (dm(i, j - 1) - dm(i, j))
        enddo
      enddo

      if (jord == 8) then
        do j = js1, je1
          do i = ifirst, ilast
            xt = 2. * dm(i, j)
            bl(i, j) = -sign(min(abs(xt), abs(al(i, j) - q(i, j))), xt)
            br(i, j) = sign(min(abs(xt), abs(al(i, j + 1) - q(i, j))), xt)
          enddo
        enddo
      elseif (jord == 11) then
        do j = js1, je1
          do i = ifirst, ilast
            xt = ppm_fac * dm(i, j)
            bl(i, j) = -sign(min(abs(xt), abs(al(i, j) - q(i, j))), xt)
            br(i, j) = sign(min(abs(xt), abs(al(i, j + 1) - q(i, j))), xt)
          enddo
        enddo
      else
        do j = js1 - 2, je1 + 1
          do i = ifirst, ilast
            dq(i, j) = 2. * (q(i, j + 1) - q(i, j))
          enddo
        enddo
        do j = js1, je1
          do i = ifirst, ilast
            bl(i, j) = al(i, j) - q(i, j)
            br(i, j) = al(i, j + 1) - q(i, j)
            if (abs(dm(i, j - 1)) + abs(dm(i, j)) + abs(dm(i, j + 1)) < near_zero) then
              bl(i, j) = 0.
              br(i, j) = 0.
            elseif (abs(3. * (bl(i, j) + br(i, j))) > abs(bl(i, j) - br(i, j))) then
              pmp_2 = dq(i, j - 1)
              lac_2 = pmp_2 - 0.75 * dq(i, j - 2)
              br(i, j) = min(max(0., pmp_2, lac_2), max(br(i, j), min(0., pmp_2, lac_2)))
              pmp_1 = -dq(i, j)
              lac_1 = pmp_1 + 0.75 * dq(i, j + 1)
              bl(i, j) = min(max(0., pmp_1, lac_1), max(bl(i, j), min(0., pmp_1, lac_1)))
            endif
          enddo
        enddo
      endif
      if (jord == 9 .or. jord == 13) then
        ! Positive definite constraint:
        do j = js1, je1
          call pert_ppm(ilast - ifirst + 1, q(ifirst, j), bl(ifirst, j), br(ifirst, j), 0)
        enddo
      endif

      if (.not. (nested .or. regional) .and. grid_type < 3) then
        if (js == 1) then
          do i = ifirst, ilast
            bl(i, 0) = s14 * dm(i, -1) + s11 * (q(i, -1) - q(i, 0))

            xt = 0.5 * (((2. * dya(i, 0) + dya(i, -1)) * q(i, 0) - dya(i, 0) * q(i, -1)) / (dya(i, -1) + dya(i, 0)) &
                      + ((2. * dya(i, 1) + dya(i,  2)) * q(i, 1) - dya(i, 1) * q(i,  2)) / (dya(i,  1) + dya(i, 2)))
!            if ( jord==8 .or. jord==10 ) then
            xt = max(xt, min(q(i, -1), q(i, 0), q(i, 1), q(i, 2)))
            xt = min(xt, max(q(i, -1), q(i, 0), q(i, 1), q(i, 2)))
!            endif
            br(i, 0) = xt - q(i, 0)
            bl(i, 1) = xt - q(i, 1)

            xt = s15 * q(i, 1) + s11 * q(i, 2) - s14 * dm(i, 2)
            br(i, 1) = xt - q(i, 1)
            bl(i, 2) = xt - q(i, 2)

            br(i, 2) = al(i, 3) - q(i, 2)
          enddo
          call pert_ppm(3 * (ilast - ifirst + 1), q(ifirst, 0), bl(ifirst, 0), br(ifirst, 0), 1)
        endif
        if ((je + 1) == npy) then
          do i = ifirst, ilast
            bl(i, npy - 2) = al(i, npy - 2) - q(i, npy - 2)

            xt = s15 * q(i, npy - 1) + s11 * q(i, npy - 2) + s14 * dm(i, npy - 2)
            br(i, npy - 2) = xt - q(i, npy - 2)
            bl(i, npy - 1) = xt - q(i, npy - 1)

            xt = 0.5 * (((2. * dya(i, npy - 1) + dya(i, npy - 2)) * q(i, npy - 1) - dya(i, npy - 1) * &
                       q(i, npy - 2)) / (dya(i, npy - 2) + dya(i, npy - 1)) + ((2. * dya(i, npy) +    &
                       dya(i, npy + 1)) * q(i, npy) - dya(i, npy) * q(i, npy + 1)) / (dya(i, npy) +   &
                       dya(i, npy + 1)))

!            if ( jord==8 .or. jord==10 ) then
            xt = max(xt, min(q(i, npy - 2), q(i, npy - 1), q(i, npy), q(i, npy + 1)))
            xt = min(xt, max(q(i, npy - 2), q(i, npy - 1), q(i, npy), q(i, npy + 1)))
!            endif

            br(i, npy - 1) = xt - q(i, npy - 1)
            bl(i, npy) = xt - q(i, npy)
            br(i, npy) = s11 * (q(i, npy + 1) - q(i, npy)) - s14 * dm(i, npy + 1)
          enddo
          call pert_ppm(3 * (ilast - ifirst + 1), q(ifirst, npy - 2), bl(ifirst, npy - 2), &
                        br(ifirst, npy - 2), 1)
        endif
      end if

    endif

    do j = js, je + 1
      do i = ifirst, ilast
        if (c(i, j) > 0.) then
          flux(i, j) = q(i, j - 1) + (1. - c(i, j)) * (br(i, j - 1) - c(i, j) * (bl(i, j - 1) + br(i, j - 1)))
        else
          flux(i, j) = q(i, j) + (1. + c(i, j)) * (bl(i, j) + c(i, j) * (bl(i, j) + br(i, j)))
        endif
      enddo
    enddo

#ifdef ENABLE_GPTL
 if (do_profile == 1) then
   ret = gptlstop('yppm')
 end if
#endif

  end subroutine yppm

  !------------------------------------------------------------------
  ! pert_ppm
  !
  !------------------------------------------------------------------
  subroutine pert_ppm(im, a0, al, ar, iv)

    integer, intent(   in) :: im
    integer, intent(   in) :: iv
    real,    intent(   in) :: a0(im)
    real,    intent(inout) :: al(im), ar(im)

    ! Local:
    real            :: a4, da1, da2, a6da, fmin
    integer         :: i
    real, parameter :: r12 = 1. / 12.
    integer         :: ret 

#ifdef ENABLE_GPTL
 if (do_profile == 1) then
   ret = gptlstart('pert_ppm')
 end if
#endif

    !-----------------------------------
    ! Optimized PPM in perturbation form:
    !-----------------------------------

    if (iv == 0) then
      ! Positive definite constraint
      do i = 1, im
        if (a0(i) <= 0.) then
          al(i) = 0.
          ar(i) = 0.
        else
          a4 = -3. * (ar(i) + al(i))
          da1 = ar(i) - al(i)
          if (abs(da1) < -a4) then
            fmin = a0(i) + 0.25 / a4 * da1**2 + a4 * r12
            if (fmin < 0.) then
              if (ar(i) > 0. .and. al(i) > 0.) then
                ar(i) = 0.
                al(i) = 0.
              elseif (da1 > 0.) then
                ar(i) = -2. * al(i)
              else
                al(i) = -2. * ar(i)
              endif
            endif
          endif
        endif
      enddo
    else
      ! Standard PPM constraint
      do i = 1, im
        if (al(i) * ar(i) < 0.) then
          da1 = al(i) - ar(i)
          da2 = da1**2
          a6da = 3. * (al(i) + ar(i)) * da1
          ! abs(a6da) > da2 --> 3. * abs(al + ar) > abs(al - ar)
          if (a6da < -da2) then
            ar(i) = -2. * al(i)
          elseif (a6da > da2) then
            al(i) = -2. * ar(i)
          endif
        else
          ! effect of dm=0 included here
          al(i) = 0.
          ar(i) = 0.
        endif
      enddo
    endif

#ifdef ENABLE_GPTL
 if (do_profile == 1) then
   ret = gptlstop('pert_ppm')
 end if
#endif

  end subroutine pert_ppm

  !------------------------------------------------------------------
  ! print_state
  !
  ! prints statistics for the kernel state variables
  !------------------------------------------------------------------
  subroutine print_state(msg)

    character(len=*) :: msg

    write(*,*)
    write(*,'(A5,A115)') "TEST ", repeat("=",115)
    write(*,'(A5,A30)') "TEST ", msg
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
    write(*,'(A5, A15,5ES20.5)') "TEST ", name, minval(DBLE(data)), maxval(DBLE(data)), DBLE(data(1,1)), &
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

  !------------------------------------------------------------------
  ! interpolate_state
  !
  ! Increase the database by interpFactor.
  !------------------------------------------------------------------
  subroutine interpolate_state(interpFactor)

    integer, intent(in) :: interpFactor

    ! locals
    ! Define new variables to replace the old variables.
    ! Once the new variables are defined, they will replace the
    ! old variables.

    integer new_isd
    integer new_ied
    integer new_jsd
    integer new_jed
    integer new_js
    integer new_je
    real, allocatable :: new_cry(:, :)
    real, allocatable :: new_q(:, :)
    real, allocatable :: new_fy2(:, :)
    real, allocatable :: new_dya(:, :)
    integer :: odims(4), idims(4)

    odims(1) = isd
    odims(2) = ied
    odims(3) = js
    odims(4) = je + 1
    call interpolate_CalculateSpace(odims, interpFactor, idims)
    allocate (new_cry(idims(1):idims(2), idims(3):idims(4)))
    call interpolate_array(cry, odims, new_cry, idims, interpFactor)
    allocate (new_fy2(idims(1):idims(2), idims(3):idims(4)))
    call interpolate_array(fy2, odims, new_fy2, idims, interpFactor)
    new_js = idims(3)
    new_je = idims(4)

    odims(1) = isd
    odims(2) = ied
    odims(3) = jsd
    odims(4) = jed
    call interpolate_CalculateSpace(odims, interpFactor, idims)
    allocate (new_q(idims(1):idims(2), idims(3):idims(4)))
    call interpolate_array(q, odims, new_q, idims, interpFactor)
    allocate (new_dya(idims(1):idims(2), idims(3):idims(4)))
    call interpolate_array(dya, odims, new_dya, idims, interpFactor)
    new_isd = idims(1)
    new_ied = idims(2)
    new_jsd = idims(3)
    new_jed = idims(4)

    ! Exchange the old dimensions to the new dimensions.
    ied = new_ied
    jed = new_jed
    je = new_je - 1

    ! Switch the old arrays to the new arrays.
    call allocate_state()
    cry = new_cry
    q = new_q
    dya = new_dya
    fy2 = new_fy2
    deallocate (new_cry); deallocate (new_q); deallocate (new_dya); deallocate (new_fy2)

  end subroutine interpolate_state

end module yppm_core_mod
