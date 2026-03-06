module eft_module
  !use, intrinsic :: iso_fortran_env, only: dp => real64
  !implicit none
  use amr_parameters
contains

  subroutine compute_eft_quantities()
    !use cosmo_commons
    use eft_commons
    !use amr_parameters
    !integer, parameter :: dp = kind(1.0d0)

    !implicit none

    ! Local variables
    real(dp) :: evt, evprime, D, om_ma, HdotbyH2, abdot
    real(dp) :: tmp1, tmp2
    real(dp) :: integral, h
    integer :: i, npts
    real(dp), parameter :: a_min = 1.0e-4_dp
    integer, parameter :: N = 1000
    real(dp) :: a_vals(N), integrand(N)
    !real(dp), dimension(N) :: a_vals
    !real(dp), dimension(100) :: integrand
    

    !real(dp) :: w0 = -1.0_dp
    !real(dp) :: wa = 0.0_dp

    ! Step 1: E(a)
    evt = aexp**(-3.0_dp * (1.0_dp + w0 + wa)) * exp(-3.0_dp * wa * (1.0_dp - aexp))
    Hv = sqrt(omega_m / aexp**3 + omega_l * evt)

    ! Step 2: Omega_m(a), alphaB, alphaM
    om_ma = omega_m / (omega_m + omega_l * evt * aexp**3)

    ! Step 3: Hdot/H^2
    !HdotbyH2 = -1.5_dp * om_ma
    HdotbyH2 = -1.5_dp * (1.0_dp + (w0 + wa*(1.0_dp-aexp))*(1.0_dp-om_ma))

    ! Step 4: abdot and alphas
    select case (trim(scaling))

    case ('a')
      alphaB = alphaB0 * aexp**nb
      alphaM = alphaM0 * aexp**nm
      abdot = nb * alphaB

    case ('de')

      alphaB = alphaB0 * ((1.0_dp - om_ma)/(1.0_dp - omega_m))**nb
      alphaM = alphaM0 * ((1.0_dp - om_ma)/(1.0_dp - omega_m))**nm
      abdot = -3.0_dp * nb * (w0 + wa*(1.0_dp-aexp)) * om_ma * alphaB

    end select

    ! derivative dEv/da = Ev * d(ln Ev)/da
    ! ln Ev = -3(1+w0+wa) ln a -3 wa(1-a)
    ! d(ln Ev)/da = -3(1+w0+wa)/a + 3 wa
    
    ! evprime = evt * ( -3.0_dp*(1.0_dp + w0 + wa)/aexp + 3.0_dp*wa )

    ! D = omega_m + omega_l * evt * aexp**3

    ! abdot = (alphaB0 / (1.0_dp - omega_m)) * (omega_m * omega_l) * &
    !     ( aexp**4 * evprime + 3.0_dp * aexp**3 * evt ) / (D*D)



    ! Step 5: Ia integral
    ! npts = 100
    ! h = (aexp - a_min) / real(npts - 1, dp)
    ! do i = 1, npts
    !   a_vals = [(a_min + real(i - 1, dp) * h, i = 1, N)]
    !   evt = a_vals(i)**(-3.0_dp * (1.0_dp + w0 + wa)) * exp(-3.0_dp * wa * (1.0_dp - a_vals(i)))
    !   om_ma = omega_m / (omega_m + omega_l * evt * a_vals(i)**3)
    !   integrand(i) = (alphaM0 * (1.0_dp - om_ma) / (1.0_dp - omega_m)) / a_vals(i)
    ! end do
    ! integral = 0.5_dp * (integrand(1) + integrand(npts))
    ! do i = 2, npts - 1
    !   integral = integral + integrand(i)
    ! end do
    ! integral = integral * h
    ! Ia = exp(-1.0_dp*integral)

    if (trim(scaling) == 'a') then


      if (abs(nm) > 1.0e-12_dp) then
          Ia = alphaM0/nm * (aexp**nm - a_min**nm)
      else
          Ia = alphaM0 * log(aexp/a_min)
      end if
      Ia = exp(-1.0_dp * Ia)
    else if (trim(scaling) == 'de') then

      npts = 100
      h = (aexp - a_min) / real(npts - 1, dp)
      do i = 1, npts
        a_vals(i) = a_min + real(i - 1, dp) * h
        evt = a_vals(i)**(-3.0_dp * (1.0_dp + w0 + wa)) * exp(-3.0_dp * wa * (1.0_dp - a_vals(i)))
        om_ma = omega_m / (omega_m + omega_l * evt * a_vals(i)**3)
        integrand(i) = alphaM0 * ((1.0_dp - om_ma)/(1.0_dp - omega_m))**nm / a_vals(i)
      end do
      integral = 0.5_dp * (integrand(1) + integrand(npts))
      do i = 2, npts - 1
        integral = integral + integrand(i)
      end do
      integral = integral * h
      Ia = exp(-1.0_dp*integral)

    end if

    ! Step 6: C2 and C4
    C2 = -alphaM + alphaB * (1.0_dp + alphaM) + (1.0_dp + alphaB) * HdotbyH2 + abdot + 1.5_dp * Ia * omega_m / (aexp**3 * Hv**2)
    !C2 = -alphaM + alphaB*(1.0_dp + alphaM) + (1.0_dp + alphaB)*HdotbyH2 &
    ! + (3.0_dp*aexp**3 * alphaB0 * omega_m) / (aexp**3 * (1.0_dp - omega_m) + omega_m)**2 &
    ! + aexp**(-3.0_dp) * 1.5_dp * Ia * omega_m / (Hv**2)
    C4 = -4.0_dp * alphaB + 2.0_dp * alphaM

    xiv = alphaB - alphaM
    nuv = -1*C2 -alphaB*(xiv - alphaM)

    mu_phi = 1 + xiv*xiv/nuv
    mu_psi = 1 + alphaB*xiv/nuv
    mu_chi = xiv/nuv

    pv = alphaB * (2.0d0 * alphaM - alphaB) - C2
    qv = -1.0d0*(alphaM - alphaB)
    rv = -0.25d0 * C4 / (aexp * aexp * Hv)**2

  end subroutine compute_eft_quantities

end module eft_module
