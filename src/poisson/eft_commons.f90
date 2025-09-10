module eft_commons
  !use, intrinsic :: iso_fortran_env, only: dp => real64
  !implicit none
  use amr_parameters


  real(dp) :: alphaB, alphaM, C2, C4, Ia, Hv, pv, rv, qv, mu_phi, xiv,nuv,mu_psi,mu_chi
  mu_phi = 1.0d0
  mu_psi = 1.0d0
  mu_chi = 0.0d0

end module eft_commons
