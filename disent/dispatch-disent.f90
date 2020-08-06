








program drivit
  use types; use consts_dp
  use disent_helper; use sub_defs_io
  use parser
  IMPLICIT NONE 

  DOUBLE PRECISION seeds(2)
  real(dp) :: ybj, muF
  integer  :: nf

  !-- unlike original disent there should be no need to reset seeds?
  !   just let it take off where previous one stopped
  !   (old seeds are to be found in other driver programs)

  call parse_setup_runs(run)

  nf  = int_val_opt('-nf', 5)
  ybj = dble_val_opt('-ybj', 0.001_dp)
  !-- this is actually the ratio of mu_F to Q
  !   recall that in the program scale = muF**2
  !muF = dble_val_opt('-muF', one)

  call dh_start(nf, ybj)


END program drivit
