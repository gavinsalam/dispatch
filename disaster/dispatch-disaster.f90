program dispatch_disaster
  use types; use consts_dp
  use dismc_helper; use sub_defs_io
  use parser
  IMPLICIT NONE 

  real(dp) :: ybj
  integer  :: nf

  !-- read things in from the file specified by the option -rundef
  call parse_setup_runs(run)

  !-- useful to have the options here....
  nf  = int_val_opt('-nf', 5)
  ybj = dble_val_opt('-ybj', 0.001_dp)
  call dh_start(nf, ybj)

end program dispatch_disaster
