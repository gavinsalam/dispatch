!======================================================================
! NB non standard module to avoid clashes with things defines 
!    in DISENT common blocks... Maybe do things in some other way?
module types
  implicit none
  integer, parameter  :: dp = kind(1.0d0), sp = kind(1.0)
  !-- this is useful when need to guarantee type even with double 128
  integer, parameter  :: dp15 = selected_real_kind(15)
end module types
!----------------------------------------------------------------------
module consts_dp
  use types
  implicit none
  private

!!$  real(dp), public, parameter :: pi =&
!!$       & 3.141592653589793238462643383279502884197_dp
  real(dp), public, parameter :: twopi =&
       & 6.283185307179586476925286766559005768394_dp
  real(dp), public, parameter :: half = 0.5_dp, two = 2.0_dp
  real(dp), public, parameter :: zero = 0.0_dp, one = 1.0_dp
  real(dp), public, parameter :: three = 3.0_dp, four = 4.0_dp
end module consts_dp
