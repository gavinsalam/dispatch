
!----------------------------------------------------------------------
! Interfaces to the dec f90 iargc and getarg routines, Known to work
! also with lahey lf95 and pgf90. 
!
! Routines for other compilers may be available on request 
! from gavin.salam@cern.ch 
!
! GPS 4/11/95 (CCN8 9)
! $Id: lcl_args.f90,v 1.2 2001/10/19 21:37:38 gsalam Exp $
!----------------------------------------------------------------------
integer function lcl_iargc()
      implicit none
      integer iargc

      lcl_iargc = iargc()
end function lcl_iargc

subroutine lcl_getarg(k, argk)
      implicit none
      integer,      intent(in)  :: k
      character(*), intent(out) :: argk

      call getarg(k, argk)
end subroutine lcl_getarg


subroutine lcl_flush(idev)
  implicit none
  integer, intent(in) :: idev
  
  call flush(idev)
end subroutine lcl_flush

! this causes difficulties sometimes
!subroutine lcl_system(string)
!      implicit none
!      character(*), intent(in) :: string
!      !------------------------------------------------------------
!      !integer return_val, system
!
!      call  system(string)
!end subroutine lcl_system

      
