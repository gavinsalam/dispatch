! SOME ROUTINES WHICH DO NOTHING, BUT WITH THE SAME NAME AS THEIR PDFLIB
! COUNTERPARTS: THE IDEA BEING TO AVOID HAVING TO LINK WITH PDFLIB

! THEY SHOULD OBVIOUSLY NOT BE CALLED: IF THEY ARE CALLED THEN THEY MUST 
! GIVE AN ERROR MESSAGE


!----- use the following version when we DO NOT want to link into pdflib
module fake_pdflib
  logical, parameter, public :: use_pdflib = .false.
end module fake_pdflib
subroutine pdfset()
  write(0,*) 'ERROR: pdfset is not just a place holder'
  write(0,*) 'do not link fake_pdflib.f if it is not wanted'
end subroutine pdfset
subroutine pftopdg()
  write(0,*) 'ERROR: pftopdg is not just a place holder'
  write(0,*) 'do not link fake_pdflib.f if it is not wanted'
end subroutine pftopdg
subroutine alphas2()
  write(0,*) 'ERROR: alphas2 is not just a place holder'
  write(0,*) 'do not link fake_pdflib.f if it is not wanted'
end subroutine alphas2
      

!!$!----- use the following version when we DO want to link into pdflib
!!$module fake_pdflib
!!$  logical, parameter, public :: use_pdflib = .true.
!!$end module fake_pdflib
