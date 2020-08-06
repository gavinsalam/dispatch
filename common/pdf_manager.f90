!======================================================================
!
! This module is designed to manage a set of PDFs, and a set of
! x,Q values
!
! The run_def will contain a pointer to a pdf_xQ_set structure, that 
! way several run_defs can in prinicple use the same pdf_xQ_set structure.
!
!
!
module pdf_manager
  use types; use consts_dp
  implicit none
  private

  !======================================================================
  ! THESE MUST NOT BE CHANGED. THEY ARE CODES WHICH ARE CARRIED
  ! ACROSS TO DISRESUM AND THEY WILL BE USED FOR COMMUNICATION PURPOSES?
  integer, parameter :: pdf_base = 100
  integer, parameter, public :: pdf_ZERO = 0

  integer, parameter, public :: pdf_CTEQ5 = 1 * pdf_base
  integer, parameter, public :: pdf_CTEQ5M  = 1 + pdf_CTEQ5
  integer, parameter, public :: pdf_CTEQ5D  = 2 + pdf_CTEQ5
  integer, parameter, public :: pdf_CTEQ5L  = 3 + pdf_CTEQ5
  integer, parameter, public :: pdf_CTEQ5HJ = 4 + pdf_CTEQ5
  integer, parameter, public :: pdf_CTEQ5M1 = 8 + pdf_CTEQ5

  
  integer, parameter, public :: pdf_GRV = 2 * pdf_base
  integer, parameter, public :: pdf_GRV98LO  = 1 + pdf_GRV
  integer, parameter, public :: pdf_GRV98NLM = 2 + pdf_GRV
  integer, parameter, public :: pdf_GRV98NLD = 3 + pdf_GRV
  
  integer, parameter, public :: pdf_MRS99    = 11 * pdf_base
  integer, parameter, public :: pdf_MRS99_1  = pdf_MRS99+1 
  integer, parameter, public :: pdf_MRS99_2  = pdf_MRS99+2 
  integer, parameter, public :: pdf_MRS99_3  = pdf_MRS99+3 
  integer, parameter, public :: pdf_MRS99_4  = pdf_MRS99+4 
  integer, parameter, public :: pdf_MRS99_5  = pdf_MRS99+5 
  integer, parameter, public :: pdf_MRS99_6  = pdf_MRS99+6 
  integer, parameter, public :: pdf_MRS99_7  = pdf_MRS99+7 
  integer, parameter, public :: pdf_MRS99_8  = pdf_MRS99+8 
  integer, parameter, public :: pdf_MRS99_9  = pdf_MRS99+9 
  integer, parameter, public :: pdf_MRS99_10 = pdf_MRS99+10
  integer, parameter, public :: pdf_MRS99_11 = pdf_MRS99+11
  integer, parameter, public :: pdf_MRS99_12 = pdf_MRS99+12

  integer, parameter, public :: pdf_MRS99DIS = 12 * pdf_base
  integer, parameter, public :: pdf_MRS99DIS_1  = pdf_MRS99DIS+1 
  integer, parameter, public :: pdf_MRS99DIS_2  = pdf_MRS99DIS+2 
  integer, parameter, public :: pdf_MRS99DIS_3  = pdf_MRS99DIS+3 
  integer, parameter, public :: pdf_MRS99DIS_4  = pdf_MRS99DIS+4 
  integer, parameter, public :: pdf_MRS99DIS_5  = pdf_MRS99DIS+5 
  integer, parameter, public :: pdf_MRS99DIS_6  = pdf_MRS99DIS+6 
  integer, parameter, public :: pdf_MRS99DIS_7  = pdf_MRS99DIS+7 
  integer, parameter, public :: pdf_MRS99DIS_8  = pdf_MRS99DIS+8 
  integer, parameter, public :: pdf_MRS99DIS_9  = pdf_MRS99DIS+9 
  integer, parameter, public :: pdf_MRS99DIS_10 = pdf_MRS99DIS+10
  integer, parameter, public :: pdf_MRS99DIS_11 = pdf_MRS99DIS+11
  integer, parameter, public :: pdf_MRS99DIS_12 = pdf_MRS99DIS+12

  integer, parameter, public :: pdf_MRS98    = 13 * pdf_base
  integer, parameter, public :: pdf_MRS98_1  = pdf_MRS98+1 
  integer, parameter, public :: pdf_MRS98_2  = pdf_MRS98+2 
  integer, parameter, public :: pdf_MRS98_3  = pdf_MRS98+3 
  integer, parameter, public :: pdf_MRS98_4  = pdf_MRS98+4 
  integer, parameter, public :: pdf_MRS98_5  = pdf_MRS98+5 

  integer, parameter, public :: pdf_MRS98LO  = 10 * pdf_base
  integer, parameter, public :: pdf_MRS98LO_1  = pdf_MRS98LO+1 
  integer, parameter, public :: pdf_MRS98LO_2  = pdf_MRS98LO+2 
  integer, parameter, public :: pdf_MRS98LO_3  = pdf_MRS98LO+3 
  integer, parameter, public :: pdf_MRS98LO_4  = pdf_MRS98LO+4 
  integer, parameter, public :: pdf_MRS98LO_5  = pdf_MRS98LO+5 

  integer, parameter, public :: pdf_MRST2001 = 14 * pdf_base
  integer, parameter, public :: pdf_MRST2001_1 = pdf_MRST2001+1
  integer, parameter, public :: pdf_MRST2001_2 = pdf_MRST2001+2
  integer, parameter, public :: pdf_MRST2001_3 = pdf_MRST2001+3
  integer, parameter, public :: pdf_MRST2001_4 = pdf_MRST2001+4
  
  !-- x^{-om), with om=1.0 in the g and q channels
  !   NB: do not change the 10g and 10q identifiers. They are
  !       used (exploiting their numerical values) later.
  integer, parameter, public :: pdf_POWER     = 30 * pdf_base
  integer, parameter, public :: pdf_POWER10g  = 1 + pdf_POWER
  integer, parameter, public :: pdf_POWER10q  = 2 + pdf_POWER
  integer, parameter, public :: pdf_POWER04g  = 3 + pdf_POWER
  integer, parameter, public :: pdf_POWER04q  = 4 + pdf_POWER
  real(dp), parameter :: pdf_omlist(2) = (/ one, 0.4_dp /)
  !======================================================================

  type xQ_pair
     real(dp) :: x, Q, muF_Q
  end type xQ_pair
  public :: xQ_pair
  

  type pdf_xQ_set
     integer,       pointer :: pdf_ids(:)
     type(xQ_pair), pointer :: xQ(:)
     integer                :: npdf, nxQ
     real(dp),      pointer :: pdf_vals(:,:,:)
     real(dp)               :: last_xi
     !--- some things which are more internal
     ! for the time being there is nothing here
  end type pdf_xQ_set
  public :: pdf_xQ_set
  
  public :: pm_setup_pdfs, pm_eval_pdfs, pm_isgluon
  
  interface alloc_copy
     module procedure alloc_copy_1d_xQ,&
          &           alloc_copy_1d_int
  end interface

  real(dp), public :: pm_smallest_x = one

contains



  !======================================================================
  ! Set things up...
  subroutine pm_setup_pdfs(pdfs, pdf_ids, xQ)
    type(pdf_xQ_set)             :: pdfs
    integer      , intent(in)    :: pdf_ids(:)
    type(xQ_pair), intent(in)    :: xQ(:)
    integer :: istat

    !call alloc_copy(pdf_ids, pdfs%pdf_ids)
    !-- this seems to be needed as a workaround for a VERY strange pgf90 bug
    allocate(pdfs%pdf_ids(size(pdf_ids,dim=1)))
    pdfs%pdf_ids = pdf_ids
    
    call alloc_copy(xQ     , pdfs%xQ     )
    pdfs%npdf = size(pdf_ids)
    pdfs%nxQ  = size(xQ     )

    !deallocate(pdfs%pdf_vals, stat=istat)
    allocate(pdfs%pdf_vals(pdfs%npdf, pdfs%nxQ, -6:6))

    !-- flag to allow generated x to be appropriate
    pm_smallest_x = min(pm_smallest_x, minval(xQ(:)%x))

    !-- this is designed to help caching
    pdfs%last_xi = -one
  end subroutine pm_setup_pdfs



  !======================================================================
  ! Routine to evaluate all the pdfs
  ! Right now, support just MRST?
  !
  ! Put the PDF [ x f(x)] results into the pdf_vals array
  subroutine pm_eval_pdfs(pdfs, xi)
    use fake_pdflib
    type(pdf_xQ_set), intent(inout) :: pdfs
    real(dp),         intent(inout) :: xi
    !--------------------------------------
    integer  :: ipdf, ixQ, mode, pdf_id, i
    real(dp) :: x, omega, muF
    real(dp) :: upv,dnv,usea,dsea,str,chm,bot,glu
    !-- the stuff for calling pdflib instead
    character(len=20) :: parm(20)
    real(dp)          :: value(20)
    !logical :: use_pdflib = .true.
    ! external routine
    real(dp) :: Ctq5Pdf_orig
    integer, save :: last_cteq_id = -1
    logical, save :: first_cteq = .true.


    !-- save excessive reevaluation --
    if (xi == pdfs%last_xi) return
    pdfs%last_xi = xi
    
    !----------------------------------------------------------------------
    ! Sanity check...
    if (xi < one) then
       if (xi > one - 1e-10_dp) then
          write(0,*) 'WARNING IN pm_eval_pdfs: xi < 1:',xi
          write(0,*) '           it will be reset to 1'
          xi = one
       else
          write(0,*) 'ERROR: In pm_eval_pdfs, xi is < 1:', xi
          stop
       end if
    end if

    do ipdf = 1, pdfs%npdf
       pdf_id = pdfs%pdf_ids(ipdf)
       select case(pdf_base*(pdf_id/pdf_base))
       case(pdf_MRST2001)
          mode = pdf_id - pdf_MRST2001
          if (use_pdflib) then
             write(0,*) 'MRST2001 not yet supported with pdflib'
             !parm(1) = 'NPTYPE';value(1) = 1
             !parm(2) = 'NGROUP';value(2) = 3
             !parm(3) = 'NSET'  ;value(3) = 89 + mode - 1
             !call pdfset(parm,value)
          end if
          
          do ixQ = 1, pdfs%nxQ
             x = pdfs%xQ(ixQ)%x * xi
             muF = pdfs%xQ(ixQ)%Q * pdfs%xQ(ixQ)%muF_Q
             if (x >= one) then
                pdfs%pdf_vals(ipdf,ixQ,:) = zero
                cycle
             end if
             if (use_pdflib) then
             call pftopdg(x,muF, pdfs%pdf_vals(ipdf,ixQ,:))
             else
                call mrst2001(x,muF,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
                pdfs%pdf_vals(ipdf,ixQ,:) = &
                     & (/ zero,bot,chm,str,usea,dsea,glu,&
                     &       dsea+dnv, usea+upv, str, chm, bot, zero /)
             end if
          end do
       case(pdf_MRS99)
          mode = pdf_id - pdf_MRS99
          if (use_pdflib) then
             parm(1) = 'NPTYPE';value(1) = 1
             parm(2) = 'NGROUP';value(2) = 3
             parm(3) = 'NSET'  ;value(3) = 89 + mode - 1
             call pdfset(parm,value)
          end if
          
          do ixQ = 1, pdfs%nxQ
             x = pdfs%xQ(ixQ)%x * xi
             muF = pdfs%xQ(ixQ)%Q * pdfs%xQ(ixQ)%muF_Q
             if (x >= one) then
                pdfs%pdf_vals(ipdf,ixQ,:) = zero
                cycle
             end if
             if (use_pdflib) then
             call pftopdg(x,muF, pdfs%pdf_vals(ipdf,ixQ,:))
             else
                call mrs99(x,muF,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
                pdfs%pdf_vals(ipdf,ixQ,:) = &
                     & (/ zero,bot,chm,str,usea,dsea,glu,&
                     &       dsea+dnv, usea+upv, str, chm, bot, zero /)
             end if
          end do
       case(pdf_CTEQ5)
          mode = pdf_id - pdf_CTEQ5
          if (use_pdflib) then
             parm(1) = 'NPTYPE';value(1) = 1
             parm(2) = 'NGROUP';value(2) = 4
             parm(3) = 'NSET'  ;
             select case (pdf_id)
             case(pdf_CTEQ5M)
                value(3) = 48
             case(pdf_CTEQ5D)
                value(3) = 47
             case(pdf_CTEQ5L)
                value(3) = 46
             case(pdf_CTEQ5HJ)
                value(3) = 49
             case(pdf_CTEQ5M1)
                value(3) = 53
             case default
                write(0,*) 'ERROR. Unsupported CTEQ5 pdf set :',mode   
             end select
             call pdfset(parm,value)
          else
             !-- if several cteq distributions are used this will be incredibly
             !   slow because it has to reread the appropriate file each time
             if (first_cteq) then
                first_cteq = .false.
                write(0,*) 'WARNING: cteq is very slow at the moment'
                write(0,*) 'ALSO: DISAGREEMENT BETWEEN PDFLIB804 and DIRECT CALL'
                write(0,*) 'THIS MUST, MUST BE INVESTIGATED BEFORE USE'
             end if
             
             if (pdf_id /= last_cteq_id) then
                if (last_cteq_id /= -1) then
                   write(0,*) 'ERROR: CURRENTLY only support CTEQ pdf at a time'
                   stop
                end if
                call SetCtq5(pdf_id - pdf_CTEQ5)
                last_cteq_id = pdf_id
             end if
             
          end if
          
          do ixQ = 1, pdfs%nxQ
             x = pdfs%xQ(ixQ)%x * xi
             muF = pdfs%xQ(ixQ)%Q * pdfs%xQ(ixQ)%muF_Q
             if (x >= one) then
                pdfs%pdf_vals(ipdf,ixQ,:) = zero
                cycle
             end if
             if (use_pdflib) then
                call pftopdg(x,muF, pdfs%pdf_vals(ipdf,ixQ,:))
             else
                pdfs%pdf_vals(ipdf,ixQ,6) = zero
                pdfs%pdf_vals(ipdf,ixQ,5) = x * Ctq5Pdf_orig(5,x,muF)
                pdfs%pdf_vals(ipdf,ixQ,4) = x * Ctq5Pdf_orig(4,x,muF)
                pdfs%pdf_vals(ipdf,ixQ,3) = x * Ctq5Pdf_orig(3,x,muF)
                pdfs%pdf_vals(ipdf,ixQ,-6) = zero
                pdfs%pdf_vals(ipdf,ixQ,-5) = pdfs%pdf_vals(ipdf,ixQ,5)
                pdfs%pdf_vals(ipdf,ixQ,-4) = pdfs%pdf_vals(ipdf,ixQ,4)
                pdfs%pdf_vals(ipdf,ixQ,-3) = pdfs%pdf_vals(ipdf,ixQ,3)
                
                do i = -1,1,2
                   pdfs%pdf_vals(ipdf,ixQ,i*2) = x * Ctq5Pdf_orig(i*1,x,muF)
                   pdfs%pdf_vals(ipdf,ixQ,i*1) = x * Ctq5Pdf_orig(i*2,x,muF)
                end do
                pdfs%pdf_vals(ipdf,ixQ,0) = x * Ctq5Pdf_orig(0,x,muF)
             end if
          end do
          
       case(pdf_POWER)
          !-- some nasty hackery (as is the rest of this routine)
          ! put a simple distribution into either gluon or quark
          ! if it goes into the gluon then arrange for the quark
          ! also to have a small component, so that normalisation
          ! still works.
          pdfs%pdf_vals(ipdf,:,:) = zero
          mode = mod(pdf_id - pdf_POWER - 1, 2)
          omega = pdf_omlist((pdf_id - pdf_POWER+1)/2)
          if (any(pdfs%xQ(:)%x*xi < 0)) write(0,*) pdfs%xQ(:)%x*xi
          where(pdfs%xQ(:)%x*xi < 1)
!!$             pdfs%pdf_vals(ipdf,:,mode) = (pdfs%xQ(:)%x*xi/(1e-6_dp))**(-omega)
             pdfs%pdf_vals(ipdf,:,mode) = (pdfs%xQ(:)%x*xi)**(-omega)*&
                  & (1-(pdfs%xQ(:)%x*xi)**2)**2
          end where
          if (mode == 0)&
               & pdfs%pdf_vals(ipdf,:,1) = pdfs%pdf_vals(ipdf,:,0)*1e-40_dp
       case default
          write(0,*) 'ERROR. Unsupported PDF set:',pdf_id
          stop
       end select
    end do
    
  end subroutine pm_eval_pdfs
  

  !-------------------------------------------------------------
  ! Indicates whether a given pdf_id corresponds to a pure gluon
  function pm_isgluon(pdf_id)
    integer, intent(in) :: pdf_id
    logical :: pm_isgluon
    pm_isgluon =  (pdf_base*(pdf_id/pdf_base) &
         &== pdf_POWER .and. mod(pdf_id,2) == 1)
  end function pm_isgluon
  


  !======================================================================
  ! Allocate a pointer array and copy array into it
  subroutine alloc_copy_1d_xQ(source, dest)
    type(xQ_pair), pointer    :: dest(:)
    type(xQ_pair), intent(in) :: source(:)
    integer :: istat
    !deallocate(dest, stat=istat)
    allocate(dest(size(source, dim=1)))
    dest = source
  end subroutine alloc_copy_1d_xQ
  subroutine alloc_copy_1d_int(source, dest)
    integer, pointer    :: dest(:)
    integer, intent(in) :: source(:)
    integer :: istat
    !deallocate(dest, stat=istat)
    allocate(dest(size(source, dim=1)))
    dest = source
  end subroutine alloc_copy_1d_int
  
end module pdf_manager
