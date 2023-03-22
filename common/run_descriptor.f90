!======================================================================
!
! This area is completely disorganised.
!
! It contains three sets of things:
!         hgrp type and routines 
!                 (a set of histograms as needed by NLO programs)
!
!         esbinned: a type which describes everything about how
!                   we want to study a variable (id, binning, upper limit)
!
!
!         run_def: a type which contains all information about the
!                  set of variables to be studied (esbinned), the
!                  structure functions to be used, the x,Q points,
!                  and additionally contains all the histograms needed
!                  to carry out the study. Also contains Elim 
!                  (should this also be set variable by variable).
!
!
!        rd_setup_run: to be called for initialisation
!
!======================================================================
module run_descriptor
  use types; use consts_dp
  use g90bk; use pdf_manager
  implicit none
  private

  !-- FORMAT VERSION 1
  !   this one did not even specify format_version
  !integer, parameter, public :: rd_format_version = 1
  !-- FORMAT VERSION 2 -------------------------------------
  !   this one specifies format version 
  !   (requires a change in disent/disaster specific routines)
  !   it also allows xq-by-xq setting of factorisation scale
  integer, parameter, public :: rd_format_version = 2

  !-- group of histograms, containing errors, means, etc...
  !   do not necessarily have to use every piece
  type hgrp
     type(histogram) :: hmn, herr, htmp
     real(dp)        :: mn, err,tmp
     character(len=24) :: name
  end type hgrp

  !--------------------------------------------------------
  type esbinned
     integer  :: id      ! esid
     real(dp) :: lo,hi   ! lo and hi values for binning
     integer  :: nbins   ! # of bins
     logical  :: log     ! logarithmic binning?
     real(dp) :: lim     ! limit above which want a fully integrated X-sct
  end type esbinned
  public  :: esbinned
  logical, public, parameter :: linear = .false., logarithmic = .true.
  character(len=12), parameter :: label_linear = ' linear     '
  character(len=12), parameter :: label_logrmc = ' logarithmic'
  logical, public, parameter :: logrmc = logarithmic
  real(dp), parameter, public :: rd_fake_limit = 1e80_dp
  real(dp), parameter, public :: esb_nolim = rd_fake_limit


  !--------------------------------------------------------
  ! This is the type which we will use to describe a run
  ! There will be a routine to set it up
  type run_def
     !-- below are the pieces which will be set by the user
     type(esbinned),   pointer :: es_list(:)
     type(pdf_xQ_set), pointer :: pdfs
     !-- the bits which are then set up by the routine itself --
     integer             :: nes, npdf, nxQ
     integer             :: nxsct, nord
     !-- distributions ----- nes
     !   first index starts from zero [npdf,nxQ,nes,nord]
     type(hgrp), pointer :: h(:,:,:,:)
     !   contains amount above a certain limit
     !   (by convention, this should mean amount above a certain limit, 
     !    given that energy in the current hemisphere is above Elim)
     type(hgrp), pointer :: hlims(:,:,:,:)
     !-- cross sections
     !   first index starts from zero and is actual cross section
     !   others are limited quantities
     !      [npdf,nxQ, {1=xsct 2=xsct E>Elim} ,nord]
     type(hgrp), pointer :: x(:,:,:,:)
     !-- temporary arrays (avoid allocation) ----------
     !   es_vals_proc may contain either the log or the linear value
     !   us it for binning, while es_vals gets used for mean value
     !   and limits on the variable(?)
     real(dp), pointer :: es_vals(:), es_vals_proc(:)
     real(dp), pointer :: weights(:,:)
     !                          ibin, ies,nord
     integer, pointer :: bin_list(:,   :,:)   ! indices of used bins
     integer, pointer :: nbin         (:,:)   ! # of bins used
     logical, pointer :: binit(:)             ! whether to bin it this time
     !-- that horrible parameter
     real(dp)          :: Elim
  end type run_def

  integer, public, parameter :: id_xsct = 1
  integer, public, parameter :: id_elim = 2


  interface alloc_copy
     module procedure alloc_copy_2d_dp, alloc_copy_1d_dp,&
          &           alloc_copy_2d_int, alloc_copy_1d_int,&
          &           alloc_copy_1d_esb
  end interface
  

  public :: hgrp, run_def, rd_setup_run, rd_write_res, rd_collate_tmp
  public :: setup_hgrp, hgrp_collate_tmp, rd_collate_tmp_list
  public :: rd_zero_tmp_list

  !-- these are needed to take into account the fact that disaster 
  !   does the onejet and the twojet runs separately
  logical, public :: rd_onejet = .true.
  logical, public :: rd_twojet = .true.
  
  public :: rd_parseiseed, rd_message

contains

  
  !--------------------------------------------------------------
  ! setup the details of a run...
  subroutine rd_setup_run (run, es_list, pdfs, Elim) 
    use event_shapes
    type(run_def)                :: run
    type(esbinned),   intent(in) :: es_list(:)
    type(pdf_xQ_set), pointer    :: pdfs
    real(dp),         intent(in) :: Elim
    !------------
    integer :: istat, ies, ipdf, ixQ, ixsct, iord


    run%nes   = size(es_list)
    run%npdf  = pdfs%npdf
    run%nxQ   = pdfs%nxQ
    run%nxsct = 2
    run%nord  = 2
    
    !call alloc_copy(es_list, run%es_list ) 
    !-- this seems to be needed as a workaround for a VERY strange pgf90 bug
    allocate(run%es_list(size(es_list)))
    run%es_list = es_list

    run%pdfs => pdfs
    run%Elim =  Elim

    !-- make sure that things are clean
    !deallocate(run%h, run%x, run%es_vals, run%es_vals_proc, run%hlims, run&
    !     & %weights, run%bin_list, run%nbin, run%binit, stat=istat)
    
    !-- allocate the necessary bits and pieces
    allocate(run%h(run%npdf, run%nxQ, run%nes, run%nord))
    allocate(run%x(run%npdf, run%nxQ, run%nxsct, run%nord))
    allocate(run%hlims(run%npdf, run%nxQ, run%nes, run%nord))
    allocate(run%es_vals(run%nes))
    allocate(run%es_vals_proc(run%nes))
    allocate(run%weights(run%npdf, run%nxQ))
    allocate(run%bin_list(20,run%nes,run%nord))
    allocate(run%nbin(run%nes,run%nord))
    allocate(run%binit(run%nes))
    run%nbin = 0


    !--- GET TO WORK HERE ON SETTING UP THE LIMITS...
    do ipdf = 1, run%npdf
       do ixQ = 1, run%nxQ
          do iord = 1, run%nord
             do ies = 1, run%nes
                !-- event shapes --
                call setup_hgrp(run%h(ipdf,ixQ,ies,iord), &
                     & 'Dist '//es_names(run%es_list(ies)%id), &
                     & run%es_list(ies)%lo, run%es_list(ies)%hi, run&
                     & %es_list(ies)%nbins )  
                !-- event shapes above a certain limit
                call setup_hgrp(run%hlims(ipdf,ixQ,ies,iord), &
                     & 'Limited '//es_names(run%es_list(ies)%id) ) 
             end do
             !-- cross sections ---------------------------------
             call setup_hgrp(run%x(ipdf,ixQ,id_xsct,iord), 'X-sect')
             call setup_hgrp(run%x(ipdf,ixQ,id_elim,iord), 'Elim')
          end do
       end do
    end do
    

  end subroutine rd_setup_run
  
  !======================================================================
  ! To be done at the end of every event-set: transfers sums of things
  ! from the temporary arrays and variables into the other vars.
  subroutine rd_collate_tmp(run)
    type(run_def), intent(inout), target :: run
    integer :: ies, ipdf, ixQ, iord, ixsct
    
    !-- at some stage we will want to think about introducing globally
    !   allocated masks -- to be dealt with in some as yet to be specified
    !   manner
    do ipdf = 1, run%npdf
       do ixQ  = 1, run%nxQ
          do iord = 1, run%nord
             if (rd_twojet) then
                do ies = 1, run%nes
                   call hgrp_collate_tmp(run%h(ipdf,ixQ, ies, iord))
                   call hgrp_collate_tmp(run%hlims(ipdf,ixQ, ies, iord))
                end do
                call hgrp_collate_tmp(run%x(ipdf,ixQ, 2, iord))  !Elim
             end if
             if (rd_onejet) &
                  &call hgrp_collate_tmp(run%x(ipdf,ixQ, 1, iord))  !X-sct
          end do
       end do
    end do
  end subroutine rd_collate_tmp

  !----------------------------------------------------------------------
  ! The same, but hopefully gets speeded up by a list
  subroutine rd_collate_tmp_list(run)
    type(run_def), intent(inout), target :: run
    integer :: ies, ipdf, ixQ, iord, ixsct
    
    !-- at some stage we will want to think about introducing globally
    !   allocated masks -- to be dealt with in some as yet to be specified
    !   manner
    do ipdf = 1, run%npdf
       do ixQ  = 1, run%nxQ
          do iord = 1, run%nord
             if (rd_twojet) then
                do ies = 1, run%nes
                   call hgrp_collate_tmp_list(run%h(ipdf,ixQ, ies, iord),&
                        & run%bin_list(:run%nbin(ies,iord),ies,iord))
                   call hgrp_collate_tmp_nohist(run%hlims(ipdf,ixQ, ies, iord))
                end do
                call hgrp_collate_tmp_nohist(run%x(ipdf,ixQ, 2, iord))  !Elim
             end if
             if (rd_onejet) &
               &call hgrp_collate_tmp_nohist(run%x(ipdf,ixQ, 1, iord))  !X-sct
          end do
       end do
    end do
    run%nbin = 0
  end subroutine rd_collate_tmp_list

  !--------------------------------------------------------------------
  ! Need to be able to zero things easily?
  subroutine rd_zero_tmp_list(run)
    type(run_def), intent(inout), target :: run
    integer :: ies, ipdf, ixQ, iord, ixsct
    
    !-- at some stage we will want to think about introducing globally
    !   allocated masks -- to be dealt with in some as yet to be specified
    !   manner
    do ipdf = 1, run%npdf
       do ixQ  = 1, run%nxQ
          do iord = 1, run%nord
             if (rd_twojet) then
                do ies = 1, run%nes
                   call hgrp_zero_tmp_list(run%h(ipdf,ixQ, ies, iord),&
                        & run%bin_list(:run%nbin(ies,iord),ies,iord))
                   call hgrp_zero_tmp_nohist(run%hlims(ipdf,ixQ, ies, iord))
                end do
                call hgrp_zero_tmp_nohist(run%x(ipdf,ixQ, 2, iord))  !Elim
             end if
             if (rd_onejet) &
               &call hgrp_zero_tmp_nohist(run%x(ipdf,ixQ, 1, iord))  !X-sct
          end do
       end do
    end do
    run%nbin = 0
  end subroutine rd_zero_tmp_list
  

  !======================================================================
  ! subroutine rd_write_res(run, nev, niter)
  ! nev corresponds to the actual number of events used
  ! RelIter corresponds to the ratio of the number of events actually used 
  ! to the number of events that would be necessary to get the proper
  ! normalisation
  subroutine rd_write_res(run, nev_xsct, RelIter_xsct, nev, RelIter)
    use sub_defs_io; use event_shapes; use pdf_names
    type(run_def), intent(in), target :: run
    integer,       intent(in) :: nev_xsct, nev
    real(dp),      intent(in) :: RelIter_xsct, RelIter
    integer :: ies, ipdf, ixQ, iord, ixsct
    type(histogram)     :: mn,err
    type(hgrp), pointer :: hh
    character(len=100) :: line, quant
    real(dp) :: norm, sigres, sigerr
    character(len=*), parameter :: format_mean = &
         &'(a,i1," = ",es22.12," +- ",es22.12)'
    

    do ipdf = 1, run%npdf
    do ixQ  = 1, run%nxQ
       call hist_tdtext('==============================================')
!!$       call hist_tdtext('PDF identifier = '&
!!$            &//trim(num2char(run%pdfs%pdf_ids(ipdf))))
       call hist_tdtext('PDF identifier = '&
            &//trim(pn_ID2Name(run%pdfs%pdf_ids(ipdf))))
       call hist_tdtext('x = '//trim(num2char(run%pdfs%xQ(ixQ)%x)))
       call hist_tdtext('Q = '//trim(num2char(run%pdfs%xQ(ixQ)%Q)))
       call hist_tdtext('muF_Q = '//trim(num2char(run%pdfs%xQ(ixQ)%muF_Q)))

       !-- write LO and NLO cross sections
       !   LO is absolute, NLO is relative to LO.
       !   write also the cross sections for E < Elim
       norm = RelIter_xsct
       do iord = 1, run%nord
          hh => run%x(ipdf,ixQ, 1, iord )
          sigres  = hh%mn / norm
          sigerr  = sqrt(abs(hh%err - hh%mn**2/nev_xsct)) / norm
          write(line,format_mean) "sigma",iord-1, sigres, sigerr
          call hist_tdtext(line)
          if (iord == 1) then
             norm = sigres*RelIter
             !-- a very nasty hack to get right normalisation
             !   do something a bit nicer some other time...
             if (pm_isgluon(run%pdfs%pdf_ids(ipdf))) norm = norm * 1e40_dp
             !if (run%pdfs%pdf_ids(ipdf) == pdf_POWER10g) norm = norm * 1e40_dp
          end if
          
       end do

       
       
       do ies = 1, run%nes
          call hist_tdtext('--------------------------------------------')
          call hist_tdtext('Event-shape identifier = '&
               &//es_names(run%es_list(ies)%id))
          !-- repeat this information for easier human reading.........
!!$          call hist_tdtext('PDF identifier = '&
!!$               &//trim(num2char(run%pdfs%pdf_ids(ipdf))))
          call hist_tdtext('PDF identifier = '&
               &//trim(pn_ID2Name(run%pdfs%pdf_ids(ipdf))))
          call hist_tdtext('x = '//trim(num2char(run%pdfs%xQ(ixQ)%x)))
          call hist_tdtext('Q = '//trim(num2char(run%pdfs%xQ(ixQ)%Q)))
          call hist_tdtext('muF_Q = '//trim(num2char(run%pdfs%xQ(ixQ)%muF_Q)))
          !-- .................................

          do ixsct = 1, 3
             do iord = 1, run%nord
                select case(ixsct)
                case(1)   !-- Elim
                   hh =>  run%x(ipdf,ixQ, 2, iord )
                   quant = "sigElim"
                   if (iord == 1) call hist_tdtext ('-- Cross sections &
                        &for E < Elim = '// num2char(run%Elim)) 
                case(2)   !-- amount above eslim
                   !write(0,*) num2char(run%es_list(ies)%lim)
                   hh =>  run%hlims(ipdf,ixQ, ies, iord )
                   quant = "sigESlim"
                   if (iord == 1) call hist_tdtext ('-- Cross sections &
                        &for shape > ESlim = '//&
                        &num2char(run%es_list(ies)%lim))
                case(3)   !-- mean values
                   hh =>  run%h(ipdf,ixQ, ies, iord )
                   if (iord == 1) call hist_tdtext ('-- Mean values') 
                   quant = "mean"
                case default
                   write(0,*) 'WHAT THE HELL IS GOING ON?'
                end select

                sigres  = hh%mn / norm
                sigerr  = sqrt(abs(hh%err - hh%mn**2/nev)) / norm
                write(line,format_mean) trim(quant),iord, sigres, sigerr
                call hist_tdtext(line)
             end do
          end do
          
          do iord = 1, run%nord
             hh =>  run%h(ipdf,ixQ, ies, iord )
             if (run%es_list(ies)%log .eqv. linear) then
                write(line,'(a,i1,a)') "Distribution-",iord,label_linear
             else
                write(line,'(a,i1,a)') "Distribution-",iord,label_logrmc
             end if
             call hist_tdtext(line)
             write(line,'(a,i6)') "NBins = ",hh%hmn%n
             call hist_tdtext(line)
             !-- manipulate copies to allow continuation --
             call hist_copy(hh%hmn,  mn)
             call hist_copy(hh%herr, err)
             !-- normalise to LL quark cross section -----------------
             mn%h = mn%h / norm
             err%h = err%h / norm**2
             err%h = sqrt(abs(err%h - mn%h**2/nev))
             !-- write out the histogram -----------------------------
             call hist_tdwrite(mn,new=.true., log=.false.,join=.false.,&
                  & histerr=err)
          end do
       end do
          
    end do
    end do

  end subroutine rd_write_res
  



  !======================================================================
  !----------------------------------------------------------------------
  ! Setup a group. 
  subroutine setup_hgrp(grp,name,min,max,n)
    type(hgrp) :: grp
    character(len=*) :: name
    real(dp), optional :: min, max
    integer,  optional :: n
    if (present(min)) then
       call hist_init(grp%hmn, name,min,max,n)
       call hist_init(grp%herr,name//' errors',min,max,n)
       call hist_init(grp%htmp,name//' temp',min,max,n)
    end if
    grp%mn  = zero
    grp%err = zero
    grp%tmp = zero
    grp%name = name
  end subroutine setup_hgrp

  !----------------------------------------------------------------------
  ! Does the job of transfering entries from the tmp vars to the final 
  ! ones
  subroutine hgrp_collate_tmp(grp,nzero_mask)
    type(hgrp)          :: grp
    logical,optional    :: nzero_mask(:)
    integer :: i,j
    grp%mn  = grp%mn  + grp%tmp
    grp%err = grp%err + grp%tmp**2
    grp%tmp = zero

    if (associated(grp%hmn%h)) then
       if (present(nzero_mask)) then
!!$          if (any(nzero_mask /= (grp%htmp%h /= zero)&
!!$               & .and. (grp%htmp%h /= zero))) then
!!$             j = lbound(grp%htmp%h,dim=1) - 1
!!$             do i = 1,size(nzero_mask)
!!$                write(0,*) i,nzero_mask(i),  grp%htmp%h(i+j)
!!$             end do
!!$             stop
!!$          end if
          where(nzero_mask)
             grp%hmn%h  = grp%hmn%h  + grp%htmp%h
             grp%herr%h = grp%herr%h + grp%htmp%h**2
             grp%htmp%h = zero
          end where
       else
          where(grp%htmp%h/= zero)
             grp%hmn%h  = grp%hmn%h  + grp%htmp%h
             grp%herr%h = grp%herr%h + grp%htmp%h**2
             grp%htmp%h = zero
          end where
       end if
    end if
  end subroutine hgrp_collate_tmp


  !----------------------------------------------------------------------
  ! Does the job of transfering entries from the tmp vars to the final 
  ! ones -- but uses a list of bins
  ! Assumes that hist is allocated
  subroutine hgrp_collate_tmp_list(grp,list)
    type(hgrp)          :: grp
    integer             :: list(:)
    integer :: i,j
    grp%mn  = grp%mn  + grp%tmp
    grp%err = grp%err + grp%tmp**2
    grp%tmp = zero

    do i = 1, size(list)
       j = list(i)
       grp%hmn%h(j)  = grp%hmn%h(j)  + grp%htmp%h(j)
       grp%herr%h(j) = grp%herr%h(j) + grp%htmp%h(j)**2
       grp%htmp%h(j) = zero
    end do
  end subroutine hgrp_collate_tmp_list

  !----------------------------------------------------------------------
  ! needed sometimes for a reset after problems
  subroutine hgrp_zero_tmp_list(grp,list)
    type(hgrp)          :: grp
    integer             :: list(:)
    integer :: i,j
    grp%mn  = grp%mn  + grp%tmp
    grp%err = grp%err + grp%tmp**2
    grp%tmp = zero

    do i = 1, size(list)
       j = list(i)
       grp%htmp%h(j) = zero
    end do
  end subroutine hgrp_zero_tmp_list

  !----------------------------------------------------------------------
  ! For use when there is no hist
  subroutine hgrp_collate_tmp_nohist(grp)
    type(hgrp)          :: grp
    grp%mn  = grp%mn  + grp%tmp
    grp%err = grp%err + grp%tmp**2
    grp%tmp = zero
  end subroutine hgrp_collate_tmp_nohist

  !----------------------------------------------------------------------
  ! For use when there is no hist
  subroutine hgrp_zero_tmp_nohist(grp)
    type(hgrp)          :: grp
    grp%tmp = zero
  end subroutine hgrp_zero_tmp_nohist

  !======================================================================
  ! Allocate a pointer array and copy array into it
  subroutine alloc_copy_2d_dp(source, dest)
    real(dp), pointer    :: dest(:,:)
    real(dp), intent(in) :: source(:,:)
    integer :: istat
    !deallocate(dest, stat=istat)
    allocate(dest(size(source, dim=1), size(source,dim=2)))
    dest = source
  end subroutine alloc_copy_2d_dp
  subroutine alloc_copy_2d_int(source, dest)
    integer, pointer    :: dest(:,:)
    integer, intent(in) :: source(:,:)
    integer :: istat
    !deallocate(dest, stat=istat)
    allocate(dest(size(source, dim=1), size(source,dim=2)))
    dest = source
  end subroutine alloc_copy_2d_int
  subroutine alloc_copy_1d_dp(source, dest)
    real(dp), pointer    :: dest(:)
    real(dp), intent(in) :: source(:)
    integer :: istat
    !deallocate(dest, stat=istat)
    allocate(dest(size(source, dim=1)))
    dest = source
  end subroutine alloc_copy_1d_dp
  subroutine alloc_copy_1d_int(source, dest)
    integer, pointer    :: dest(:)
    integer, intent(in) :: source(:)
    integer :: istat
    !deallocate(dest, stat=istat)
    allocate(dest(size(source, dim=1)))
    dest = source
  end subroutine alloc_copy_1d_int
  subroutine alloc_copy_1d_esb(source, dest)
    type(esbinned), pointer    :: dest(:)
    type(esbinned), intent(in) :: source(:)
    integer :: istat
    !deallocate(dest, stat=istat)
    allocate(dest(size(source, dim=1)))
    dest = source
  end subroutine alloc_copy_1d_esb
  


  !=========================================================================
  ! not strictly the right place
  ! But it is useful to have it in a common location
  ! It looks for the Current ISEED string and reads it in
  ! It is up to the caller to ensure that the length of iseed is appropriate
  subroutine rd_parseiseed(iseed)
    integer, intent(out) :: iseed(:)
    !-----------------------------
    integer :: dev, iostat, ix
    character(len=200) :: line
    character(len=*), parameter :: string='Current ISEED is'

    iseed(:) = -1

    dev = hist_tdopenlast()
    if (dev < 0) return

    do
       read(dev,'(a)', iostat=iostat) line
       if (iostat /= 0) exit
       ix = index(line,string)
       if (ix > 0) then
          !write(0,*) trim(line)
          read(line(ix+len(string):),*,iostat=iostat) iseed
          !-- make sure caller got what he/she wanted
          if (iostat /= 0) iseed = -1
          exit
       end if
    end do
    
    !-- clean up...
    close(dev)
  end subroutine rd_parseiseed
  

  subroutine rd_message()
    write(0,*) '-----------------------------------------------------------------'
    write(0,*) 'DISPATCH. An optimised DISENT/DISASTER driver.'
    write(0,*) 'Version 1.0.5'
    write(0,*) 'Written by Gavin Salam, June-August 2000'
    write(0,*) 'with updates in October 2001 & 2020-2023'
    write(0,*) 'If you use this program please refer to'
    write(0,*) 'M. Dasgupta and G.P. Salam, JHEP 0208 (2002) 032 [hep-ph/0208073]'
    write(0,*) 'and to the original DISENT/DISASTER publications and arXiv:2010.07354'
    write(0,*) '-----------------------------------------------------------------'
  end subroutine rd_message
  


end module run_descriptor
