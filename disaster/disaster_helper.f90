!======================================================================
! Module for checking to see whether a number is a NaN
! Will fail in some very particular cases (if it is very close to huge) 
! Not at all portable, since NaNs behave differently everywhere
!
! Tested only with the absoft compiler
!
! GPS 30/08/2000
module check_nan
  implicit none
  integer, parameter :: dp = kind(1d0)
  private

  real(dp), parameter :: huge_pos =  (huge(1.0_dp)*0.9999999999999_dp)
  real(dp), parameter :: huge_neg = -(huge(1.0_dp)*0.9999999999999_dp)
  
  interface is_a_nan
     module procedure is_a_nan_sclr, is_a_nan_1d, is_a_nan_2d
  end interface
  
  public :: is_a_nan

contains
  
  function is_a_nan_sclr(x) result(is_a_nan)
    real(dp), intent(in) :: x
    logical              :: is_a_nan!(size(x))
    is_a_nan = (x > huge_pos) .or. (x < huge_neg)
  end function is_a_nan_sclr

  function is_a_nan_1d(x) result(is_a_nan)
    real(dp), intent(in) :: x(:)
    logical              :: is_a_nan(size(x))
    is_a_nan = (x > huge_pos) .or. (x < huge_neg)
  end function is_a_nan_1d
  function is_a_nan_2d(x) result(is_a_nan)
    real(dp), intent(in) :: x(:,:)
    logical              :: is_a_nan(size(x,dim=1),size(x,dim=2))
    is_a_nan = (x > huge_pos) .or. (x < huge_neg)
  end function is_a_nan_2d
end module check_nan


!======================================================================
! The common routines that are of use to disaster
module disaster_helper
  use types; use consts_dp
  use run_descriptor
  use event_shapes!, except => dot
  use pdf_manager
  use sub_defs_io
  implicit none

  private


  !== lower limit on number of events to have in an interation
  integer, public :: nev_iter_limit = 10000

  !-- constants for disaster --------------------------------------
  real(dp), parameter  :: alphaem = 1/137.0_dp
  integer,  parameter  :: scheme = 0
  real(dp)             :: nf = 5, scale = one
  real(dp)             :: cutoff
  integer              :: nfint = 5
  !-- rho2pdf: contains all the stuff relating to nf-dependence
  real(dp), target :: rho2pdf(-2:2,11)
  
  !--- need a pointer for parser to be allocate it without referring
  !    to disaster_helper explicitly
  !type(run_def), allocatable, target :: run(:)
  type(run_def), pointer :: run(:)
  
  !---- various things for communicating the stage which has been
  !     reached
  integer, public :: nevt_xsct, nevt_wanted, nevt_jet, nevt
  logical, public :: final, jetrun
  !-- how to deal with output file
  logical         :: noclose
  integer         :: writefreq

  integer, parameter :: adapt_one  = 1
  integer, parameter :: adapt_sqrt = 2
  integer, parameter :: adapt_lin  = 3
  integer, parameter :: adapt_abs  = 4
  integer, parameter :: adapt_abs2 = 5
  integer            :: adapt = adapt_one
  !-- the index (within the sequence of specified vars) of the var 
  !   and pdf on which to adapt
  integer            :: ivar_adapt = 1
  integer            :: ipdf_adapt = 1


  public :: dh_start, dh_user2!, dh_final, dh_event
  public :: run

  !-- the following are public to allow compatibility with older programs
  real(dp), public :: dh_ybj = -one

  integer :: seed_init(3)

  !-- help deal with NaNs -------------
  logical :: detected_nan = .false.
  ! might need nevt_nan to reestablish correct normalisation
  ! at least if it is a significant fraction of total events
  integer :: nevt_nan = zero

contains

  !----------------------------------------------------------------------
  ! Once run(:) has been set up, then call here
  subroutine dh_start(nf_in, ybj)
    integer,  intent(in) :: nf_in
    real(dp), intent(in) :: ybj
    INTEGER  :: NEV, NEV_INT, IN
    integer  :: ndiv = 20,  ORDER, iseq, niter, irun
    integer  :: points
    real(dp) :: S, fact_prep, fact_final
    character(len=70) :: tdfl


    ! --- redirect all output to files with names similar to tdfile.
    tdfl = string_val_opt('-tdfl','')
    if (tdfl /= '' ) then
       open(FILE=trim(tdfl)//'.log', UNIT=42, STATUS='UNKNOWN') 
       open(FILE=trim(tdfl)//'.cot', UNIT=41, STATUS='UNKNOWN') 
       open(FILE=trim(tdfl)//'.err', UNIT=0, STATUS='UNKNOWN') 
       open(FILE=trim(tdfl)//'.out', UNIT=6, STATUS='UNKNOWN') 
    end if

    !-- tell user what is happening...
    call rd_message()

    nf     = nf_in
    nfint  = nint(nf)
    call dh_setup_rho2pdf
    dh_ybj = ybj


    !-- details of how to write the output...
    noclose = log_val_opt('-noclose',.false.)
    writefreq = nint(dble_val_opt('-wf',5e5_dp))

    S=4*30*820 
    
    !-- initialise the random number generator 
    call rd_parseiseed(seed_init)
    iseq = int_val_opt('-iseq',1)
    if (seed_init(1) > 0 .and. iseq == seed_init(1)) then
       write(0,*) 'Using last seed from previous run'
       write(0,*) 'Seed:', seed_init
       call vginit_ran3(seed_init)
       !call dh_SetSeed(seed_init)
    else
       call vginit_ran(iseq)
       write(0,*) 'Using default seed with iseq =', iseq
    end if
    

    adapt = int_val_opt('-adapt', adapt_one)
    ivar_adapt = int_val_opt('-ivar', 1)
    ipdf_adapt = int_val_opt('-ipdf', 1)
    write(0,*) 'Adaptation scheme will be of type:', adapt
    write(0,*) 'and  based on var #:', ivar_adapt
    write(0,*) 'and  based on pdf #:', ipdf_adapt

    if (ivar_adapt > ubound(run(1)%h,dim=3)) then
       stop 'ERROR: ivar_adapt is larger than number of event shapes&
            & being consided' 
    end if
    if (ipdf_adapt > ubound(run(1)%h,dim=1)) then
       stop 'ERROR: ipdf_adapt is larger than number of PDFs&
            & being consided' 
    end if
    
    
    do irun = 1, size(run)
       if (any(run(irun)%pdfs%xQ(:)%muF_Q /= one)) then
          write(0,*) 'Disaster driver does not currently support muF_Q /= 1'
          stop
       end if
    end do

    write(0,*) '--------------------------------------------------------'

    ! modified this because of problems with core dumps with pgf90
    !seed_init = dh_GetSeed()
    call rm48ut(seed_init(1), seed_init(2), seed_init(3))
    !-- make this default -- life will be simpler




    ! --- start DISASTER++
    call disaster_ca() 
    ! --- set some parameters
    call disaster_cd('ECM', sqrt(S)) 
    call disaster_ci('LEPTON_INTEGRATION', 1) 
    call disaster_ci('PROCESS_INDEX', 2) 
    call disaster_cd('FACT_PREP',  0.25_dp) 
    call disaster_cd('FACT_FINAL', 0.75_dp) 

    !-- have a routine to do this, as was the case for disent
    !   or too ugly?
    call dh_cuts

    cutoff = dble_val_opt('-cutoff',1e-12_dp)
    call disaster_cd('TECHNICAL_CUT_PARAMETER', cutoff)

    !-- info for distributions: get it straight away so that
    !   we can do all checks about arguments
    points     = nint(dble_val_opt('-nev',1e6_dp))
    niter      = int_val_opt('-niter', 5)
    if (.not. CheckAllOptsUsed(ErrDev=0)) stop


    !-- run it: X-sctn. Do no adaptation; this ensures correct
    !                   LO cross section. NLO piece must be
    !                   calculated by hand
    jetrun = .false.
    call disaster_ci('NUMBER_OF_FINAL_STATE_PARTONS_IN_BORN_TERM', 1)
    call disaster_ci('ITERATIONS', 0) 
    call disaster_ci('POINTS', 10000) 
!!$    ! TMP TMP TMP
!!$    call disaster_ci('POINTS', 10) 
    call disaster_cc('RUN_MC') 
    
    nevt_xsct = nevt

    !-- run-it: distributions
    fact_prep  = 0.1_dp
    fact_final = one - fact_prep
    nevt_wanted = nint(points * fact_final)
    jetrun = .true.
    call disaster_cd('FACT_PREP',  fact_prep) 
    call disaster_cd('FACT_FINAL', fact_final) 
    call disaster_ci('NUMBER_OF_FINAL_STATE_PARTONS_IN_BORN_TERM', 2)
    call disaster_ci('ITERATIONS', niter) 
    call disaster_ci('POINTS', points) 
    call disaster_cc('RUN_MC') 
    
    nevt_jet = nevt

    call dh_final(nevt_jet, one, interim=.false.)
    write(0,*) "Just after dh_final..."
    
    ! --- end DISASTER++
    call disaster_cb() 
    write(0,*) "Just after disaster_cb"
    
    ! --- close the logfile
    close(UNIT=42) 
    write(0,*) "Just after close 42"
    call lcl_flush(0)
  end subroutine dh_start



  !----------------------------------------------------------------------
  ! A peculiarity of disaster...
  ! only construct the first two flavours, since the others are identical
  subroutine dh_setup_rho2pdf
    real(dp), parameter :: eq(6) = (/ -1/three, 2/three, -1/three, 2/three,&
         & -1/three, 2/three /)
    integer :: i, j

    rho2pdf = zero

    do i = 1, min(2,nfint)
       rho2pdf( i, 1) = eq(i)**2
       rho2pdf(-i, 2) = eq(i)**2
       rho2pdf( i, 4) = eq(i)**2
       rho2pdf(-i, 5) = eq(i)**2
       rho2pdf( i, 6) = eq(i)**2 * (nfint-1)
       rho2pdf(-i, 7) = eq(i)**2 * (nfint-1)
    end do

    do i = 1, nfint
       rho2pdf( 0, 3) = rho2pdf(0, 3) + eq(i)**2
    end do
    
    do i = 1, min(2,nfint)
       do j = 1, nfint
          if (j == i) cycle
          rho2pdf( i, 8) = rho2pdf( i, 8) + eq(j)**2
          rho2pdf(-i, 9) = rho2pdf(-i, 9) + eq(j)**2
          rho2pdf( i,10) = rho2pdf( i,10) + eq(i)*eq(j)
          rho2pdf(-i,11) = rho2pdf(-i,11) + eq(i)*eq(j)
       end do
    end do
    
    
  end subroutine dh_setup_rho2pdf
  


  !----------------------------------------------------------------------
  ! Outputs the results
  ! nev is total number of events so far; in is the number of divisions so
  ! far (each division sums to total cross section)
  ! RelIter is defined as in rd_write_res
  subroutine dh_final(nev,RelIter, interim)
    use g90bk
    integer,  intent(in) :: nev
    real(dp), intent(in) :: RelIter
    logical,  intent(in), optional :: interim
    integer :: irun
    if (noclose) call hist_tdrewind   ! reset histo out file
    write(0,*) "Just before dh_describe details...."
    call dh_DescribeDetails(interim)
    write(0,*) "Just after dh_describe details...."
    call hist_tdtext(' NEV = '//num2char(nev),'(X ')
    call hist_tdtext('******************************************************')
    call hist_tdtext('ybj = '//num2char(dh_ybj))
    do irun = 1, size(run)
       call rd_write_res(run(irun), nevt_xsct, one, nev, RelIter)
    end do
    if (noclose) then
       call hist_tdflush
    else
       call hist_tdclose   ! ensure that everything is cleanly output
    end if
    
  end subroutine dh_final


  !------- say what sort of run this was: DISENT specific details?
  subroutine dh_DescribeDetails(interim)
    use g90bk
    logical,  intent(in), optional :: interim
    character(len=*), parameter :: cc='(X '
    character(len=100) :: line

    integer :: seed(3)
    !-- first some things related to the technical details of the run
    write(line,'(a,i2)') 'DISASTER RUN -- format version ', rd_format_version
    call hist_tdtext(trim(line),cc)
    !call hist_tdtext('DISASTER RUN',cc)
    if (present(interim)) then
       if (interim) call hist_tdtext(' interim result',cc)
    end if
    
    write(line,*) 'Initial ISEED was ',seed_init
    call hist_tdtext(line,cc)
    ! modified this because of problems with core dumps with pgf90
    !seed = dh_GetSeed()
    call rm48ut(seed(1), seed(2), seed(3))
    write(line,*) 'Current ISEED is  ',seed
    call hist_tdtext(line,cc)
    call hist_tdtext(' adapt  = '//num2char(adapt), cc)
    call hist_tdtext(' ivar_adapt  = '//num2char(ivar_adapt), cc)
    call hist_tdtext(' ipdf_adapt  = '//num2char(ipdf_adapt), cc)
    call hist_tdtext(' cutoff = '//num2char(cutoff), cc)
    call hist_tdtext(' scheme = '//num2char(scheme), cc)
    !-- not in format version >= 2
    !call hist_tdtext(' scale  = '//num2char(scale ), cc)
    write(line,&
         & '(" CA =",f10.6,"  CF =",f10.6,"  TR =",f10.6,"  nf =",i3)')&
         & three, four/three, half, nfint
    call hist_tdtext(line,cc)
  end subroutine dh_DescribeDetails
  

  !======================================================================
  ! nr is number of phase-spaces
  ! nl is number of weights (each associated with phase space = irps(il)
  ! fvect(0:3,-10:10,1:*) where * = npartons(ir)
  ! Presume xb(:), q2(:), xi(:) all go up to nr
  ! 
  ! ialphas(1:*), ialphaem(1:*), lognf(1:*) all go up to nl
  !
  ! weight goes up to nl
  !
  !
  function dh_user2( nr, nl, fvect, npartons, xb, q2, xi, weight, irps,&
       & ialphas, ialphaem, lognf ) 
    use types
    implicit none
    real(dp) dh_user2
    integer  nr, nl 
    real(dp) fvect(0:3, -10:10, 1:*) 
    integer  npartons(1:*) 
    real(dp) xb(1:*), q2(1:*), xi(1:*), weight(1:11, 1:*) 
    integer  irps(1:*), ialphas(1:*), ialphaem(1:*), lognf(1:*)    
    !---- up to three partons allowed -----------
    real(dp) :: parr(4,3), photon(4), aobs
    real(dp) :: weight_pdf(-6:6)
    integer  :: np
    integer :: ir, il, i, na, irun, irho
    logical :: ir_already_used
    real(dp) :: ffactin(4), ffactout(13), dxpdf(-6:6), rhoxx(-6:6,1:11)

    dh_user2 = zero
    nevt = nevt + 1

    detected_nan = .false.

!!$    write(0,*) "REMEMBER TO SWITCH TO USING f(x), rather than xf(x)"

    do ir = 1, nr
       photon(1:3) = fvect(1:3,-4,ir) - fvect(1:3,-5,ir)
       photon(4  ) = fvect(0  ,-4,ir) - fvect(0  ,-5,ir)
       np = npartons(ir)
       parr(1:3,:np) = fvect(1:3, 0:np-1, ir)
       parr(4  ,:np) = fvect(0  , 0:np-1, ir)
       !-- now find the weight sets which correspond to this ir
       ir_already_used = .false.
       do il = 1, nl
          if (irps(il) /= ir) cycle
          !-- for now do not support variable scales
          if (lognf(il) /= 0) cycle
          !-- get weights into a simple-to-use form
          weight_pdf = zero
          do irho = 1, 11
             weight_pdf(-2:2) = weight_pdf(-2:2) +&
                  & weight(irho, il) * rho2pdf(:,irho)
          end do
          na = ialphas(il)
          !-- 1/xi -> use f(xi), not xi*f(xi)
          weight_pdf(-2:2) = weight_pdf(-2:2) *&
               & ( (twopi)**na * alphaem**ialphaem(il)) / xi(ir)
          do i = 3, nfint
             weight_pdf( i) = weight_pdf( i-2)
             weight_pdf(-i) = weight_pdf(-i+2)
          end do

          !-- for testing purposes...
!!$          !-- use MRST99 == 1,3,89
!!$          rhoxx = zero
!!$          rhoxx(-2:2,:) = rho2pdf(:,:)
!!$          do i = 3, nfint
!!$             rhoxx( i,:) = rhoxx( i-2,:)
!!$             rhoxx(-i,:) = rhoxx(-i+2,:)
!!$          end do
!!$          
!!$          ffactin = (/xi(ir), 20.0_dp, 20.0_dp, 20.0_dp/)
!!$          call disaster_cff(1,3,89,1,2,0.3_dp,2,nfint,ffactin, ffactout)
!!$          call pftopdg(xi(ir), 20.0_dp, dxpdf)
!!$          dxpdf = dxpdf
!!$          do irho = 1, 11
!!$             write(0,*) ir, irho, ffactout(irho),&
!!$                  & sum(rhoxx(:,irho)*dxpdf)/xi(ir)
!!$          end do
!!$          write(0,*) '------',sum(ffactout(1:11)*weight(:,il))&
!!$               & ,sum(weight_pdf*dxpdf)/( (twopi)**na * alphaem&
!!$               & **ialphaem(il))
          
          !-- call a routine similar to the usual disent one
          call dh_event_list(np, na, parr, photon, xb(ir), xi(ir), q2(ir),&
               & weight_pdf, ir_already_used, aobs)
          dh_user2 = dh_user2 + aobs
          !-- this way the above routine can avoid recalculating
          !   event-shapes and x-dependent things
          ir_already_used = .true.
       end do
    end do
    
    if (detected_nan) then
       write(0,*) 'WARNING from dh_user2: a NaN &
            &has been detected and the event will be thrown away'
       write(0,*) 'Event number is:', nevt
       nevt_nan = nevt_nan + 1
       write(0,*) 'Number of events thrown away so far is', nevt_nan
       dh_user2 = zero
       nevt = nevt - 1
       do irun = 1, size(run)
          call rd_zero_tmp_list(run(irun))
       end do
    end if
    

    if (final) then
       !-- this is soooo ugly, but it is how run_descriptor was "designed"
       rd_onejet = .not. jetrun
       rd_twojet =       jetrun
       do irun = 1, size(run)
          call rd_collate_tmp_list(run(irun))
       end do
    end if
    


    !-- arrange for the occasional output of the results
    if (final .and. jetrun .and. mod(nevt,writefreq)==0) then
!    if (final .and. jetrun .and. mod(nevt,1000)==0) then
       call dh_final(nevt, real(nevt,dp)/nevt_wanted, interim=.true.)
    end if
    
  end function dh_user2
  




  !----------------------------------------------------------------------
  ! This part does the hard work of receiving a disent event and doing
  ! something with it. This version uses a list of bins
  subroutine dh_event_list(n, na, p, photon, xb, xidstr, q2, weight,&
       & config_already_used, aobs) 
    use check_nan
    use event_shapes, ldot => dot; use pdf_manager
    use g90bk
    integer,  intent(in)    :: n, na
    real(dp), intent(in)    :: p(4,3), photon(4), xb, xidstr, q2, weight(-6:6)
    logical,  intent(in)    :: config_already_used
    real(dp), intent(inout) :: aobs
    !------ work vars ------------------------------------
    real(dp) :: EoverQ
    integer  :: i, irun, ibin, naxx 
    integer, pointer :: nbin
    integer :: istat, ies, ipdf, ixQ, ixsct, iord
    type(hgrp), pointer :: hh, hlims
    real(dp) :: esvp, xi
    type(run_def), pointer :: runi


    xi = xidstr/xb

    !-- if have multiple runs, then the loop starts here
    do irun = 1, size(run)
       runi => run(irun)

       if (.not. config_already_used)  then
          !====== start of semi-independent stuff =========================
          if (jetrun) then
             !-- ensure that event-shape cache is started off properly --
             call es_cache_new
             !-- get the event shapes ------------------------------------
             call es_manyvals(p(:,1:n), photon, runi%es_list%id, runi&
                  & %es_vals)

             !-- check for NANs ---------------------------------------------
             ! should cause event to be thrown away?
             detected_nan = detected_nan .or. any(is_a_nan(runi%es_vals))

             EoverQ = es_val(p(:,1:n), photon, esid_CrntEQ)
             !-- take logs where necessary (would like nested where, but
             !   not in f90)
             do ies = 1, runi%nes
                if (runi%es_list(ies)%log) then
                   if (runi%es_vals(ies) > zero) then
                      runi%es_vals_proc(ies) = log(runi%es_vals(ies))
                   else
                      runi%es_vals_proc(ies) = es_undefined_value
                   end if
                else
                   runi%es_vals_proc(ies) = runi%es_vals(ies)
                end if
             end do
          end if

          !--- get the pdfs ----------------
          call pm_eval_pdfs(runi%pdfs, xi)
       end if
    
       !-- the rest always happens even for configs that have already been
       !   analysed.
       runi%weights = zero
       do i = -6, 6
          runi%weights = runi%weights&
               &           + weight(i)*runi%pdfs%pdf_vals(: ,:,i)
       end do

       !-- do a safety check here --------------------------
       !   not necessarily the most efficient way of doing things?
       ! Maybe should give up if any of them is bad? Yes
       detected_nan = detected_nan .or. any(is_a_nan(runi%weights))

       !-- just need to calculate adaption observable -- do it elsewhere
       !   uses just run(1)
       if (.not. final) exit

       !--- total cross sections -----------------------------------------
       if (na < 2 .and. (.not. jetrun) ) then
          runi%x(:,:,1,na+1)%tmp = runi%x(:,:,1,na+1)%tmp + runi%weights
       end if

       !--- stuff related to final states ------------
       if (na >= 1  .and. jetrun) then
          !-- NB use >= so that 0.0 means no limit...
          if (EoverQ >= runi%Elim) then
             !-- first determine which bins will be used
             ! following line allows change from na=fixed, to variable na
             ! remember to make corresponding change in run_descriptor
             ! (i.e. determines we keep a single list of bins used (naxx=1)
             !  or a separate list for each na value (naxx=na))
             naxx = na
             do ies = 1, runi%nes
                ibin = hist_ibin(runi%h(1,1,ies,naxx)%htmp,&
                     & runi%es_vals_proc(ies))
                !-- following is signal that a valid bin has been chosen
                if (ibin >  0) then
                   nbin => runi%nbin(ies, naxx)
                   if (nbin < size(runi%bin_list,dim=1)) then
                      nbin = nbin + 1
                      runi%bin_list(nbin, ies, naxx) = ibin
                      runi%binit(ies) = .true.
                   else
                      write(0,*) 'OVERRAN MAX # OF BINS IN LIST'
                      stop
                   end if
                else
                   runi%binit(ies) = .false.
                end if
!!$                !-- start of alt version (turned out to be marginally slower
!!$                do ipdf = 1, runi%npdf
!!$                   do ixQ = 1, runi%nxQ
!!$                      !-- for brevity (hope it does not take too much time? --
!!$                      hh    => runi%h(ipdf,ixQ, ies, na )
!!$                      hlims => runi%hlims(ipdf,ixQ, ies, na )
!!$                      !-- distribution --
!!$                      if (runi%binit(ies)) then
!!$                         hh%htmp%h(ibin) = hh%htmp%h(ibin) +&
!!$                           & runi%weights(ipdf,ixQ) / hh%htmp%spc
!!$                      end if
!!$                      !-- mean value ----
!!$                      hh%tmp = hh%tmp +&
!!$                           & runi%weights(ipdf,ixQ)*runi%es_vals(ies)
!!$                      !-- part above certain value --
!!$                      if (runi%es_vals(ies) > runi%es_list(ies)%lim)&
!!$                           &hlims%tmp = hlims%tmp + runi%weights(ipdf,ixQ)
!!$                   end do
!!$                end do
!!$                !runi%h(:,:, ies, na )%tmp = runi%h(:,:, ies, na )%tmp + &
!!$                !     &runi%weights(:,:)*runi%es_vals(ies)
!!$                !if (runi%es_vals(ies) > runi%es_list(ies)%lim) then
!!$                !   runi%hlims(:,:, ies, na )%tmp &
!!$                !        &= runi%hlims(:,:, ies, na )%tmp + runi%weights(:,:)
!!$                !end if
!!$                !-- end of alt version
             end do
             

             do ipdf = 1, runi%npdf
                do ixQ = 1, runi%nxQ
                   do ies = 1, runi%nes
                      !-- for brevity (hope it does not take too much time? --
                      hh    => runi%h(ipdf,ixQ, ies, na )
                      hlims => runi%hlims(ipdf,ixQ, ies, na )
                      !-- distribution --
                      if (runi%binit(ies)) then
                         ibin = runi%bin_list(runi%nbin(ies, naxx), ies,naxx)
                         hh%htmp%h(ibin) = hh%htmp%h(ibin) +&
                           & runi%weights(ipdf,ixQ) / hh%htmp%spc
                      end if
                      !-- mean value ----
                      hh%tmp = hh%tmp +&
                           & runi%weights(ipdf,ixQ)*runi%es_vals(ies)
                      !-- part above certain value --
                      if (runi%es_vals(ies) > runi%es_list(ies)%lim)&
                           &hlims%tmp = hlims%tmp + runi%weights(ipdf,ixQ)
                   end do
                end do
             end do
          else
             !-- the probability of being below Elim
             runi%x(:,:,2,na)%tmp = runi%x(:,:,2,na)%tmp + runi%weights
          end if
       end if
       !-- loop over runs
    end do


    !------------------------------------------------------------------
    ! Things for adaptation: do them always since they do not cost much
    runi => run(1)
    if (.not. jetrun) then
       if (na == 1) aobs = runi%weights(1,1)
    else
       !-- use adaption observable which consists of probability of
       !   being binned or of there being a low-energy current
       !   hemisphere + twice the mean value
       !   (HOW UGLY DOES ONE GET?)
       hh   => runi%h(ipdf_adapt,1,ivar_adapt,na)
       esvp =  runi%es_vals_proc(ivar_adapt)
       if (esvp >= hh%htmp%min .and. esvp <= hh%htmp%max&
            & .and. EoverQ >= runi%Elim) then
          !-- put in some effective value for na
          select case(adapt)
          case(adapt_one)
             aobs = one
          case(adapt_sqrt)
             aobs = two*sqrt(abs(esvp))
          case(adapt_lin)
             aobs = four * esvp
          !-- the following two are designed for fixed-order tests
          case(adapt_abs)
             aobs = abs(esvp) / (0.4d0/twopi)**min(1,na-1)
          case(adapt_abs2)
             !-- multiply by 0.2 so that at LH edge it is not
             !   to large compared to some other contributions?
             aobs = 0.2_dp * abs(esvp)**2 / (0.4d0/twopi)**min(1,na-1)
          case default
             write(0,*) 'ERROR in dh_event_list: unrecognized adapt val:',adapt
             stop
          end select
       else
          aobs = zero
       end if
       if (EoverQ < runi%Elim) aobs = aobs + 6.0_dp
       if (EoverQ >= runi%Elim) aobs = aobs + 6.0_dp * runi%es_vals(ivar_adapt)
       aobs = aobs * runi%weights(ipdf_adapt,1) * (0.4d0/twopi)**na
       !-- try to avoid problems with NaNs? ----------
       if (is_a_nan(aobs)) detected_nan = .true.
    end if
    

  end subroutine dh_event_list


  !----------------------------------------------------------------------
  ! Used by disent to establish the cuts that will be used.
  subroutine dh_cuts
    if (pm_smallest_x >= one) then
       stop 'BAD pm_smallest_x in dh_cuts'
    end if
    
    call disaster_cd('XBMIN', pm_smallest_x * (one - 1e-4_dp)) 
    call disaster_cd('XBMAX', pm_smallest_x) 

    
    if (dh_ybj < zero) then
       stop 'BAD dh_ybj in dh_cuts'
    end if
    
    call disaster_cd('YMIN', dh_ybj * (one - 5e-5_dp) ) 
    call disaster_cd('YMAX', dh_ybj * (one + 5e-5_dp) ) 
  end subroutine dh_cuts
  
  !------- return the disent seed ------------------------------
  !        requires a modified RANGEN which when called with N=0
  !        and NINT(R(1)=0) return the seed rather than the usual
  !        behaviour of setting it. The code is given here:
  !      ELSEIF (N.EQ.0) THEN
  !        IF(NINT(R(1)) .eq. 0) then
  !           !-- retrieve seed --
  !           R(1) = ISEED(1)
  !           R(2) = ISEED(2)
  !        else
  !           !-- set seed -------
  !           ISEED(1)=NINT(R(1))
  !           ISEED(2)=NINT(R(2))
  !        end if
  !      ENDIF
  function dh_GetSeed() result(seed)
    integer :: seed(3)

    call rm48ut(seed(1), seed(2), seed(3))
  end function dh_GetSeed
  


end module disaster_helper




!======================================================================
! Standard name by which to access either the disaster_helper or
! disent_helper 
module dismc_helper
  use disaster_helper
end module dismc_helper






! --------------------------------------------------------------------  
!                                                                       
! ALL THE BITS AND PIECES NEEDED BY DISASTER
!
! THESE MUST BE KEPT OUTSIDE OF A MODULE, SINCE THEY MUST BE
! ACCESSIBLE FROM C++.
!                                                                       
! --------------------------------------------------------------------  
                                                                        
subroutine user1(iaction) 
  use types; use consts_dp; use disaster_helper
  implicit none 
  
  integer iaction 
  
  integer i 
  
  ! write(6,*) 'iaction=', iaction                                   
  
  !-- end of run ---
  if (iaction .eq. 2) then 
     
     
  !-- vegas configuration run 
  elseif (iaction .eq. 3) then 
     write(0,*) 'Doing a vegas configuration run'
     final = .false.
     nevt = 0

  !-- normal run (after vegas adaptation)
  elseif (iaction .eq. 4) then 
     write(0,*) 'Doing a vegas final run'
     final = .true.
     nevt = 0 
  endif 
      
  return 
end subroutine user1


!---------------------------------------------
function user2( nr, nl, fvect, npartons, xb, q2, xi, weight, irps,&
     & ialphas, ialphaem, lognf )
  use types; use disaster_helper
  implicit none
  real(dp) user2 
  
  integer nr, nl 
  real(dp) fvect(0:3, -10:10, 1:*) 
  integer npartons(1:*) 
  real(dp) xb(1:*), q2(1:*), xi(1:*), weight(1:11, 1:*) 
  integer irps(1:*), ialphas(1:*), ialphaem(1:*), lognf(1:*)    

  user2 = dh_user2(nr, nl, fvect, npartons, xb, q2, xi, weight, irps,&
     & ialphas, ialphaem, lognf )
end function user2


! --------------------------------------------------------------------  
! needed for the adaptive integrator
subroutine user3(average, error) 
  use types
  implicit none 
  real(dp) average, error 
  
100 format('Result of Integration =', 1e16.6, ' +- ', 1e16.6) 
  write(6, 100) average, error 
  return 
end subroutine user3


! --------------------------------------------------------------------  
!                                                                       
! --> printing routine (char* from C++)                                 
!                            
!
! GPS: modified in order to be able to send C++ info to stderr
! --------------------------------------------------------------------  
subroutine print_f(string, length, unit) 
  
  implicit none 

  character*1000 string 
  integer length, unit 

  character*7 fmts 

  integer i1, i2, i3 

  integer :: lunit

  character*1 number(0:9) 
  data                                                              &
       &   number(0) /'0'/,                                               &
       &   number(1) /'1'/,                                               &
       &   number(2) /'2'/,                                               &
       &   number(3) /'3'/,                                               &
       &   number(4) /'4'/,                                               &
       &   number(5) /'5'/,                                               &
       &   number(6) /'6'/,                                               &
       &   number(7) /'7'/,                                               &
       &   number(8) /'8'/,                                               &
       &   number(9) /'9'/                                                


  if (unit == 6) then
     lunit = 41
  else
     lunit = unit
  end if

  i3 = length / 100 
  i2 = (length - i3 * 100) / 10 
  i1 = (length - i3 * 100 - i2 * 10) 

  if (length .eq. 0) then 

     write(UNIT=lunit, FMT=*) 

  else 

     fmts(1:1) = '(' 
     fmts(2:2) = '1' 
     fmts(3:3) = 'A' 
     fmts(4:4) = number(i3) 
     fmts(5:5) = number(i2) 
     fmts(6:6) = number(i1) 
     fmts(7:7) = ')' 

     write(UNIT=lunit, FMT=fmts) string(1:length) 

  endif

  return 
end subroutine print_f                                  
                                                                        

! --------------------------------------------------------------------  
!                                                                       
! --> provoke a floating point exception                                
!                                                                       
! --------------------------------------------------------------------  
subroutine pfpx_f() 

  implicit none 

  double precision d 

  write(6,*) 'provoked FPE (FORTRAN)' 
  d = dsqrt(log(0.5d0)) 
  write(6,*) 'FPE (FORTRAN) value = ', d 

  return 
END subroutine pfpx_f
