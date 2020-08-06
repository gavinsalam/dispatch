!======================================================================
! The common routines that are of use to disent
module disent_helper
  use types; use consts_dp
  use run_descriptor
  use event_shapes!, except => dot
  use pdf_manager
  use sub_defs_io
  implicit none

  !== lower limit on number of events to have in an interation
  integer, public :: nev_iter_limit = 10000

  private
  !type(run_def), allocatable, target :: run(:)
  type(run_def), pointer :: run(:)
  !type(run_def), target :: run(1)



  public :: dh_start!, dh_final, dh_event
  public :: run

  logical, allocatable :: nzmask(:,:)

  !-- the following are public to allow compatibility with older programs
  real(dp), public :: dh_ybj = -one
  public :: dh_cuts

  integer, parameter :: nseed = 3
  integer :: seed_init(nseed)


contains

  !----------------------------------------------------------------------
  ! Once run(:) has been set up, then call here
  subroutine dh_start(nf, ybj)
    integer,  intent(in) :: nf
    real(dp), intent(in) :: ybj
    INTEGER  :: NEV, NEV_INT, IN, iseq, irun, writefreq
    integer  :: ndiv = 20,  ORDER
    real(dp) :: S, NPOW1,NPOW2,CUTOFF, cf, ca, tr, SCALE, muF_Q
    character(len=70) :: tdfl

    
    ! --- redirect all output to files with names similar to tdfile.
    tdfl = trim(string_val_opt('-tdfl',''))
    if (tdfl /= '' ) then
       open(FILE=trim(tdfl)//'.err', UNIT=0, STATUS='UNKNOWN') 
       open(FILE=trim(tdfl)//'.out', UNIT=6, STATUS='UNKNOWN') 
    end if


    !-- tell user what is happening...
    call rd_message()

    call rd_parseiseed(seed_init)
    iseq = int_val_opt('-iseq',1)
    if (seed_init(1) > 0 .and. iseq == seed_init(1)) then
       !if (seed_init(1) > 0) then
       write(0,*) 'Using last seed from previous run'
       write(0,*) 'Seed:', seed_init
       call dh_SetSeed(seed_init)
    else
       call dh_SetSeed((/iseq,0,0/))
       write(0,*) 'Using default seed with iseq =', iseq
    end if
    write(0,*) '--------------------------------------------------------'
    
    !-- now see what the actual seed is...
    seed_init = dh_GetSeed()

    S=4*30*820 
    
    NEV=nint(dble_val_opt('-nev',1e6_dp15))
    NPOW1 = int_val_opt('-npow1',2)
    NPOW2 = int_val_opt('-npow2',4)
    CUTOFF = dble_val_opt('-cutoff',1e-12_dp15)
    ORDER = int_val_opt('-order',2)
    !-- test that we always have the same renorm scale for DISENT.
    muF_Q = run(1)%pdfs%xQ(1)%muF_Q
    do irun = 1, size(run)
       if (any(run(irun)%pdfs%xQ(:)%muF_Q /= muF_Q)) then
          write(0,*) 'ERROR: Disent driver requires all muF_Q values to be the same'
          write(0,*) 'This is not satisfied in run', irun
          stop
       end if
    end do
    SCALE = muF_Q**2

    cf = dble_val_opt('-cf',4.0_dp/3.0_dp)
    ca = dble_val_opt('-ca',3.0_dp)
    tr = dble_val_opt('-tr',half)
    
    !-- ensure that things get written about once an hour by default.
    writefreq = nint(dble_val_opt('-wf',5e6_dp))
    ndiv = max(ndiv, nev / writefreq)
    NEV_INT = NEV / ndiv
    if (NEV_INT < nev_iter_limit) then
       ndiv = ceiling(real(nev)/10000)
       nev_int = nev/ndiv
       write(0,*) 'NDIV changed to',ndiv
    end if

    dh_ybj = ybj
    
    if (nf /= 5) write(0,*) 'WARNING: NF .NE. 5'

    if (.not. CheckAllOptsUsed(ErrDev=0)) stop

    do in = 1, ndiv
       write(6,'(A,I2,A,I4)') 'STARTING EVOLUTION SET ',IN,' OF ',NDIV
       !CALL DISENT(NEV_INT,S,nf,dh_event,dh_cuts,NPOW1,NPOW2,CUTOFF,ORDER) 
       CALL DISENT(NEV_INT,S,nf,dh_event_list,dh_cuts,&
            &NPOW1,NPOW2,CUTOFF,SCALE,ORDER,cf,ca,tr) 
       CALL dh_final(NEV_INT*in,in) 
    end do
    
    
  end subroutine dh_start
  


  !----------------------------------------------------------------------
  ! Outputs the results
  ! nev is total number of events so far; in is the number of divisions so
  ! far (each division sums to total cross section)
  subroutine dh_final(nev,niter)
    use g90bk
    integer, intent(in) :: nev,niter
    integer :: irun, dev
    
    !dev = hist_tdopen() ! no need to call this: it gets done automatically

    !call hist_tdrewind   ! not needed since we have tdclose at the end

    call dh_DescribeDetails
    call hist_tdtext(' NEV = '//num2char(nev),'(X ')
    call hist_tdtext('******************************************************')
    call hist_tdtext('ybj = '//num2char(dh_ybj))
    do irun = 1, size(run)
       call rd_write_res(run(irun), nev, real(niter,dp), nev, real(niter,dp))
    end do

    call hist_tdclose   ! try to ensure that results will be available

  end subroutine dh_final


  !------- say what sort of run this was: DISENT specific details?
  subroutine dh_DescribeDetails
    use g90bk
    !-- common blocks from disent
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
    !---
    character(len=*), parameter :: cc='(X '
    character(len=100) :: line
    integer :: seed(nseed)
    !-- first some things related to the technical details of the run
    write(line,'(a,i2)') 'DISENT RUN -- format version ', rd_format_version
    call hist_tdtext(trim(line),cc)
    write(line,*) 'Initial ISEED was ',seed_init
    call hist_tdtext(line,cc)
    seed = dh_GetSeed()
    write(line,*) 'Current ISEED is  ',seed
    call hist_tdtext(line,cc)
    call hist_tdtext(' cutoff = '//num2char(cutoff), cc)
    call hist_tdtext(' scheme = '//num2char(scheme), cc)
    !-- not in format version >= 2
    !call hist_tdtext(' scale  = '//num2char(scale ), cc)
    write(line,*) 'NPOW = ', NPOW
    call hist_tdtext(line,cc)
    write(line,&
         & '(" CA =",f10.6,"  CF =",f10.6,"  TR =",f10.6,"  nf =",i3)')&
         & ca,cf,tr,nf
    call hist_tdtext(line,cc)

  end subroutine dh_DescribeDetails
  

  !----------------------------------------------------------------------
  ! This part does the hard work of receiving a disent event and doing
  ! something with it.
  subroutine dh_event(n, na, nt, p, s, weight)
    use event_shapes, ldot => dot; use pdf_manager
    use g90bk
    integer,  intent(in) :: n, na, nt
    real(dp), intent(inout) :: p(4,7), s, weight(-6:6)
    INTEGER  :: SCHEME,NF 
    real(dp) :: CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE 
    COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF 
    !------ work vars ------------------------------------
    real(dp) :: eta,Q2,x,y,a2pi, rtQ2
    real(dp) :: xi, EoverQ
    integer  :: i, irun
    integer :: istat, ies, ipdf, ixQ, ixsct, iord
    type(hgrp), pointer :: hh, hlims
    real(dp) :: dot
    type(run_def), pointer :: runi

    !-- to be done at the end of a bunch of events ----------
    !   to ensure cancellation of large positive and negative weights in
    !   the calculation of the errors
    if (n == 0) then
       do irun = 1, size(run)
          call rd_collate_tmp(run(irun))
       end do
       return
    end if

    !--- calculate momentum fractions
    eta=2*ldot(p,1,6)/s      ! momentum fraction of incoming parton
    q2=abs(ldot(p,5,5))      ! Bjorken Q2
    x=eta*q2/(2*ldot(p,1,5)) ! Bjorken x
    y=ldot(p,1,5)/ldot(p,1,6) ! Bjorken y
    !write(0,*) n-1,na, eta

    !-- for the pdf weights ------------------
    xi = eta / x            ! x_incoming/ x_Bj
    !write(0,*) xi

    !====== start of semi-independent stuff =========================
    !-- ensure that event-shape cache is started off properly --
    call es_cache_new

    !-- if have multiple runs, then the loop starts here
    do irun = 1, size(run)
       runi => run(irun)  
       !-- get the event shapes ------------------------------------------
       call es_manyvals(p(:,2:n), p(:,5), runi%es_list%id, runi%es_vals)
       EoverQ = es_val(p(:,2:n), p(:,5), esid_CrntEQ)
       !-- take logs where necessary (would like nested where, but not in
       ! f90)
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

       !--- get the weights: DISENT/DISASTER dependent -------------------
       call pm_eval_pdfs(runi%pdfs, xi)
       runi%weights = zero
       do i = -6, 6
          runi%weights = runi%weights + weight(i)*runi%pdfs%pdf_vals(:,:,i)
       end do



       !--- total cross sections -----------------------------------------
       if (na < 2) then
          runi%x(:,:,1,na+1)%tmp = runi%x(:,:,1,na+1)%tmp + runi%weights
       end if

       !--- stuff related to final states ------------
       if (na >= 1) then
          !-- NB use >= so that 0.0 means no limit...
          if (EoverQ >= runi%Elim) then
             do ipdf = 1, runi%npdf
                do ixQ = 1, runi%nxQ
                   do ies = 1, runi%nes
                      !-- for brevity (hope it does not take too much time? --
                      hh    => runi%h(ipdf,ixQ, ies, na )
                      hlims => runi%hlims(ipdf,ixQ, ies, na )
                      !-- distribution --
                      call hist_adddat(hh%htmp,runi%es_vals_proc(ies), &
                           &runi%weights(ipdf,ixQ) / hh%htmp%spc)
                      !-- mean value ----
                      hh%tmp = hh%tmp +&
                           & runi%weights(ipdf,ixQ)*runi%es_vals(ies)
!!$                      !-- A TEMPORARY HACK --
!!$                      ! try to get probability of having just gluon
!!$                      if (na == 1 .and. n==3 .and. p(3,2) > zero .and. p(3&
!!$                           & ,3) < zero) then
!!$                         hlims%tmp = hlims%tmp + runi%weights(ipdf,ixQ)
!!$                      end if
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

  end subroutine dh_event



  !----------------------------------------------------------------------
  ! This part does the hard work of receiving a disent event and doing
  ! something with it. This version uses a list of bins
  subroutine dh_event_list(n, na, nt, p, s, weight)
    use event_shapes, ldot => dot; use pdf_manager
    use g90bk
    integer,  intent(in) :: n, na, nt
    real(dp), intent(inout) :: p(4,7), s, weight(-6:6)
    INTEGER  :: SCHEME,NF 
    real(dp) :: CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE 
    COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF 
    !------ work vars ------------------------------------
    real(dp) :: eta,Q2,x,y,a2pi, rtQ2
    real(dp) :: xi, EoverQ
    integer  :: i, irun, ibin, naxx 
    integer, pointer :: nbin
    integer :: istat, ies, ipdf, ixQ, ixsct, iord
    type(hgrp), pointer :: hh, hlims
    real(dp) :: dot
    type(run_def), pointer :: runi

    !-- to be done at the end of a bunch of events ----------
    !   to ensure cancellation of large positive and negative weights in
    !   the calculation of the errors
    if (n == 0) then
       do irun = 1, size(run)
          call rd_collate_tmp_list(run(irun))
       end do
       return
    end if

    !--- calculate momentum fractions
    eta=2*ldot(p,1,6)/s      ! momentum fraction of incoming parton
    q2=abs(ldot(p,5,5))      ! Bjorken Q2
    x=eta*q2/(2*ldot(p,1,5)) ! Bjorken x
    y=ldot(p,1,5)/ldot(p,1,6) ! Bjorken y
    !write(0,*) n-1,na, eta

    !-- for the pdf weights ------------------
    xi = eta / x            ! x_incoming/ x_Bj
    !write(0,*) xi

    !====== start of semi-independent stuff =========================
    !-- ensure that event-shape cache is started off properly --
    call es_cache_new

    !-- if have multiple runs, then the loop starts here
    do irun = 1, size(run)
       runi => run(irun)  
       !-- get the event shapes ------------------------------------------
       call es_manyvals(p(:,2:n), p(:,5), runi%es_list%id, runi%es_vals)
       EoverQ = es_val(p(:,2:n), p(:,5), esid_CrntEQ)
       !-- take logs where necessary (would like nested where, but not in
       ! f90)
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

       !--- get the weights: DISENT/DISASTER dependent -------------------
       call pm_eval_pdfs(runi%pdfs, xi)
       runi%weights = zero
       do i = -6, 6
          runi%weights = runi%weights + weight(i)*runi%pdfs%pdf_vals(:,:,i)
       end do



       !--- total cross sections -----------------------------------------
       if (na < 2) then
          runi%x(:,:,1,na+1)%tmp = runi%x(:,:,1,na+1)%tmp + runi%weights
       end if

       !--- stuff related to final states ------------
       if (na >= 1) then
          !-- NB use >= so that 0.0 means no limit...
          if (EoverQ >= runi%Elim) then
             !-- first determine which bins will be used
             ! following line allows change from na=fixed, to variable na
             ! remember to make corresponding change in run_descriptor
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
             end do
             
             do ipdf = 1, runi%npdf
                do ixQ = 1, runi%nxQ
                   do ies = 1, runi%nes
                      !-- for brevity (hope it does not take too much time? --
                      hh    => runi%h(ipdf,ixQ, ies, na )
                      hlims => runi%hlims(ipdf,ixQ, ies, na )
                      !-- distribution --
                      if (runi%binit(ies)) then
                         ibin = runi%bin_list(runi%nbin(ies, naxx),ies,naxx)
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

  end subroutine dh_event_list


  !----------------------------------------------------------------------
  ! Used by disent to establish the cuts that will be used.
  subroutine dh_cuts(S,XMIN,XMAX,QMIN,QMAX,YMIN,YMAX)
    real(dp) :: S,XMIN,XMAX,QMIN,QMAX,YMIN,YMAX 
    !-- ensure that this is always suitable.
  
    if (pm_smallest_x >= one) then
       stop 'BAD pm_smallest_x in dh_cuts'
    end if
    XMIN = pm_smallest_x
    XMAX = pm_smallest_x
    
    if (dh_ybj < zero) then
       stop 'BAD dh_ybj in dh_cuts'
    end if
    
    YMIN = dh_ybj
    YMAX = dh_ybj
    
    !-- any Q value is good
    QMIN=0
    QMAX=1d40
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
    integer :: seed(nseed)
    real(dp) :: rseed(nseed)

    rseed = zero
    call RANGEN(0,rseed)
    seed = nint(rseed)
    if (seed(1) == 0) then
       stop 'ERROR: dh_GetSeed assumes modified disent RANGEN; see program&
            & for details'
    end if
!!$    write(0,*) 'WARNING, using modified dh_GetSeed for compatibility &
!!$         &with old rangen'
!!$    seed = 0
  end function dh_GetSeed
  

  subroutine dh_SetSeed(iseed)
    integer, intent(in) :: iseed(nseed)
    call RANGEN(0,real(iseed,kind=dp))
  end subroutine dh_SetSeed
  
end module disent_helper
