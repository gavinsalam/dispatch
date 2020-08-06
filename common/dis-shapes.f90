!======================================================================
!
! Nomenclature:
!      Tax = thrust axis
!      Gax = Gamma axis
!
!      E = normalised to E_vis
!      Q = normalised to Q
!
! Need to think a little about how to introduce Elim: as a global
! parameter with the option of modifying with an optional argument 
!--------------------------------------------------------------------
module event_shapes
  use types; use consts_dp; implicit none
  private

  !------- the ids for the different event shapes  ------
  integer, parameter, public :: esid_illegal = -1
  integer, parameter, public :: esid_tau_zQ = 1
  integer, parameter, public :: esid_tau_zE = 2
  integer, parameter, public :: esid_tau_tQ = 3
  integer, parameter, public :: esid_tau_tE = 4
  integer, parameter, public :: esid_brd_zE = 5
  integer, parameter, public :: esid_brd_tE = 6
  integer, parameter, public :: esid_rho_E  = 7
  integer, parameter, public :: esid_Cpr_E  = 8
  integer, parameter, public :: esid_CrntEQ = 9  ! current energy norm Q
  integer, parameter, public :: esid_glninC = 10
  integer, parameter, public :: esid_rho23E = 11
  integer, parameter, public :: n_esid = 11

  real(dp), parameter, public :: es_undefined_value = -1e90_dp

  !------- the names for the different event shapes ------
  character(len=6), parameter, public :: es_names(n_esid) = (/&
       &'tau_zQ',&
       &'tau_zE',&
       &'tau_tQ',&
       &'tau_tE',&
       &'brd_zE',&
       &'brd_tE',&
       &'rho_E ',&
       &'Cpr_E ',&
       &'CrntEQ',&
       &'glninC',&
       &'rho23E'/)     

  real(dp), public :: store(2,2)


  !--- allow caching ----------
  logical :: event_cache_on = .false.
  integer :: event_id       = 0, id_save = -1
  real(dp) :: es_res_cache(n_esid)


  !-- keep these private?
  interface dot
     module procedure dot_m, dot_2v
  end interface
  !-- not really what I want: should be thought about a little...
  public :: dot

  ! not used?
  !real(dp) :: Elim = 1e-19_dp

  logical, public :: es_verbose = .false.


  public :: es_manyvals, es_val

  public :: es_cache_new, es_cache_off


contains

  !======================================================================
  function es_val(parr, Qvec, esid) result(res)
    real(dp), intent(in)  :: parr(:,:), Qvec(:)
    integer,  intent(in)  :: esid
    real(dp)              :: res

    if ((.not. (event_cache_on)) .or. id_save /= event_id)&
         & call es_buildcache(parr, Qvec)
    if (esid <= n_esid .and. esid > 0) then
       res = es_res_cache(esid)
    else
       write(0,*) 'Illegal esid:', esid
    end if
  end function es_val


  subroutine es_manyvals(parr, Qvec, esid_list, esvals_list)
    real(dp), intent(in)  :: parr(:,:), Qvec(:)
    integer,  intent(in)  :: esid_list(:)
    real(dp), intent(out) :: esvals_list(:)
    !--------------
    integer :: i, esid

    if ((.not. (event_cache_on)) .or. id_save /= event_id)&
         & call es_buildcache(parr, Qvec)
    do i = 1, size(esid_list)
       esid = esid_list(i)
       if (esid <= n_esid .and. esid > 0) then
          esvals_list(i) = es_res_cache(esid)
       else
          write(0,*) 'Illegal esid:', esid
       end if
    end do
  end subroutine es_manyvals
  
  !======================================================================
  ! paradoxically it seems to be quicker to caculate all event shapes
  ! and then access them as needed, rather than to have separate routines
  ! which get called and do caching
  subroutine es_buildcache(parr, Qvec)
    real(dp), intent(in)  :: parr(:,:), Qvec(:)
    !------------------------------------------
    !real(dp), pointer :: pcur(:,:)
    real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
    real(dp)          :: Ecur
    integer           :: n, i, j
    !----------------------------------------
    real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4)
    real(dp) :: cpr, costheta, glninC


    !-- tell cache that this event will be cached.
    id_save = event_id

    !-- start off cleanly: if anything is forgotten then it will be undefined
    es_res_cache = es_undefined_value

    !- if uncommented, following line needs pcur to be a pointer
    !call es_select_current(parr, Qvec, pcur, Ecur, n)
    n = 0; Ecur = zero
    do i = 1, size(parr,dim=2)
       if (parr(3,i)*Qvec(3) > zero) then
          n = n+1
          pcur(:,n) = parr(:,i)
          Ecur = Ecur + parr(4,i)
       end if
    end do
    

    !--- assume that photon is along 3-direction
    EoverhalfQ = Ecur / (half*abs(Qvec(3)) )
    !------ tau_z ----
    if (n == 0) then
       thr = 0
    else
       thr = sum(pcur(3,:n)) / (half * Qvec(3))
    end if
    es_res_cache(esid_tau_zQ) = one - thr 
    if (Ecur /= 0) es_res_cache(esid_tau_zE) = one - thr / EoverhalfQ
    !------ tau_t ----
    call es_thr_tQ_selctd(pcur, Ecur, Qvec, n, thr, a)
    es_res_cache(esid_tau_tQ) = one - thr
    if (Ecur /= 0) es_res_cache(esid_tau_tE) = one - thr / EoverhalfQ
    !------ brd_tE ---
    select case(n)
    case(0)
       brd = es_undefined_value
    case(1)
       brd = zero
    case default
       brd = zero
       do i = 1, n
          brd = brd + sqrt(mod3sq(pcur(:3,i) - a*dot_product(a,pcur(:3,i))))
       end do
       brd = half*brd / Ecur
    end select
    es_res_cache(esid_brd_tE) = brd
    !------ brd_zE ---
    select case(n)
    case(0)
       brd = es_undefined_value
    case default
       brd = zero
       do i = 1, n
          !brd = brd + sqrt(mod3sq(pcur(:3,i) - a*dot_product(a,pcur(:3,i))))
          !-- assume z axis is along 3 direction!
          brd = brd + sqrt(pcur(1,i)**2 + pcur(2,i)**2)
       end do
       brd = half*brd / Ecur
    end select
    es_res_cache(esid_brd_zE) = brd
    !------ rho_E ----
    select case(n)
    case(0)
       rho = es_undefined_value
    case default
       ptot = sum(pcur(:4,:n),dim=2)
       rho  = max(zero,dot(ptot,ptot))
       !write(0,'(4f15.8)') real(pcur(:4,:n))
       rho = rho / (two*Ecur)**2
    end select
    es_res_cache(esid_rho_E) = rho
    !------ cpr_E ----
    select case(n)
    case(0)
       Cpr = es_undefined_value
    case(1)
       Cpr = zero
    case default
       Cpr = zero
       do i = 1, n-1
          do j = i+1, n
             costheta = dot_product(pcur(:3,i),pcur(:3,j)) /&
                                         & (pcur(4,i)*pcur(4,j))
             Cpr = Cpr + (pcur(4,i)*pcur(4,j))*(one - costheta**2)
          end do
       end do
       !--- factor of half cancelled by not doing symmetrized sum
       Cpr = three * Cpr / (Ecur**2)
    end select
    es_res_cache(esid_cpr_E ) = cpr
    !------ crntEQ ----
    es_res_cache(esid_CrntEq) = Ecur/(abs(Qvec(3)))
    !------ glninC ----
    glninC = zero
    if (size(parr,dim=2) == 2) then
       if (parr(3,2)*Qvec(3) > 0 .and. parr(3,1)*Qvec(3) < 0) then
          glninC = one
       end if
    end if
    es_res_cache(esid_glninC) = glninC
    !------ rho23E ----
    es_res_cache(esid_rho23E) = zero
    if (size(parr,dim=2) == 3) then
       if (parr(3,3)*Qvec(3) > 0 .and. parr(3,2)*Qvec(3) > 0&
            & .and. parr(3 ,1)*Qvec(3) < 0) then
          es_res_cache(esid_rho23E) = rho
       end if
    end if
  end subroutine es_buildcache


  !======================================================================
  subroutine es_thr_tQ_selctd(pcur, Ecur, Qvec, n, thr, a)
    real(dp), intent(in)  :: pcur(:,:), Ecur, Qvec(:)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: thr, a(3)
    !------------------------------------------------
    real(dp) :: thrtmp, atmp(3)
    integer  :: i, j, maxiter

    thr = zero
    select case (n)
    case (0)
       !-- give a value for safety
       a = (/ zero, zero, Qvec(3)/abs(Qvec(3)) /)
    case (1)
       thr = pcur(4,1)
       a   = pcur(:3,1) / thr
    case default
       maxiter = 2**(n-1) - 1
       if (n > 3) write(0,*) 'VERY INEFFICIENT ALGORTHM FOR tau_zQ&
            & for large n!' 
       !-- algorithm tries all sign combinations of the momenta (except for
       ! overall sign which is adjusted at the end) so as to get all potential
       ! thrust axes, and then calculates the thrust for each one.
       ! 
       thr = zero
       do i = 0, maxiter
          atmp = pcur(:3,1)
          do j = 2, n
             if (.not. btest(i,j-2)) then
                atmp = atmp + pcur(:3,j)
             else
                atmp = atmp - pcur(:3,j)
             end if
          end do
          !-- on 10^4 events, this seems to give identical results to below.
          thrtmp = mod3sq(atmp)
          !a = a / sqrt(mod3sq(a))
          !thrtmp = zero
          !do j = 1, n
          !   thrtmp = thrtmp + abs(dot_product(a,pcur(:3,j)))
          !end do
          if (thrtmp > thr) then
             thr = thrtmp
             a   = atmp
          end if
       end do
       !-- take square roots, normalise, and get axis in right direction
       !-- thr should always be > 0?
       if (thr <= zero) then
          write(0,*) 'WARNING: thr^2 had illegal value:',thr
          thr = zero
          a = (/zero,zero,one/)
       else
          thr = sqrt(thr)
          a = a / thr
       end if
       
       if (dot_product(a,Qvec(:3)) < 0) a = -a
    end select
    thr = thr / (half*abs(Qvec(3)))
  end subroutine es_thr_tQ_selctd
  


    
  !--- start caching ---------
  subroutine es_cache_new()
    event_cache_on = .true.
    event_id     = event_id + 1
  end subroutine es_cache_new

  !--- stop caching & clear ---------
  subroutine es_cache_off()
    event_cache_on = .false.
  end subroutine es_cache_off
  


  !-- used by Mike ----------------------------
  function dot_m(p,i,j) result(dot)
    real(dp), intent(in) :: p(:,:)
    integer,  intent(in) :: i,j
    real(dp) :: dot
    dot = p(4,i)*p(4,j) - p(3,i)*p(3,j) - p(2,i)*p(2,j) - p(1,i)*p(1,j)
  end function dot_m

  !-- used instead of Mike --------------------
  function dot_2v(p,q) result(dot)
    real(dp), intent(in) :: p(:),q(:)
    real(dp) :: dot
    dot = p(4)*q(4) - p(3)*q(3) - p(2)*q(2) - p(1)*q(1)
  end function dot_2v

  !---------------------------------------------------
  ! converts a logical to an integer with value 0 or 1
  function log2int(l) result(i)
    logical, intent(in) :: l
    integer :: i
    if (l) then
       i = 1
    else
       i = 0
    end if
  end function log2int
  
  !----------------------------------------------------
  ! squared modulus of 3 vector
  function mod3sq(p) result(sqr)
    real(dp), intent(in) :: p(:)
    real(dp) :: sqr
    sqr = p(1)**2 + p(2)**2 + p(3)**2
  end function mod3sq
  

end module event_shapes
