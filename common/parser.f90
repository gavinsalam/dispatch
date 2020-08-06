module parser
  use types; use consts_dp
  implicit none
  private
  
  integer, parameter :: cid_eof       = -100
  integer, parameter :: cid_data      = -50
  integer, parameter :: cid_newrun    = 1
  integer, parameter :: cid_eslist    = 2
  integer, parameter :: cid_pdflist   = 3
  integer, parameter :: cid_xqlist    = 4
  integer, parameter :: cid_Elim    = 5
  integer, parameter :: ncid = 5
  character(len=20), parameter :: command_strings(ncid) = (/&
       &'new run             ',&
       &'event-shape list    ',&
       &'pdf list            ',&
       &'xQ list             ',&
       &'Elim                ' /)
  

  integer, parameter :: maxlen = 100
  logical               :: pushback_done = .false. 
  character(len=maxlen) :: pushback_line
  integer               :: pushback_iostat = -1

  public :: parse_setup_runs

!!$  interface label_val
!!$     module procedure label_val_gl!, label_val_ln
!!$  end interface
!!$
!!$  interface which_key
!!$     module procedure which_key_1, which_key_2
!!$  end interface
  

contains
  

  !-------------------------------------------------------------
  ! sets up the run on the basis of the contents of the file
  ! specified by the option -rundef
  subroutine parse_setup_runs(run)
    use sub_defs_io; use run_descriptor; use pdf_manager; use pdf_names
    type(run_def), pointer :: run(:)
    !-------------------------------------------
    character(len=maxlen) :: filename, line
    integer :: idev, nruns, ix
    integer  :: nxQ=0, npdf=0, nes=0, irun=0, command_id
    real(dp) :: Elim=zero, muF_Q_default
    integer, parameter :: maxn = 50
    integer            :: pdf_ids(maxn), iostat, ihash
    type(esbinned)     :: eslist(maxn)
    type(xQ_pair)      :: xQlist(maxn)
    character(len=3) :: loglin


    filename = trim(string_val_opt('-rundef'))
    idev     = get_new_device()
    open(unit=idev, file=trim(filename),status='OLD')

    muF_Q_default = dble_val_opt('-muF',one)

    ! determine number of runs
    nruns = 0
    do 
       call get_line(idev, line, command_id)
       if (command_id == cid_eof) exit
       if (command_id == cid_newrun) nruns = nruns + 1
    end do

    allocate(run(nruns))
    rewind(idev)


    !-- now read in the details of each run
    irun = 0
    do 
       call get_line(idev, line, command_id)
       select case(command_id)
       case(cid_eof)
          call setup_single_run(run(irun), pdf_ids(1:npdf), &
               &                xQlist(1:nxQ), &
               &                eslist(1:nes), Elim)
          exit
       case(cid_newrun)
          if (irun /= 0) call setup_single_run(run(irun), &
               &                pdf_ids(1:npdf), &
               &                xQlist(1:nxQ), &
               &                eslist(1:nes), Elim)
          irun = irun + 1
          nxQ  = 0 
          npdf = 0
          nes  = 0
          Elim = zero
       case(cid_pdflist)
          do
             call get_line(idev, line, command_id)
             if (command_id /= cid_data) then
                call pushback(line)
                exit
             end if
             npdf = npdf + 1
             if (npdf > maxn) stop 'Error in parse_setup_runs: &
                  &npdf exceeded maxn'
             !-- allow both numeric and symbolic PDF identifiers
             read(line,*,iostat=iostat) pdf_ids(npdf)
             if (iostat /= 0) then
                !-- this will not work if there are comments...
                ihash = index(line,'#') - 1
                if (ihash < 1) ihash = len(line)
                pdf_ids(npdf) = pn_Name2ID(trim(adjustl(line(1:ihash))))
             end if
          end do
       case(cid_xQlist)
          do
             call get_line(idev, line, command_id)
             if (command_id /= cid_data) then
                call pushback(line)
                exit
             end if
             nxQ = nxQ + 1
             if (nxQ > maxn) stop 'Error in parse_setup_runs: &
                  &nxQ exceeded maxn'
             !-- try reading muF_Q as well -- if it is not available
             !   then go back to reading just x and Q.
             read(line,*,iostat=iostat) xQlist(nxQ)%x,xQlist(nxQ)%Q,&
                  &xQlist(nxQ)%muF_Q
             if (iostat /= 0) then
                read(line,*) xQlist(nxQ)%x,xQlist(nxQ)%Q
                xQlist(nxQ)%muF_Q = muF_Q_default
             end if
             
          end do
       case(cid_eslist)
          do
             call get_line(idev, line, command_id)
             if (command_id /= cid_data) then
                call pushback(line)
                exit
             end if
             nes = nes + 1
             if (nes > maxn) stop 'Error in parse_setup_runs: &
                  &nes exceeded maxn'
             line = adjustl(line)
             !-- could make use of name, but save work
             ix = index(line,' ')
             eslist(nes)%id = es_idOfName(line(1:ix))
             read(line(ix+1:),*) eslist(nes)%lo, eslist(nes)%hi, &
                  &              eslist(nes)%nbins, loglin, &
                  &              eslist(nes)%lim
             if (loglin == "log") then
                eslist(nes)%log = .true.
             else if (loglin == "lin") then
                eslist(nes)%log = .false.
             else
                write(0,*) 'Error in parse_setup_run: &
                     &loglin should be either log or lin' 
                write(0,*) 'instead it is: ', loglin
                stop
             end if
          end do
       case(cid_Elim)
          call get_line(idev, line, command_id)
          read(line,*) Elim
       case default
          write(0,*) 'Error in parse_setup_run: unrecognized command_id'
          stop
       end select
       
    end do
    
    
    close(idev)
  end subroutine parse_setup_runs
  
  
  !-----------------------------------------------------------------
  ! interfaces to the various routines for setting up a run
  subroutine setup_single_run(run, pdf_ids, xQlist, eslist, Elim)
    use run_descriptor; use pdf_manager
    type(run_def),  intent(inout) :: run
    type(xQ_pair),  intent(in)    :: xQlist(:)
    integer,        intent(in)    :: pdf_ids(:)
    type(esbinned), intent(in)    :: eslist(:)
    real(dp),       intent(in)    :: Elim
    !------------------------------------------
    type(pdf_xQ_set), pointer :: pdfs

    !-- each time want to allocate a new piece of memory,
    !   because this will be pointed to, not copied
    nullify(pdfs)
    allocate(pdfs)

    call pm_setup_pdfs(pdfs, pdf_ids(:), xQlist(:))
    call rd_setup_run(run, eslist(:), pdfs, Elim)
  end subroutine setup_single_run
  


  !--------------------------------------------------
  ! gets the next line and prcesses it as follows:
  ! lines starting with # are comments
  ! lines starting with ( are command
  ! anything else could be anything else
  ! end-of-file is considered a "command"
  subroutine get_line(idev, line, command_id)
    integer,          intent(in)   :: idev
    character(len=*), intent(out)  :: line
    integer,          intent(out)  :: command_id
    !-----------------
    integer :: iostat, i

    do 
       !-- figure out whether to do a true read or not
       if (pushback_done) then
          line = pushback_line
          pushback_done = .false.
          iostat = pushback_iostat
       else
          read(idev, '(a)', iostat=iostat) line
          pushback_iostat = iostat
       end if
       if (iostat /= 0) then
          command_id = cid_eof
          return
       end if
       if (line(1:1) == '#') cycle ! skip comments
       if (len(trim(line)) == 0) cycle 
       exit
    end do
    
    if (line(1:1) /= '(') then
       command_id = cid_data
       return
    end if
    
    do i = 1, ncid
       if (index(line, trim(command_strings(i))) >= 1) then
          command_id = i
          return
       end if
    end do

    write(0,*) 'Error in parse:get_line. Unrecognized command: ', trim(line)
    stop
  end subroutine get_line
  

  !------------------------------------------------------
  ! pretend that the line is put back onto the buffer...
  subroutine pushback(line)
    character(len=*), intent(in)  :: line
    pushback_done = .true.
    pushback_line = line
  end subroutine pushback
  

  !--------------------------------------------------------------
  ! returns the id corresponding to a given name, or esid_illegal
  ! if the name is not supported.
  function es_idOfName(name) result(id)
    use event_shapes
    character(len=*) :: name
    integer :: id
    do id = 1, n_esid
       if (trim(name) == trim(es_names(id))) then
          return
       end if
    end do
    write(0,*) 'ERROR in es_idOfName. Unrecognized esname: ',name
    stop
  end function es_idOfName


!!$  !----------------------------------------------------------------
!!$  ! searches for the label until the end of the file 
!!$  ! (using this presumes that the order of labels is well-defined
!!$  function label_val_gl(label) result(res)
!!$    character(len=*), intent(in) :: label
!!$    real(dp)                     :: res
!!$    integer :: iostat, ix
!!$    
!!$    do
!!$       ix = index(gl_line, label)
!!$       if (ix > 1) then
!!$          read(gl_line(ix+len(label):), *) res
!!$          exit
!!$       end if
!!$
!!$       read(gl_idev, '(a)', iostat=iostat) gl_line
!!$       if (iostat /= 0) then
!!$          write(0,*) 'ERROR in label_val: requested label "'//label//'" not&
!!$               & present'
!!$          stop
!!$       end if
!!$    end do
!!$  end function label_val_gl
!!$  
!!$  
!!$
!!$  !----------------------------------------------------------------------
!!$  ! gives the index of the key that is present
!!$  ! whichever key it comes across first is good.
!!$  function which_key_2(key1, key2) result(res)
!!$    character(len=*), intent(in) :: key1, key2
!!$    integer :: res
!!$    integer :: iostat, nl1, nl2
!!$    
!!$    nl1 = len(key1)
!!$    nl2 = len(key2)
!!$    do 
!!$       if (gl_line(1:nl1) == key1) then
!!$          res = 1; exit
!!$       else if (gl_line(1:nl2) == key2) then
!!$          res = 2; exit
!!$       end if
!!$       
!!$       read(gl_idev, '(a)', iostat=iostat) gl_line
!!$       if (iostat /= 0) then
!!$          res = 0; exit
!!$       end if
!!$    end do
!!$  end function which_key_2
!!$  
!!$  function which_key_1(key1) result(res)
!!$    character(len=*), intent(in) :: key1
!!$    integer :: res
!!$    integer :: iostat, nl1
!!$    
!!$    nl1 = len(key1)
!!$    do 
!!$       if (gl_line(1:nl1) == key1) then
!!$          res = 1; exit
!!$       end if
!!$       read(gl_idev, '(a)', iostat=iostat) gl_line
!!$       if (iostat /= 0) then
!!$          res = 0; exit
!!$       end if
!!$    end do
!!$  end function which_key_1

end module parser

