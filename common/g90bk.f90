!======================================================================
! g90bk: A very basic histogramming package;
!        basic functionality modelled on gbook.
!        let user do most of the hard work, since that is his/her job
!
! It will help itself to devs 55 and 56 for the output of 2d histograms
! in PAW format
!======================================================================
module g90bk
  use types; use consts_dp; implicit none
  private
!!$  type, public ::  histogram_vw
!!$     real(dp)          :: min, max, spc
!!$     integer           :: n
!!$     real(dp), pointer :: vertices(:),h(:)
!!$     character(len=16) :: name
!!$  end type histogram_vw
!!$  public :: hist_vw_init, hist_vw_adddat!, hist_vw_write
!!$  public :: hist_vw_mean, hist_vw_close
  type, public ::  histogram
     real(dp)          :: min, max, spc
     integer           :: n
     real(dp), pointer :: h(:)
     character(len=32) :: name
  end type histogram
  type, public :: histogram2d
     real(dp)          :: min1, max1, min2, max2, spc1, spc2, tot_wgt
     integer           :: n1, n2
     real(dp), pointer :: h(:,:)
     character(len=32) :: name
  end type histogram2d

  integer, parameter, private :: devd=55, devk = 56
  integer :: glb_dev = -1
  integer, parameter :: glb_maxnn = 999, glb_biglen = 200
  logical :: glb_started = .false.
  character(len=160) :: glb_file 

  public :: hist_init, hist_adddat, hist_addtri, hist_write, hist_copy
  public :: hist_ibin
  public :: hist_close, hist_clear
  public :: hist_tdwrite, hist_tdrewind, hist_tdtext
  public :: hist_tdclose, hist_tdopen, hist_tdopenlast, hist_tdflush


  public :: hist2d_init, hist2d_adddat, hist2d_close
  public :: hist2d_paw

contains


  !======================================================================
  ! All the dp 1-d histogram routines
  subroutine hist_init(hist, name, min, max, n)
    use consts_dp; implicit none
    type(histogram)                :: hist
    character(len=*),   intent(in) :: name
    real(dp),           intent(in) :: min, max
    integer,            intent(in) :: n
    !------------------------------------------------------------
    integer :: istat

    !-- initialise -----------------
    !deallocate(hist%h, stat=istat)  ! remove any old rubbish
    allocate(hist%h(1:n))
    hist%min = min
    hist%max = max
    hist%n   = n
    hist%spc = (max-min)/n
    hist%h   = zero
    hist%name = trim(name)
    !write(0,*) 'N is:',hist%n,hist%spc
  end subroutine hist_init

  !----------------------------------------------------------------------
  ! Make a copy of a histogram...
  subroutine hist_copy(hist_orig,hist_cp)
    type(histogram), intent(in) :: hist_orig
    type(histogram)             :: hist_cp
    integer :: istat

    !-- initialise -----------------
    !deallocate(hist_cp%h, stat=istat)  ! remove any old rubbish
    allocate(hist_cp%h(1:size(hist_orig%h)))
    hist_cp%min  = hist_orig%min 
    hist_cp%max  = hist_orig%max 
    hist_cp%n    = hist_orig%n   
    hist_cp%spc  = hist_orig%spc 
    hist_cp%h    = hist_orig%h   
    hist_cp%name = hist_orig%name 
  end subroutine hist_copy
  

  subroutine hist_adddat(hist,dat,weight)
    use consts_dp; implicit none
    type(histogram), intent(inout)    :: hist
    real(dp),  intent(in) :: dat
    real(dp), intent(in), optional :: weight
    !------------------------------------------------------------
    integer  :: i
    real(dp) :: add_to_bin

    !------------------------------------------------------------
    ! first decide what needs to be added
    if (present(weight)) then
       add_to_bin = weight
    else
       add_to_bin = 1.0_dp
    end if
    !------------------------------------------------------------
    ! now decide whether to add it
    i = ceiling((dat-hist%min)/hist%spc)
    if (i>=1 .and. i<=hist%n) then
       hist%h(i) = hist%h(i) + add_to_bin
    end if
  end subroutine hist_adddat


  !----------------------------------------------------------------------
  ! Returns bin corresponding to the data, or 
  ! -1 if the bin is outside the limits of the data
  ! Not very nice, since it encourages direct access to components
  ! But we do that quite often anyway...
  function hist_ibin(hist, dat) result(ibin)
    use consts_dp; implicit none
    type(histogram), intent(in)    :: hist
    real(dp),        intent(in)    :: dat
    integer                        :: ibin
    if (dat >= hist%min .and. dat <= hist%max) then
       ibin = ceiling((dat-hist%min)/hist%spc)
    else
       ibin = -1
    end if
  end function hist_ibin
  

  !----------------------------------------------------------------------
  ! Add data to bin using a triangular sampling. Triangles have their
  ! maximum at the center of each bin, extending to the center of the bins
  ! on either side.
  subroutine hist_addtri(hist,dat,weight)
    use consts_dp; implicit none
    type(histogram), intent(inout)    :: hist
    real(dp),  intent(in) :: dat
    real(dp),  intent(in) :: weight
    !------------------------------------------------------------
    integer  :: i, j
    real(dp) :: add_to_bin, posn, dm
    
    add_to_bin = weight
    !------------------------------------------------------------
    ! decide whether, where & how to add it
    dm = (dat-hist%min)/hist%spc
    i = ceiling(dm)
    !-- relative position --
    posn = (dm - (i-half))
    if (i>=1 .and. i<=hist%n) then
       if (posn > zero) then
          !-- take care of boundaries so that they are full half rectangle,
          !   rather than a triangle -- will ensure that integrated
          !   distribution comes out correctly --- though differential
          !   distribution will look strange there
          j = min(hist%n,i+1)
       else
          posn = - posn
          j = max(1,i-1)
       end if
       
       hist%h(i) = hist%h(i) + add_to_bin*(one-posn)
       hist%h(j) = hist%h(j) + add_to_bin*(posn)
    end if
  end subroutine hist_addtri

  
  

  !----------------------------------------------------------------------
  ! Explicitly add comment text to the topdraw file.
  ! Intended to allow the user to add comments or other useful things
  ! (e.g. also formatting instructions)
  subroutine hist_tdtext(text,comment_string)
    character(len=*) :: text
    character(len=*), optional :: comment_string
    character(len=30) :: comment
    integer           :: lencom
    if (present(comment_string)) then
       comment = comment_string
       lencom  = len(comment_string)
    else
       comment = '( '
       lencom  = 2
    end if
    
    if (glb_dev < 0) glb_dev = hist_tdopen()
    !-- trim it, to save hassle for the calling routine.
    write(glb_dev, '("'//comment(1:lencom)//'",a)') trim(text)
  end subroutine hist_tdtext
  

  !----------------------------------------------------------------------
  subroutine hist_tdwrite(hist,new,log,join, histerr)
    use sub_defs_io
    type(histogram), intent(in) :: hist
    logical,         intent(in) :: new, log, join
    type(histogram), intent(in), optional :: histerr
    integer :: dev
    integer :: i
    real(dp) :: ymin, x
    real(dp15) :: xmin, xmax
    character(len=100) :: line
    if (glb_dev < 0) glb_dev = hist_tdopen()
    dev = glb_dev
    if (new) then
       write(dev,*) 'NEW FRAME'
       write(dev,*) 'SET WINDOW X 4.9 9.5 Y 2 9'
       write(dev,*) 'SET FONT DUPLEX'
       write(dev,*) "TITLE TOP '"//trim(hist%name)//"'"
       !-- was supposed to be a workaround on the pgf90 compiler
       xmin = hist%min
       xmax = hist%max
       line = "SET LIMITS X "// trim(num2char(xmin))//" "&
            & //trim(num2char(xmax))
       write(dev,*) trim(line)
       !-- below was the original form which caused problems:
       !   recursive write is possibly a valid cause for complaint
       !   But real(hist%min,dp15) was also a cause of problems,
       !   giving non-sensical answers...
       !write(dev,*) "SET LIMITS X "//&
       !     & trim(num2char(real(hist%min,dp15)))//" "&
       !     & //trim(num2char(real(hist%max,dp15)))
       !- CARRY ON MAKING QP SAFE2
       write(dev,*) 'SET PATT 0.02 0.08'
    end if
    if (log) then
       write(0,*) 'LOG PLOTS NOT YET SUPPORTED'
    end if

    if (present(histerr)) then
       write(dev,*) 'SET ORDER X Y DY'
    else
       write(dev,*) 'SET ORDER X Y'
    end if
    
    do i = 1, hist%n
       x = hist%min+(i-0.5_dp)*hist%spc
       if (present(histerr)) then
          write(dev,'(3es25.15)') x, hist%h(i), histerr%h(i)
       else
          write(dev,'(3es25.15)') x, hist%h(i)
       end if
    end do
    
    if (join) then
       write(dev,*) 'JOIN DASH ; PLOT'
    else
       write(dev,*) 'HIST ; PLOT'
    end if
    
  end subroutine hist_tdwrite


  !----------------------------------------------------------------------
  ! allows the name to be set
  subroutine hist_tdsetname(basename,nn,filename)
    use sub_defs_io
    character(len=*), intent(in)  :: basename
    integer,          intent(in)  :: nn
    character(len=*), intent(out) :: filename
    character(len=30) :: numb
    if (nn == 0) then
       numb = ''
    else if (nn <= 99) then
       numb = trim(num2char(nn,'(i2.2)'))
    else if (nn <= glb_maxnn) then
       numb = trim(num2char(nn,'(i3.3)'))
    else
       write(0,*) 'ERROR: hist_tdsetname should not be called with nn >',nn
    end if
    filename = trim(basename)//trim(numb)//'.top'
  end subroutine hist_tdsetname
  

  !----------------------------------------------------------------------
  ! Opens (for read only) the last td file in the sequence and return its
  ! device number
  function hist_tdopenlast() result(dev)
    use sub_defs_io
    integer :: dev
    integer :: i, iostat
    character(len=glb_biglen) :: name, lastname
    character(len=glb_biglen) :: tdbase
    logical :: exist

    tdbase=string_val_opt('-tdfl','gtopdraw')
    name = ''
    do i = 0, glb_maxnn
       lastname = name
       call hist_tdsetname(tdbase, i, name)
       inquire (FILE = trim(name), EXIST = exist, NUMBER=dev)
       !-- opened file is current one: want previous one?
       if ((.not. exist) .or. dev >= 0) exit
    end do
    
    name = lastname
    if (trim(name) /= '') then
       dev = get_new_device()
       open(dev,file=trim(name),status='OLD',iostat=iostat)
       if (iostat == 0) then
          write(0,*) 'Opened '//trim(name)//' as file from previous run'
       else
          write(0,*) 'Failed to open '&
               & //trim(name)//' as file from previous run'
          dev = -1
       end if
    else
       dev = -1
    end if
  end function hist_tdopenlast
  


  !----------------------------------------------------------------------
  ! Returns device number of next available td file (first time around:
  ! subsequently returns device number of the same file)
  function hist_tdopen() result(dev)
    use sub_defs_io
    integer :: dev
    integer :: i, iostat
    character(len=glb_biglen) :: name
!!$    character(len=2)   :: nn
    character(len=glb_biglen) :: tdbase

    !-- signal for nothing to be done (hope user used hist_tdclose?)
    if (glb_dev > 0) then
       dev = glb_dev
       return
    end if
    
    if (.not. glb_started) then
       tdbase=string_val_opt('-tdfl','gtopdraw')
       dev = get_new_device()
       do i = 0, glb_maxnn
          call hist_tdsetname(tdbase, i, name)
!!$          if (i==0) then
!!$             nn = ''
!!$          else
!!$             nn = trim(num2char(i,'(i2.2)'))
!!$          end if
!!$          
!!$          name=trim(tdbase)//trim(nn)//'.top'
          
          open(dev,file=trim(name),status='NEW',iostat=iostat)
          if (iostat == 0) exit
       end do
       glb_started = .true.
       glb_file = trim(name)
    else
       dev = get_new_device()
       open(dev,file=trim(glb_file),iostat=iostat)
    end if
    

    if (iostat /= 0) then
       dev = 6
       write(0,*) 'was not able to open topdrawer file'
       write(0,*) 'will send info to stdout'
    else
       write(0,*) 'OUTPUT WILL BE SENT TO FILE ',trim(glb_file)
    end if
  end function hist_tdopen
  

  !----------------------------------------------------------------------
  ! allow user to start from scratch
  subroutine hist_tdrewind
    if (glb_dev >= 0) then
       rewind(glb_dev)
       write(0,*) 'REWOUND FILE:',trim(glb_file)
    end if
  end subroutine hist_tdrewind
  
  !----------------------------------------------------------------------
  subroutine hist_tdclose
    if (glb_dev >= 0) then
       close(glb_dev)
       glb_dev = -1
    end if
  end subroutine hist_tdclose

  !----------------------------------------------------------------------
  subroutine hist_tdflush
    use sub_defs_io
    if (glb_dev >= 0) then
       call lcl_flush(glb_dev)
    end if
  end subroutine hist_tdflush


  !----------------------------------------------------------------------
  subroutine hist_close(hist)
    use consts_dp; implicit none
    type(histogram)  :: hist
    deallocate(hist%h)
  end subroutine hist_close


  !------------------------------------------------------------
  ! output the histogram to the requested device, optionally to 
  ! a writing routine specified by the user (which is useful e.g.
  ! if one wants to change the measure, normalisation, etc...)
  subroutine hist_write(hist,dev,user_write)
    use consts_dp; implicit none
    type(histogram), intent(in) :: hist
    integer, intent(in) :: dev
    interface 
       subroutine user_write(dev,x,h)
         use types; implicit none
         integer,  intent(in) :: dev
         real(dp), intent(in) :: x, h
       end subroutine user_write
    end interface
    optional :: user_write
    !------------------------------------------------------------
    integer  :: i
    real(dp) :: x, h

    do i = 1, hist%n
       x = hist%min+(i-0.5_dp)*hist%spc
       h = hist%h(i)
       if (present(user_write)) then
          call user_write(dev,x,h)
       else
          write(dev,*) x, hist%h(i)
       end if
    end do
  end subroutine hist_write




  subroutine hist_clear(hist)
    use consts_dp; implicit none
    type(histogram)  :: hist
    hist%h = zero
  end subroutine hist_clear






  !======================================================================
  ! All the dp 2-d histogram routines
  subroutine hist2d_init(hist, min1, max1, n1, min2, max2, n2)
    use consts_dp; implicit none
    type(histogram2d), pointer    :: hist
    real(dp),             intent(in) :: min1,min2,max1,max2
    integer,              intent(in) :: n1, n2

    !-- initialise -----------------
    allocate(hist)
    allocate(hist%h(1:n1, 1:n2))
    hist%min1 = min1
    hist%min2 = min2
    hist%max1 = max1
    hist%max2 = max2
    hist%n1   = n1
    hist%n2   = n2
    hist%spc1 = (max1-min1)/n1
    hist%spc2 = (max2-min2)/n2
    hist%h = 0.0_dp
    hist%tot_wgt = 0.0_dp
  end subroutine hist2d_init

  subroutine hist2d_adddat(hist,dat1,dat2,weight,tot_wgt)
    use consts_dp; implicit none
    type(histogram2d), intent(inout)    :: hist
    real(dp),  intent(in) :: dat1,dat2
    real(dp), intent(in), optional :: weight,tot_wgt
    !------------------------------------------------------------
    integer  :: i1, i2
    real(dp) :: add_to_tot_wgt, add_to_bin

    !------------------------------------------------------------
    ! first decide what needs to be added
    if (present(weight)) then
       add_to_bin = weight
    else
       add_to_bin = 1.0_dp
    end if
    if (present(tot_wgt)) then
       add_to_tot_wgt = tot_wgt
    else
       add_to_tot_wgt = add_to_bin
    end if
    hist%tot_wgt = hist%tot_wgt  + add_to_tot_wgt

    !------------------------------------------------------------
    ! now decide whether to add it
    i1 = ceiling((dat1-hist%min1)/hist%spc1)
    i2 = ceiling((dat2-hist%min2)/hist%spc2)
    if (i1>=1 .and. i1<=hist%n1 .and. i2>=1 .and. i2<=hist%n2) then
       hist%h(i1,i2) = hist%h(i1,i2) + add_to_bin
    end if
  end subroutine hist2d_adddat

  subroutine hist2d_close(hist)
    use consts_dp; implicit none
    type(histogram2d), pointer :: hist
    deallocate(hist%h)
    deallocate(hist)
  end subroutine hist2d_close

  subroutine hist2d_paw(hist,filebase,macroname,nm1,nm2,nm3)
    use consts_dp; use sub_defs_io; implicit none
    type(histogram2d), intent(in)   :: hist
    character(*),         intent(in)   :: filebase
    character(*), optional, intent(in) :: macroname,nm1,nm2,nm3
    !------------------------------------------------------------
    integer        :: i1, i2, histn = 123
    character(300) :: line
    character(20)  :: macnm, vecnm, histnm


    !-- produce a data file ------------------------------
    open(devd,file=trim(filebase)//'.dat')
    do i2 = 1, hist%n2
       do i1 = 1, hist%n1
          write(devd,*) hist%h(i1,i2)/hist%tot_wgt
       end do
    end do
    close(devd)
    !-- produce kumac file ------------------------------
    if (present(macroname)) then
       macnm = trim(macroname)
       if (trim(macnm) == 'none') return
    else
       macnm = trim(filebase)
    end if
    vecnm = trim(macnm)//'v'; histnm = trim(num2char(histn))

    open(devk,file=trim(filebase)//'.kumac')
    write(devk,11) 'MACRO '//trim(macnm)
    write(devk,11) ' ve/cre '//trim(vecnm)//&
         & '('//trim(num2char(hist%n1))//','//&
         &      trim(num2char(hist%n2))//')'
    write(devk,11) ' ve/re '//trim(vecnm)//&
         & " '"//trim(filebase)//".dat'"
    write(devk,11) ' cre/2dhisto '//trim(histnm)// ' title '//&
         & trim(num2char(hist%n1))//' '//&
         & trim(num2char(real(hist%min1)))//' '//&
         & trim(num2char(real(hist%max1)))//' '//&
         & trim(num2char(hist%n2))//' '//&
         & trim(num2char(real(hist%min2)))//' '//&
         & trim(num2char(real(hist%max2)))
    write(devk,11) ' put/cont '//trim(histnm)//' '//&
         & trim(vecnm)//'(1:'//trim(num2char(hist%n1))//&
         & ',1:'//trim(num2char(hist%n2))//')'
!!$    write(devk,11) ' 2d_plot/lego '//trim(histnm)//' 30 30 1'
!!$    write(devk,11) ' histo/plot '//trim(histnm)//' colz'
    write(devk,11) ' Set NCOl 28;Pal 1;Conto '//trim(histnm)//' 20 3'
    if (present(nm1) .and. present(nm2) .and. present(nm3)) then
       write(devk,11) ' atitle "'//trim(nm1)//'" "'//trim(nm2)//&
            & '" '//trim(nm3)
    end if
    write(devk,11) 'RETURN'
    close(devk)
11  format(a)
  end subroutine hist2d_paw
end module g90bk





!O   !======================================================================
!O   ! All the variable-width 1-d histogram routines
!O   subroutine hist_vw_init(hist, vertices)
!O     use consts_dp; implicit none
!O     type(histogram_vw), pointer    :: hist
!O     real(dp),           intent(in) :: vertices(0:)
!O     !----------------------------------------------------------------------
!O     integer :: n
!O 
!O     n = ubound(vertices,dim=1)
!O 
!O     !-- initialise -----------------
!O     allocate(hist)
!O     allocate(hist%h(1:n),hist%vertices(0:n))
!O     hist%vertices = vertices
!O     hist%min = hist%vertices(0)
!O     hist%max = hist%vertices(n)
!O     hist%n   = n
!O     hist%h = 0.0_dp
!O     hist%tot_wgt = 0.0_dp
!O     hist%mean    = 0.0_dp
!O   end subroutine hist_vw_init
!O 
!O   subroutine hist_vw_adddat(hist,dat,weight,tot_wgt)
!O     use consts_dp; implicit none
!O     type(histogram_vw), intent(inout)    :: hist
!O     real(dp),  intent(in) :: dat
!O     real(dp), intent(in), optional :: weight,tot_wgt
!O     !------------------------------------------------------------
!O     integer  :: i, ih, il
!O     real(dp) :: add_to_tot_wgt, add_to_bin
!O 
!O     !------------------------------------------------------------
!O     ! first decide what needs to be added
!O     if (present(weight)) then
!O        add_to_bin = weight
!O     else
!O        add_to_bin = 1.0_dp
!O     end if
!O     if (present(tot_wgt)) then
!O        add_to_tot_wgt = tot_wgt
!O     else
!O        add_to_tot_wgt = add_to_bin
!O     end if
!O     hist%tot_wgt = hist%tot_wgt  + add_to_tot_wgt
!O     hist%mean = hist%mean + add_to_bin*dat
!O 
!O     !------------------------------------------------------------
!O     ! now decide whether to add it
!O     if (dat < hist%min .or. dat > hist%max) return
!O     ih = hist%n; il = 0
!O     do 
!O        !-- ensure that difference between ih and il always gets reduced,
!O        ! but never goes below 1
!O        i = (ih - il) / 2 + il
!O        if (dat < hist%vertices(i)) then
!O           ih = i
!O        else
!O           il = i
!O        end if
!O        if (ih -il == 1) exit
!O     end do
!O     hist%h(ih) = hist%h(ih) + add_to_bin
!O   end subroutine hist_vw_adddat
!O 
!O   !------------------------------------------------------------
!O   ! NO NEED FOR OUTPUT ROUTINE -- LET THE USER DO THE WORK, SINCE WHAT
!O   ! S/HE WANTS MAY BE QUITE DIFFERENT FROM WHAT I AM EXPECTING
!O   ! output the histogram to the requested device, optionally to 
!O   ! a writing routine specified by the user (which is useful e.g.
!O   ! if one wants to change the measure, normalisation, etc...)
!O !O?   subroutine hist_vw_write(hist,dev,user_write)
!O !O?     use consts_dp; implicit none
!O !O?     type(histogram_vw), intent(in) :: hist
!O !O?     integer, intent(in) :: dev
!O !O?     interface 
!O !O?        subroutine user_write(dev,x,h)
!O !O?          use types; implicit none
!O !O?          integer,  intent(in) :: dev
!O !O?          real(dp), intent(in) :: x, h
!O !O?        end subroutine user_write
!O !O?     end interface
!O !O?     optional :: user_write
!O !O?     !------------------------------------------------------------
!O !O?     integer  :: i
!O !O?     real(dp) :: x, h
!O !O? 
!O !O?     do i = 1, hist%n
!O !O?        x = hist%min+(i-0.5_dp)*hist%spc
!O !O?        h = hist%h(i) / hist%tot_wgt / hist%spc
!O !O?        if (present(user_write)) then
!O !O?           call user_write(dev,x,h)
!O !O?        else
!O !O?           write(dev,*) x, hist%h(i) / hist%tot_wgt / hist%spc
!O !O?        end if
!O !O?     end do
!O !O?   end subroutine hist_vw_write
!O 
!O 
!O   function hist_vw_mean(hist)
!O     use consts_dp; implicit none
!O     real(dp)                        :: hist_vw_mean
!O     type(histogram_vw), intent(in) :: hist
!O     integer :: i
!O     hist_vw_mean = hist%mean / hist%tot_wgt
!O   end function hist_vw_mean
!O 
!O 
!O 
!O   subroutine hist_vw_close(hist)
!O     use consts_dp; implicit none
!O     type(histogram_vw), pointer :: hist
!O     deallocate(hist%h,hist%vertices)
!O     deallocate(hist)
!O   end subroutine hist_vw_close



!!$program test
!!$  use types; use consts_dp; use g90bk
!!$  implicit none
!!$  
!!$  integer :: i
!!$  type(histogram) :: h1, h2
!!$  call hist_init(h1,'thrust LL',zero,one,20)
!!$  call hist_init(h2,'thrust NLL',zero,one,20)
!!$  h1%h = one
!!$  h2%h = 0.1_dp
!!$  call hist_tdwrite(h1,new=.true.,log=.false.,join=.false.)
!!$  call hist_tdwrite(h2,new=.false.,log=.false.,join=.true.)
!!$end program test
