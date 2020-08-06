!--------------------------------------------------------------
! File generated automatically by GenNames.pl
! Contains routines for converting between pdf names and IDs
!
module pdf_names
  use warnings_and_errors
  use pdf_manager
  implicit none
  private

  public :: pn_ID2Name, pn_Name2ID

contains

  function pn_ID2Name(pdfID) result(name)
	character(len=30)   :: name
	integer, intent(in) :: pdfID

	select case(pdfID)
	case (pdf_ZERO)
	  name = 'pdf_ZERO'
	case (pdf_CTEQ5)
	  name = 'pdf_CTEQ5'
	case (pdf_CTEQ5M)
	  name = 'pdf_CTEQ5M'
	case (pdf_CTEQ5D)
	  name = 'pdf_CTEQ5D'
	case (pdf_CTEQ5L)
	  name = 'pdf_CTEQ5L'
	case (pdf_CTEQ5HJ)
	  name = 'pdf_CTEQ5HJ'
	case (pdf_CTEQ5M1)
	  name = 'pdf_CTEQ5M1'
	case (pdf_GRV)
	  name = 'pdf_GRV'
	case (pdf_GRV98LO)
	  name = 'pdf_GRV98LO'
	case (pdf_GRV98NLM)
	  name = 'pdf_GRV98NLM'
	case (pdf_GRV98NLD)
	  name = 'pdf_GRV98NLD'
	case (pdf_MRS99)
	  name = 'pdf_MRS99'
	case (pdf_MRS99_1)
	  name = 'pdf_MRS99_1'
	case (pdf_MRS99_2)
	  name = 'pdf_MRS99_2'
	case (pdf_MRS99_3)
	  name = 'pdf_MRS99_3'
	case (pdf_MRS99_4)
	  name = 'pdf_MRS99_4'
	case (pdf_MRS99_5)
	  name = 'pdf_MRS99_5'
	case (pdf_MRS99_6)
	  name = 'pdf_MRS99_6'
	case (pdf_MRS99_7)
	  name = 'pdf_MRS99_7'
	case (pdf_MRS99_8)
	  name = 'pdf_MRS99_8'
	case (pdf_MRS99_9)
	  name = 'pdf_MRS99_9'
	case (pdf_MRS99_10)
	  name = 'pdf_MRS99_10'
	case (pdf_MRS99_11)
	  name = 'pdf_MRS99_11'
	case (pdf_MRS99_12)
	  name = 'pdf_MRS99_12'
	case (pdf_MRS99DIS)
	  name = 'pdf_MRS99DIS'
	case (pdf_MRS99DIS_1)
	  name = 'pdf_MRS99DIS_1'
	case (pdf_MRS99DIS_2)
	  name = 'pdf_MRS99DIS_2'
	case (pdf_MRS99DIS_3)
	  name = 'pdf_MRS99DIS_3'
	case (pdf_MRS99DIS_4)
	  name = 'pdf_MRS99DIS_4'
	case (pdf_MRS99DIS_5)
	  name = 'pdf_MRS99DIS_5'
	case (pdf_MRS99DIS_6)
	  name = 'pdf_MRS99DIS_6'
	case (pdf_MRS99DIS_7)
	  name = 'pdf_MRS99DIS_7'
	case (pdf_MRS99DIS_8)
	  name = 'pdf_MRS99DIS_8'
	case (pdf_MRS99DIS_9)
	  name = 'pdf_MRS99DIS_9'
	case (pdf_MRS99DIS_10)
	  name = 'pdf_MRS99DIS_10'
	case (pdf_MRS99DIS_11)
	  name = 'pdf_MRS99DIS_11'
	case (pdf_MRS99DIS_12)
	  name = 'pdf_MRS99DIS_12'
	case (pdf_MRS98)
	  name = 'pdf_MRS98'
	case (pdf_MRS98_1)
	  name = 'pdf_MRS98_1'
	case (pdf_MRS98_2)
	  name = 'pdf_MRS98_2'
	case (pdf_MRS98_3)
	  name = 'pdf_MRS98_3'
	case (pdf_MRS98_4)
	  name = 'pdf_MRS98_4'
	case (pdf_MRS98_5)
	  name = 'pdf_MRS98_5'
	case (pdf_MRS98LO)
	  name = 'pdf_MRS98LO'
	case (pdf_MRS98LO_1)
	  name = 'pdf_MRS98LO_1'
	case (pdf_MRS98LO_2)
	  name = 'pdf_MRS98LO_2'
	case (pdf_MRS98LO_3)
	  name = 'pdf_MRS98LO_3'
	case (pdf_MRS98LO_4)
	  name = 'pdf_MRS98LO_4'
	case (pdf_MRS98LO_5)
	  name = 'pdf_MRS98LO_5'
	case (pdf_MRST2001)
	  name = 'pdf_MRST2001'
	case (pdf_MRST2001_1)
	  name = 'pdf_MRST2001_1'
	case (pdf_MRST2001_2)
	  name = 'pdf_MRST2001_2'
	case (pdf_MRST2001_3)
	  name = 'pdf_MRST2001_3'
	case (pdf_MRST2001_4)
	  name = 'pdf_MRST2001_4'
	case (pdf_POWER)
	  name = 'pdf_POWER'
	case (pdf_POWER10g)
	  name = 'pdf_POWER10g'
	case (pdf_POWER10q)
	  name = 'pdf_POWER10q'
	case (pdf_POWER04g)
	  name = 'pdf_POWER04g'
	case (pdf_POWER04q)
	  name = 'pdf_POWER04q'
	case default
	  call wae_error('pn_ID2Name','Unrecognized PDF ID')
	end select
	
  end function pn_ID2Name

  function pn_Name2ID(name) result(pdfID)
	character(len=*), intent(in) :: name
	integer                      :: pdfID

	select case(trim(name))
	case ('pdf_ZERO')
	  pdfID = pdf_ZERO
	case ('pdf_CTEQ5')
	  pdfID = pdf_CTEQ5
	case ('pdf_CTEQ5M')
	  pdfID = pdf_CTEQ5M
	case ('pdf_CTEQ5D')
	  pdfID = pdf_CTEQ5D
	case ('pdf_CTEQ5L')
	  pdfID = pdf_CTEQ5L
	case ('pdf_CTEQ5HJ')
	  pdfID = pdf_CTEQ5HJ
	case ('pdf_CTEQ5M1')
	  pdfID = pdf_CTEQ5M1
	case ('pdf_GRV')
	  pdfID = pdf_GRV
	case ('pdf_GRV98LO')
	  pdfID = pdf_GRV98LO
	case ('pdf_GRV98NLM')
	  pdfID = pdf_GRV98NLM
	case ('pdf_GRV98NLD')
	  pdfID = pdf_GRV98NLD
	case ('pdf_MRS99')
	  pdfID = pdf_MRS99
	case ('pdf_MRS99_1')
	  pdfID = pdf_MRS99_1
	case ('pdf_MRS99_2')
	  pdfID = pdf_MRS99_2
	case ('pdf_MRS99_3')
	  pdfID = pdf_MRS99_3
	case ('pdf_MRS99_4')
	  pdfID = pdf_MRS99_4
	case ('pdf_MRS99_5')
	  pdfID = pdf_MRS99_5
	case ('pdf_MRS99_6')
	  pdfID = pdf_MRS99_6
	case ('pdf_MRS99_7')
	  pdfID = pdf_MRS99_7
	case ('pdf_MRS99_8')
	  pdfID = pdf_MRS99_8
	case ('pdf_MRS99_9')
	  pdfID = pdf_MRS99_9
	case ('pdf_MRS99_10')
	  pdfID = pdf_MRS99_10
	case ('pdf_MRS99_11')
	  pdfID = pdf_MRS99_11
	case ('pdf_MRS99_12')
	  pdfID = pdf_MRS99_12
	case ('pdf_MRS99DIS')
	  pdfID = pdf_MRS99DIS
	case ('pdf_MRS99DIS_1')
	  pdfID = pdf_MRS99DIS_1
	case ('pdf_MRS99DIS_2')
	  pdfID = pdf_MRS99DIS_2
	case ('pdf_MRS99DIS_3')
	  pdfID = pdf_MRS99DIS_3
	case ('pdf_MRS99DIS_4')
	  pdfID = pdf_MRS99DIS_4
	case ('pdf_MRS99DIS_5')
	  pdfID = pdf_MRS99DIS_5
	case ('pdf_MRS99DIS_6')
	  pdfID = pdf_MRS99DIS_6
	case ('pdf_MRS99DIS_7')
	  pdfID = pdf_MRS99DIS_7
	case ('pdf_MRS99DIS_8')
	  pdfID = pdf_MRS99DIS_8
	case ('pdf_MRS99DIS_9')
	  pdfID = pdf_MRS99DIS_9
	case ('pdf_MRS99DIS_10')
	  pdfID = pdf_MRS99DIS_10
	case ('pdf_MRS99DIS_11')
	  pdfID = pdf_MRS99DIS_11
	case ('pdf_MRS99DIS_12')
	  pdfID = pdf_MRS99DIS_12
	case ('pdf_MRS98')
	  pdfID = pdf_MRS98
	case ('pdf_MRS98_1')
	  pdfID = pdf_MRS98_1
	case ('pdf_MRS98_2')
	  pdfID = pdf_MRS98_2
	case ('pdf_MRS98_3')
	  pdfID = pdf_MRS98_3
	case ('pdf_MRS98_4')
	  pdfID = pdf_MRS98_4
	case ('pdf_MRS98_5')
	  pdfID = pdf_MRS98_5
	case ('pdf_MRS98LO')
	  pdfID = pdf_MRS98LO
	case ('pdf_MRS98LO_1')
	  pdfID = pdf_MRS98LO_1
	case ('pdf_MRS98LO_2')
	  pdfID = pdf_MRS98LO_2
	case ('pdf_MRS98LO_3')
	  pdfID = pdf_MRS98LO_3
	case ('pdf_MRS98LO_4')
	  pdfID = pdf_MRS98LO_4
	case ('pdf_MRS98LO_5')
	  pdfID = pdf_MRS98LO_5
	case ('pdf_MRST2001')
	  pdfID = pdf_MRST2001
	case ('pdf_MRST2001_1')
	  pdfID = pdf_MRST2001_1
	case ('pdf_MRST2001_2')
	  pdfID = pdf_MRST2001_2
	case ('pdf_MRST2001_3')
	  pdfID = pdf_MRST2001_3
	case ('pdf_MRST2001_4')
	  pdfID = pdf_MRST2001_4
	case ('pdf_POWER')
	  pdfID = pdf_POWER
	case ('pdf_POWER10g')
	  pdfID = pdf_POWER10g
	case ('pdf_POWER10q')
	  pdfID = pdf_POWER10q
	case ('pdf_POWER04g')
	  pdfID = pdf_POWER04g
	case ('pdf_POWER04q')
	  pdfID = pdf_POWER04q
	case default
	  call wae_error('pn_Name2ID','Unrecognized PDF name '//trim(name))
	end select
	
  end function pn_Name2ID
  
end module pdf_names
