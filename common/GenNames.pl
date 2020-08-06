#!/usr/local/bin/perl -w
#----------------------------------------------------------------------
# Perl script for generating the pdf_name routine, which converts
# between PDF ids and pdf names
#
# $Id: GenNames.pl,v 1.1 2001/10/19 22:10:57 gsalam Exp $
open (DEST,"> pdf_names.f90") || die "Could not open pdf_names.f90 for output";

open (PDF,  "grep 'parameter.*pdf_[A-Z]' pdf_manager.f90 | sed 's/.*: *//' | sed 's/ *=.*//'|") || die "Could not extract list of PDF IDs";

# write this fortran style...
$n=0;
while (<PDF>){
  $n=$n+1;
  $IDS[$n]=$_;
  chop($IDS[$n]);
  #print $IDS[$n];
}

#exit;


print DEST "!--------------------------------------------------------------\n";
print DEST "! File generated automatically by GenNames.pl\n";
print DEST "! Contains routines for converting between pdf names and IDs\n";
print DEST "!\n";
print DEST "module pdf_names\n";
print DEST "  use warnings_and_errors\n";
print DEST "  use pdf_manager\n";
print DEST "  implicit none\n";
print DEST "  private\n";
print DEST "\n";
print DEST "  public :: pn_ID2Name, pn_Name2ID\n";
print DEST "\n";
print DEST "contains\n";
print DEST "\n";
print DEST "  function pn_ID2Name(pdfID) result(name)\n";
print DEST "	character(len=30)   :: name\n";
print DEST "	integer, intent(in) :: pdfID\n";
print DEST "\n";
print DEST "	select case(pdfID)\n";

# now put in some meat!
for($i=1; $i<=$n; $i++) {
  print DEST "	case (".$IDS[$i].")\n";
  print DEST "	  name = '".$IDS[$i]."'\n";
}
print DEST "	case default\n";
print DEST "	  call wae_error('pn_ID2Name','Unrecognized PDF ID')\n";
print DEST "	end select\n";
print DEST "	\n";
print DEST "  end function pn_ID2Name\n";
print DEST "\n";
print DEST "  function pn_Name2ID(name) result(pdfID)\n";
print DEST "	character(len=*), intent(in) :: name\n";
print DEST "	integer                      :: pdfID\n";
print DEST "\n";
print DEST "	select case(trim(name))\n";
# now put in the rest of the meat!
for($i=1; $i<=$n; $i++) {
  print DEST "	case ('".$IDS[$i]."')\n";
  print DEST "	  pdfID = ".$IDS[$i]."\n";
}
print DEST "	case default\n";
print DEST "	  call wae_error('pn_Name2ID','Unrecognized PDF name '//trim(name))\n";
print DEST "	end select\n";
print DEST "	\n";
print DEST "  end function pn_Name2ID\n";
print DEST "  \n";
print DEST "end module pdf_names\n";
