// ============================================================================
// 
// --> macros to use cfortran.h
//
// file:              fortran.h
// created:           14.11.1997
// last modification: 16.11.1997
//
// ============================================================================

// load only once...
#ifndef FORTRAN_H
#define FORTRAN_H

// ----------------------------------------------------------------------------

#define USE_NEW_DELETE_DUMMY_DUMMY

#define USE_NEW_DELETE 1
// GPS: try more recent version
//#include "cfortran4.0.h"
#include "cfortran4.3.h"

// ============================================================================
//
// --> macros
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> names for external FORTRAN routines to be used
//
// ----------------------------------------------------------------------------

#ifndef Underscore_FORTRAN_Name
#define Underscore_FORTRAN_Name(name) name
#endif

// ============================================================================

#endif // of FORTRAN_H

// ============================================================================
//
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

#if 0
// ====================== GARBAGE =============================================
// ====================== END OF GARBAGE ======================================
#endif
