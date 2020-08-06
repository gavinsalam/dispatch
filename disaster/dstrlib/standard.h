// ============================================================================
// 
// --> standard definitions for C / C++
//
// file:              standard.h
// created:           14.11.1997
// last modification: 11.12.1997
//
// ============================================================================

// load only once...
#ifndef STANDARD_H
#define STANDARD_H

// ----------------------------------------------------------------------------

// debugging macros

#if 0
#endif

// ----------------------------------------------------------------------------

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

// ----------------------------------------------------------------------------

enum Boolean {False = FALSE, True = TRUE};

// ----------------------------------------------------------------------------

#ifndef NULL
#define NULL 0
#endif

// ----------------------------------------------------------------------------

#ifdef __GNUC__
#define templateX  template
#endif
#ifdef __DECCXX
#define templateX  template
#endif
#ifndef templateX
#define templateX
#endif

// ----------------------------------------------------------------------------

#ifndef instantiate_Template_1
#define instantiate_Template_1(C, T1) \
           templateX class C < T1 >;
#endif

// ----------------------------------------------------------------------------

#ifndef instantiate_Template_2
#define instantiate_Template_2(C, T1, T2) \
           templateX class C < T1, T2 >;
#endif

// ----------------------------------------------------------------------------

#ifndef instantiate_Template_3
#define instantiate_Template_3(C, T1, T2, T3) \
           templateX class C < T1, T2, T3 >;
#endif

// ============================================================================

#endif // of STANDARD_H

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
