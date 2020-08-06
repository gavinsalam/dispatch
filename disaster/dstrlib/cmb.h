// ============================================================================
// 
// --> include file for combinatorical functions
//
// file:              cmb.h
// created:           29.03.1997
// last modification: 04.09.1997
//
// ============================================================================

// load only once...
#ifndef CMBH
#define CMBH

// ============================================================================

#include "switch.h"

// ============================================================================
//
// --> classes
//
// ============================================================================

class Permute {

private:
protected:   
public:

   // public for simplicity
   int number,*length;
   int ***list;
   int **startAt;
   double *dataS;
   int nS;
 
//   Permute();
   Permute(int);
   ~Permute();
   INLINE void Define(int);
   INLINE void CreateList(int);
   void PrintList(FILE*);
   INLINE void CreateListOfLists(int nmin,int nmax);
   INLINE int **GetList(int);
   INLINE int *GetStartAt(int);
   void PrintPartial(FILE *);
   INLINE double SumOverPermutations(double *data, int n, int ifixed);
};

// ----------------------------------------------------------------------------

class MapSquareLinear {

private:

public:

   int number, in;
   int ***mapij;
   int **imap;
   int **jmap;
   int *ijmaplength;

   MapSquareLinear(int,int);
   ~MapSquareLinear();
   INLINE void Define(int,int);
   INLINE void Delete();
   void Print(FILE *);
   INLINE void CreateList(int);
   INLINE void CreateListOfLists(int nmin,int nmax);
   int **GetSquareToLinearSafe(int);
   inline int *GetSquareToLinearProjection1Unsafe(int in_L) 
                  { return imap[in_L]; }
   inline int *GetSquareToLinearProjection2Unsafe(int in_L)
                  { return jmap[in_L]; }
   inline int  GetSquareToLinearLengthUnsafe(int in_L)
                  { return ijmaplength[in_L]; }
};

// ============================================================================

#endif // of CMBH

// ============================================================================
//
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
