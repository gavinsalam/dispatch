// ============================================================================
// 
// --> Disaster++ matrix elements
//     class definitions
//
// file:              me.h
// created:           05.04.1997
// last modification: 15.09.1997
//
// ============================================================================

// load only once...
#ifndef ME_H
#define ME_H

// ============================================================================

#include <stdio.h>

#include "switch.h"

#include "enum.h"
#include "proc.h"
#include "strng.h"

// ============================================================================
//      
// symbols
//
// ============================================================================

// ============================================================================
//      
// enumerations
//
// ============================================================================

// ============================================================================
//
// --> dummy declarations for later reference
//
// ============================================================================

// ============================================================================
//
// --> global variables
//
// ============================================================================

// ============================================================================
//
// --> functions
//
// ============================================================================

// ============================================================================
//
// --> classes
//
// ============================================================================

// DIS matrix elements

class DISMatrixElement : public MatrixElement {

   
private:

   TwoD_Array <double> termArray;
   double **term;
   int nGraphs;

   double bornTerm;

   double *softCoefficient;
   double *softAngle; 
   int *softContributes;
   int nCoefficient;
   int nSoftLength;
   int *index1;  // :MOD: these are here for test purposes only
   int *index2;

   double collinearCoefficient;
   int collinearContributes;
   SplittingFunction splittingFunction; 

   DISSubProcess mel_index;

   double virtq, virtg;

   static long lastVirtual;
   static double virtqSave, virtgSave;

//   double dummy;

protected:

public:

   DISMatrixElement(DISSubProcess);
   ~DISMatrixElement();

   void Define();
   INLINE int GetIn();
   INLINE int GetOut();
   INLINE int GetSpecialCut();
   void Evaluate(
             Event *,
             Contribution *con,
             double factor
        );
   void FillLorentzLqXq();
   void FillLorentzLqXqg();
   void FillLorentzLgXqq();
   void FillLorentzLqXqgg();
   void FillLorentzLgXqqg();
   void FillLorentzLqXqqq();
   double GetBornTerm(enum DISSubProcess);
//   double FillLorentzLqXqgV(Event*);
//   double FillLorentzLgXqqV(Event*);
   void FillLorentzLpXppV(Event*, double& xq, double& xg);
   double GetLorentzSubtraction(int list, LimitType, int i, int j);
   void   GetLorentzAddedSubtraction(int list, LimitType, int i, int j);
   double EvaluateLorentzAddedSubtraction(
             Event *, int list, LimitType, int i, int j);
   double EvaluateLorentzCollinearInitial(
             Event *, SplittingFunction, LimitType, double);
   void EvaluateLorentzAndColour(Event*, double*);
};

// ----------------------------------------------------------------------------

class DISProcessNoSubtraction_Test : public Process {
public:
   DISProcessNoSubtraction_Test();
   virtual ~DISProcessNoSubtraction_Test();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class DISProcessSubtraction_Test : public Process {
public:
   DISProcessSubtraction_Test();
   virtual ~DISProcessSubtraction_Test();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class DISProcessAddedSubtraction_Test : public Process {
public:
   DISProcessAddedSubtraction_Test();
   virtual ~DISProcessAddedSubtraction_Test();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class DISProcessCollinearInitial_Test : public Process {
public:
   DISProcessCollinearInitial_Test();
   virtual ~DISProcessCollinearInitial_Test();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class DISProcessLO : public Process {
public:
   DISProcessLO();
   virtual ~DISProcessLO();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class DISProcessNLO : public Process {
public:
   DISProcessNLO();
   virtual ~DISProcessNLO();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class DISProcessNLO_Light : public Process {
public:
   DISProcessNLO_Light();
   virtual ~DISProcessNLO_Light();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class DISProcessLO_Test : public Process {
public:
   DISProcessLO_Test();
   virtual ~DISProcessLO_Test();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class DISProcessNLO_Test : public Process {
public:
   DISProcessNLO_Test();
   virtual ~DISProcessNLO_Test();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ============================================================================

#endif // of ME_H

// ============================================================================
//
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
