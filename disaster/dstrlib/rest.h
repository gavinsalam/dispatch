// ============================================================================
// 
// --> Disaster++ main test routine include file: 
//     class definitions etc.
//
// file:              rest.h
// created:           09.12.1996
// last modification: 12.05.1997
//
// ============================================================================

// load only once...
#ifndef RESTH
#define RESTH

// ============================================================================

#include "switch.h"
#include "proc.h"

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

class Event;
class MatrixElement;
class Process;
class Global;

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

double rescaled(Global *global,Event *ein,
                LimitType ltype,double frac,
                double (*pmel)(Global*,Event*)
               );

double rescalednew(Global *global,Event *ein,
                LimitType ltype,double frac,
                int variant, Process *process, int ic,
                double jacobian
               );

double testmel(Global *global, Event *event);
double testmeli(Global *global, Event *event);

// ============================================================================
//
// --> classes
//
// ============================================================================

// matrix elements for test purposes

class TestMatrixElement : public MatrixElement {

private:

   int mel_index;

protected:

public:

   TestMatrixElement(int);
   ~TestMatrixElement();

   INLINE int GetIn();
   INLINE int GetOut();
   void Evaluate(
             Event *,
             Contribution *con,    
             double factor
          );
};

// ----------------------------------------------------------------------------

class TestProcess_1 : public Process {

private:

protected:

public:

   TestProcess_1();
   virtual ~TestProcess_1();

   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class TestProcess_2 : public Process {
public:
   TestProcess_2();
   virtual ~TestProcess_2();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ----------------------------------------------------------------------------

class TestProcess_3 : public Process {
public:
   TestProcess_3();
   virtual ~TestProcess_3();
   INLINE void Define(Global*);
   INLINE void AssignMatrixElements();
};

// ============================================================================

#endif // of RESTH

// ============================================================================
//
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
