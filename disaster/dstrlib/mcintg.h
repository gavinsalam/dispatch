// ============================================================================
// 
// --> include file for MC integration
//
// file:              mcintg.h
// created:           29.03.1997
// last modification: 31.03.1997
//
// ============================================================================

// load only once...
#ifndef MCINTGH
#define MCINTGH

// ============================================================================

#include "switch.h"
#include "global.h"

// ============================================================================
//
// --> types
//
// ============================================================================

// ... to be used for casts
typedef double (*MCIntegrandPointer)(double*,double,int,int);

// ============================================================================
//
// --> classes
//
// ============================================================================

class MCFunction {

private:

protected:

public:

   MCFunction();
   virtual ~MCFunction();

   virtual double Integrand(double* ,double, int, int) = 0;

   virtual void TestMCFunction() = 0;
};

// ----------------------------------------------------------------------------

// the general framework for MC integration

class MCframe 
   : public MCFunction 
   {

private:

   Global* global;

protected:

   int dimension;
   double* unitCube;

public:

   int Get_dimension() { return dimension; }
   double* Get_unitCube() { return unitCube; }

   MCframe();
   ~MCframe();

   void AssignGlobal(Global*);   
   void Integrate(int iterations,int npoints1,int npoints2,
                  double *value,double *error);
};

// ----------------------------------------------------------------------------

extern "C"{
double fvegas1_c_from_f(double* x, double* wgt, double* calls, int* itmx);
}

// ============================================================================

#endif // of MCINTGH

// ============================================================================
//
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
