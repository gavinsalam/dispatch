// ============================================================================
//
// --> include file for mathematical routines
//
// file:              mth.h
// created:           29.03.1997
// last modification: 28.11.1997
//
// ============================================================================

// load only once...
#ifndef MTHH
#define MTHH

// ============================================================================

#include "switch.h"
#include "fold.h"

// ============================================================================
//
// global variables
//
// ============================================================================

// mathematical constants

extern double Pi;
extern double PiSquared;
extern double TwoPi;
extern double FourPi;
    
// ============================================================================
//
// --> class-unrelated functions
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> set mathematical constants
//
// ----------------------------------------------------------------------------

INLINE void SetMathConst();

// ----------------------------------------------------------------------------
//
// --> spherical geometry
//
// ----------------------------------------------------------------------------

double ThreePolarVWCos(double v, double w, double cosphi);
double ThreePolar(double v, double w, double phi);
double GetCosPhi(double e, double v, double w);
double GetPhi(double e, double v, double w);
double AbsSinFromCos(double cosPhi);
double atan2ZeroTwoPi(double, double);

long lfactorial(int n);
double dfactorial(int n);
int min(int a,int b);
double dmin(double a,double b);
double dmin3(double a,double b,double c);
double dmin4(double a,double b,double c,double d);
double dmin5(double a,double b,double c,double d,double e);
int max(int a,int b);
double dmax(double a,double b);
double dmax3(double a,double b,double c);
double dmax4(double a,double b,double c,double d);
double dmax5(double a,double b,double c,double d,double e);
int sign(double x);
int lsign(long x);
void correct_range01(double *x);
long rclv(double x);
long mapDoubleToLong(double x, long m);
void gmae(double, double&, int&);
void mfffp(double in, char *str);
int inRange(int in, int imin, int imax);
int inRange(double in, double imin, double imax);
int compareAbsRel(double a, double b, double epsAbs, double epsRel);
int compareDouble(
       double a,
       double b,
       double eps,
       int method
    );
double ReturnZeroZero();
double naiveRandomNumber(); 

double variableTransformation(
          double  x,
          int     index,
          double  alpha,
          double  beta,
          double  lambda,
          double  mu,   
          double* jacobian
       );

double variableTransformationInverse(
          double  y,
          int     index,
          double  alpha,
          double  beta,
          double  lambda,
          double  mu
       );

double variableTransformationLower(
          double  x,
          int     index,
          double  alpha,
          double* jacobian
       );

// ============================================================================ 

#endif // of MTHH

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
