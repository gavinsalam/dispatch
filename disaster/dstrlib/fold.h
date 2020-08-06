// ============================================================================
//
// --> Interface of old FORTRAN routines to C++
//
// file:              fold.h.h
// created:           04.05.1997
// last modification: 17.10.1997
//
// ============================================================================

void evalvirtC(
      double *par,
      int i0loop, int i1loop, int ihelicity,
      double& xq, double& xg
     );

// ----------------------------------------------------------------------------

double RealSpence(double x);

// ----------------------------------------------------------------------------

void ciarrayC(
        double *parray, int npa, double dcscale2,
        int& njetsout, double etain, int& icut, 
        int iflags
     );

// ----------------------------------------------------------------------------

int ktclus_C(
       int imode,
       double* pp,
       int n,
       float ecut,
       float* y
    );
           
int ktreco_C(
       int reco,
       double* pp,
       int nn,
       float ecut,
       float ycut,
       float ymac,
       double* pjet,
       int* jet,
       int njet, 
       int nsub
    );
           
// ----------------------------------------------------------------------------

void user1_C(int iaction);

double user2_C(
          int nr,
          int nl,
          double* fvect,
          int* npartons,
          double* xb,
          double* q2,
          double* xi,
          double* weight,
          int* irps,
          int* ialphas,
          int* ialphaem,
          int* lognf
       );

void user3_C(double value, double error);

void print_f_C(
        const char* string,
        int length,
        int unit
     );

// ----------------------------------------------------------------------------

void pfpx_f_C();

// ----------------------------------------------------------------------------

void strFORT(const char* str, int length, int i);

void strFORT_x(const char* str, int length, int i);

// ============================================================================
//        
// --> End of file.
//      
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

