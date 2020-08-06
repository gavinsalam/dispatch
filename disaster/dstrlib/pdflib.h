// ============================================================================
//
// --> Interface of PDFLIB to C++
//
// file:              pdflib.h
// created:           16.04.1997
// last modification: 15.10.1997
//
// ============================================================================

void pdfsetparC(
        int ipd, int ipdset,
        int& nfl, int& norder,
        double& dlambda4,
        double& xmin, double& xmax, double& scale_min, double& scale_max
     );

void pdfgetC(
        double x, double scale, double *pdf
     );

void pdfassetparC(
        int ilo,
        double* dlambda4
     );

double pdfasC(
        double scale
     );

// ============================================================================
//        
// --> End of file.
//      
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

