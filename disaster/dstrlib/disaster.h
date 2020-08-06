// ============================================================================
// 
// --> DISASTER++ include file: 
//     class definitions etc.
//
// file:              disaster.h
// created:           31.03.1997
// last modification: 02.12.1997
//
// ============================================================================

// load only once...
#ifndef DISASTER_H
#define DISASTER_H

// ============================================================================

#include "switch.h"

#include "qcd.h"
#include "mcintg.h"

// ============================================================================
//
// --> declarations
//
// ============================================================================

//class QCD;

// ============================================================================
//
// --> classes
//
// ============================================================================

class Disaster 
   : public TheoryEvaluation,
     public MCframe 
   {

private:

   // --> debugging parameters

   int print_debug;
   double *x_variable;
   int dimension_debug;
   int i_dimension_debug;
   double wgt_debug;
   int icalls_debug;
   int iterations_debug;

   // --> counters to keep track of things

   int number_of_integrations_done;
   int event_counter;

   int fpError_counter;
   Counter fpError_print_counter;
   Counter fpError_log_counter;

   // --> variables to guide the parsing of an input file

   int jump_to;
   int execute_commands;

   // --> the return value of disaster
   int disaster_return_value;

   // --> parameters for the technical cut

   int scaleChoiceForTechnicalCut;
   double technicalCutParameter;
   Boolean technicalCut_flag;
   int technicalCut_counter;
   Counter technicalCut_print_counter;
   Counter technicalCut_log_counter;

   // --> parameters for MC integration

   int mc_iterations;     // number of VEGAS iterations for grid preparation
   int mc_points;         // number of MC points (reference)
   int mc_points_prep;    // number of MC points for grid preparation 
   int mc_points_final;   // number of MC points for final integration
   double mc_fact_prep;   // fraction...
   double mc_fact_final;  //     "
   double mc_eps;         // cuts of boundaries of the box

   // --> parameters which control the output
   
   int print_progress;    // do not print progress information 

   // --> explicit choice of the reference frame

   int explicitFrameChoice; 

   // --> phase space parameters

   int lepton_integration; // 0: xB, y fixed, 1: integrate

   double xB_fixed;
   double y_fixed;

   int xi_integration;     // 0: xi fixed, 1: integrate
   double xi_fixed;   

   int ps_permute_max;     // maximum number of permuted partons
   int *ps_permute_n;      // number of indices required for PS permutation
   int ps_permute_c;       // counter
   int **ps_permute;       // array containing the PS permutation

   // kinematical ranges/cuts
   double ECMFULL;         // total CM energy
   double Q2min, Q2max, W2min, W2max, xBmin, xBmax, ymin, ymax;  // lepton cuts

   // actual kinematical ranges
   double Q2min_actual, Q2max_actual, W2min_actual, W2max_actual;
   double xBmin_actual, xBmax_actual;

   // CM energy squared
   double SH;

   // --> parameters for variable transformations
   
   double e_gamma;        // first energy
   double e_gamma_upper;  // first energy, upper limit of the unit interval
   double e_lambda;       // first energy, separation of upper and lower region
   double e_mu;           // first energy, value at the separation point

   double e_gamma1;       // the other energies

   double v_delta;        // first angle
   double v_delta_upper;
   double v_delta_lambda;
   double v_delta_mu;
   
   double v_delta1;       // the other angles
   double v_delta1_upper;
   double v_delta1_lambda;
   double v_delta1_mu;

   int    xB_y_parametrization;
   double xB_y_trho;
   double xB_y_tsigma;

   double u_alpha;
   double xi_alpha;

   double shifted_parton_variables;

   // --> for test purposes:

   int limitTestFlag;
   int limittypetest;
   int limitQuantity;
   double fracstart;
   double fracmin;
   double fracstep;

   int min_points;        // minimum number of points: skipped
   int max_points;        // maximum number of points, then: exit

   ElectroWeak* electroWeak;
   
protected:

public:

   inline int Get_print_debug() { return print_debug; }

   inline int Get_event_counter()        { return event_counter; }
   inline int Get_fpError_counter()      { return fpError_counter; }
   inline int Get_technicalCut_counter() { return technicalCut_counter; }

   inline int Get_disaster_return_value() { return disaster_return_value; }

   inline void Set_technicalCut_flag() { technicalCut_flag = True; }

   inline double Get_xBmin_actual() { return xBmin_actual; }
   inline double Get_xBmax_actual() { return xBmax_actual; }
   inline double Get_ymin()         { return ymin; }
   inline double Get_ymax()         { return ymax; }
   inline double Get_Q2min_actual() { return Q2min_actual; }
   inline double Get_Q2max_actual() { return Q2max_actual; }
   inline double Get_W2min_actual() { return W2min_actual; }
   inline double Get_W2max_actual() { return W2max_actual; }
   inline double Get_xB_fixed()     { return xB_fixed; }
   inline double Get_y_fixed()      { return y_fixed; }
   inline int    Get_lepton_integration() { return lepton_integration; }
   inline int    Get_xi_integration()     { return xi_integration; }
   inline double Get_xi_fixed()           { return xi_fixed; }
   inline int*   Get_ps_permute_n() { return ps_permute_n; }
   inline int**  Get_ps_permute()   { return ps_permute; }

   inline double Get_ECMFULL() { return ECMFULL; }

   inline double Get_e_gamma()       { return e_gamma; }
   inline double Get_e_gamma_upper() { return e_gamma_upper; }
   inline double Get_e_lambda()      { return e_lambda; }
   inline double Get_e_mu()          { return e_mu; }
   inline double Get_e_gamma1()      { return e_gamma1; }

   inline double Get_v_delta()        { return v_delta; }
   inline double Get_v_delta_upper()  { return v_delta_upper; }
   inline double Get_v_delta_lambda() { return v_delta_lambda; }
   inline double Get_v_delta_mu()     { return v_delta_mu; }

   inline double Get_v_delta1()        { return v_delta1; }
   inline double Get_v_delta1_upper()  { return v_delta1_upper; }
   inline double Get_v_delta1_lambda() { return v_delta1_lambda; }
   inline double Get_v_delta1_mu()     { return v_delta1_mu; }
    
   inline int    Get_xB_y_parametrization() { return xB_y_parametrization; }
   inline double Get_xB_y_trho()            { return xB_y_trho; }
   inline double Get_xB_y_tsigma()          { return xB_y_tsigma; }

   inline double Get_u_alpha()  { return u_alpha; }
   inline double Get_xi_alpha() { return xi_alpha; }

   inline double Get_shifted_parton_variables() 
                    { return shifted_parton_variables; }

   QCD* qcd;
   
   QCD*         Get_qcd()         { return qcd; }
   ElectroWeak* Get_electroWeak() { return electroWeak; }

   int Get_process_index() { return process_index; }   
   int Get_numberOfFinalStatePartonsInBornTerm() 
          { return numberOfFinalStatePartonsInBornTerm; }   

   Disaster();
   ~Disaster();   

   inline void SetScaleChoiceForTechnicalCut(int in) 
                  { scaleChoiceForTechnicalCut = in; }
   inline int  GetScaleChoiceForTechnicalCut() 
                  { return scaleChoiceForTechnicalCut; }

   inline void   SetTechnicalCutParameter(double in) 
                    { technicalCutParameter = in; }
   inline double GetTechnicalCutParameter() 
                    { return technicalCutParameter; }

   void Define();
   void Delete();

   void SetProcess(Process*);
   void SetUser(MCUser*);
   void RunMC();
   void RunOneEvent();
   int PPCommand(char*, char*);
   int PPCommandInt(char*, int);
   int PPCommandDouble(char*, double);
   void ParseFile(FILE*);
   void Start();
   void Terminate();
   double ScalingCheck(double*, double);
   double Integrand(double*, double, int, int);
   void PrepareForMC();
   void CleanUpAfterMC();

   void TestMCFunction();
};

// ============================================================================

#endif // of DISASTER_H

// ============================================================================
//
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
