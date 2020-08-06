// ============================================================================
//
// DISASTER++: Deeply Inelastic Scattering:
//             All Subtractions Through Evaluated Residues
//
// ============================================================================
//
// --> Disaster++ classes
//
// file:              disaster.cc
// created:           31.03.1997
// last modification: 18.12.1997
//
// ============================================================================

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "fortran.h"

#include "switch.h"
#include "global.h"
#include "mth.h"
#include "strng.h"
#include "user.h"
#include "mcintg.h"
#include "cmb.h"
#include "rest.h"
#include "me.h"
#include "qcd.h"
#include "fold.h"
#include "container.h"

#include "disaster.h"

// ============================================================================
//
// definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

// ============================================================================
//
// --> classes
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> The Disaster class
//
// ----------------------------------------------------------------------------

Disaster :: Disaster() {

   MCframe :: AssignGlobal(this);

   fpnan = NULL;
   number_of_integrations_done = 0;

   qcd         = this;
   electroWeak = this;

   Global :: disaster = this;
}

// ----------------------------------------------------------------------------

Disaster :: ~Disaster() {

   Delete();
}

// ----------------------------------------------------------------------------

void Disaster :: Define() {

// ----------------------------------------------------------------------------

   // print the header message

   To_status() << "\n\n\n";
   To_status() <<
      "  +-------------------------------------------------------------+\n";
   To_status() <<
      "  |                                                             |\n";
   To_status() <<
      "  | DISASTER++ 1.0.1                                            |\n";
   To_status() <<
      "  |                                                             |\n";
   To_status() <<
      "  | Deeply Inelastic Scattering:                                |\n";
   To_status() <<
      "  | All Subtractions Through Evaluated Residues                 |\n";
   To_status() <<
      "  |                                                             |\n";
   To_status() <<
      "  | Dirk Graudenz                                               |\n";
   To_status() <<
      "  | Paul Scherrer Institute                                     |\n";
   To_status() <<
      "  | Department of Nuclear and Particle Physics                  |\n";
   To_status() <<
      "  | -- Theory Group --                                          |\n";
   To_status() <<
      "  | CH-5232 Villigen PSI                                        |\n";
   To_status() <<
      "  | Switzerland                                                 |\n";
   To_status() <<
      "  |                                                             |\n";
   To_status() <<
      "  | Electronic Mail: Dirk.Graudenz@psi.ch                       |\n";
   To_status() <<
      "  | WWW URL:         http://www.hep.psi.ch/graudenz/index.html  |\n";
   To_status() <<
      "  | Phone:           +41-56-310-3661                            |\n";
   To_status() <<
      "  | FAX:             +41-56-310-3294                            |\n";
   To_status() <<
      "  |                                                             |\n";
   To_status() <<
      "  +-------------------------------------------------------------+\n";
   To_status() << "\n\n\n";

// ----------------------------------------------------------------------------

   TheoryEvaluation :: Define();
   mc_function = this;

   // permutations of up to 5 objects
   pPermute -> Define(5);
   pPermute -> CreateListOfLists(1, 5);
//   pPermute->PrintList(stdout);

   // floating point errors: get an error handler
   fpnan = new FPNaN(this, if_fpHandler);

#if 0

   // test the floating point handler
   gl -> To_status() << "test the floating point handler\n";
   gl -> To_status() << "\nstatus 1:\n";
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 2:\n";
   double a = 1.;
   double b = 0.;
   double c = a / b;
   gl -> To_status() << FS("number is %16.6e\n", c);
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 2a:\n";
   provoke_floating_point_exception();
   fpnan -> PrintStatus(gl -> To_status());

//   fpnan -> fpsh -> Reinstall_new();
//   fpnan -> fpsh -> Reinstall_old();

   gl -> To_status() << "\nstatus 3:\n";
   provoke_floating_point_exception();
   double d = 5.0;
   gl -> To_status() << FS("number is %16.6e\n", d);
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 4\n";
//   fpnan -> fpsh -> Reinstall_old(); 
   fpnan -> Hold(); 
   fpnan -> Continue(); 
   fpnan -> Hold(); 
   pfpx_f_C();
//   double e = a / b;
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 5\n";
   fpnan -> Continue(); 
//   fpnan -> fpsh -> Reinstall_new();
   provoke_floating_point_exception();
   double f = 2.0;
   gl -> To_status() << FS("number is %16.6e\n", f);
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 6\n"; 
   double x1 = sqrt(-1.0);
   gl -> To_status() << FS("number is %16.6e\n", x1);
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 6a\n";
   errno = 1;
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 7\n";
   fpnan -> Set_flag_mth(20, "twenty");
   fpnan -> Set_flag_mth(21, "twentyone");
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 8\n"; 
   provoke_floating_point_exception();
   double g = 3.0;
   gl -> To_status() << FS("number is %16.6e\n", g);
   fpnan -> LogDoubleStatus(g);
   fpnan -> LogDoubleStatus(g);
   fpnan -> LogDoubleStatus(g);
   gl -> To_status() << FS("number is %16.6e\n", x1);
   fpnan -> LogDoubleStatus(x1);   
   fpnan -> LogDoubleStatus(x1);
   fpnan -> LogDoubleStatus(x1);
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 9\n"; 
   fpnan -> Hold();
   fpnan -> Set_flag_mth(22, "twentytwo");
//   provoke_floating_point_exception();
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 10\n"; 
   fpnan -> Continue();
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 11\n"; 
   provoke_floating_point_exception();
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "\nstatus 99\n";
   fpnan -> PrintStatus(gl -> To_status());

   fpnan -> Reset();
   gl -> To_status() << "\nstatus 100\n";
   fpnan -> PrintStatus(gl -> To_status());

   gl -> To_status() << "End of test.\n";

   exit(-1);

#endif

   // initialize mathematical constants such as Pi, etc.
   SetMathConst();

   jump_to = 0;
   execute_commands = TRUE;

// ----------------------------------------------------------------------------

   // default parameters

   print_debug = 0;
   x_variable = NULL;
   dimension_debug = -1;

   mc_iterations   = 10;
   mc_points       = 1000;
   mc_points_prep  = 1000;
   mc_points_final = 1000;
   mc_fact_prep    = 1.;
   mc_fact_final   = 1.;
   mc_eps          = 1.e-12;

   print_progress = 0;

   lepton_integration = 1;
   xB_fixed           = 0.5;
   y_fixed            = 0.5;
   xi_integration     = 1;
   xi_fixed           = 1.0;
 
   ECMFULL            = 100.;
   Q2min              = 0.;
   Q2max              = 1.e15;
   W2min              = 0.;
   W2max              = 1.e15;
   xBmin              = 0.;
   xBmax              = 1.;
   ymin               = 0.;
   ymax               = 1.;
   
   register int i;

   ps_permute_c   = -1;
   ps_permute_max = 4;
   ps_permute_n   = new int [ps_permute_max+1];
   ps_permute     = new int* [ps_permute_max+1];
   for (i = 1; i <= ps_permute_max; ++i) {
       ps_permute_n[i] = 0;
       ps_permute[i]   = new int [i];
   }

#if 0
   for (i = 1; i <= ps_permute_max; ++i) {
       delete [] ps_permute[i];
   }  
   delete [] ps_permute;
   delete [] ps_permute_n;

   ps_permute_n   = new int [ps_permute_max+1];
   ps_permute     = new int* [ps_permute_max+1];
   for (i = 1; i <= ps_permute_max; ++i) {
       ps_permute_n[i] = 0;
       ps_permute[i]   = new int [i];
   }
#endif   

   disaster_return_value = 0; 

   SetScaleChoiceForTechnicalCut(0); // :_CAUTION_: these two should be 
   SetTechnicalCutParameter(1.e-8);  // :_CAUTION_: process dependent

   explicitFrameChoice = -1;

   // variable transformations
   e_gamma        = 0.5; // 0
   e_gamma_upper  = 0.5;
   e_lambda       = 1.0;
   e_mu           = 1.0;

   e_gamma1       = 0.5; // 2

   // perhaps related to imom=0 branching (== 1->2?)
   v_delta        = 0.5; //  1
   v_delta_upper  = 0.5; // 13                
   v_delta_lambda = 1.0; // 14
   v_delta_mu     = 1.0; // 15
   // perhaps related to imom=1 branching (== 2->3?) other way around
   v_delta1        = 0.5; // 3
   v_delta1_upper  = 0.5;
   v_delta1_lambda = 1.0;
   v_delta1_mu     = 1.0;

   xB_y_parametrization = 1; // 0: logarithmic, 1: power law
   xB_y_trho   = 0.01;        
   xB_y_tsigma = 2.0;

   // mapParameter( 9)
   u_alpha  = 1.;
   // larger xi_alpha gives xi closer to xB? 
   // (normally log distribution). mapParameter(16)
   xi_alpha = 1.;

   shifted_parton_variables = 0;

   // for test purposes

   limitTestFlag = FALSE;
   limittypetest = SOFT;
   limitQuantity = 0;
   fracstart     = 1.0;
   fracmin       = 1.e-6;
   fracstep      = 0.1;

   min_points =  0;
   max_points = -1;
}

// ----------------------------------------------------------------------------

void Disaster :: Delete() {

   if (fpnan != NULL) {
      delete fpnan;
      fpnan = NULL;
   }   

   delete x_variable;

   register int i;
   for (i = 1; i <= ps_permute_max; ++i) {
       delete [] ps_permute[i];
   }
   delete [] ps_permute;
   delete [] ps_permute_n;

   TheoryEvaluation :: Delete();
}

// ----------------------------------------------------------------------------

// set the process under consideration

void Disaster :: SetProcess(Process* processI) {

   process = processI;
}

// ----------------------------------------------------------------------------

// set the pointer to the user routines

void Disaster :: SetUser(MCUser* userIn) {

   user = userIn;
}

// ----------------------------------------------------------------------------

// run the Monte Carlo with predefined parameters

void Disaster :: RunMC() {

   PrepareForMC();

   double value, error;
   Integrate(mc_iterations, mc_points_prep, mc_points_final, &value, &error);

   ++number_of_integrations_done;

   To_status() << "\n\n";
   To_status() << 
      "===========================================================\n";
   To_status() << FS("Integration no. %3d\n", number_of_integrations_done);
   To_status() << FS("value=%e\n", value);
   To_status() << FS("error=%e\n", error);
   if (value != 0.0)
      To_status() << FS("error in per cent=%e\n", 100. * error / value);
   else
      To_status() << "value is 0.0, no relative error!\n";
   To_status() << 
      "===========================================================\n";
   To_status() << "\n\n";

   CleanUpAfterMC();
}

// ----------------------------------------------------------------------------

// --> evaluate the integrand for one specific set of parameters

void Disaster :: RunOneEvent() {

   PrepareForMC();

   if (x_variable == NULL)
      errf(-1, "Disaster :: RunOneEvent: error (1)");
   if (i_dimension_debug != dimension_debug)
      errf(-1, "Disaster :: RunOneEvent: error (2)");
   if (dimension_debug != dimension)
      errf(-1, "Disaster :: RunOneEvent: error (3)");

   user -> StartIntegration1();
   user -> StartIntegration2();

   double value = Integrand(
                     x_variable, 
                     wgt_debug,
                     icalls_debug,
                     iterations_debug 
                  );
   
   user -> EndIntegration(value, 0.e0);

   To_status() << "\n\n";
   To_status() << 
      "===========================================================\n";
   To_status() << FS("value=%e\n", value);
   To_status() << 
      "===========================================================\n";
   To_status() << "\n\n";

   CleanUpAfterMC();
}

// ----------------------------------------------------------------------------
//
// --> read and execute parameters:
//
//     input format is: ``pname'' is a string with the parameter name,
//                      ``parameter'' is a string with the value.
// 
//     returns: 0 if ok, 1 for termination, -1 for error
//
// ----------------------------------------------------------------------------

int Disaster :: PPCommand (char* pname, char* parameter) {

   int ret = 0;
   
   long trap;

// ----------------------------------------------------------------------------

   // navigation in the parameter file

//   To_status() << FS("pname=:%s: ", pname) 
//               << FS("parameter=:%s:\n", parameter);
//   if (!execute_commands)
//      To_Status() << "(-) ";

   if (pname != NULL && pname[0] == '!') {
      // do not use this parameter
      return 0;
   }

   if (!strcmp(pname, "LABEL")) {
      // this is a label
      int label = ReadInt(this, pname, parameter);
      if (jump_to != 0 && label == jump_to) {
         jump_to = 0;
         execute_commands = TRUE;
      }
      return 0;
   }

   if (!execute_commands)
      {} // do nothing
   else

   if (!strcmp(pname, "JUMP")) {
      // this is a jump location
      jump_to = ReadInt(this, pname, parameter);
      if (jump_to != 0)
         execute_commands = FALSE;
   }
   else

   if (!strcmp(pname, "END")) {
      if (execute_commands) {
         To_status() << "END parameter read.\n";
         To_log() << "END __dummy__\n";
         To_log() << "# --------------------------------------"
                  << "----------------------------------------\n";
         ret = 1;
      }
   }
   else

   if (!strcmp(pname, "ECHO")) {
      String echo(100);
      ReadString(this, pname, parameter, &echo);
      echo.BackslashToNewline();
      To_status() << FS("%s\n", echo.data);
      To_log()    << "ECHO __dummy__\n";
   }
   else

   if (!strcmp(pname,"ECHO_")) {
      String echo(100);
      ReadString(this, pname, parameter, &echo);
      echo.BackslashToNewline();
//    print to TestOut...
      To_log() << "ECHO_ __dummy__\n";
   }
   else

   if (!strcmp(pname,"ECHO_BOTH")) {
      String echo(100);
      ReadString(this, pname, parameter, &echo);
      echo.BackslashToNewline();
      To_status() << FS("%s\n", echo.data);
      To_log() << "ECHO_BOTH __dummy__\n";
   }
   else

// ----------------------------------------------------------------------------

   // --> miscellaneous parameters

   if (!strcmp(pname,"DISASTER_RETURN_VALUE")) {
      disaster_return_value = ReadInt(this,pname,parameter);
   }
   else

// ----------------------------------------------------------------------------

   // Monte Carlo parameters

   if (!strcmp(pname,"ITERATIONS")) {
      mc_iterations=ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"POINTS")) {
      mc_points=ReadInt(this,pname,parameter);
      mc_points_prep = (int) (mc_points * mc_fact_prep);
      mc_points_final = (int) (mc_points * mc_fact_final);
   }
   else

   if (!strcmp(pname,"POINTS_PREP")) {
      mc_points_prep=ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"POINTS_FINAL")) {
      mc_points_final=ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"FACT_PREP")) {
      mc_fact_prep=ReadDouble(this,pname,parameter);
      mc_points_prep = (int) (mc_points * mc_fact_prep);
   }
   else

   if (!strcmp(pname,"FACT_FINAL")) {
      mc_fact_final=ReadDouble(this,pname,parameter);
      mc_points_final = (int) (mc_points * mc_fact_final);
   }
   else

// ----------------------------------------------------------------------------

   // --> technical cuts

   if (!strcmp(pname,"MC_EPS")) {
      mc_eps
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"SCALE_CHOICE_FOR_TECHNICAL_CUT")) {
      SetScaleChoiceForTechnicalCut(
         ReadInt(this,pname,parameter)
      );
   }
   else

   if (!strcmp(pname,"TECHNICAL_CUT_PARAMETER")) {
      SetTechnicalCutParameter(
         ReadDouble(this,pname,parameter)
      );
   }
   else

   if (!strcmp(pname,"EXPLICIT_FRAME_CHOICE")) {
      explicitFrameChoice=ReadInt(this,pname,parameter);
   }
   else

   // -------------------------------------------------------------------------

   if (!strcmp(pname,"PRINT_PROGRESS")) {
      print_progress = ReadInt(this, pname, parameter);
   }
   else

   // -------------------------------------------------------------------------

   // variable transformations
   
   if (!strcmp(pname,"E_GAMMA")) {
      e_gamma
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"E_GAMMA_UPPER")) {
      e_gamma_upper
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"E_LAMBDA")) {
      e_lambda
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"E_MU")) {
      e_mu
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"E_GAMMA1")) {
      e_gamma1
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"V_DELTA")) {
      v_delta
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"V_DELTA_UPPER")) {
      v_delta_upper
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"V_DELTA_LAMBDA")) {
      v_delta_lambda
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"V_DELTA_MU")) {
      v_delta_mu
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"V_DELTA1")) {
      v_delta1
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"V_DELTA1_UPPER")) {
      v_delta1_upper
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"V_DELTA1_LAMBDA")) {
      v_delta1_lambda
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"V_DELTA1_MU")) {
      v_delta1_mu
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"XB_Y_PARAMETRIZATION")) {
      xB_y_parametrization = ReadInt(this, pname, parameter);
   }   
   else

   if (!strcmp(pname,"XB_Y_TRHO")) {
      xB_y_trho
         = ReadDouble(this, pname, parameter);
   }
   else

   if (!strcmp(pname,"XB_Y_TSIGMA")) {
      xB_y_tsigma
         = ReadDouble(this, pname, parameter);
   }
   else

   if (!strcmp(pname,"U_ALPHA")) {
      u_alpha
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"XI_ALPHA")) {
      xi_alpha
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"SHIFTED_PARTON_VARIABLES")) {
      shifted_parton_variables
         = ReadInt(this, pname, parameter);
   }
   else

   // -------------------------------------------------------------------------

   if (!strcmp(pname,"MIN_POINTS")) {
      min_points
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"MAX_POINTS")) {
      max_points
         = ReadInt(this,pname,parameter);
   }
   else

// ----------------------------------------------------------------------------

   // Lepton phase space, CM energy and lepton cuts

   if (!strcmp(pname,"LEPTON_INTEGRATION")) {
      lepton_integration
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"XB_FIXED")) {
      xB_fixed
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"Y_FIXED")) {
      y_fixed
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"XI_INTEGRATION")) {
      xi_integration
         = ReadBool(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"XI_FIXED")) {
      xi_fixed
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"ECM")) {
      ECMFULL
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"QMIN")) {
      Q2min
         = pow(ReadDouble(this,pname,parameter),2);
   }
   else

   if (!strcmp(pname,"QMAX")) {
      Q2max
         = pow(ReadDouble(this,pname,parameter),2);
   }
   else

   if (!strcmp(pname,"WMIN")) {
      W2min
         = pow(ReadDouble(this,pname,parameter),2);
   }
   else

   if (!strcmp(pname,"WMAX")) {
      W2max
         = pow(ReadDouble(this,pname,parameter),2);
   }
   else

   if (!strcmp(pname,"XBMIN")) {
      xBmin
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"XBMAX")) {
      xBmax
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"YMIN")) {
      ymin
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"YMAX")) {
      ymax
         = ReadDouble(this,pname,parameter);
   }
   else

// ----------------------------------------------------------------------------

   // process selection, permutation and subtraction types

   if (!strcmp(pname,"PROCESS_INDEX")) {
      process_index
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"NUMBER_OF_FINAL_STATE_PARTONS_IN_BORN_TERM")) {
      numberOfFinalStatePartonsInBornTerm
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"SUBTRACTION_TYPE")) {
      subtraction_type
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"PS_PERMUTE_N")) {
      ps_permute_c=ReadInt(this,pname,parameter);
      if (ps_permute_c > ps_permute_max)
         errf(-1, " Disaster::PPCommand: too many particles to permute");
      if (ps_permute_c <= 0)
         errf(-1, " Disaster::PPCommand: 0 particles to permute");
      ps_permute_n[ps_permute_c]=0;
   } 
   else

   if (!strcmp(pname,"PS_PERMUTE")) {
      if (ps_permute_c < 0)
         errf(-1, " Disaster::PPCommand: ps_permute_c < 0");
      if (ps_permute_n[ps_permute_c] >= ps_permute_c)
         errf(-1,"Permutation too long.");
      ps_permute[ps_permute_c][ps_permute_n[ps_permute_c]++]
         = ReadInt(this,pname,parameter);
      int co;
      if (ps_permute_n[ps_permute_c]==ps_permute_c) {
         To_status() << "Permutation:\n";
         for (co=0;co<ps_permute_n[ps_permute_c];++co)
             To_status() << FS("%2d -> ", co)
                         << FS("%2d\n", ps_permute[ps_permute_c][co]);
      }
   }
   else

// ----------------------------------------------------------------------------

   // test parameters

   if (!strcmp(pname,"FRACSTART")) {
      fracstart
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"FRACMIN")) {
      fracmin
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"FRACSTEP")) {
      fracstep
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"ICATCH")) {
      icatch
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"JCATCH")) {
      jcatch
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"TCATCH")) {
      tcatch
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"IA")) {
      ia
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"I0")) {
      i0
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"I1")) {
      i1
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"I2")) {
      i2
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"I3")) {
      i3
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"I2P")) {
      i2p
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"I3P")) {
      i3p
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"IB")) {
      ib
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"ITEST1")) {
      itest1
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"ITEST2")) {
      itest2
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"ITEST3")) {
      itest3
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"MEL_INDEX")) {
      mel_index
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"RHO_LIM")) {
      rho_lim
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"ALPHA_LIM")) {
      alpha_lim
         = ReadDouble(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"LIMIT_TEST_FLAG")) {
      limitTestFlag
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"LIMIT_TYPE_TEST")) {
      limittypetest
         = ReadInt(this,pname,parameter);
   }
   else

   if (!strcmp(pname,"LIMIT_QUANTITY")) {
      limitQuantity
         = ReadInt(this,pname,parameter);
   }
   else

// ----------------------------------------------------------------------------

   // debugging parameters

   if (!strcmp(pname, "FP_HANDLER")) {
      if_fpHandler
         = (Boolean) ReadBool(this, pname, parameter);
   }
   else

   if (!strcmp(pname, "PRINT_DEBUG")) {
      print_debug
         = ReadInt(this, pname, parameter);
   }
   else

   if (!strcmp(pname, "DIMENSION_DEBUG")) {
      dimension_debug
         = ReadInt(this, pname, parameter);
      i_dimension_debug = 0;
      delete x_variable;
      x_variable = new double[dimension_debug];
   }
   else

   if (!strcmp(pname, "X_VARIABLE")) {
      if (i_dimension_debug >= dimension_debug) {
         To_error() << FS("%d ", dimension_debug)
                    << FS("%d\n", i_dimension_debug);
         errf(-1, "Disaster :: PPCommand: X_VARIABLE");
      }
      x_variable[i_dimension_debug++]
         = ReadDouble(this, pname, parameter);
   }
   else

   if (!strcmp(pname, "WGT_DEBUG")) {
      wgt_debug
         = ReadDouble(this, pname, parameter);
   }
   else

   if (!strcmp(pname, "ICALLS_DEBUG")) {
      icalls_debug
         = ReadInt(this, pname, parameter);
   }
   else

   if (!strcmp(pname, "ITERATIONS_DEBUG")) {
      iterations_debug
         = ReadInt(this, pname, parameter);
   }
   else

   if (!strcmp(pname, "TRAP_SET")) {
      trap_set
         = ReadInt(this, pname, parameter);
      is_soft_trap = -1;
      trap_check_count = -1;
   }
   else

   if (!strcmp(pname, "IS_SOFT_TRAP")) {
      is_soft_trap
         = ReadInt(this, pname, parameter);
   }
   else

   if (!strcmp(pname, "TRAP_CHECK_COUNT")) {
      trap_check_count
         = ReadInt(this, pname, parameter);
   }
   else

   if (!strcmp(pname, "TRAP")) {
      trap
         = ReadInt(this, pname, parameter);
      trap = trap;  // :_TAWM_:
#if CHECKDESTRUCT
         gl -> tlist[trap_set] -> SetTrap(trap, is_soft_trap, trap_check_count);
         gl -> tlist[trap_set] -> Print(stdout);
#endif
   }
   else

// ----------------------------------------------------------------------------

   // function triggers

   if (!strcmp(pname, "RUN_MC")) {
      To_status() << "\n\n\nRunning MC:\n\n";
      To_log() << "#\n# RUN_MC __dummy__\n#\n";
      To_log().Flush();
      RunMC();
      To_log().Flush();
   }
   else

   if (!strcmp(pname, "RUN_ONE_EVENT")) {
      To_status() << "\n\n\nRunning one event:\n\n";
      To_log() << "# RUN_ONE_EVENT __dummy__\n";
      To_log().Flush();
      RunOneEvent();
      To_log().Flush();
   }
   else

   if (!strcmp(pname, "RUN_VEGAS")) {
      To_status() << "\n\n\nRunning VEGAS:\n\n";
      To_log() << "# RUN_VEGAS __dummy__";
//      RunVEGAS(this);
      errf(-1,"RunVEGAS currently not implemented");
   }
   else

// ----------------------------------------------------------------------------

   // user parameters

   {
      ret = user -> SetParameter(pname, parameter);
   }

   return ret;
}

// ----------------------------------------------------------------------------
//
// --> integer parameter
//
// ----------------------------------------------------------------------------

int Disaster :: PPCommandInt(char* pName, int in) {

   char parameter[1000];
   sprintf(parameter, "%20d", in);
   return PPCommand(pName, parameter);
}

// ----------------------------------------------------------------------------
//
// --> double parameter
//
// ----------------------------------------------------------------------------

int Disaster :: PPCommandDouble(char* pName, double in) {

   char parameter[1000];
   sprintf(parameter, "%30.22e", in);
   return PPCommand(pName, parameter);
}

// ----------------------------------------------------------------------------

// parse a complete input file

void Disaster :: ParseFile(FILE* in) {

   FILE* ParametersIn = in;

#if 0
   // read ``data_path'' and ``job_name''
   data_path -> ReadBoundedFromFile(ParametersIn);
   To_status() << FS("data path = :%s:\n", data_path -> data);
   job_name -> ReadBoundedFromFile(ParametersIn);
   To_status() << FS("job name = :%s:\n", job_name -> data);
   To_status() << "\n\n\n";

   // open the file for the output of plot data
   String PlotOutName(1000);
   PlotOutName.CopyFromString(data_path);
   PlotOutName.AppendFromChar(".plot");
   PlotOut = fopen(PlotOutName.data, "w");

   // open the file for the output of the log file
   String logOutName(1000);
   logOutName.CopyFromString(data_path);
   logOutName.AppendFromChar(".log");
   FILE* logOut = fopen(logOutName.data, "w");
   Mstream* old_logMstreamOut = To_log().Get_mstreamOut();
   To_log().Associate(logOut);   
#endif

   // read parameters from the input file

   String str1(100), str2(100);

   int eof = FALSE;
   int ret = 0;

   jump_to = 0;
   execute_commands = TRUE;

   while (!eof && ret == 0) {
      eof = SkipWhiteSpace(ParametersIn);
      if (!eof) {
         str1.ReadNonWhiteFromFile(ParametersIn);
         eof = SkipWhiteSpace(ParametersIn);
         if (!eof) {
            str2.ReadNonWhiteFromFile(ParametersIn);
            ret = PPCommand (str1.data, str2.data);
            if (ret == -1) {
               To_status() << FS("str1=:%s:\n",str1.data);
               To_status() << FS("str2=:%s:\n",str2.data);
               errf(-1,"DISASTER++: unknown command.");
            }
         } else {
            To_status() << FS("str1=:%s:\n",str1.data);
            errf(-1,"parameter missing.");
         }
      }
   }

#if 0
   // close the plot file
   if (PlotOut != NULL)
      fclose(PlotOut);

   // restore log Mstream to old value and close the log file
   To_log().Associate(old_logMstreamOut);
   fclose(logOut);
#endif
}

// ----------------------------------------------------------------------------

// starts a Disaster run (opens files, signals the user, etc.)

void Disaster :: Start() {

   To_status() << "\n\n\n";
   To_status() <<
      "  +-------------------------------------------------------------+";
   To_status() <<
      "\n... Disaster :: Start.\n\n";

   To_status() << "\n\n\n";

   user -> Start();
}

// ----------------------------------------------------------------------------

// terminates a Disaster run (closes files, etc.)

void Disaster :: Terminate() {

   user -> End();

   To_status() << "\n\n\n";

//   delete fpnan;
//   fpnan = NULL;

   To_status() << 
      "\n... Disaster++ terminated successfully.\n\n";
   To_status() << 
      "  +-------------------------------------------------------------+";
   To_status() << "\n\n\n";
}

// ----------------------------------------------------------------------------

// --> scaling check

double Disaster :: ScalingCheck(double* y, double trueweight) {

      double frac;

      double ysave[10];
      int a;
      for (a = 0; a < 10; ++a) {
          ysave[a] = y[a];
          To_status() << FS("y[%d] = ", a)
                      << FS("%16.6e\n", y[a]);
      }

      int offset = 0;

      if (lepton_integration == 1)
         offset += 2;

      if (process -> GetFrame() == pCMS && xi_integration)
         ++offset;

      // offset now points to the first parton phase space variable
      // for (2+1) and (3+1) parton Born term-like events
      // or to the u-variable for collinear terms

      To_status() << FS("offset = %d\n", offset);

      for (frac = fracstart; frac > fracmin; frac *= fracstep) {

          double compensate = pow(frac, rho_lim);     
          double log1 = - log(frac);
          double log2 =   pow(log(frac), 2);

          for (a = 0; a < 10; ++a)                                
              y[a] = ysave[a];

          switch (limitQuantity) {

             case 0:

                // limits of 3+1 subtracted, leading singularities

                switch (limittypetest) {
                   case SOFT:
                      y[offset] *= frac;
                      break;
                   case COLLINEAR:
                      y[offset + 1] *= frac;
                      break;
                   case SOFT_AND_COLLINEAR:
                      y[offset]     *= frac;
                      y[offset + 1] *= frac;
                      break;
                   default:
                      errf(-1, "(0) not known");
                }

                break;

             case 1:

                // collinear limit of 3+1 subtracted, subleading singularity

                switch (limittypetest) {
                   case COLLINEAR:
                      y[offset + 3] *= frac;
                      break;
                   default:
                      errf(-1, "(1) not known");
                }

                break;

             case 2:

                // collinear limit of 2+1 Born, leading singularity, pCMS

                switch (limittypetest) {
                   case COLLINEAR:
                      y[offset] *= frac;
                      break;
                   default:
                      errf(-1, "(2) not known");
                }

                break;

             case 3:

                // limits of 2+1 Born, leading singularity, hCMS

                switch (limittypetest) {
                   case SOFT:
                      y[offset] *= frac;
                      break;
                   case COLLINEAR:
                      y[offset + 1] *= frac;
                      break;
                   case SOFT_AND_COLLINEAR:
                      y[offset]     *= frac;
                      y[offset + 1] *= frac;
                      break;
                   default:
                      errf(-1, "(3) not known");
                }

                break;

             case 4:

                // limits of added subtraction for 2+1 jets

                switch (limittypetest) {
                   case -10:
                      y[offset] = 1.- frac*(1.-y[offset]);
                      break;
                   case COLLINEAR:
                      y[offset + 1] *= frac;
                      break;
                   default:
                      errf(-1, "(4) not known");
                }

                break;

             case 10:

                // limits for variable -> 0

                y[offset + limittypetest] 
                   = frac * y[offset + limittypetest];

                break;

             case 11:

                // limits for variable -> 1

                y[offset + limittypetest] 
                   = 1. - frac * (1. - y[offset + limittypetest]);

                break;

             case -10:

                // limits for variable -> 0

                y[offset - 1]
                   = frac * y[offset - 1];

                y[offset + limittypetest]
                   = 0.00001 * y[offset + limittypetest];
                
                break;

             case -11:

                // limits for variable -> 1

                y[offset - 1]
                   = 1. - frac * (1. - y[offset - 1]);

                y[offset + limittypetest] 
                   = 1. - 0.00001 * (1. - y[offset + limittypetest]);

                break;

             default:
                errf(-1, "limitQuantity not known");
          }

          ContributionArray* ca;
          ca = caFree -> GetDefinedObject(this, 100, 11);

          // check the list with the events
          if (process -> GetEventList() -> IsNotEmpty()) 
             errf(-1, "Disaster::Integrand: event list is not empty");

          user -> BeginEvent(process);

          process -> CalculateWeight(y, ca);

          double return_weight = user 
                                 -> EndEvent( 
                                       process, 
                                       ca, 
                                       finalrun
                                    );

          user -> EndOfEvent(process);

          if (fabs(trueweight) > 1.e-50)
              return_weight /= trueweight;
          else
              return_weight = 0.;

          caFree -> ReturnObject(ca);

         // clean up the list with the events
         GetEFree() -> ReturnObjectsFromUseList(
                          &Event :: GetH,
                          process -> GetEventList()
                       );
     
         if (log1 > 0.00000000001 && log2 > 0.00000000001)
            To_status() << FS("f=%13.6e ", frac)
                        << FS("w=%13.6e ", return_weight)
                        << FS("p=%10.3e ", return_weight * compensate)
                        << FS("1=%10.3e ",  return_weight / log1)
                        << FS("2=%10.3e\n", return_weight / log2);
         else
            To_status() << FS("f=%13.6e ", frac) 
                        << FS("w=%13.6e ", return_weight)
                        << FS("p=%10.3e\n", return_weight * compensate); 
      }

   return 1.;
}

// ----------------------------------------------------------------------------

// --> the Monte Carlo Integrand

double Disaster :: Integrand(
                      double* x, 
                      double wgt, 
                      int icalls, 
                      int iterations
                   ) {

   ++ event_counter;

   if (print_progress > 0 && event_counter % print_progress == 0) {
      To_status() << FS("# = %9d\n", event_counter);
      To_status().Flush();
   }

   if (event_counter < min_points)
      return 0.;

   technicalCut_flag = False;
   Boolean fpError   = False;

   if (fpnan -> Get_Error()) {
      fpError = True;
      ++ fpError_counter;
      if (fpError_print_counter.IncrementAndCheck()) {
         To_error() 
            << "Floating point error in Disaster :: Integrand (1)\n"
            << FS("event #   %8d\n", event_counter);
         fpnan -> PrintStatus(To_error());
      }
      if (fpError_print_counter.BorderlineCase())
         To_error() << "fpError: no more warning messages\n";
      fpError_log_counter.Increment();
      fpnan -> Reset();
   }

   // store this for use in the user routine
   unitCube = x;

   double trueweight = wgt / (icalls * iterations);

   double return_weight;

   // transform to avoid phase space boundaries
   double y[10];
   double eps = mc_eps;
   register int i;
   for (i = 0; i < dimension; ++i)
       y[i] = eps + (1. - 2. * eps) * x[i];

   VEGAS_weight_current = trueweight;
   ++VEGAS_pts_in_last_iteration;

   if ( ! limitTestFlag) {

      ContributionArray* ca;

      ca = caFree -> GetDefinedObject(this, 100, 11);

      // check the list with the events
      if (process -> GetEventList() -> IsNotEmpty()) 
         errf(-1, "Disaster :: Integrand: event list is not empty");

      user -> BeginEvent(process);

      process -> CalculateWeight(y, ca);

      // check for NaN --- temporarily disabled
      // fpnan -> LogDoubleStatus(return_weight);

      if (fpnan -> Get_Error()) {
         fpError = True;
         ++ fpError_counter;
         if (fpError_print_counter.IncrementAndCheck()) {
            To_error() 
               << "Floating point error in Disaster :: Integrand (2)\n"
               << FS("event #   %8d\n", event_counter);
            fpnan -> PrintStatus(To_error());
         }
         if (fpError_print_counter.BorderlineCase())
            To_error() << "fpError: no more warning messages\n";
         fpError_log_counter.Increment();
         fpnan -> Reset();
      }

      if (technicalCut_flag) {
         ++ technicalCut_counter;
         if (technicalCut_print_counter.IncrementAndCheck())
            To_error() 
               << "technical cut applied in Disaster :: Integrand\n"
               << FS("event #   %8d\n", event_counter);
         if (technicalCut_print_counter.BorderlineCase())
            To_error() << "technical cut applied: no more warning messages\n";
         technicalCut_log_counter.Increment();
      }

      if (fpError || technicalCut_flag) {
         user -> DropEvent(process);
         return_weight = 0.0;
      } else {
         user -> AcceptEvent(process);
         return_weight = user -> EndEvent(process, ca, finalrun);
         // VEGAS may return a ``trueweight'' of 0. Since we multiply, we have
         // to divide here. If the ``trueweight'' is zero, then return_weight
         // ought to be zero.
         if (fabs(trueweight) > 1.e-50)
            return_weight /= trueweight;
         else
            return_weight = 0.;
      }

      user -> EndOfEvent(process);

      caFree -> ReturnObject(ca);

      // clean up the list with the events
      GetEFree() -> ReturnObjectsFromUseList(
                       &Event :: GetH,
                       process -> GetEventList()
                    );

   } else {  // ! limitTestFlag

      // check the scaling in the limit of the integration variables
      return_weight = ScalingCheck(y, trueweight);

   }  // ! limitTestFlag

   if (max_points > 0 && event_counter >= max_points) {
      errf(-1, "forced stop by max_points.");
   }

   if (fpnan -> Get_Error()) {
      fpError = True;
      ++ fpError_counter;
      if (fpError_print_counter.IncrementAndCheck()) {
         To_error() 
            << "Floating point error in Disaster :: Integrand (3)\n"
            << FS("event #   %8d\n", event_counter);
         fpnan -> PrintStatus(To_error());
      }
      if (fpError_print_counter.BorderlineCase())
         To_error() << "fpError: no more warning messages\n";
      fpError_log_counter.Increment();
      fpnan -> Reset();
   }

   if (fpError || technicalCut_flag) {

      if (   fpError           && fpError_log_counter.Check()
          || technicalCut_flag && technicalCut_log_counter.Check()
         ) {
         // dump event to log file
         To_log() << "# --------------------------------------"
                  << "----------------------------------------\n";
         To_log() << FS("# event #   %8d\n", event_counter);
         if (fpError)
            To_log() << "# fpError\n";
         if (technicalCut_flag)
            To_log() << "# technicalCut_flag\n";
         To_log() << "# --------------------------------------"
                  << "----------------------------------------\n";
         To_log() << FS("LABEL             %30d\n", - event_counter);
         To_log() << FS("DIMENSION_DEBUG   %30d\n", dimension);
         for (int co = 0; co < dimension; ++co)
             To_log() << FS("X_VARIABLE        %30.22e", x[co])
                      << FS("   # [%2d]\n", co);
         To_log() << FS("WGT_DEBUG         %30.22e\n", wgt);
         To_log() << FS("ICALLS_DEBUG      %30d\n", icalls);
         To_log() << FS("ITERATIONS_DEBUG  %30d\n", iterations);
         To_log() << "RUN_ONE_EVENT __dummy__\n";
         To_log().Flush();
      }

      if (technicalCut_log_counter.BorderlineCase())
         To_log() << "# technical cut applied: no more debug records\n";

      if (fpError_log_counter.BorderlineCase())
         To_log() << "# fpError applied: no more debug records\n";
   }

   return return_weight;
}

// ----------------------------------------------------------------------------

// prepare for Monte Carlo:
// parton densities, coupling constant, integration boundaries, etc.

void Disaster :: PrepareForMC() {

   fpnan -> Modify_fpHandler(if_fpHandler);

   if (fpnan -> Get_Error()) {
      To_error() << "Floating point error in Disaster :: PrepareForMC (1)\n";
      To_log() << "#\n# --------------------------------------"
               << "----------------------------------------\n";
      To_log() << "# Floating point error in Disaster :: PrepareForMC (1)\n";
      To_log() << "# --------------------------------------"
               << "----------------------------------------\n#\n";
      fpnan -> PrintStatus(To_error());
      fpnan -> Reset();
   }

   To_status() << FS("mc_points_prep = %8d\n", mc_points_prep);
   To_status() << FS("mc_points_final = %8d\n", mc_points_final);

   switch (process_index) {
      case -7: 
         process = new DISProcessCollinearInitial_Test;
         break;
      case -6: 
         process = new DISProcessAddedSubtraction_Test;
         break;
      case -5: 
         process = new DISProcessSubtraction_Test;
         break; 
      case -4: 
         process = new DISProcessNoSubtraction_Test;
         break;
      case -3: 
         process = new TestProcess_3;
         break;
      case -2: 
         process = new TestProcess_2;
         break;
      case -1: 
         process = new TestProcess_1;
         break;
      case  1: 
         process = new DISProcessLO;
         break;
      case  2: 
         process = new DISProcessNLO;
         break;
      case  3: 
         process = new DISProcessLO_Test;
         break;
      case  4: 
         process = new DISProcessNLO_Test;
         break;
      case  5: 
         process = new DISProcessNLO_Light;
         break;
      default:
         errf(-1,"Disaster :: RunMC: process unknown");
   }

   process -> Define(this);
   process -> AssignMatrixElements();

   To_status() << FS("n_in_partons=        %d\n\n",process -> GetIn());

   // prepare for up to 5 final state partons
   pMapSquareLinear -> Define(5, process -> GetIn()); 
   pMapSquareLinear -> CreateListOfLists(1, 5); 
//   pMapSquareLinear->Print(stdout);

   To_status() << FS("min. final state partons=%d\n\n", process -> GetMinOut());
   To_status() << FS("max. final state partons=%d\n\n", process -> GetMaxOut());

   dimension = process -> Get_Max_nIntVar();

   switch(process -> GetIn()) {
      case 0: 
         break;
      case 1: 
         if (   lepton_integration && !xi_integration
             && xi_fixed <0.999999)
            errf(-1, "xi_fixed must be 1. if xB, y are integrated over!");
         if (lepton_integration == 1)
            dimension += 2; 
         if (xi_integration == 1)
            dimension += 1; 
         break;
      case 2: 
         if (xi_integration == 1)
            dimension += 2; 
         break;
      default: 
         errf(-1, "too many incident particles");
   }

   if (dimension > 10)
      errf(-1, "too many integration variables");

   To_status() << FS("Number of integration variables: %d\n\n", dimension);

   if (explicitFrameChoice >= 0)
      process -> SetFrame((Frame) explicitFrameChoice);
   To_status() << FS("Frame: %s\n", FrameName[process -> GetFrame()]);

   event_counter = 0;

   fpError_counter = 0;
   fpError_print_counter.Set_maximum(5);
   fpError_print_counter.Reset();
   fpError_log_counter.Set_maximum(100);
   fpError_log_counter.Reset();

   technicalCut_counter = 0;
   technicalCut_print_counter.Set_maximum(5);
   technicalCut_print_counter.Reset();
   technicalCut_log_counter.Set_maximum(10);
   technicalCut_log_counter.Reset();

   // check momenta permutations (if any)
   register int i;
   for (i = 1; i < ps_permute_max; ++i)
       if (ps_permute_n[i] > 0 && ps_permute_n[i] != i)
          errf(-1, "Disaster :: PrepareForMC: permutation incomplete");

   double dlbd, dubd;

   SH = pow(ECMFULL, 2);

   switch (process -> GetIn()) {

      case 1:

         Q2min_actual = Q2min;
         Q2max_actual = Q2max;
         W2min_actual = W2min;
         W2max_actual = W2max;

         if (Q2max_actual >= SH)
            Q2max_actual = SH;

         if (W2max_actual >= SH)
            W2max_actual = SH;

         if (ymin < 1.e-10) {
            dlbd = 0.;
            dubd = 1.;
         } else {
            dlbd = 1. - W2max_actual / SH / ymin;
            dubd = Q2max_actual / SH / ymin;
         }
         xBmin_actual 
            = dmax5(xBmin,
                    Q2min_actual / SH,
                    Q2min_actual / (Q2min_actual + W2max_actual),
                    dlbd,
                    Q2min_actual / SH / ymax
                   );
         xBmax_actual 
            = dmin5(xBmax,
                    1. - W2min_actual / SH,
                    Q2max_actual / (Q2max_actual + W2min_actual),
                    dubd,
                    1. - W2min_actual / SH / ymax
                   );

         To_status() << "      Range of Lepton variables:\n";
         To_status() << FS("%16.6e <= x_B <= ", xBmin_actual)
                     << FS("%16.6e\n",          xBmax_actual);
         To_status() << FS("%16.6e <=  y  <= ", ymin)
                     << FS("%16.6e\n",          ymax);
         To_status() << FS("%16.6e <=  Q  <= ", sqrt(Q2min_actual))
                     << FS("%16.6e\n",          sqrt(Q2max_actual));
         To_status() << FS("%16.6e <=  W  <= ", sqrt(W2min_actual))
                     << FS("%16.6e\n",          sqrt(W2max_actual));
         To_status() << "\n";

         if (     
                xBmin_actual >= xBmax_actual
             || ymin         >= ymax
             || Q2min_actual >= Q2max_actual
             || W2min_actual >= W2max_actual
            )
            errf(-1,"Lepton phase space is empty");

         break;
   }

   if (fpnan -> Get_Error()) {
      To_error() << "Floating point error in Disaster :: PrepareForMC (2)\n";
      To_log() << "#\n# --------------------------------------"
               << "----------------------------------------\n";
      To_log() << "# Floating point error in Disaster :: PrepareForMC (2)\n";
      To_log() << "# --------------------------------------"
               << "----------------------------------------\n#\n";
      fpnan -> PrintStatus(To_error());
      fpnan -> Reset();
   }
}

// ----------------------------------------------------------------------------

// clean up after running the Monte Carlo

void Disaster :: CleanUpAfterMC() {

   // delete the mess
   pMapSquareLinear -> Delete();

   // delete the process
   delete process;

   if (fpnan -> Get_Error()) {
      To_error() 
         << "Floating point error in Disaster :: CleanupUpAfterMC (1)\n";
      To_log() << "#\n# --------------------------------------"
               << "----------------------------------------\n";
      To_log() 
         << "# Floating point error in Disaster :: CleanupUpAfterMC (1)\n";
      To_log() << "# --------------------------------------"
               << "----------------------------------------\n#\n";
      fpnan -> PrintStatus(To_error());
      fpnan -> Reset();
   }
}

// ----------------------------------------------------------------------------

void Disaster :: TestMCFunction() {
   
   To_status() << "Disaster :: TestMCFunction called\n";
}

// ----------------------------------------------------------------------------
//
// --> the main routine
//
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------

extern "C" { 
int ppcmain_c_from_f() {

   // -------------------------------------------------------------------------

#if PSW_MACHINE

   // at PSI: need to reconfigure memory allocation
   // otherwise the allocation of long dynamical objects
   // uses up lots of system time!
   // (due to incompetence of the system manangers?)

#error To_status() not yet defined

   // dynamic memory allocation parameters on the local PSI machine
   To_status() << "dynamic memory allocation on the local PSI machine:\n";
   To_status() << FS("__noshrink=       %ld\n",__noshrink);
   To_status() << FS("__minshrink=      %ld\n",(unsigned long)__minshrink);
   To_status() << FS("__minshrinkfactor=%e\n",__minshrinkfactor);
   To_status() << FS("__mingrow=        %ld\n",(unsigned long)__mingrow);
   To_status() << FS("__mingrowfactor=  %e\n",__mingrowfactor);
   To_status() << FS("__madvisor=       %ld\n",__madvisor);
   To_status() << FS("__small_buff=     %ld\n",__small_buff);
   To_status() << FS("__fast_free_mac=  %d\n",__fast_free_max);
   To_status() << FS("__sbrk_override=  %ld\n",__sbrk_override);
   To_status() << FS("__taso_mode=      %ld\n",__taso_mode);
   To_status() << "\n\n\n";

   // this is required to avoid problems with Event Copy(100); etc.
   __noshrink = 1;

#endif

   // -------------------------------------------------------------------------

   Disaster* disaster = new Disaster;

   // -------------------------------------------------------------------------

   FILE* ParametersIn;
   ParametersIn = stdin;

   // --> job and file names
   
   String job_name(1000);
   String data_path(1000);

   // read ``data_path'' and ``job_name''
   data_path.ReadBoundedFromFile(ParametersIn);
   job_name.ReadBoundedFromFile(ParametersIn);

   // -------------------------------------------------------------------------

   // open files and streams

   // open the file for the output of plot data
   String PlotOutName(1000);
   PlotOutName.CopyFromString(&data_path);
   PlotOutName.AppendFromChar(".plot");
   FILE* PlotOut = fopen(PlotOutName.data, "w");

   // open the file for the output of the log file
   String logOutName(1000);
   logOutName.CopyFromString(&data_path);
   logOutName.AppendFromChar(".log");
   FILE* logOut = fopen(logOutName.data, "w");

   disaster -> To_status().Associate(print_f_C, 6);
   disaster -> To_error().Associate(print_f_C, 6);
   disaster -> To_log().Associate(logOut);   

//   To_status().Associate(strFORT_x, 6);
//   To_error().Associate(strFORT_x, 6);

   disaster -> To_status() << FS("data path = :%s:\n", data_path.data);
   disaster -> To_status() << FS("job name = :%s:\n", job_name.data);
   disaster -> To_status() << "\n\n\n";

   // -------------------------------------------------------------------------

   disaster -> Define();

   // -------------------------------------------------------------------------

   // define a user interface and associate it with disaster

   MyUserClass* user;
   user = new MyUserClass;
   user -> Define(disaster);
   disaster -> SetUser(user);

   // -------------------------------------------------------------------------

   // now run the MC: read parameters from the file ``ParametersIn''

   disaster -> Start();
   disaster -> ParseFile(ParametersIn);
   int disaster_return_value = disaster -> Get_disaster_return_value();
   disaster -> Terminate();

// ----------------------------------------------------------------------------

   // close files and streams

   disaster -> To_status().Flush();
   disaster -> To_error().Flush();
   disaster -> To_log().Flush();

   // just in case...
   disaster -> To_status().Associate(&cout);
   disaster -> To_error().Associate(&(disaster -> To_status()));
   disaster -> To_log().Associate(&(disaster -> To_status()));

   if (PlotOut != NULL)
      fclose(PlotOut);

   if (logOut != NULL)
      fclose(logOut);

   // -------------------------------------------------------------------------

   delete user;
   delete disaster;

   // -------------------------------------------------------------------------

   return disaster_return_value;
}}

// ----------------------------------------------------------------------------
//
// wrapper for C functions
//
// ----------------------------------------------------------------------------
       
   FCALLSCFUN0(                           \
      INT, ppcmain_c_from_f, PPCMAIN, ppcmain  \
   )

// ----------------------------------------------------------------------------

// dummy to avoid warning messages (static functions...)
   
void dummyfunctionDisaster() {
   static char c[]="a";
   c2fstrv(c,c,0,0);
   f2cstrv(c,c,0,0);
   kill_trailing(c,'a');
   vkill_trailing(c,0,0,'a');
   num_elem(c,0,0,0);
}
    
// ----------------------------------------------------------------------------

// template instantiations

I_FreeList_1(Particle, TheoryEvaluation*)

I_FreeListArray_1(Event, TheoryEvaluation*)

I_Array(Contribution)
I_FreeListArray_1(Contribution, Global*)

I_FreeListTwoD_Array_1(ContributionArray, Global*)

I_Array(PartonList)
I_FreeListArray_1(PartonList, Global*)

I_FreeListArray_1(PartonDensity, Global*)
I_Cache(PartonDensity)

I_ArrayOfPointers(Book*)

I_FreeList_1(PDFC, Disaster*)
I_Cache(PDFC)

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
// ============================================================================
//
// --> MC integration classes
//
// file:              mcintg.cc
// created:           29.03.1996
// last modification: 14.11.1997
//
// ============================================================================

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "fortran.h"

#include "global.h"
#include "mcintg.h"
#include "mth.h"
#include "user.h"

// ============================================================================
//
// definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

// ============================================================================
//
// --> class-unrelated functions
//
// ============================================================================

// ============================================================================
//
// --> FORTRAN interface routines
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> interface to VEGAS
//
// ----------------------------------------------------------------------------

// wrapper for FORTRAN function to call VEGAS

   PROTOCCALLSFSUB2(
      Underscore_FORTRAN_Name(VEGAS1_F), Underscore_FORTRAN_Name(vegas1_f),
      DOUBLEV, 
      INTV
   )
        
#define vegas_f_Macro(arg1, arg2) \
           CCALLSFSUB2( \
              Underscore_FORTRAN_Name(VEGAS1_F), \
              Underscore_FORTRAN_Name(vegas1_f), \
              DOUBLEV, \
              INTV, \
              arg1, arg2 \
           )

// the following is the function (with C linkage!) to be called.
                    
void vegas1_C(
        double* acc,
        int* ndim,
        int* ncall,
        int* itmxin,
        int* nprn,
        int* igraph,
        int* ientry,
        int* ninpin,
        int* noutpin,
        double* alphin, 
        int* mdsin,
        double* avgtout,
        double* errtout,  
        double* avgiout,
        double* erriout,
        double* callsout,
        int* itmxout
     ) {

   double dpar[20];
   int    ipar[20];

   dpar[0] = *acc;
   dpar[6] = *alphin;
      
   ipar[0] = *ndim;
   ipar[1] = *ncall;
   ipar[2] = *itmxin;
   ipar[3] = *nprn;
   ipar[4] = *igraph;
   ipar[5] = *ientry;
   ipar[6] = *ninpin;
   ipar[7] = *noutpin;
   ipar[9] = *mdsin;

   vegas_f_Macro(dpar, ipar);

   *avgtout  = dpar[1];
   *errtout  = dpar[2];
   *avgiout  = dpar[3];
   *erriout  = dpar[4];
   *callsout = dpar[5];
      
   *itmxout  = ipar[8];
}

// ----------------------------------------------------------------------------

// dummy to avoid warning messages (static functions...)

void dummyfunctionFold_mcintg() {
   static char c[]="a";
   c2fstrv(c,c,0,0);
   f2cstrv(c,c,0,0);
   kill_trailing(c,'a');
   vkill_trailing(c,0,0,'a');
   num_elem(c,0,0,0);
}

// ----------------------------------------------------------------------------
//
// --> function that is integrated by VEGAS; 
//     reference to 'mcfunction'
//
// ----------------------------------------------------------------------------

extern "C"{

    double fvegas1_c_from_f(double* x, double* wgt, double* calls, int* itmx) {
       
       Global* global = gl;  
       
       return global 
              -> mc_function
              -> Integrand(x, *wgt, (int) rclv(*calls), *itmx);
    }
}

// ----------------------------------------------------------------------------
//
// wrapper for C functions
//
// ----------------------------------------------------------------------------

   FCALLSCFUN4(                           \
      DOUBLE,                             \
      fvegas1_c_from_f, FVEGAS1, fvegas1, \
      DOUBLEV,                            \
      PDOUBLE,                            \
      PDOUBLE,                            \
      PINT                                \
   )     

// ============================================================================
//
// --> classes
//
// ============================================================================
   
// ----------------------------------------------------------------------------
//
// --> Functions to be integrated
//
// ----------------------------------------------------------------------------
   
MCFunction::MCFunction(){
}

// ----------------------------------------------------------------------------
   
MCFunction::~MCFunction(){
}

// ----------------------------------------------------------------------------
//
// --> Monte Carlo framework
//
// ----------------------------------------------------------------------------

MCframe :: MCframe() {

   global = NULL;
}

// ----------------------------------------------------------------------------

MCframe :: ~MCframe() {
}

// ----------------------------------------------------------------------------

void MCframe :: Integrate(
                   int iterations,
                   int npoints1,
                   int npoints2,  
                   double* value,
                   double* error
                ) {

   double avgtOut, errtOut, callsOut;

   if (dimension < 0) {
      errf(-1, "MCframe :: Integrate: dimension < 0");
   } else
   if (dimension == 0) {

      // 0-dimensional integration: evaluate integrand just once!
      global -> To_status() 
         << "No integration; integrand is evaluated once ...\n";
      global -> finalrun = TRUE;
      global -> user -> StartIntegration1();
      global -> user -> StartIntegration2();

      callsOut = 1.;
      errtOut  = 0.;
      avgtOut  = Integrand(NULL, 1.0, 1, 1);

   } else {
 
      // integrate 

      int ncall, itmxIn, ientry;
      double erriOut, avgiOut;
      int itmxOut;

      // VEGAS parameters (default)
      double acc     = 1.e-10;
      int    ndim    = dimension;
      int    nprn    = 1;
      int    igraph  = 0;
      int    ninpIn  = 5;
      int    noutpIn = 6;
      double alphIn  = 1.5;
      int    mdsIn   = 0;   // use importance sampling only!

      global -> user -> StartIntegration1();
      if (iterations != 0) {
          global -> To_status() << "Grid definition run ...\n";
          global -> finalrun = FALSE;
          itmxIn = iterations;
          ncall = npoints1 / iterations;
          ientry = 0;
          vegas1_C(
             &acc, &ndim, &ncall, &itmxIn,
             &nprn, &igraph, &ientry,
             &ninpIn, &noutpIn,
             &alphIn, &mdsIn,
             &avgtOut, &errtOut, &avgiOut,
             &erriOut,
             &callsOut, &itmxOut
          );
      }

      global -> To_status() << "Final run ...\n";
      global -> finalrun = TRUE;
      global -> VEGAS_pts_in_last_iteration=0;
      itmxIn = 1;
      ncall  = npoints2;
      ientry = (iterations == 0) ? 0 : 1;
      global -> user -> StartIntegration2();
      vegas1_C(
         &acc, &ndim, &ncall, &itmxIn,
         &nprn, &igraph, &ientry,
         &ninpIn, &noutpIn,
         &alphIn, &mdsIn,
         &avgtOut, &errtOut, &avgiOut,
         &erriOut,
         &callsOut, &itmxOut
      );

   }

   global -> To_status() << FS("avgt=%e\n", avgtOut);
   global -> To_status() << FS("errt=%e\n", errtOut);
   global -> To_status() << FS("calls=%e\n", callsOut);

   global -> user -> EndIntegration(avgtOut, errtOut);

   *value = avgtOut;
   *error = errtOut;
}

// ----------------------------------------------------------------------------

void MCframe :: AssignGlobal(Global* globalIn) {

   global = globalIn;
}

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
// ============================================================================
//
// --> mathematical routines
//
// file:              mth.cc
// created:           29.03.1997
// last modification: 30.11.1997
//
// ============================================================================

#include <stdio.h>
#include <math.h>

#include "global.h"
#include "mth.h"

// ============================================================================
//
// definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

// mathematical constants

double Pi;
double PiSquared;
double TwoPi;
double FourPi;
    
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

void SetMathConst() {
   Pi        = 4. * atan(1.);
   PiSquared = Pi * Pi;
   TwoPi     = 2. * Pi;
   FourPi    = 4. * Pi;
}

// ----------------------------------------------------------------------------
//
// --> spherical geometry
//
// ----------------------------------------------------------------------------

// calculate polar angle

double ThreePolarVWCos(double v, double w, double cosphi) {

   double t = v * (1. - v) * w * (1. - w);
   if (t < 0. || t > 1.) {
      gl -> To_error() << FS("t=%30.22e\n", t);
      const char msg1[] = "ThreePolar out of range!";
      fpnan -> Set_flag_mth(-1, msg1);
      errf(0, msg1);
      if (t < 0.)
         t = 0.;
      else
         t = 1.;
   }   

   return v + w - 2. * v * w -2. * sqrt(t) * cosphi;
}

// ----------------------------------------------------------------------------

// calculate polar angle

double ThreePolar(double v, double w, double phi) {

   return ThreePolarVWCos(v, w, cos(phi));
}

// ----------------------------------------------------------------------------

// calculate cosine of the azimuthal angle

double GetCosPhi(double e, double v, double w) {

   double cosphi;

   double t = v * (1. - v) * w * (1. - w);

   if (t < 0.) 
      if (t < -1.e-8) {
         gl -> To_error() << FS("e=%30.22e\n", e);
         gl -> To_error() << FS("v=%30.22e\n", v);
         gl -> To_error() << FS("w=%30.22e\n", w);
         gl -> To_error() << FS("t=%30.22e\n", t);
         const char msg1[] = "GetCosPhi out of range (1)!";
         fpnan -> Set_flag_mth(-1, msg1);
         errf(0, msg1);
         t = 0.;
      } else {
         t = 0.;
      }
   if (t > 1.) 
      if (t > 1 + 1.e-8) {
         gl -> To_error() << FS("e=%30.22e\n", e);
         gl -> To_error() << FS("v=%30.22e\n", v);
         gl -> To_error() << FS("w=%30.22e\n", w);
         gl -> To_error() << FS("t=%30.22e\n", t);
         const char msg3[] = "GetCosPhi out of range (3)!";
         fpnan -> Set_flag_mth(-1, msg3);
         errf(0, msg3);
         t = 1.;
      } else {
         t = 1.;
      }

   if (t == 0.) {

      cosphi = 1.;

   } else {

      double ratio = (e - (v + w - 2. * v * w)) / (-2. * sqrt(t));

      // do this silently (--> e == 0, v == w happens frequently
      //                       in collinear subtraction terms)
      
      if (ratio >  1.0 && ratio < 1.0000001)
         ratio = 1.0;
      else
      if (ratio < -1.0 && ratio > -1.0000001)
         ratio = -1.0;
      else
      if (ratio > 1.0 || ratio < -1.0) {
         gl -> To_error() << FS("e=%30.22e\n", e);
         gl -> To_error() << FS("v=%30.22e\n", v);
         gl -> To_error() << FS("w=%30.22e\n", w);
         gl -> To_error() << FS("ratio=%30.22e\n", ratio);
         const char msg5[] = "GetCosPhi out of range (5)!";
         fpnan -> Set_flag_mth(-1, msg5);
         errf(0, msg5);
         if (ratio > 1.0)
            ratio = 1.0;
         else
            ratio = -1.0;
      }

      cosphi = ratio;

   }

   return cosphi;
}

// ----------------------------------------------------------------------------

// calculate azimuthal angle

double GetPhi(double e, double v, double w) {

   return acos(GetCosPhi(e, v, w));
}

// ----------------------------------------------------------------------------

// absolute sin(phi) from cos(phi)

double AbsSinFromCos(double cosPhi) {

   return sqrt(1. - cosPhi * cosPhi);
}

// ----------------------------------------------------------------------------

double atan2ZeroTwoPi(double y, double x) {

   double phi = atan2(y, x);
   if (phi < 0)
      phi += TwoPi;
   return phi;
}

// ----------------------------------------------------------------------------
//
// --> mathematics
//
// ----------------------------------------------------------------------------

// factorial (long)

long lfactorial(int n) {
   int i;
   long fact=1;
   for (i=1;i<=n;++i) {
       fact*=i;
   }
   return fact;
}

// ----------------------------------------------------------------------------

// factorial

double dfactorial(int n) {
   return (double) lfactorial(n);
}

// ----------------------------------------------------------------------------

// minimum

int min(int a,int b) {
   return a < b ? a : b;
}

// ----------------------------------------------------------------------------

double dmin(double a,double b) {
   return a < b ? a : b;
}

// ----------------------------------------------------------------------------

double dmin3(double a,double b,double c) {
   return dmin(dmin(a,b),c);
}

// ----------------------------------------------------------------------------

double dmin4(double a,double b,double c,double d) {
   return dmin(dmin(dmin(a,b),c),d);
}

// ----------------------------------------------------------------------------

double dmin5(double a,double b,double c,double d,double e) {
   return dmin(dmin(dmin(dmin(a,b),c),d),e);
}

// ----------------------------------------------------------------------------

// maximum

int max(int a,int b) {
   return a > b ? a : b;
}

// ----------------------------------------------------------------------------

double dmax(double a,double b) {
   return a > b ? a : b;
}

// ----------------------------------------------------------------------------

double dmax3(double a,double b,double c) {
   return dmax(dmax(a,b),c);
}

// ----------------------------------------------------------------------------

double dmax4(double a,double b,double c,double d) {
   return dmax(dmax(dmax(a,b),c),d);
}

// ----------------------------------------------------------------------------

double dmax5(double a,double b,double c,double d,double e) {
   return dmax(dmax(dmax(dmax(a,b),c),d),e);
}

// ----------------------------------------------------------------------------

// sign

int sign(double x) {
   if (x>0.)
      return 1;
   else 
   if (x<0.)
      return -1;
   else 
      return 0;
}

// ----------------------------------------------------------------------------

// long sign

int lsign(long x) {
   if (x>0)
      return 1;
   else 
   if (x<0)
      return -1;
   else 
      return 0;
}

// ----------------------------------------------------------------------------

// correct range [0,1]

void correct_range01(double *x) {
   if (*x < 0.)
      if (*x > -1.e-8) {
         *x = 0.;
      } else {
         gl -> To_error() << FS("x=%30.22e\n", *x);
         const char msg2[] = " correct_range01 (2)!";
         fpnan -> Set_flag_mth(-1, msg2);
         errf(0, msg2);
         *x = 0.;
      }

   if (*x > 1.)
      if (*x < 1. + 1.e-8) {
         *x = 1.;
      } else {
         gl -> To_error() << FS("x=%30.22e\n", *x);
         const char msg4[] = " correct_range01 (4)!";
         fpnan -> Set_flag_mth(-1, msg4);
         errf(0, msg4);
         *x = 1.;
      }
}

// ----------------------------------------------------------------------------

// return closest long value

INLINE long rclv(double x) {
   return x>=0?(long)(floor(x+0.5)+0.1):(long)(floor(x+0.5)-0.1);
}

// ----------------------------------------------------------------------------

// maps the interval [0,1] to the long interval [0,m] (inclusive)

long mapDoubleToLong(double x, long m) {

   long ret = rclv(floor(x * (m+1)));
   if (ret < 0)
      ret = 0;
   else if (ret > m)
      ret = m;

   return ret;
}

// ----------------------------------------------------------------------------

// get mantissa and exponent

void gmae(double in, double &mantissa, int &exponent) {
   int expo;
   mantissa=frexp(in, &expo);
   exponent = expo;
}

// ----------------------------------------------------------------------------

// Mathematica format for floating point

void mfffp(double in, char *str) {
   double mantissa;
   int expo;
   gmae(in, mantissa, expo);
   sprintf(str, "%26.18f * 2^%d", mantissa, expo);
}

// ----------------------------------------------------------------------------

// check whether i0 <= i <= i1.
// returns 0 or 1.

   // -------------------------------------------------------------------------

int inRange(int i, int i0, int i1) {
   return (i0 <= i && i <= i1) ? 1 : 0;
}

   // -------------------------------------------------------------------------

int inRange(double i, double i0, double i1) {
   return (i0 <= i && i <= i1) ? 1 : 0;
}

// ----------------------------------------------------------------------------
//
// check whether a and b are compatible.
// epsAbs is absolute measure (for small values), 
// epsRel is relative measure (for large values)
//
// problem: if a, b are equal with oposite sign, // :_CAUTION_:
//          this routine crashes 
//
// ----------------------------------------------------------------------------

int compareAbsRel(double a, double b, double epsAbs, double epsRel) {
    if (   fabs(a) < epsAbs && fabs(b) < epsAbs
        || fabs(a/(a+b)-0.5) < epsRel           )
       return 1;
    else
       return 0;
}

// ----------------------------------------------------------------------------
//
// compare two ``double''s 
//
// ----------------------------------------------------------------------------

int compareDouble(
       double a, 
       double b,
       double eps, 
       int method
    ) {

   int equal = 1;

   switch (method) {

      case 0:
         if (  fabs(a - b)
             > fabs(b * eps))
            equal = 0;
         break;

     default:
        gl -> To_error() << FS("method=%d\n", method);
        errf(-1, "compareDouble: method not known");
   }

   return equal;
}

// ----------------------------------------------------------------------------
//
// returns 0.0; useful to provoke floating point exceptions
// (compilers optimize!)
//
// ----------------------------------------------------------------------------

double ReturnZeroZero() {
   return 0.0;
}

// ----------------------------------------------------------------------------
//
// returns a random number (naive) between 0 and 1
//
// ----------------------------------------------------------------------------

double naiveRandomNumber() {
   return (double) rand() / (double) (RAND_MAX+1);
}

// ----------------------------------------------------------------------------
//
// variable transformation f : [0,1] -> [0,1]
//
// index = 1:
//    for x <= lambda: ~ x^alpha
//    for x >  lambda: ~ 1 - b * (1-x)^beta
//    at  x == lambda" mu
//
// also calculates the jacobian for the measure df = jacobian dx
//
// ----------------------------------------------------------------------------

double variableTransformation(
          double  x, 
          int     index,
          double  alpha,
          double  beta, 
          double  lambda,
          double  mu,
          double* jacobian
       ) {

   double term;

   switch (index) {

      case 1:
         if (x <= lambda) {
            if (x < 1.e-10) {
               term = 0.;
               *jacobian = 0.;
               const char msg1[] = "variableTransformation (1)";
               fpnan -> Set_flag_mth(-1, msg1);
            } else {
               term = mu * pow(x / lambda, alpha);
               *jacobian = alpha / x * term;
            }
         } else {
            if (x > 1.-1.e-10) {
               term = 1.;
               *jacobian = 0.;
               const char msg2[] = "variableTransformation (2)";
               fpnan -> Set_flag_mth(-1, msg2);
            } else {
               term = (1. - mu) * pow((1. - x) / (1. - lambda), beta);
               *jacobian = beta / (1. - x) * term;
               term = 1. - term;
            }
         }
         break;
         
      default:
         errf(-1, "variableTransformation: no such index");
         *jacobian = 0.;
         term = 0.;
   }

#if 0
   gl -> To_status() << "transformation:\n";
   gl -> To_status() << FS("%30.22e\n", x);
   gl -> To_status() << FS("%30.22e\n", alpha);
   gl -> To_status() << FS("%30.22e\n", beta);
   gl -> To_status() << FS("%30.22e\n", lambda);
   gl -> To_status() << FS("%30.22e\n", mu);
   gl -> To_status() << FS("%30.22e\n", term);
   gl -> To_status() << FS("%30.22e\n", *jacobian);
#endif

   return term;
}

// ----------------------------------------------------------------------------
//
// inverse variable transformation f^(-1): [0,1] -> [0,1]
//
// ----------------------------------------------------------------------------

double variableTransformationInverse(
          double  y, 
          int     index,
          double  alpha,
          double  beta, 
          double  lambda,
          double  mu
       ) {

   double term;

   switch (index) {

      case 1:
         if (y <= mu) {
            term = lambda * pow(y / mu, 1. / alpha);
         } else {
            term = 1. - (1. - lambda) * pow((1. - y) / (1. - mu), 1. / beta);
         }
         break;
         
      default:
         errf(-1, "variableTransformationInverse: no such index");
         term = 0.;
   }

#if 0
   gl -> To_status() << "transformation:\n";
   gl -> To_status() << FS("%30.22e\n", x);
   gl -> To_status() << FS("%30.22e\n", alpha);
   gl -> To_status() << FS("%30.22e\n", beta);
   gl -> To_status() << FS("%30.22e\n", lambda);
   gl -> To_status() << FS("%30.22e\n", mu);
   gl -> To_status() << FS("%30.22e\n", term);
   gl -> To_status() << FS("%30.22e\n", *jacobian);
#endif

   return term;
}

// ----------------------------------------------------------------------------
//
// variable transformation f : [0,1] -> [0,1]
//
// index = 1:
//    x^alpha
//
// also calculates the jacobian for the measure df = jacobian dx
//
// ----------------------------------------------------------------------------

double variableTransformationLower(
          double  x, 
          int     index,
          double  alpha,
          double* jacobian
       ) {

   double term;

   switch (index) {

      case 1:
         term = pow(x, alpha);
         *jacobian = alpha / x * term;
         break;
         
      default:
         errf(-1, "variableTransformationLower: no such index");
         *jacobian = 0.;
         term = 0.;
   }

   return term;
}

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
// ============================================================================
//
// --> process classes
//
// file:              proc.cc
// created:           01.03.1997
// last modification: 15.12.1997
//
// ============================================================================

#include <stdlib.h>

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "global.h"
#include "mth.h"
#include "cmb.h"
#include "user.h"
#include "strng.h"
#include "rest.h"
#include "qcd.h"
#include "disaster.h"
#include "container.h"

#include "proc.h"

// ============================================================================
//
// definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

const char* FrameName[] = {
               "none",
               "original",
               "pCMS",
               "hCMS",
               "breit",
               "lab"
            };  

const char* P0_ReferenceName[] = {
               "IS_REFERENCE",
               "FS_REFERENCE"
            };

const char* LimitTypeName[] = {
               "NO_LIMIT",
               "SOFT",
               "COLLINEAR",
               "SOFT_AND_COLLINEAR",
               "SOFT_ADDED",
               "COLLINEAR_ADDED",
               "SOFT_AND_COLLINEAR_ADDED"
            };

const char* SubtractionTypeName[] = {
               "NoSubtraction",
               "Subtraction",
               "AddedSubtraction",
               "CollinearInitial",
               "Collinear"
            };

const char* EventTypeName[] = {
               "EventTypeNone",
               "Additional",
               "Final"
            };

// ============================================================================
//
// --> class-unrelated functions
//
// ============================================================================

// ============================================================================
//
// --> classes
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> Particles: 
//     the default representation is in polar coordinates
//     `energy', `momentum', `v', `phi'
//
// ----------------------------------------------------------------------------

Particle :: Particle() {

#if CHECKDESTRUCT
      ident = gl -> tParticle -> Mark();
#else
      ident = gl -> tParticle -> Unique();
#endif

   Define(NULL);
}

// ----------------------------------------------------------------------------

Particle :: ~Particle() {

#if CHECKDESTRUCT
      gl -> tParticle.Remove(ident);
#endif

   Delete();
}

// ----------------------------------------------------------------------------

// the overloaded assignment operator

Particle& Particle :: operator=(const Particle &a) {

   register int i;

   for (i = 0; i < 4; ++i) {
       fvect[i] = a.fvect[i];
       nvect[i] = a.nvect[i];
   }

   energy = a.energy;
   momentum = a.momentum;
   v = a.v;
   sv = a.sv;
   sinPhi = a.sinPhi;
   cosPhi = a.cosPhi;

   sinPhiRef = a.sinPhiRef;
   cosPhiRef = a.cosPhiRef;

   tEval = a.tEval;
   label = a.label;
   type  = a.type;

   return *this;
}

// ----------------------------------------------------------------------------

void Particle :: Define(TheoryEvaluation* tIn) {

   tEval = tIn;
   type  = 'X';
   label = -99;
}

// ----------------------------------------------------------------------------

void Particle :: Delete() {
}

// ----------------------------------------------------------------------------

void Particle :: Reset(TheoryEvaluation *tIn) {

   tEval = tIn;

#if CHECKDESTRUCT
      (void) gl->tParticle.Check(ident, 1);
#endif
}

// ----------------------------------------------------------------------------

void Particle :: ReturnToFreeList() {

#if CHECKDESTRUCT
      if (gl -> tParticle.Check(ident, 2) == 1) {
         tEval -> To_status() << "the event:\n";
         GetEvent() -> Print(stdout);
      }
#endif
}

// ----------------------------------------------------------------------------

// calculate the cross product p x q
// only cartesian components! Energy set to 0.

void Particle :: CrossFromCartesian(Particle* p, Particle* q) {

   fvect[0] = 0.;
   fvect[1] = p -> fvect[2] * q -> fvect[3] - p -> fvect[3] * q -> fvect[2];
   fvect[2] = p -> fvect[3] * q -> fvect[1] - p -> fvect[1] * q -> fvect[3];
   fvect[3] = p -> fvect[1] * q -> fvect[2] - p -> fvect[2] * q -> fvect[1];
}

// ----------------------------------------------------------------------------

// calculate the expression ``this . (p x q)''

double Particle :: DotIntoCrossFromCartesian(Particle* p, Particle* q) {

   return
       fvect[1] * (  p -> fvect[2] * q -> fvect[3] 
                   - p -> fvect[3] * q -> fvect[2])
     + fvect[2] * (  p -> fvect[3] * q -> fvect[1] 
                   - p -> fvect[1] * q -> fvect[3])
     + fvect[3] * (  p -> fvect[1] * q -> fvect[2] 
                   - p -> fvect[2] * q -> fvect[1]);
}

// ----------------------------------------------------------------------------

// Set to z-axis; energy component is 0.

void Particle :: SetToZ() {

   energy   = 0.;
   momentum = 1.;
   v        = 0.;
   sv       = 0.;
   sinPhi   = 0.;
   cosPhi   = 1.;
}

// ----------------------------------------------------------------------------

// Set to x-axis; energy component is 0.

void Particle :: SetToX() {

   energy   = 0.;
   momentum = 1.;
   SetV(0.5);
   sinPhi   = 0.;
   cosPhi   = 0.;
}

// ----------------------------------------------------------------------------

// Set to y-axis; energy component is 0.

void Particle :: SetToY() {

   energy   = 0.;
   momentum = 1.;
   SetV(0.5);
   SetPhi(0.5 * Pi);
}

// ----------------------------------------------------------------------------

// Set to energy axis; space component is 0

void Particle :: SetToE() {

   energy   = 1.;
   momentum = 0.;
   v        = 0;
   sv       = 0.;
   sinPhi   = 0.;
   cosPhi   = 0.;
}

// ----------------------------------------------------------------------------

// Calculate relative azimuth of a 3-vector.
// Relative to ``zaxis'' as the z-axis, and the projection 
// of ``xaxis'' as the x-axis.
// Cartesian components and nvects of all three particles are assumed to be 
// defined.

void Particle :: CalculateRelativeAzimuth(
                    Particle* zaxis,  
                    Particle* xaxis,
                    double* sinPhi_L,
                    double* cosPhi_L
                 ) {

   xaxis -> CalculateCartesianAndNormalFromPolar();
   zaxis -> CalculateCartesianAndNormalFromPolar();
   CalculateCartesianAndNormalFromPolar();

   double vxz   = VNormal(xaxis, zaxis);
   double vx    = VNormal(this, xaxis);
   double vz_L  = VNormal(this, zaxis);
   *cosPhi_L = ::GetCosPhi(vx, vz_L, vxz);
   *sinPhi_L = AbsSinFromCos(*cosPhi_L);
   if (zaxis -> DotIntoCrossFromCartesian(xaxis, this) < 0.)
      *sinPhi_L = - *sinPhi_L;
}

// ----------------------------------------------------------------------------

// :_MOD_: define by operator overloading!

void Particle :: CopyInto(Particle* into) {

   into -> energy   = energy;
   into -> momentum = momentum;
   into -> v        = v;
   into -> sv       = sv;
   into -> sinPhi   = sinPhi;
   into -> cosPhi   = cosPhi;

   into -> sinPhiRef = sinPhiRef;
   into -> cosPhiRef = cosPhiRef;

   into -> tEval = tEval;
   into -> label = label;
   into -> type  = type;
}

// ----------------------------------------------------------------------------

// ``adds'' particles

void Particle :: Add(Particle* a, Particle* b) {

   a -> CalculateCartesianFromPolar();
   b -> CalculateCartesianFromPolar();

   AddCartesian(a, b);

   CalculatePolarFromCartesian();
}

// ----------------------------------------------------------------------------

// ``adds'' particles

void Particle :: AddCartesian(Particle* a, Particle* b) {

   fvect[0] = a -> fvect[0] + b -> fvect[0];
   fvect[1] = a -> fvect[1] + b -> fvect[1];
   fvect[2] = a -> fvect[2] + b -> fvect[2];
   fvect[3] = a -> fvect[3] + b -> fvect[3];
}

// ----------------------------------------------------------------------------

// ``adds'' particles according to the E0 scheme

void Particle :: Add_E0Cartesian(Particle* a, Particle* b) {

   AddCartesian(a, b);

   double alpha =   fvect[0]
                  / sqrt(  fvect[1] * fvect[1]
                         + fvect[2] * fvect[2]  
                         + fvect[3] * fvect[3]
                    );

   fvect[1] *= alpha;
   fvect[2] *= alpha;
   fvect[3] *= alpha;
}

// ----------------------------------------------------------------------------

// ``adds'' particles according to the P scheme

void Particle :: Add_PCartesian(Particle* a, Particle* b) {

   AddCartesian(a, b);

   fvect[0] = sqrt(  fvect[1] * fvect[1]
                   + fvect[2] * fvect[2]  
                   + fvect[3] * fvect[3]
              );
}

// ----------------------------------------------------------------------------

// ``subtracts'' particles (a - b)

void Particle :: Subtract(Particle* a, Particle* b) {

   a -> CalculateCartesianFromPolar();
   b -> CalculateCartesianFromPolar();

   for (register int i = 0; i < 4; ++i) {
       fvect[i] = a -> fvect[i] - b -> fvect[i];
   }

   CalculatePolarFromCartesian();
}

// ----------------------------------------------------------------------------

// ``scales'' particle

void Particle :: ScaleBy(double scale) {

   energy   *= scale;
   momentum *= scale;
}

// ----------------------------------------------------------------------------

// ``divides'' particles (for test purposes only)

void Particle :: Divide(Particle* a, Particle* b) {

   a -> CalculateCartesianFromPolar();
   b -> CalculateCartesianFromPolar();

   for (register int i = 0; i < 4; ++i)
       fvect[i] = a -> fvect[i] / b -> fvect[i];

   CalculatePolarFromCartesian();

   label = -3;
   errf(-1, "Particle :: Divide: inconsistent!");
}

// ----------------------------------------------------------------------------

// ``scales'' particle

void Particle :: ScaleByCartesian(double scale) {

   fvect[0] *= scale;
   fvect[1] *= scale;
   fvect[2] *= scale;
   fvect[3] *= scale;
}

// ----------------------------------------------------------------------------

// linear combination of four-vectors

void Particle :: LinearCombinationCartesian(
                    double alpha, 
                    Particle* a, 
                    double beta,  
                    Particle* b
                 ) {

   fvect[0] = alpha * a -> fvect[0] + beta * b -> fvect[0];
   fvect[1] = alpha * a -> fvect[1] + beta * b -> fvect[1];
   fvect[2] = alpha * a -> fvect[2] + beta * b -> fvect[2];
   fvect[3] = alpha * a -> fvect[3] + beta * b -> fvect[3];
}

// ----------------------------------------------------------------------------

// linear combination of four-vectors

void Particle :: LinearCombination(
                    double alpha, 
                    Particle* a, 
                    double beta,  
                    Particle* b
                 ) {

   a -> CalculateCartesianFromPolar();
   b -> CalculateCartesianFromPolar();

   LinearCombinationCartesian(alpha, a, beta, b);

   CalculatePolarFromCartesian();
}

// ----------------------------------------------------------------------------

// transfers labels

void Particle :: TransferLabelsInto(Particle* p) {
   p -> label = label;
   p -> type  = type;
}

// ----------------------------------------------------------------------------

// print contents

void Particle :: Print(FILE* fout) {

   Mstream s_out;
   s_out.Associate(fout);

   s_out << FS("ident=%ld\n", ident);
   s_out << FS("N %c ", type)
         << FS("%2d ", label)
         << FS("%14.6e", energy)
         << FS("%14.6e", momentum)
         << FS("%14.6e", nvect[1])
         << FS("%14.6e", nvect[2])
         << FS("%14.6e\n", nvect[3]);
   s_out << FS("C %c ", type)
         << FS("%2d ", label)
         << FS("%14.6e", fvect[0])
         << FS("%14.6e", fvect[1])
         << FS("%14.6e", fvect[2])
         << FS("%14.6e\n", fvect[3]);
   s_out << FS("P %c ", type)
         << FS("%2d ", label)
         << FS("%14.6e", energy)
         << FS("%14.6e", momentum)
         << FS("%14.6e", v)
         << FS("%14.6e\n", CalculatePhi());
}

// ----------------------------------------------------------------------------

void Particle :: SetToZeroMomentum() {

   energy   = 0.;
   momentum = 0.;
   v        = 0.;
   sv       = 0.;
   sinPhi   = 0.;
   cosPhi   = 1.;
}

// ----------------------------------------------------------------------------

double Particle :: Get_pT() {

   return momentum * 2 * sqrt(v * (1. - v));
}

// ----------------------------------------------------------------------------

void Particle :: EpsilonCartesian(Particle* a, Particle* b, Particle* c) {

   double* afvect = a -> fvect;
   double* bfvect = b -> fvect;
   double* cfvect = c -> fvect;

   fvect[0] =
       afvect[1] * bfvect[2] * cfvect[3]
     + afvect[2] * bfvect[3] * cfvect[1]
     + afvect[3] * bfvect[1] * cfvect[2]
     - afvect[2] * bfvect[1] * cfvect[3]
     - afvect[1] * bfvect[3] * cfvect[2]
     - afvect[3] * bfvect[2] * cfvect[1]
   ;

   fvect[1] =
       afvect[2] * bfvect[3] * cfvect[0]
     + afvect[3] * bfvect[0] * cfvect[2]
     + afvect[0] * bfvect[2] * cfvect[3]
     - afvect[3] * bfvect[2] * cfvect[0]
     - afvect[2] * bfvect[0] * cfvect[3]
     - afvect[0] * bfvect[3] * cfvect[2]
   ;

   fvect[2] =
     - afvect[3] * bfvect[0] * cfvect[1]
     - afvect[0] * bfvect[1] * cfvect[3]
     - afvect[1] * bfvect[3] * cfvect[0]
     + afvect[0] * bfvect[3] * cfvect[1]
     + afvect[3] * bfvect[1] * cfvect[0]
     + afvect[1] * bfvect[0] * cfvect[3]
   ;

   fvect[3] =
       afvect[0] * bfvect[1] * cfvect[2]
     + afvect[1] * bfvect[2] * cfvect[0]
     + afvect[2] * bfvect[0] * cfvect[1]
     - afvect[1] * bfvect[0] * cfvect[2]
     - afvect[0] * bfvect[2] * cfvect[1]
     - afvect[2] * bfvect[1] * cfvect[0]
   ;
}

// ----------------------------------------------------------------------------

double Particle :: CalculateMass() {

   return sqrt(energy * energy - momentum * momentum);
}

// ----------------------------------------------------------------------------

double Particle :: CalculateImaginaryMass() {

   return sqrt(momentum * momentum - energy * energy);
}

// ----------------------------------------------------------------------------

double Particle :: CalculateMassCartesian() {

   double square =   fvect[0] * fvect[0]
                   - fvect[1] * fvect[1]
                   - fvect[2] * fvect[2]
                   - fvect[3] * fvect[3];

   if (square < 0)
      errf(-1, "Particle :: CalculateMassCartesian: negative argument");

   return sqrt(square);
}

// ----------------------------------------------------------------------------

double Particle :: CalculateImaginaryMassCartesian() {

   double square = - fvect[0] * fvect[0]
                   + fvect[1] * fvect[1]
                   + fvect[2] * fvect[2]
                   + fvect[3] * fvect[3];

   if (square < 0)
      errf(-1, "Particle :: CalculateImaginaryMassCartesian: "
               "negative argument");

   return sqrt(square);
}

// ----------------------------------------------------------------------------

void Particle :: CalculateSV() {
   sv = 2. * sqrt( v * (1. - v) );
}

// ----------------------------------------------------------------------------

// calculate components in a new coordinate system
// this == old is possible

void Particle :: NewComponentsCartesian(
                    Particle* old, 
                    Particle* nE,
                    Particle* nX, 
                    Particle* nY, 
                    Particle* nZ) {

   double nfvect[4];

   nfvect[0] =   DotCartesian(nE, old); 
   nfvect[1] = - DotCartesian(nX, old); 
   nfvect[2] = - DotCartesian(nY, old); 
   nfvect[3] = - DotCartesian(nZ, old); 

   fvect[0] = nfvect[0];
   fvect[1] = nfvect[1];
   fvect[2] = nfvect[2];
   fvect[3] = nfvect[3];
}

// ----------------------------------------------------------------------------

double Particle :: Calculate_pT() {

   return momentum * sv;
}

// ----------------------------------------------------------------------------

double Particle :: Calculate_PseudoRapidity() {

   if (v < 1.e-20)
      return 100.0;
   else 
   if (v > 1.-1e-10)
      return - 100.0;
   else
      return - log(2.*v/sv);
}

// ----------------------------------------------------------------------------

int Particle :: UseListIsSame(Particle* other) {

   // :_TAWM_:
   other = other;

   errf(-1, "Particle :: UseListIsSame: not yet implemented");
   return 0;
}

// ----------------------------------------------------------------------------

void Particle :: CalculateNormalFromPolar() {

   nvect[1] = sv * cosPhi;
   nvect[2] = sv * sinPhi;
   nvect[3] = 1. - 2. * v;
}

// ----------------------------------------------------------------------------

void Particle :: CalculateCartesianFromPolar() {

   fvect[0] = energy;
   fvect[1] = momentum * sv * cosPhi;
   fvect[2] = momentum * sv * sinPhi;
   fvect[3] = momentum * (1. - 2. * v);
}

// ----------------------------------------------------------------------------

// directly, not by involing other routines to make it more efficient

void Particle :: CalculateCartesianAndNormalFromPolar() {

   fvect[0] = energy;
   fvect[1] = momentum * (nvect[1] = sv * cosPhi);
   fvect[2] = momentum * (nvect[2] = sv * sinPhi);
   fvect[3] = momentum * (nvect[3] = 1. - 2. * v);
}

// ----------------------------------------------------------------------------

// calculate polar from Cartesian coordinates
// :_MOD_: optimize 1/(momentum * sv)

void Particle :: CalculatePolarFromCartesian() {

   energy = fvect[0];

   momentum = sqrt(  fvect[1] * fvect[1]
                   + fvect[2] * fvect[2]
                   + fvect[3] * fvect[3]);

   if (momentum < 1.e-40) {

      v      = 0.0;
      sv     = 0.0;
      sinPhi = 0.0;
      cosPhi = 0.0;

   } else {

      v = 0.5 * (1. - fvect[3] / momentum);
      sv = 2. * sqrt(v * (1.0 - v));

      if (sv < 1.e-40) {
         sinPhi = 0.;
         cosPhi = 1.;
      } else {
         if (fabs(fvect[1]) < 1.e-40)
            cosPhi = 0.;
         else
            cosPhi = fvect[1] / (momentum * sv);
         if (fabs(fvect[2]) < 1.e-40)
            sinPhi = 0.;
         else
            sinPhi = fvect[2] / (momentum * sv);
      }
   }
}

// ----------------------------------------------------------------------------

double Dot(Particle* p, Particle* q) {

   return   p -> energy * q -> energy
          - p -> momentum * q -> momentum
            * ( (1. - 2 * p -> v) * ( 1. - 2 * q -> v)
               + p -> sv * q -> sv
                 * (  p -> sinPhi * q -> sinPhi 
                    + p -> cosPhi * q -> cosPhi)
              );
}

// ----------------------------------------------------------------------------

double DotCartesian(Particle* p, Particle* q) {

   return     p -> fvect[0] * q -> fvect[0]
         - (  p -> fvect[1] * q -> fvect[1]
            + p -> fvect[2] * q -> fvect[2]
            + p -> fvect[3] * q -> fvect[3]
           );
}

// ----------------------------------------------------------------------------

double DotSpaceCartesian(Particle* p, Particle* q) {

   return   p -> fvect[1] * q -> fvect[1]
          + p -> fvect[2] * q -> fvect[2]
          + p -> fvect[3] * q -> fvect[3]
   ;
}

// ----------------------------------------------------------------------------

double InvariantMassCartesian(Particle *p, Particle *q) {

   double a0 = p -> fvect[0] + q -> fvect[0];
   double a1 = p -> fvect[1] + q -> fvect[1];
   double a2 = p -> fvect[2] + q -> fvect[2];
   double a3 = p -> fvect[3] + q -> fvect[3];

   return a0 * a0 - (a1 * a1 + a2 * a2 + a3 * a3);
}

// ----------------------------------------------------------------------------

double V(Particle* p, Particle* q) {

   return 0.5 * (   
                   2. * (p -> v + q -> v) 
                 - 4. *  p -> v * q -> v
                 - p -> sv * q -> sv
                     * (  p -> sinPhi * q -> sinPhi 
                        + p -> cosPhi * q -> cosPhi)
                );
}

// ----------------------------------------------------------------------------

double VCartesian(Particle* p, Particle* q) {

   double squareP = DotSpaceCartesian(p, p);
   double squareQ = DotSpaceCartesian(q, q);

   if (squareP < 1.e-40 || squareQ < 1.e-40)
      return 1.0;
   else
      return 0.5 * ( 1. - DotSpaceCartesian(p, q) / sqrt(squareP * squareQ));
}

// ----------------------------------------------------------------------------

double VNormal(Particle* p, Particle* q) {

   return 0.5 * (  1. 
                 - (  p -> nvect[1] * q -> nvect[1]
                    + p -> nvect[2] * q -> nvect[2]
                    + p -> nvect[3] * q -> nvect[3]
                   )
                );
}

// ----------------------------------------------------------------------------

void DotAndVMasslessNormal(Particle* p, Particle* q, double* d, double* v) {

   *v = 0.5 * (  1. 
               - (  p -> nvect[1] * q -> nvect[1]
                  + p -> nvect[2] * q -> nvect[2]
                  + p -> nvect[3] * q -> nvect[3]
                 )
              );

   *d =  2. * p -> energy * q -> energy * *v;
}

// ----------------------------------------------------------------------------

void InvariantAndVMasslessNormal(
        Particle* p,
        Particle* q, 
        double* d, 
        double* v) {

   *v = 0.5 * (  1. 
               - (  p -> nvect[1] * q -> nvect[1]
                  + p -> nvect[2] * q -> nvect[2]
                  + p -> nvect[3] * q -> nvect[3]
                 )
              );

   *d =  4. * p -> energy * q -> energy * *v;
}

// ----------------------------------------------------------------------------
//
// --> Sets of invariants: 
//
// ----------------------------------------------------------------------------

Invariant :: Invariant() {
}

// ----------------------------------------------------------------------------

Invariant :: ~Invariant() {
}

// ----------------------------------------------------------------------------

// calculate the invariants of the phase space variables

void Invariant :: CalculatePartonInvariants(Event *e) {
 
   n_in_partons_local = e -> n_in_partons;
   npartons_local     = e -> npartons;

   Particle** mappingL = e -> GetMapping();

   Q2inv = e -> Q2;

   double dummy;

   e -> CalculateNormalFromPolarRangeMap( - n_in_partons_local, 
                                            npartons_local - 1);

   e0 = mappingL[0] -> GetEnergy();
   f0 = e0 * e0;

   if (npartons_local >= 2) {
      e1 = mappingL[1] -> GetEnergy();
      f1 = e1 * e1;
      e -> InvariantAndVMasslessNormalMap(0, 1, &s01, &v01);
   }

   if (npartons_local >= 3) {
      e2 = mappingL[2] -> GetEnergy();
      f2 = e2 * e2;
      e -> InvariantAndVMasslessNormalMap(0, 2, &s02, &v02);
      e -> InvariantAndVMasslessNormalMap(1, 2, &s12, &v12);
   }

   if (npartons_local >= 4) {
      e3 = mappingL[3] -> GetEnergy();
      f3 = e3 * e3;
      v03 = e -> VNormalMap(0, 3);
      v13 = e -> VNormalMap(1, 3);
      v23 = e -> VNormalMap(2, 3);
   }

   if (n_in_partons_local >= 1) {   

      mappingL[-4] -> CalculateNormalFromPolar();
      mappingL[-5] -> CalculateNormalFromPolar();
      ei = mappingL[-1] -> GetEnergy();
      fi = ei * ei;
      e -> InvariantAndVMasslessNormalMap(-1,  0, &si0, &vi0);
      e -> InvariantAndVMasslessNormalMap(-4, -5, &skl, &vkl);
      e -> InvariantAndVMasslessNormalMap(-1, -4, &sik, &vik);
      e -> InvariantAndVMasslessNormalMap(-1, -5, &sil, &vil);
      e -> InvariantAndVMasslessNormalMap( 0, -4, &s0k, &v0k);
      e -> InvariantAndVMasslessNormalMap( 0, -5, &s0l, &v0l);

      if (npartons_local >= 2) {
         e -> InvariantAndVMasslessNormalMap(-1,  1, &si1, &vi1);
         e -> InvariantAndVMasslessNormalMap( 1, -4, &s1k, &v1k);
         e -> InvariantAndVMasslessNormalMap( 1, -5, &s1l, &v1l);
      }

      if (npartons_local >= 3) {
         e -> InvariantAndVMasslessNormalMap(-1,  2, &si2, &vi2);
         e -> InvariantAndVMasslessNormalMap( 2, -4, &s2k, &v2k);
         e -> InvariantAndVMasslessNormalMap( 2, -5, &s2l, &v2l);
      }

      if (npartons_local >= 4) {
         vi3 = e -> VNormalMap(-1, 3);
      }
   }

   LimitType lt = e -> GetLimitType();

   // calculate invariants involving the ``combined cluster''
   if (   lt == COLLINEAR 
       || lt == SOFT_AND_COLLINEAR
       || lt == SOFT_ADDED          // :_MOD_: this is too general!
       || lt == COLLINEAR_ADDED
       || lt == SOFT_AND_COLLINEAR_ADDED
      ) {

      clusterDefined = TRUE;

      mappingL[-3] -> CalculateNormalFromPolar();

      eh = mappingL[-3] -> GetEnergy();
      fh = eh * eh;

      e -> InvariantAndVMasslessNormalMap(0, -3, &s0h, &dummy);

      if (npartons_local >= 2)
         e -> InvariantAndVMasslessNormalMap(1, -3, &s1h, &dummy);

      if (npartons_local >= 3)
         e -> InvariantAndVMasslessNormalMap(2, -3, &s2h, &dummy);

      if (n_in_partons_local >= 1) {   
         e -> InvariantAndVMasslessNormalMap(-1, -3, &sih, &dummy);
         e -> InvariantAndVMasslessNormalMap(-4, -3, &skh, &dummy);
         e -> InvariantAndVMasslessNormalMap(-5, -3, &slh, &dummy);
      }

   } else

      clusterDefined = FALSE;

   // calculate azimuthal subtraction terms
   if (lt == COLLINEAR && npartons_local >= 3) {

      dtTermsDefined = TRUE;

      e -> CalculateWList();

      dt0 = e -> dtlist[0];
      dt1 = e -> dtlist[1];
      dt2 = e -> dtlist[2];

      if (npartons_local >= 4)
         dt3 = e -> dtlist[3];

      if (n_in_partons_local == 1) {   
         dti = e -> dtlist[-1];
         dtk = e -> dtlist[-4];
         dtl = e -> dtlist[-5];
      } else
      if (n_in_partons_local > 1)
         errf(-1, "Invariant :: CalculatePartonInvariants: "
                  "not yet implemented");
   } else

      dtTermsDefined = FALSE;
}

// ----------------------------------------------------------------------------

void Invariant :: CopyPartonInvariantsFrom(Invariant *inv) {

   n_in_partons_local = inv -> n_in_partons_local;
   npartons_local     = inv -> npartons_local;
   clusterDefined     = inv -> clusterDefined;
   dtTermsDefined     = inv -> dtTermsDefined;

   Q2inv = inv -> Q2inv;

   e0 = inv -> e0;
   f0 = inv -> f0;

   if (npartons_local >= 2) {
      e1  = inv -> e1;
      f1  = inv -> f1;
      v01 = inv -> v01;
      s01 = inv -> s01;
   }

   if (npartons_local >= 3) {
      e2  = inv -> e2;
      f2  = inv -> f2;
      v02 = inv -> v02;
      v12 = inv -> v12;
      s02 = inv -> s02;
      s12 = inv -> s12;
   }

   if (npartons_local >= 4) {
      e3  = inv -> e3;
      f3  = inv -> f3;
      v03 = inv -> v03;
      v13 = inv -> v13;
      v23 = inv -> v23;
   }

   if (n_in_partons_local >= 1) {   

      ei  = inv -> ei;
      fi  = inv -> fi;
      vi0 = inv -> vi0;
      vkl = inv -> vkl;
      vik = inv -> vik;
      vil = inv -> vil;
      v0k = inv -> v0k;
      v0l = inv -> v0l;
      si0 = inv -> si0;
      skl = inv -> skl;
      sik = inv -> sik;
      sil = inv -> sil;
      s0k = inv -> s0k;
      s0l = inv -> s0l;

      if (npartons_local >= 2) {
         vi1 = inv -> vi1;
         v1k = inv -> v1k;
         v1l = inv -> v1l;
         si1 = inv -> si1;
         s1k = inv -> s1k;
         s1l = inv -> s1l;
      }

      if (npartons_local >= 3) {
         vi2 = inv -> vi2;
         v2k = inv -> v2k;
         v2l = inv -> v2l;
         si2 = inv -> si2;
         s2k = inv -> s2k;
         s2l = inv -> s2l;
      }

      if (npartons_local >= 4) {
         vi3 = inv -> vi3;
      }
   }

   if (clusterDefined) {

      lambda = inv -> lambda;
      eta    = inv -> eta; 

      eh = inv -> eh;
      fh = inv -> fh;

      s0h = inv -> s0h;

      if (npartons_local >= 2)
         s1h = inv -> s1h;

      if (npartons_local >= 3)
         s2h = inv -> s2h;

      if (n_in_partons_local >= 1) {   
         sih = inv -> sih;
         skh = inv -> skh;
         slh = inv -> slh;
      }
   }

   if (dtTermsDefined) {

      dt0 = inv -> dt0;
      dt1 = inv -> dt1;
      dt2 = inv -> dt2;

      if (npartons_local>=4)
         dt3 = inv -> dt3;

      if (n_in_partons_local == 1) {   
         dti = inv -> dti;
         dtk = inv -> dtk;
         dtl = inv -> dtl;
      }
   }
}

// ----------------------------------------------------------------------------

void Invariant :: Print(FILE* out) {

   Mstream s_out;
   s_out.Associate(out);

   s_out << FS("------------- INVARIANTS ------- %d", n_in_partons_local)
         << FS(" %d \n", npartons_local);

   s_out << FS("e0 =%30.33e\n",e0);
   s_out << FS("e1 =%30.33e\n",e1);
   s_out << FS("e2 =%30.33e\n",e2);
   s_out << FS("eh =%30.33e\n",eh);
   s_out << FS("f0 =%30.33e\n",f0);
   s_out << FS("f1 =%30.33e\n",f1);
   s_out << FS("f2 =%30.33e\n",f2);
   s_out << FS("fh =%30.33e\n",fh);
   s_out << FS("v01=%30.33e\n",v01);
   s_out << FS("v02=%30.33e\n",v02);
   s_out << FS("v12=%30.33e\n",v12);
   s_out << FS("v0k=%30.33e\n",v0k);
   s_out << FS("v1k=%30.33e\n",v1k);
   s_out << FS("v2k=%30.33e\n",v2k);
   s_out << FS("v0l=%30.33e\n",v0l);
   s_out << FS("v1l=%30.33e\n",v1l);
   s_out << FS("v2l=%30.33e\n",v2l);
   s_out << FS("vkl=%30.33e\n",vkl);
   s_out << FS("s01=%30.33e\n",s01);
   s_out << FS("s02=%30.33e\n",s02);
   s_out << FS("s12=%30.33e\n",s12);
   s_out << FS("s0k=%30.33e\n",s0k);
   s_out << FS("s1k=%30.33e\n",s1k);
   s_out << FS("s2k=%30.33e\n",s2k);
   s_out << FS("s0l=%30.33e\n",s0l);
   s_out << FS("s1l=%30.33e\n",s1l);
   s_out << FS("s2l=%30.33e\n",s2l);
   s_out << FS("skl=%30.33e\n",skl);
   s_out << FS("s0h=%30.33e\n",s0h);
   s_out << FS("s1h=%30.33e\n",s1h);
   s_out << FS("s2h=%30.33e\n",s2h);
   s_out << FS("skh=%30.33e\n",skh);
   s_out << FS("slh=%30.33e\n",slh);
   s_out << FS("Q2inv=%30.33e\n",Q2inv);
   s_out << FS("lambda=%30.33e\n",lambda);
   s_out << FS("eta   =%30.33e\n",eta);
   if (n_in_partons_local >= 1) {   
      s_out << FS("ei =%30.33e\n",ei);
      s_out << FS("fi =%30.33e\n",fi);
      s_out << FS("vi0=%30.33e\n",vi0);
      s_out << FS("vi1=%30.33e\n",vi1);
      s_out << FS("vi2=%30.33e\n",vi2);
      s_out << FS("vik=%30.33e\n",vik);
      s_out << FS("vil=%30.33e\n",vil);
      s_out << FS("si0=%30.33e\n",vi0);
      s_out << FS("si1=%30.33e\n",vi1);
      s_out << FS("si2=%30.33e\n",vi2);
      s_out << FS("sik=%30.33e\n",sik);
      s_out << FS("sil=%30.33e\n",sil);
      s_out << FS("sih=%30.33e\n",sih);
   }
   if (npartons_local >= 4) {
      s_out << FS("e3 =%30.33e\n",e3);
      s_out << FS("f3 =%30.33e\n",f3);
      s_out << FS("v03=%30.33e\n",v03);
      s_out << FS("v13=%30.33e\n",v13);
      s_out << FS("v23=%30.33e\n",v23);
      if (n_in_partons_local >= 1) {
         s_out << FS("vi3=%30.33e\n",vi3);
      }
   }
}

// ----------------------------------------------------------------------------
//
// --> Events: 
//     given as a pointer array to a set of `Particle's
//
// ----------------------------------------------------------------------------

Event :: Event() {

   hl = new HorizontalList <Event> (this, &Event :: GetH);   

#if CHECKDESTRUCT
      ident=gl->tEvent -> Mark();
#else
      ident=gl->tEvent -> Unique();
#endif

   global = NULL;
   tEval  = NULL;

   addentries = 9; // number of additional entries with negative `label'
                   // -1, -2: incident partons
   
   maxparticles = -1;

   ecm     = -1.;
   ecmfull = -1.;

//   EventWithAlphaEM = NULL;

   limitType = NO_LIMIT;
   userFrame = none;
}

// ----------------------------------------------------------------------------

Event :: ~Event() {

   delete hl;

#if CHECKDESTRUCT
      gl->tEvent.Remove(ident);
#endif

   register int i;

   if (maxparticles >= 0) {

      for (i = 0; i < maxparticles + addentries; ++i)
          tEval -> GetPFree() -> ReturnObject(particle[i]);

      delete [] particle;
      delete [] (mapping     - addentries);
      delete [] (mappingSave - addentries);
      delete [] vlist;
      delete [] flist;
      delete [] iwlist;
      for (i = -addentries; i < maxparticles; ++i)
          delete [] (wlist[i] - addentries);
      delete [] (wlist  - addentries);
      delete [] (dtlist - addentries);
   }
}

// ----------------------------------------------------------------------------

void Event :: Define(int maxparticlesIn, Global *globalIn) {

   DefineX(maxparticlesIn, globalIn);
   errf(-1, "Event :: Define: direct call...");
}

// ----------------------------------------------------------------------------

void Event :: DefineX(int maxparticlesIn, Global* globalIn) {

   if (IsDefined())
      errf(-1, "Event::DefineX: already defined");

   // fudge here! If tEval is not defined, then use the global object...
   if (tEval == NULL)
      tEval = gl -> disaster;

   register int i;

   global = globalIn;

   disaster = global -> disaster;

   if (maxparticles != -1) {
      global -> To_error()
         << FS("maxparticles=%d", maxparticles)
         << FS(", nparticlesIn=%d\n", maxparticlesIn);
      errf(-1,"Event :: Define: error.");
   }

   maxparticles = maxparticlesIn;

   particle = new Particle* [maxparticles+addentries];
   mapping = (new Particle* [maxparticles+addentries]) + addentries;
   mappingSave = (new Particle* [maxparticles+addentries]) +  addentries;
   vlist = new double [(   maxparticles + addentries) 
                       * (maxparticles + addentries)];
   flist = new double [maxparticles + addentries];

   // index list 
   iwlist = new int [maxparticles + addentries];
   // list of limits
   wlist = (new double* [maxparticles + addentries]) + addentries;
   for (i = -addentries; i < maxparticles; ++i)
       wlist[i] = (new double [maxparticles + addentries] + addentries);
   dtlist = (new double [maxparticles + addentries]) + addentries;

   for (i = 0; i < maxparticles + addentries; ++i) {
       particle[i] = tEval -> GetPFree() -> GetDefinedObject(tEval);
       particle[i] -> SetEvent(this);
       particle[i] -> SetLocation(i);
   }

   Reset();
}

// ----------------------------------------------------------------------------

void Event :: Define(TheoryEvaluation *tIn, int iIn) {

   tEval = tIn;
   DefineX(iIn, tEval);
}

// ----------------------------------------------------------------------------

void Event :: Reset() {

   hl -> Reset();

#if CHECKDESTRUCT
      (void) gl->tEvent.Check(ident, 1);
#endif

   nparticles = 0;
   npartons   = 0;

   SetMappingToNULL();

   id = gl -> idEvent -> Unique(); // set identification

   userFrame = original;

   isConstructed = FALSE;
   eventType     = EventTypeNone;

   contribution = NULL;

   alreadyInUseList = FALSE;
}

// ----------------------------------------------------------------------------

void Event :: Reset(TheoryEvaluation* tIn) {

   tEval = tIn;
   Reset();
}

// ----------------------------------------------------------------------------

void Event :: ReturnToFreeList() {

#if CHECKDESTRUCT
      if (gl->tEvent.Check(ident, 2) == 1) {
         Print(stdout);
      }
#endif
}

// ----------------------------------------------------------------------------

int Event :: GetMaxNumber() {

   return maxparticles;
}

// ----------------------------------------------------------------------------

void Event :: SetPartonNumber(int npartonsIn) {

   if (npartonsIn > maxparticles || nparticles != 0)
      errf(-1,"Event :: SetNumber: error.");
   npartons = npartonsIn;
}

// ----------------------------------------------------------------------------

int Event :: GetPartonNumber() {

   return npartons;
}

// ----------------------------------------------------------------------------

void Event :: SetIn(int inI) {

   n_in_partons = inI;
}

// ----------------------------------------------------------------------------

int Event :: GetIn() {

   return n_in_partons;
}

// ----------------------------------------------------------------------------

Particle** Event :: GetParticle() {

   return particle;
}

// ----------------------------------------------------------------------------

void Event :: Print(Mstream& s_out) {

   int i;

//   for (i = 0; i < nparticles; ++i)
//       particle[i] -> CalculateCartesianAndNormalFromPolar();

   s_out << FS("limitType=%s", LimitTypeName[(int) limitType])
         << FS("  limit1=%d", limit1)
         << FS(" limit2=%d", limit2)
         << FS(",    maxparticles=%d\n", maxparticles);
   s_out << FS("npartons=%d", npartons)
         << FS(", n_in_partons=%d\n", n_in_partons);
   s_out << FS("ident=%ld", GetIdent())
         << FS(" userFrame=%s\n", FrameName[GetUserFrame()]);
   s_out << FS("frame=%s", FrameName[GetFrame()])
         << FS(" eventType=%s\n", EventTypeName[GetEventType()]);
   s_out << FS("toBeSubtracted=%s\n", 
               SubtractionTypeName[(int) GetToBeSubtracted()]
              );
   s_out << FS("xi = %16.6e", xi)
         << FS("  xi_hard   = %16.6e\n", xi_hard);
   s_out << FS("xi_int = %16.6e", xi_int)
         << FS("  xi_pden   = %16.6e\n", xi_pden);
   s_out << FS("u  = %16.6e", u)
         << FS("  ujacobian = %16.6e\n", ujacobian);
   s_out << FS("xB = %16.6e", xB)
         << FS("  y = %16.6e\n", y);
   s_out << FS("p0ref = %s\n", P0_ReferenceName[(int) p0ref]);
   s_out << FS("Event in Normalized Coordinates (nparticles=%d", 
               nparticles
              )
         << FS(",npartons=%d):\n", npartons);
   for (i = 0; i < nparticles; ++i) {
       s_out << FS("%2d ", i)
             << FS("%c ", particle[i] -> GetType())
             << FS("%2d ", particle[i] -> GetLabel())
             << FS("%14.6e", particle[i] -> GetEnergy())
             << FS("%14.6e", particle[i] -> GetMomentum())
             << FS("%14.6e", particle[i] -> GetNvect(1))
             << FS("%14.6e", particle[i] -> GetNvect(2))
             << FS("%14.6e\n", particle[i] -> GetNvect(3));
   }
#if 1
   s_out << "Event in Cartesian Coordinates:\n";
   for (i = 0; i < nparticles; ++i) {
       s_out << FS("%2d ", i)
             << FS("%c ", particle[i]->GetType())
             << FS("%2d ", particle[i]->GetLabel())
             << FS("%16.6e ", particle[i]->GetFvect(0))
             << FS("%16.6e ", particle[i]->GetFvect(1))
             << FS("%16.6e ", particle[i]->GetFvect(2))
             << FS("%16.6e\n", particle[i]->GetFvect(3));
   }
#endif
#if 1
   s_out << "Event in Polar Coordinates:\n";

   s_out << "                 energy         momentum                v"
         << "              phi\n";
   for (i = 0; i < nparticles; ++i) {
       s_out << FS("%2d ", i)
             << FS("%c ", particle[i] -> GetType())
             << FS("%2d ", particle[i] -> GetLabel())
             << FS("%16.6e ", particle[i] -> GetEnergy())
             << FS("%16.6e ", particle[i] -> GetMomentum())
             << FS("%16.6e ", particle[i] -> GetV())
             << FS("%16.6e\n", particle[i] -> CalculatePhi());
   }
#endif
#if 1
   s_out << "sv, sinPhi, cosPhi:\n";

   for (i = 0; i < nparticles; ++i) {
       s_out << FS("%2d ", i)
             << FS("%c ", particle[i] -> GetType())
             << FS("%2d ", particle[i] -> GetLabel())
             << FS("%16.6e ", particle[i] -> GetSV())
             << FS("%16.6e ", particle[i] -> GetSinPhi())
             << FS("%16.6e\n", particle[i] -> GetCosPhi());
   }
#endif
#if 1
   s_out << "Identities\n";
   for (i = 0; i < maxparticles + addentries; ++i) {
       s_out << FS("%2d ", i)
             << FS("%c ", particle[i] -> GetType())
             << FS("%2d ", particle[i] -> GetLabel())
             << FS("%10ld\n", particle[i] -> GetIdent());
   }
#endif
#if 1
   s_out << "mapping\n";
   for (i=-addentries;i<maxparticles;++i) {
       if (mapping[i] != NULL) {
          s_out << FS("%2d", i)
                << FS(" --> %2d", mapping[i] -> GetLocation())
                << FS(" %c", mapping[i] -> GetType())
                << FS(" %2d", mapping[i] -> GetLabel())
                << FS(" %10ld\n", mapping[i] -> GetIdent()); 
       } else {
          s_out << FS("%2d --> NULL\n", i);
          }
   }
#endif
#if 0
   int j;
   if (wlist!=NULL) {
      s_out << FS("collinear limits: (%d)\n", iwentries);
      for (i = 0; i < iwentries; ++i) {
          j = iwlist[i];
          s_out << FS("%2d", j)
                << FS(" %c", mapping[j]->GetType())
                << FS(" %2d", mapping[j]->GetLabel())
                << FS(" %2d", mapping[j]->GetLocation())
                << fs("%16.6e", mapping[j]->cdp)
                << FS(" %16.6e\n",mapping[j]->sq);
      }
      s_out << "collinear limits:\n";
      for (i = 0;i < iwentries - 1; ++i)
          for (j = i+1; j < iwentries; ++j) {
          s_out << FS("%2d", iwlist[i])
                << FS(" %2d", iwlist[j])
                << FS(" %16.6e\n", wlist[iwlist[i]][iwlist[j]]);
      }
   }
#endif
   s_out << "--------------------------------------"
         << "--------------------------------------\n";
   s_out.Flush();
}

// ----------------------------------------------------------------------------

void Event :: Print(FILE* fout) {

   Mstream s_out;
   s_out.Associate(fout);

   Print(s_out);
}

// ----------------------------------------------------------------------------

void Event :: PrintMapping(FILE* fout) {

   Mstream s_out;
   s_out.Associate(fout);

   for (register int i = -addentries; i < maxparticles; ++i) {
       if (mapping[i]==NULL) {
          s_out << "i --> NULL\n";
       } else {
          s_out << FS("i --> %3d", mapping[i] -> GetLabel())
                << FS("  %3d\n",   mapping[i] -> GetLocation()); 
       }
   }
}

// ----------------------------------------------------------------------------

void Event :: CopyPSInformationInto(Event* into) {

   into -> frame = frame;
   into -> p0ref = p0ref;

   into -> limit1    = limit1;
   into -> limit2    = limit2;
   into -> limitType = limitType;

   into -> ecm = ecm;
   into -> ecmfull = ecmfull;
   into -> SH = SH;
   into -> proton_energy = proton_energy;

   into -> xB = xB;
   into -> y  = y;
   into -> xi = xi;
   into -> Q  = Q;
   into -> W  = W;
   into -> Q2 = Q2;
   into -> W2 = W2;

   into -> u = u;
   into -> ujacobian = ujacobian;
   into -> xi_hard = xi_hard;

   into -> xi_int = xi_int;
   into -> xi_pden = xi_pden;

//   into -> EventWithAlphaEM = EventWithAlphaEM;

   into -> n_in_partons = n_in_partons;

   // copy this for the time being. Should in principle be done by Define().
   into -> tEval  = tEval;
   into -> global = global;

   into -> userFrame = userFrame;

   into -> isConstructed = isConstructed;
   into -> toBeSubtracted = toBeSubtracted;
   into -> eventType = eventType;

   into -> contribution = contribution;

   into -> vtp = vtp;
}

// ----------------------------------------------------------------------------

void Event :: CopyInto(Event* into) {

   if (   into -> maxparticles < npartons
       || into -> addentries != addentries) {
      Print(global -> To_status());
      tEval -> To_error()
         << FS("%d", into -> maxparticles)
         << FS(" %d", nparticles)
         << FS(" %d", into -> addentries)
         << FS(" %d\n", addentries);
      errf(-1, "Event :: CopyInto: array too small or mismatch.");
   }

   CHECK_EVENT_RECORD(this)

   into -> Reset();

   Particle* a;
   Particle** here = particle;
   Particle** stop = here + nparticles;
   while (here != stop) {
      a = into -> AddParticle();
      (*here) -> CopyInto(a);
      into -> mapping[a -> GetLabel()] = a;
      ++here;
   }

   into -> npartons   = npartons;
   into -> nparticles = nparticles;

   into -> mapij = mapij;
   into -> imap  = imap;
   into -> jmap  = jmap;

   CopyPSInformationInto(into);

   // set identification (contents may have changed)
   into -> id = gl -> idEvent -> Unique(); 

   CHECK_EVENT_RECORD(into)
}

// ----------------------------------------------------------------------------

void Event :: CopyAdditionalParticlesInto(Event* into) {

   CHECK_EVENT_RECORD(this)
   CHECK_EVENT_RECORD(into)

   Particle* a;
   Particle** map = GetMapping() - addentries;
   Particle** stop = GetMapping(); // this is the address of mapping[0];
   while (map != stop) {
      if (*map != NULL) {
         a = into -> AddParticle();
         (*map) -> CopyInto(a);
         into -> mapping[a -> GetLabel()] = a;
      }
      ++map;
   }

   CHECK_EVENT_RECORD(into)
}

// ----------------------------------------------------------------------------

// recombine the 4-vectors labelled by i and j into i
// according to a particular scheme; 
// rename such that negative labels do not disappear
// rearrange pointers!

void Event :: SchemeCombineCartesianMap(RecombinationType rt, int i, int j) {

   CHECK_EVENT_RECORD(this)

   if (i == j || i < 0 && j < 0)
      errf(-1, "Event :: SchemeCombineCartesianMap: error.");

   switch (rt) {
      case JADE_Scheme:
      case E_Scheme:
         mapping[i] -> AddCartesian(mapping[i], 
                                    mapping[j]);
         break;
      case E0_Scheme:
         mapping[i] -> Add_E0Cartesian(mapping[i], 
                                       mapping[j]);
         break;
      case P_Scheme:
         mapping[i] -> Add_PCartesian(mapping[i], 
                                      mapping[j]);
         break;
      case P0_Scheme:
         errf(-1, "Event :: SchemeCombineCartesianMap: P0 not implemented");
         break; 
      default:
         errf(-1, "Event :: SchemeCombineCartesianMap: rt not known");
   }

   Particle* other = mapping[j];

   if (i >= 0 && j < 0) {

      mapping[i] -> SetLabel(j);
      mapping[j] = mapping[i];
      mapping[i] = NULL;

   } else

      mapping[j] = NULL;

   // now the particle pointed to by ``other'' is obsolete

   --npartons;
   --nparticles;

   particle[other -> GetLocation()] = particle[nparticles];
   particle[nparticles]             = other;

   particle[other -> GetLocation()] -> SetLocation(other -> GetLocation());
   particle[nparticles]             -> SetLocation(nparticles);

   CHECK_EVENT_RECORD(this)
}

// ----------------------------------------------------------------------------

// renames all entries with a label >= 0 in a consecutive order
// :_MOD_: this should be beased on the mapping[] variable (!= NULL)

void Event :: RenamePartons() {

   register int i;
   register int co = 0;

   for (i = 0; i < nparticles; ++i)
       if (particle[i] -> GetLabel() >= 0) {
          if (particle[i] -> GetLabel() > co)
             mapping[particle[i] -> GetLabel()] = NULL;
          mapping[co] = particle[i];
          particle[i] -> SetLabel(co++);
       }

   CHECK_EVENT_RECORD(this)
}

// ----------------------------------------------------------------------------

void Event :: Divide(Event* a, Event* b) {
   register int i;
   for (i=0;i<nparticles;++i) {
       particle[i]->Divide(a->particle[i],b->particle[i]);
   }
}

// ----------------------------------------------------------------------------

void Event :: Permute(int value1,...) {

   va_list ap;
   int nread, value, value2;
   Particle* temp;

   va_start(ap, value1);
   nread = 1;
   value = value1;
   while (value != -1) {
      ++nread;
      value = va_arg(ap, int);
      if (nread == 2) {
         value2 = value;
         nread = 0;
         if (value == -1 || value1 >= nparticles || value2 >= nparticles)
            errf(-1, "Event :: Permute: illegal permutation.");      
         temp = particle[value2];
         particle[value2] = particle[value1];
         particle[value1] = temp;
         particle[value1] -> SetLocation(value1);
         particle[value2] -> SetLocation(value2);
      } else {
         value1 = value; 
      }
   }
   va_end(ap);
}

// ----------------------------------------------------------------------------

// permute two entries into the 0th and first position, respectively.
// take care of special cases!
// assume i2 and i3 are distinct

void Event :: PermuteTwoTo01(int i2, int i3) {
   if (i3 == 0)
      if (i2 == 1)
         Permute(0, 1, -1);
      else
         Permute(0, 1, i2, 0, -1);
   else
      Permute(i2, 0, i3, 1, -1);
}

// ----------------------------------------------------------------------------

// permute one entry into the 0th position.

void Event :: PermuteOneTo0(int i1) {
   if (i1 != 0)
      Permute(i1, 0, -1);
}

// ----------------------------------------------------------------------------

// undo PermuteTwoTo01

void Event :: PermuteBackTwoTo01(int i2, int i3) {
   if (i3 == 0)
      if (i2 == 1)
         Permute(0, 1, -1);
      else
         Permute(i2, 0, 0, 1, -1);
   else
      Permute(i3, 1, i2, 0, -1);
}

// ----------------------------------------------------------------------------

// creates a mapping of labels to pointers --> OBSOLETE

void Event :: CreateMapping() {

   CHECK_EVENT_RECORD(this)

   return;
}

// ----------------------------------------------------------------------------

// Sets the array ``mapping'' to NULL.

void Event :: SetMappingToNULL() {

   register int i;

   for (i = -addentries; i < maxparticles; ++i)
       mapping[i] = NULL;
}

// ----------------------------------------------------------------------------

// checks the consistency of an event record

void Event :: CER() {

//   tEval -> To_status() << "CHECKING...\n";

   register int i, j;

   for (i = 0; i < maxparticles+addentries; ++i) {
       if (particle[i] == NULL) {
          tEval -> To_error() << FS("i=%d\n", i);
          errf(-1, "Event::CER: particle is NULL");
       }
       if (particle[i] -> GetEvent() != this) {
          tEval -> To_error() << FS("i=%d\n", i);
          errf(-1, "Event::CER: particle is not in event");
       }  
       if (particle[i] -> GetLocation() != i) {
          tEval -> To_error()
             << FS("i=%d", i)
             << FS(" location=%d\n", particle[i] -> GetLocation());
          errf(-1, "Event::CER: particle has wrong location");
       }  
       for (j = 0; j < i; ++j) {
           if (particle[i] == particle[j]) {
              tEval -> To_error()
                 << FS("i=%d", i)
                 << FS(" j=%d\n", j);
              errf(-1, "Event::CER: particles identical");
           }
       }
       if (i < nparticles) {
          j = -addentries;
          while (j < maxparticles && mapping[j] != particle[i])
             ++j;
          if (j == maxparticles) { 
             Print();
             tEval -> To_error() 
                << FS("i=%d", i)
                << FS(" j=%d\n", j);
             errf(-1, "Event::CER: mapping not found");
          }
       }
   }

   int npco = 0;

   for (i = -addentries; i < maxparticles; ++i) {
       if (mapping[i] != NULL) {
          // search for particle
          j = 0;
          while (j < maxparticles+addentries && mapping[i] != particle[j])
                ++j;
          if (j == maxparticles+addentries) { 
             tEval -> To_error()
                << FS("i=%d", i)
                << FS(" j=%d\n", j);
             errf(-1, "Event::CER: particle not found");
          }
          if (j >= nparticles) { 
             Print(stdout);
             tEval -> To_error()
                << FS("i=%d", i)
                << FS(" j=%d\n", j);
             errf(-1, "Event::CER: particle found in wrong range");
          }
          if (particle[j] -> GetLabel() != i) { 
             Print();
             tEval -> To_error()
                << FS("i=%d", i)
                << FS(" j=%d", j)
                << FS(" label=%d\n", particle[j] -> GetLabel());
             errf(-1, "Event::CER: particle has wrong label");
          }
          if (i >= 0)
             ++npco;
       }
   }

   if (npco != npartons) { 
      Print();
      tEval -> To_error()
         << FS("npco=%d", npco)
         << FS(" npartons=%d\n", npartons);
      errf(-1, "Event::CER: npartons is inconsistent");
   }
   
}

// ----------------------------------------------------------------------------

// transfers particle labels
// works correctly only of we have only partons in the event record

void Event :: TransferPartonLabelsInto(Event* e) {

   CHECK_EVENT_RECORD(this)
   CHECK_EVENT_RECORD(e)

   register int i;
   for (i = 0; i < npartons; ++i) {
       particle[i] -> TransferLabelsInto(e -> particle[i]);
   }

   // have to deal with indices >= 0 only
   for (i = 0; i < maxparticles; ++i)
       e -> mapping[i] = NULL;
   for (i = 0; i < npartons; ++i)
       e -> mapping[e -> particle[i] -> GetLabel()] = e -> particle[i];
   
   CHECK_EVENT_RECORD(e)
}

// ----------------------------------------------------------------------------

// modifies particle labels
// *** ineffective, for test purposes only ***

void Event :: Rename(int oldlabel, int newlabel) {

   register int i;
   for (i = 0; i < nparticles; ++i)
       if (particle[i] -> GetLabel() == oldlabel)
          particle[i]  -> SetLabel(newlabel);
}

// ----------------------------------------------------------------------------

// swaps particle labels and types, 
// leaves momenta unmodified

void Event :: SwapLabelAndType(int label1, int label2) {

   Particle* p1 = mapping[label1];
   Particle* p2 = mapping[label2];

   if (p1 == NULL || p2 == NULL) {
      Print();
      tEval -> To_error()
         << FS("label1 = %d", label1)
         << FS(", label2 = %d\n", label2);
      errf(-1, "Event :: SwapLabelAndType: label1 or label2 out of range.");
   }

   p1 -> SetLabel(label2);
   p2 -> SetLabel(label1);

   char typeT;
   typeT = p1 -> GetType();
   p1 -> SetType(p2 -> GetType());
   p2 -> SetType(typeT);

   mapping[label1] = p2;
   mapping[label2] = p1;
}

// ----------------------------------------------------------------------------

// swaps two particle labels and types to positions 0 and 1

void Event :: SwapLabelAndTypeTwoTo01(int label1, int label2) {

   SwapLabelAndType(particle[0] -> GetLabel(), label1);
   SwapLabelAndType(particle[1] -> GetLabel(), label2);
   limit1 = label1;
   limit2 = label2;
}

// ----------------------------------------------------------------------------

// swaps one particle label and type to position 0

void Event :: SwapLabelAndTypeOneTo0(int label1) {

   SwapLabelAndType(particle[0] -> GetLabel(), label1);
   limit1 = label1;
}

// ----------------------------------------------------------------------------

// note: imin may be incident!

void Event :: SmallestInvariantCartesian(
                 ClusterAlgorithmType cat,
                 RecombinationType rt,
                 double &minscale2,
                 int &imin, 
                 int &jmin
              ) {

   register int i;
   register int j;
   Particle* pi;
   Particle* pj;
   int defd;
   double scale2;

   // :_TAWM_:
   imin = -99;
   jmin = -99;

   // find smallest invariant

   defd = FALSE;

   for (i = - n_in_partons; i < npartons - 1; ++i)
       for (j = max(0, i + 1); j < npartons; ++j) {
           pi = mapping[i];
           pj = mapping[j];
           switch (cat) {
              case JADE_algorithm:
                 switch (rt) {
                    case JADE_Scheme:
                       scale2 = 4 * pi -> GetFvect(0)
                                  * pj -> GetFvect(0)
                                  * VCartesian(pi, pj);
                       break;
                    case E_Scheme:
                    case E0_Scheme:
                    case P_Scheme:
                    case P0_Scheme:
                       scale2 = InvariantMassCartesian(pi, pj);
                       break;
                    default:
                       errf(-1, "Event :: SmallestInvariant: rt not known");
                       scale2=0; //TAWM
                 }
                 break;
              case kT_algorithm:
                 if (i == -1) {
                    scale2 = 4 * pj -> GetFvect(0)
                               * pj -> GetFvect(0)
                               * VCartesian(pi, pj);
                 } else {
                    scale2 = 4 * dmin(  pi -> GetFvect(0)
                                      * pi -> GetFvect(0), 
                                        pj -> GetFvect(0)
                                      * pj -> GetFvect(0)) 
                               * VCartesian(pi, pj);
                 }
                 break;
              default:
                 errf(-1, "Event :: SmallestInvariantCartesian: "
                          "cat not known");
                 scale2 = 0; // :_TAWM_:
           }
           if (!defd || minscale2 > scale2) {
              defd = TRUE;
              minscale2 = scale2;
              imin = i;
              jmin = j;
           }
       }           
}

// ----------------------------------------------------------------------------

void Event :: ClusterGeneral( 
                 ClusterAlgorithmType cat,
                 RecombinationType rt,
                 double scale2,
                 int &startCluster,         
                 int &finalCluster
              ) {
   
   int imin;
   int jmin;
   double mininvariant;

   // rename remnant(s) for the time being (to -1 (and -2))
   if (n_in_partons >= 1) {
      SwapLabelAndType(-1, -8);
      if (n_in_partons >= 2) {
         errf(-1, "Event :: ClusterGeneral: not yet implemented");
      }
      CHECK_EVENT_RECORD(this);
   }

   CalculateCartesianFromPolarRangeMap(-n_in_partons, npartons-1);

   startCluster = npartons + n_in_partons;
   int nclusters = startCluster;

   int ncombined = 0;

   int changed = TRUE;

   while (changed && nclusters >= 2 ) {

      // find smallest invariant and scale
      SmallestInvariantCartesian(cat, rt, mininvariant, imin, jmin);

      // combine particles?
      if (mininvariant < scale2) {

         SchemeCombineCartesianMap(rt, imin, jmin);
         RenamePartons();
         ++ncombined;
         --nclusters;

      } else

         changed=FALSE;
   }

   CalculatePolarFromCartesianRangeMap(-n_in_partons, npartons-1);

   // re-rename remnant(s) (to -8 (and -?))
   if (n_in_partons >= 1) {
      SwapLabelAndType(-1, -8);
      if (n_in_partons >= 2) {
      }
      CHECK_EVENT_RECORD(this);
   }

   finalCluster = startCluster - ncombined;
}

// ----------------------------------------------------------------------------

// ``scale2'' and ``nCluster'' contain the information on when an 
// event switches the jet multiplicity.
// the variables of the event clustered down to 2 jets are available
// only in Cartesian format...
// Moreover, the remnant and incident particles are still exchanged.

void Event :: ClusterGeneralHistory(
                 ClusterAlgorithmType cat,
                 RecombinationType rt,
                 Array_double *scale2,
                 Array_int *nCluster,
                 int &nEntries
              ) {

   int imin;
   int jmin;
   double mininvariant;

   // rename remnant(s) for the time being (to -1 (and -2))
   if (n_in_partons >= 1) {
      SwapLabelAndType(-1, -8);
      if (n_in_partons >= 2) {
         errf(-1, "Event :: ClusterGeneralHistory: not yet implemented");
      }
      CHECK_EVENT_RECORD(this)
   }

   CalculateCartesianFromPolarRangeMap(-n_in_partons, npartons-1);

   int nclusters = npartons + n_in_partons;
   
   nEntries = 0;
   nCluster -> SetData(0, nclusters);

   while (nclusters >= 2 ) {

      // find smallest invariant and scale
      SmallestInvariantCartesian(cat, rt, mininvariant, imin, jmin);
      scale2 -> SetData(nEntries, mininvariant);

      // Combine clusters
      SchemeCombineCartesianMap(rt, imin, jmin);
      RenamePartons();
      --nclusters;
      nCluster -> SetData(++nEntries, nclusters);
   }
}

// ----------------------------------------------------------------------------

void Event :: ClusterGeneralOneStep(
                 ClusterAlgorithmType cat,
                 RecombinationType rt
              ) {
   
   int imin;
   int jmin;
   double mininvariant;

   // rename remnant(s) for the time being (to -1 (and -2 ?))
   if (n_in_partons >= 1) {
      SwapLabelAndType(-1, -8);
      if (n_in_partons >= 2) {
         errf(-1, "ClusterGeneralOneStep: not yet implemented");
      }
      CHECK_EVENT_RECORD(this)
   }

   CalculateCartesianFromPolarRangeMap(-n_in_partons, npartons-1);

   // find smallest invariant and scale
   SmallestInvariantCartesian(cat, rt, mininvariant, imin, jmin);

//   tEval -> T-Status() << FS("smallest invariants are %d", imin)
//                       << FS(" %d\n", jmin);

   // combine particles

   SchemeCombineCartesianMap(rt, imin, jmin);                
   RenamePartons();      
      
   CalculatePolarFromCartesianRangeMap(-n_in_partons, npartons-1);

   // re-rename remnant(s) (to -8 (and -2 ?))
   if (n_in_partons >= 1) {
      SwapLabelAndType(-1, -8);
      if (n_in_partons >= 2) {
      }
      CHECK_EVENT_RECORD(this)
   }
}

// ----------------------------------------------------------------------------

void Event :: ClusterGeneralWithEvent(
                 Event *copy,
                 ClusterAlgorithmType cat,
                 RecombinationType rt,
                 double scale2,
                 int &startCluster,         
                 int &finalCluster
              ) {
   
   CopyInto(copy);
   copy -> ClusterGeneral(cat, rt, scale2, startCluster, finalCluster);
}

// ----------------------------------------------------------------------------

void Event :: ClusterGeneralHistoryWithEvent(
                 Event *copy,
                 ClusterAlgorithmType cat,
                 RecombinationType rt,
                 Array_double *scale2,
                 Array_int *nCluster,
                 int &nEntries
              ) {

   CopyInto(copy);
   copy -> ClusterGeneralHistory(cat, rt, scale2, nCluster, nEntries);
}

// ----------------------------------------------------------------------------

// clusters according to a Lorentz-invariant algorithm, old version
// (i.e. routine written in FORTRAN),
// returns the number of removed clusters;
// works on final state particles
// *** inefficient, for test purposes only; uses old FORTRAN routine ***
// assumes ``fvect''s to be defined
// does _not_ destroy the original record

int Event :: ClusterInvariantOld(
                       ClusterAlgorithm cluster_algorithm,
                       double scale, int &startCluster,
                       int &etaCut, int iflags
                    ) {

   register int i;

   if (cluster_algorithm!=LORENTZ_INVARIANT_OLD_1)
      errf(-1,"Event :: ClusterInvariantOld: algorithm not known");

   startCluster = npartons + 1;

   CalculateCartesianFromPolarRangeMap(-1, npartons - 1);
   mapping[-4] -> CalculateCartesianFromPolar();
   mapping[-7] -> CalculateCartesianFromPolar();

   // copy momenta into a linear array
   double parray[200];
   parray[0] = mapping[-1] -> GetFvect(1);
   parray[1] = mapping[-1] -> GetFvect(2);
   parray[2] = mapping[-1] -> GetFvect(3);
   parray[3] = mapping[-1] -> GetFvect(0);
   int store = 4;
   for (i = 0; i < npartons; ++i) {
       parray[store++] = mapping[i] -> GetFvect(1);
       parray[store++] = mapping[i] -> GetFvect(2);
       parray[store++] = mapping[i] -> GetFvect(3);
       parray[store++] = mapping[i] -> GetFvect(0);
   }
   parray[store++] = mapping[-4] -> GetFvect(1);
   parray[store++] = mapping[-4] -> GetFvect(2);
   parray[store++] = mapping[-4] -> GetFvect(3);
   parray[store++] = mapping[-4] -> GetFvect(0);

   // and call the clustering routine
   int njetsout;
   ciarrayC(
        parray, npartons+1, scale*scale,
        njetsout, mapping[-1] -> GetFvect(0) / mapping[-7] -> GetFvect(0), 
        etaCut, iflags
   ); 

   return startCluster - njetsout;
}

// ----------------------------------------------------------------------------

// creates list of the v_ij; assumes that polar and
// normal variables are defined.

void Event :: MakeVList() {

   register int k, i, j;

   // serves also for error detection (npartons)
   // ``mapij'' is stored here and used later
   mapij  = global  
            -> pMapSquareLinear 
            -> GetSquareToLinearSafe(npartons); 

   imap   = global 
            -> pMapSquareLinear 
            -> GetSquareToLinearProjection1Unsafe(npartons); 
   jmap   = global 
            -> pMapSquareLinear 
            -> GetSquareToLinearProjection2Unsafe(npartons); 
   
   nvlist = global 
            -> pMapSquareLinear 
            -> GetSquareToLinearLengthUnsafe(npartons);
   
   // angles including also the initial state partons,
   // but not, for two incident partons, their inner product

   CalculateNormalFromPolarRangeMap(- n_in_partons, npartons - 1);

   for (k = 0; k < nvlist; ++k) {
       i = imap[k];
       j = jmap[k];
       vlist[k] = VNormalMap(i, j);
   }
}

// ----------------------------------------------------------------------------

// creates list of the energies

void Event :: MakeEList() {

   register int k;

   // energies only for the final state partons   
   double normalizeScale = 1. / (0.5*ecm);
   for (k = 0; k < npartons; ++k) {
       flist[k] = mapping[k] -> GetEnergy() * normalizeScale;
       flist[k] *= flist[k];
   }
}

// ----------------------------------------------------------------------------

void Event::PrintVELists(FILE *file) {

   Mstream s_out;
   s_out.Associate(file);

   int k;
   s_out << "----- PrintVELists -----\n";
   for (k=0;k<npartons;++k)
       s_out << FS("energy^2: %d", k)
             << FS(" -> %16.6e\n", flist[k]);
   for (k=0;k<nvlist;++k)
       s_out << FS("angle: %d", k)
             << FS(" (%d", imap[k])
             << FS(",%d", jmap[k])
             << FS(") -> %16.6e\n", vlist[k]);
   s_out << "----- End of Print -----\n";
}

// ----------------------------------------------------------------------------

// perform sum over partial fractions for angular variables; 
// assume MakeVList has already been called,
// keep v_(ifixed, jfixed) fixed.
 
double Event :: SumPartialV(int ifixed, int jfixed) {

   return global
          -> pPermute
          -> SumOverPermutations(
                vlist,
                nvlist,
                GetMapij()[ifixed][jfixed]
             );
}

// ----------------------------------------------------------------------------

// perform sum over partial fractions for energies; 
// assume MakeEList has already been called,
// keep E_ifixed fixed.
 
double Event :: SumPartialE(int ifixed) {

   return global -> pPermute -> SumOverPermutations(flist, npartons, ifixed);
}

// ----------------------------------------------------------------------------

// performs the sum over all subtraction terms
// assumes all types of phase space variables to be defined
// (i.e. Cartesian, polar and normal)

void Event :: DoSumOverSubtractionTerms(
                 Process* process, 
                 int ic,
                 double weight_factor,
                 ContributionArray* ca
              ) {

   int inloc = n_in_partons;

   int ienergy, iangle, jangle, irelative;
   int flag;
   double energy, v, jacrel, pffactor, sfactor;

   Event *full = tEval 
                 -> GetEFree()
                 -> GetDefinedObject(tEval, GetMaxNumber());

   // the full term
   global -> lstore = NO_LIMIT;
   CopyInto(full);
   if (GetFrame() == hCMS)
      full -> EcmAddRemnantAndIncident();

   // check the technical cut
   int tcpass = full -> CheckTechnicalCut(
                           process, 
                           "Event :: DoSumOverSubtractionTerms"
                        );
   if (tcpass) {

      Event* limit;

      process -> GetEventList() -> InsertObject(full);
      global -> user -> PhaseSpace(
                                   process, full,
                                   flag,
                                   ca -> NextContribution());
      process -> Evaluate(flag, ic, NO_LIMIT, -1, -2, full,
                          weight_factor,
                          ca -> GetActual());
   
      // sum over all related terms (T1)
      global -> tstore = 1;
      for (ienergy = 0; ienergy < npartons; ++ienergy)
          for (jangle = -inloc; jangle < npartons; ++jangle)
   
              if (ienergy != jangle) {
   
                 global -> istore = ienergy;
                 global -> jstore = jangle;
   
                 // soft limit
                 global -> lstore = SOFT;
                 limit = tEval 
                         -> GetEFree()
                         -> GetDefinedObject(tEval, GetMaxNumber());
                 Limit(SOFT, ienergy, jangle, limit,
                       &jacrel, &energy, &v);
                 if (GetFrame() == hCMS)
                    limit -> EcmAddRemnantAndIncident();
                 process -> GetEventList() -> InsertObject(limit);
                 global -> user -> PhaseSpace(process, limit,
                                              flag,
                                              ca -> NextContribution());
                 sfactor = 1./(energy*energy);
                 limit -> MakeVList();
                 pffactor = limit -> SumPartialV(ienergy, jangle);
                 process -> Evaluate(flag, ic, SOFT, ienergy, jangle, limit,
                                     -sfactor*jacrel*pffactor*weight_factor,
                                     ca -> GetActual());
   
                 // collinear limit
                 global -> lstore = COLLINEAR;
                 limit = tEval 
                         -> GetEFree()
                         -> GetDefinedObject(tEval, GetMaxNumber());
                 Limit(COLLINEAR, ienergy, jangle, limit,
                       &jacrel, &energy, &v);
                 if (GetFrame() == hCMS)
                    limit -> EcmAddRemnantAndIncident();
                 process -> GetEventList() -> InsertObject(limit);
                 global -> user -> PhaseSpace(
                                              process, limit,
                                              flag,
                                              ca -> NextContribution());
                 sfactor = 1./v;
                 limit -> MakeEList();
                 pffactor = limit -> SumPartialE(ienergy);
                 process -> Evaluate(flag, ic, COLLINEAR, 
                                     ienergy, jangle, limit,
                                     -sfactor*jacrel*pffactor*weight_factor,
                                     ca -> GetActual());
   
                 // soft and collinear limit
                 global -> lstore = SOFT_AND_COLLINEAR;
                 limit = tEval 
                         -> GetEFree()
                         -> GetDefinedObject(tEval, GetMaxNumber());
                 Limit(SOFT_AND_COLLINEAR, ienergy, jangle, limit,
                       &jacrel, &energy, &v);
                 if (GetFrame() == hCMS)
                    limit -> EcmAddRemnantAndIncident();
                 process -> GetEventList() -> InsertObject(limit);
                 global -> user -> PhaseSpace(process,
                                              limit,
                                              flag,
                                              ca -> NextContribution());
                 sfactor = 1./(energy*energy*v);
                 pffactor = 1.;
                 process -> Evaluate(flag, ic, SOFT_AND_COLLINEAR,
                                     ienergy, jangle, limit,
                                     -(-sfactor*jacrel*pffactor)
                                      *weight_factor,
                                     ca -> GetActual());
              }
   
      // sum over all unrelated terms, energy (T2)
      global -> tstore = 2;
      for (ienergy = 0; ienergy < npartons; ++ienergy) {
   
          global -> istore = ienergy;
   
          irelative = (ienergy == 0) ? 1 : 0;
             
          // soft limit
          global -> lstore = SOFT;
          limit = tEval 
                  -> GetEFree()
                  -> GetDefinedObject(tEval, GetMaxNumber());
          Limit(SOFT, ienergy, irelative, limit,
                &jacrel, &energy, &v);
          if (GetFrame() == hCMS)
             limit -> EcmAddRemnantAndIncident();
          process -> GetEventList() -> InsertObject(limit);
          global -> user -> PhaseSpace(
                                       process, limit,
                                       flag,
                                       ca -> NextContribution());
          sfactor = 1./(energy*energy);
          limit -> MakeVList();
          pffactor = 0.;
          for (jangle = -inloc; jangle < npartons-1; ++jangle)
              for (iangle = jangle+1; iangle < npartons; ++iangle)
                  if (   ienergy != iangle && ienergy != jangle
                      && !(iangle < 0 && jangle < 0)) {
   
                     pffactor += limit -> SumPartialV(iangle, jangle);
                  }
          process -> Evaluate(flag, ic, SOFT, ienergy, irelative, limit,
                              -sfactor*jacrel*pffactor*weight_factor,
                              ca -> GetActual());
      }
   
      // sum over all unrelated terms, angle (T3)
      global -> tstore=3;
      for (jangle = -inloc; jangle < npartons-1; ++jangle)
          for (iangle = jangle+1; iangle < npartons; ++iangle) {
   
              if (!(iangle < 0 && jangle < 0)) {
       
                 global -> istore = iangle;
                 global -> jstore = jangle;
   
                 // collinear limit
                 global -> lstore = COLLINEAR;
                 limit = tEval 
                         -> GetEFree()
                         -> GetDefinedObject(tEval, GetMaxNumber());
                 Limit(COLLINEAR, iangle, jangle, limit,
                       &jacrel, &energy, &v);
                 if (GetFrame() == hCMS)
                    limit -> EcmAddRemnantAndIncident();
                 process -> GetEventList() -> InsertObject(limit);
                 global -> user -> PhaseSpace(
                                              process, limit,
                                              flag,
                                              ca -> NextContribution());
                 sfactor = 1./v;
                 limit -> MakeEList();
                 pffactor = 0;
                 for (ienergy = 0; ienergy < npartons; ++ienergy)
                     if (ienergy != iangle && ienergy != jangle) {
                        pffactor += limit -> SumPartialE(ienergy);
                     }
                 process -> Evaluate(flag, ic, COLLINEAR, iangle, jangle, limit,
                                     -sfactor*jacrel*pffactor
                                      *weight_factor,
                                      ca -> GetActual());
              }
          }

   } else {  // of if (tcpass)

      tEval -> GetEFree() -> ReturnObject(full);

   }
}

// ----------------------------------------------------------------------------

// performs the sum over all differences of the full matrix element
// and subtraction terms
// assumes all types of phase space variables to be defined
// (i.e. Cartesian, polar and normal)

void Event :: DoSumOverDiffSubtractionTerms(
                 Process* process, 
                 int ic,
                 double weight_factor,
                 ContributionArray* ca
              ) {

   int vmin, vmax, vmaxp;
   if (p0ref == FS_REFERENCE) {
      vmin = 0;
      vmax = npartons - 1;
      vmaxp = npartons;
   } else {
      vmin = - n_in_partons;
      vmax = 0;
      vmaxp = 0;
   }

   int ienergy, iangle, jangle, irelative;
   double pffactorlimit, pffactorfull, sfactor;
   double efactorfull, vfactorfull;

   double energyS,  vS,  jacrelS;
   double energyC,  vC,  jacrelC;
   double energySC, vSC, jacrelSC;

   int iF, iS, iC, iSC;
   int flagF, flagS, flagC, flagSC;

   double sumPartialEF, sumPartialEC;
   double sumPartialVF, sumPartialVS;

   // full phase space
   Event *fullmaster = tEval -> GetEFree() 
                       -> GetDefinedObject(tEval, GetMaxNumber());
   CopyInto(fullmaster);
   if (GetFrame() == hCMS)
      fullmaster -> EcmAddRemnantAndIncident();

   // check the technical cut
   int tcpass = fullmaster -> CheckTechnicalCut(
                                 process, 
                                 "Event :: DoSumOverDiffSubtractionTerms"
                              );
   if (tcpass) {

      Event *full       = tEval -> GetEFree() 
                          -> GetDefinedObject(tEval, GetMaxNumber());
      Event *limit      = tEval -> GetEFree() 
                          -> GetDefinedObject(tEval, GetMaxNumber());
      Event *limitS     = tEval -> GetEFree() 
                          -> GetDefinedObject(tEval, GetMaxNumber());
      Event *limitC     = tEval -> GetEFree() 
                          -> GetDefinedObject(tEval, GetMaxNumber());
      Event *limitSC    = tEval -> GetEFree() 
                          -> GetDefinedObject(tEval, GetMaxNumber());

      process -> GetEventList() -> InsertObject(fullmaster);
      global -> user -> PhaseSpace(
                                   process, fullmaster,
                                   flagF,
                                   ca -> NextContribution());
      iF = ca -> GetActual();
      fullmaster -> MakeVList();
      fullmaster -> MakeEList();
      sumPartialEF = fullmaster 
                     -> SumPartialE(fullmaster -> particle[0] -> GetLabel());
      if (p0ref == FS_REFERENCE)
         sumPartialVF = fullmaster -> SumPartialV(
                           fullmaster -> particle[0] -> GetLabel(),
                           fullmaster -> particle[1] -> GetLabel()
                        );
      else
         sumPartialVF = fullmaster -> SumPartialV(
                           fullmaster -> particle[0] -> GetLabel(),
                           -1
                        );

      // soft limit
      Limit1(SOFT, p0ref, limitS, &jacrelS, &energyS, &vS);
      if (GetFrame() == hCMS)
         limitS -> EcmAddRemnantAndIncident();
      process -> GetEventList() -> InsertObject(limitS);
      global -> user -> PhaseSpace(
                           process, 
                           limitS,
                           flagS,
                           ca -> NextContribution()
                        );
      iS = ca -> GetActual();
      limitS -> MakeVList();
      if (p0ref == FS_REFERENCE)
         sumPartialVS = limitS 
                        -> SumPartialV(
                              limitS -> particle[0] -> GetLabel(),
                              limitS -> particle[1] -> GetLabel()
                           );
      else
         sumPartialVS = limitS 
                        -> SumPartialV(
                              limitS -> particle[0] -> GetLabel(),
                              -1
                           );
   
      // collinear limit
      Limit1(COLLINEAR, p0ref, limitC, &jacrelC, &energyC, &vC);
      if (GetFrame() == hCMS)
         limitC -> EcmAddRemnantAndIncident();
      process -> GetEventList() -> InsertObject(limitC);
      global -> user -> PhaseSpace(
                           process, 
                           limitC,
                           flagC,
                           ca -> NextContribution()
                        );
      iC = ca -> GetActual();
      limitC -> MakeEList();
      sumPartialEC = limitC -> SumPartialE(limitC -> particle[0] -> GetLabel());
   
      // soft and collinear limit
      Limit1(SOFT_AND_COLLINEAR, p0ref, limitSC, &jacrelSC, &energySC, &vSC);
      if (GetFrame() == hCMS)
         limitSC -> EcmAddRemnantAndIncident();
      process -> GetEventList() -> InsertObject(limitSC);
      global -> user -> PhaseSpace(
                           process, 
                           limitSC,
                           flagSC,
                           ca -> NextContribution()
                        );
      iSC = ca -> GetActual();
   
      // sum over all related terms (T1)
      global -> tstore = 1;
      for (ienergy = 0; ienergy < npartons; ++ienergy)
          for (jangle = vmin; jangle < vmaxp; ++jangle) {
              if (ienergy != jangle) {
   
                 global -> istore = ienergy;
                 global -> jstore = jangle;
   
                 // full term
   
                 global -> lstore = NO_LIMIT;
                 fullmaster -> CopyInto(full);
   
                 if (p0ref == FS_REFERENCE)
                    full -> SwapLabelAndTypeTwoTo01(ienergy, jangle);
                 else
                    full -> SwapLabelAndTypeOneTo0(ienergy);
                 
                 pffactorfull = sumPartialEF * sumPartialVF;
                 process -> Evaluate(flagF, ic, NO_LIMIT, -1, -2, full,
                                     pffactorfull*weight_factor,
                                     iF);

                 // soft limit
   
                 global -> lstore = SOFT;
                 limitS -> CopyInto(limit);
   
                 if (p0ref == FS_REFERENCE)
                    limit -> SwapLabelAndTypeTwoTo01(ienergy, jangle);
                 else
                    limit -> SwapLabelAndTypeOneTo0(ienergy);
   
                 pffactorlimit = sumPartialVS;
                 sfactor = 1./(energyS*energyS);
                 process -> Evaluate(flagS, ic,
                                     SOFT, ienergy, jangle, limit,
                                     -sfactor*jacrelS*pffactorlimit
                                      *weight_factor,
                                     iS);
   
                 // collinear limit
   
                 global -> lstore = COLLINEAR;
                 limitC -> CopyInto(limit);
   
                 if (p0ref == FS_REFERENCE)
                    limit -> SwapLabelAndTypeTwoTo01(ienergy, jangle);
                 else
                    limit -> SwapLabelAndTypeOneTo0(ienergy);
   
                 pffactorlimit = sumPartialEC;
                 sfactor = 1./vC;
                 process -> Evaluate(flagC, ic, COLLINEAR,
                                     ienergy, jangle, limit,
                                     -sfactor*jacrelC*pffactorlimit
                                      *weight_factor,
                                     iC);

                 // soft and collinear limit
   
                 global -> lstore = SOFT_AND_COLLINEAR;
                 limitSC -> CopyInto(limit);
   
                 if (p0ref == FS_REFERENCE)
                    limit -> SwapLabelAndTypeTwoTo01(ienergy, jangle);
                 else
                    limit -> SwapLabelAndTypeOneTo0(ienergy);
   
                 pffactorlimit = 1.;
                 sfactor = 1./(energySC*energySC*vSC);
                 process -> Evaluate(flagSC, ic, SOFT_AND_COLLINEAR,
                                     ienergy, jangle, limit,
                                     sfactor
                                        * jacrelSC
                                        * pffactorlimit
                                        * weight_factor,
                                     iSC);
              }
          }
      // sum over all unrelated terms, energy (T2)
      global -> tstore = 2;
      for (ienergy=0; ienergy < npartons; ++ienergy) {
   
          global -> istore = ienergy;
   
          // this is the second momentum used as a reference for
          // the first angle (in principle arbitrary)
          irelative = (ienergy==0) ? 1 : 0;
             
          // full term
   
          global -> lstore = NO_LIMIT;
          fullmaster -> CopyInto(full);
   
          if (p0ref == FS_REFERENCE)
             full -> SwapLabelAndTypeTwoTo01(ienergy, irelative);
          else
             full -> SwapLabelAndTypeOneTo0(ienergy);
    
          full -> MakeVList();
          full -> MakeEList();
          efactorfull = sumPartialEF;
   
          // soft limit
   
          global -> lstore = SOFT;
          limitS -> CopyInto(limit);
   
          if (p0ref == FS_REFERENCE)
             limit -> SwapLabelAndTypeTwoTo01(ienergy, irelative);
          else
             limit -> SwapLabelAndTypeOneTo0(ienergy);
   
          limit -> MakeVList();
          sfactor = 1./(energyS*energyS);
   
          vfactorfull = 0.;
          pffactorlimit = 0.;
          for (jangle = vmin; jangle < vmax; ++jangle)
              for (iangle = jangle+1; iangle < npartons; ++iangle)
                  if (   ienergy != iangle && ienergy != jangle
                      && !(iangle < 0 && jangle < 0)) {
                     vfactorfull
                        += 1./(1. + full -> flist[ienergy]
                                   /full -> vlist[full -> GetMapij()
                                                             [iangle][jangle]]
                              )
                         *full -> SumPartialV(iangle, jangle);
                     pffactorlimit += limit -> SumPartialV(iangle, jangle);
                  }
   
          process -> Evaluate(flagF, ic, NO_LIMIT, -1, -2, full,
                              vfactorfull*efactorfull*weight_factor,
                              iF);
   
          process -> Evaluate(flagS, ic, SOFT, ienergy, irelative, limit,
                              -sfactor*jacrelS*pffactorlimit
                               *weight_factor,
                              iS);
      }
   
      // sum over all unrelated terms, angle (T3)
      global -> tstore = 3;
      for (jangle = vmin; jangle < vmax; ++jangle)
          for (iangle = jangle+1; iangle < npartons; ++iangle) {
   
              if (!(iangle < 0 && jangle < 0)) {
   
                 global -> istore = iangle;
                 global -> jstore = jangle;
   
                 // full term
   
                 global -> lstore = NO_LIMIT;
                 fullmaster -> CopyInto(full);
   
                 if (p0ref == FS_REFERENCE)
                    full -> SwapLabelAndTypeTwoTo01(iangle, jangle);
                 else
                    full -> SwapLabelAndTypeOneTo0(iangle);
   
                 full -> MakeVList();
                 full -> MakeEList();
                 vfactorfull = sumPartialVF;
   
                 // collinear limit
   
                 global -> lstore = COLLINEAR;
                 limitC -> CopyInto(limit);
   
                 if (p0ref == FS_REFERENCE)
                    limit -> SwapLabelAndTypeTwoTo01(iangle, jangle);
                 else
                    limit -> SwapLabelAndTypeOneTo0(iangle);
       
                 limit -> MakeEList();
                 sfactor = 1./vC;
   
                 efactorfull = 0;
                 pffactorlimit = 0;
                 for (ienergy = 0; ienergy < npartons; ++ienergy)
                     if (ienergy != iangle && ienergy != jangle) {
                        efactorfull
                           += 1./(  1. 
                                  + full -> vlist[full -> GetMapij()
                                                             [iangle][jangle]]
                                   /full -> flist[ienergy]
                                 )
                            * full -> SumPartialE(ienergy);
                        pffactorlimit += limit -> SumPartialE(ienergy);
                     }
   
                 process -> Evaluate(flagF, ic, NO_LIMIT, -1, -2, full,
                                     vfactorfull*efactorfull
                                      *weight_factor,
                                     iF);
   
                 process -> Evaluate(flagC, ic, COLLINEAR, iangle, jangle,
                                     limit,
                                     -sfactor*jacrelC*pffactorlimit
                                      *weight_factor,
                                     iC);
              }
          }

      tEval -> GetEFree() -> ReturnObject(full);
      tEval -> GetEFree() -> ReturnObject(limit);

   } else {  // of if (tcpass)

      tEval -> GetEFree() -> ReturnObject(fullmaster);

   }
}

// ----------------------------------------------------------------------------

// adds the finite parts of the subtraction terms
// assumes all types of phase space variables to be defined
// (i.e. Cartesian, polar and normal)
// One more parton in events (to get the calculation of the invariants right)
// For simplicity, we only re-use the phase space. The limiting 
// event records are re-created several times.

void Event :: DoSumOverAddedSubtractionTerms(
                 Process* process,
                 int ic,
                 double weight_factor,
                 ContributionArray* ca
              ) {

   int inloc = n_in_partons;

   int iS, iCf, iCi, iSCi;
   int flagS, flagCf, flagCi, flagSCi;

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // soft limit

   Event* limitS = CreateAdded(process, SOFT_ADDED, 0, -1);
   global->user->PhaseSpace(process, limitS,
                            flagS,
                            ca -> NextContribution());
   iS = ca -> GetActual();

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // collinear limit final state

   Event* limitCf = CreateAdded(process, COLLINEAR_ADDED, 0, 1);
   global->user->PhaseSpace(process, limitCf,
                            flagCf,
                            ca -> NextContribution());
   iCf = ca -> GetActual();

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // collinear limit initial state

   Event* limitCi = CreateAdded(process, COLLINEAR_ADDED, 0, -1);
   global->user->PhaseSpace(process, limitCi,
                            flagCi,
                            ca -> NextContribution());
   iCi = ca -> GetActual();

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // soft and collinear limit final state

   Event* limitSCi = CreateAdded(process, SOFT_AND_COLLINEAR_ADDED, 0, -1);
   global->user->PhaseSpace(process, limitSCi,
                            flagSCi,
                            ca -> NextContribution());
   iSCi = ca -> GetActual();

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Event *f;

   int ienergy;
   int iangle, jangle;
   int flag, index;

      global->tstore=1;
      // collinear terms
      for (iangle=0; iangle < npartons+1; ++iangle)
          for (jangle=-inloc; jangle < iangle; ++jangle)
   
              if (iangle != jangle) {
   
                 // -----------------------------------------------------------

                 global->istore=iangle;
                 global->jstore=jangle;

                 // -----------------------------------------------------------
   
                 // collinear

                 if (jangle >= 0) {
                    flag = flagCf;
                    index = iCf;
                    f = limitCf;
                 } else {
                    flag = flagCi;
                    index = iCi;
                    f = limitCi;
                 }                    
                 f -> AddedRenumbering(iangle, jangle); 
                 process -> Evaluate(flag, ic, 
                                   COLLINEAR_ADDED, iangle, jangle, 
                                   f, weight_factor,
                                   index);

                 // -----------------------------------------------------------

                 // soft and collinear: for initial state singularity only

                 if ( jangle < 0 ) {
                    flag = flagSCi;
                    index = iSCi;
                    f = limitSCi;
                    f -> AddedRenumbering(iangle, jangle);
                    process -> Evaluate(
                                  flag, ic, 
                                  SOFT_AND_COLLINEAR_ADDED, iangle, jangle, 
                                  f, weight_factor,
                                  index);
                 }

                 // -----------------------------------------------------------
              }
   
      // soft terms
      global->tstore=2;
      for (ienergy=0; ienergy < npartons+1; ++ienergy) {
   
          global->istore=ienergy;

          int irelative = (ienergy==0) ? 1 : 0;
             
          flag = flagS;
          index = iS;
          f = limitS;
          f -> AddedRenumbering(ienergy, -1);
          process->Evaluate(flag, ic, SOFT_ADDED, ienergy, irelative, 
                            f, weight_factor,
                            index);
      }
}

// ----------------------------------------------------------------------------

// adds the finite parts of the subtraction terms
// assumes all types of phase space variables to be defined
// (i.e. Cartesian, polar and normal)
// One more parton in events (to get the calculation of the invariants right)
// For simplicity, we only re-use the phase space. The limiting 
// event records are re-created several times.
// for an ``npartons'' final state, we always add the collinear parton 
// with label ``npartons''.

void Event :: DoSumOverCollinearInitialTerms(
                 Process* process,
                 int ic,
                 double weight_factor,
                 ContributionArray* ca
              ) {

   int outloc = npartons;

   int iS, iCi, iSCi;
   int flagS, flagCi, flagSCi;

   // -------------------------------------------------------------------------

   Event* limitS = CreateAdded(process, SOFT_ADDED, outloc, -1);
   limitS -> AddedRenumbering(outloc, -1);
   global -> user -> PhaseSpace(process, limitS,
                                flagS,
                                ca -> NextContribution());
   iS = ca -> GetActual();
   process -> Evaluate(flagS, ic, SOFT_ADDED, 
                       outloc, 0, 
                       limitS, weight_factor,
                       iS);

   // -------------------------------------------------------------------------

   Event* limitCi = CreateAdded(process, COLLINEAR_ADDED, outloc, -1);
   limitCi -> AddedRenumbering(outloc, -1);
   global->user->PhaseSpace(process, limitCi,
                            flagCi,
                            ca -> NextContribution());
   iCi = ca -> GetActual();
   process -> Evaluate(flagCi, ic, 
                       COLLINEAR_ADDED, outloc, -1, 
                       limitCi, weight_factor,
                       iCi);

   // -------------------------------------------------------------------------

   Event *limitSCi = CreateAdded(process, 
                                 SOFT_AND_COLLINEAR_ADDED, outloc, -1);
   limitSCi -> AddedRenumbering(outloc, -1);
   global->user->PhaseSpace(process, limitSCi,
                            flagSCi,
                            ca -> NextContribution());
   iSCi = ca -> GetActual();
   process -> Evaluate(flagSCi, ic, 
                       SOFT_AND_COLLINEAR_ADDED, outloc, -1, 
                       limitSCi, weight_factor,
                       iSCi);
}

// ----------------------------------------------------------------------------

// permute particle labels according to ``list'' from co_min to co_max-1

void Event :: PermuteLabels(int co_max, int* list) {

   register int i;

   for (i = 0; i < co_max; ++i) {
       if (mapping[i] == NULL) {
          tEval -> To_error() << FS("i=%d\n", i);
          errf(-1, "Event :: PermuteLabels: permutation is illegal.");
       }
       mapping[i] -> SetLabel(list[i]);
   }
}

// ----------------------------------------------------------------------------

// :_MOD_: ecm is used only once -- to normalize the f_i

void Event :: DefineEcm() {

   switch (n_in_partons) {

      case 0:

         ecm = ecmfull;
         break;

      case 1:

         switch (frame) {

            case pCMS:
               double diff;
               diff = xi - xB;
               if (diff < 0.) {
                  tEval -> To_error()
                     << FS("%e", xB)
                     << FS(" %e", xi)
                     << FS(" %e\n", diff);
                  errf(-1,"Event :: DefineEcm: xi-xB < 0");
               }
               ecm = ecmfull * sqrt(y * diff);
               break;

            case hCMS:

               ecm = W;
               break;

            default:
               errf(-1, "Event :: DefineEcm: frame not known");
         }
         break;
      case 2:
         errf(-1, "not yet implemented");
         break;
   }
}

// ----------------------------------------------------------------------------

double Event :: GetEcm() {

   return ecm;
}

// ----------------------------------------------------------------------------

void Event :: SetEcmFull(double ecmIn) {

   ecmfull = ecmIn;
   SH = ecmfull * ecmfull;
}

// ----------------------------------------------------------------------------

double Event :: GetEcmFull() {

   return ecmfull;
}

// ----------------------------------------------------------------------------

void Event :: SetLeptonVariables(double xBI, double yI) {

   xB = xBI;
   y  = yI;
}

// ----------------------------------------------------------------------------

double Event :: GetxB() {

   return xB;
}

// ----------------------------------------------------------------------------

double Event :: Gety() {

   return y;
}

// ----------------------------------------------------------------------------

// this sets xi and xi_hard to the same value. Setting xi_hard is redundant
// when it is not needed; but if it is needed, xi is modified later in the
// program (initial-state splitting)
// :MOD: should be renamed Set_xi_And_xi_hard
// also sets xi_int and xi_pden

void Event :: Setxi(double xiI) {

   xi      = xiI;
   xi_hard = xiI;
   xi_int  = xiI;
   xi_pden = xiI;
}

// ----------------------------------------------------------------------------

double Event :: Getxi() {

   return xi;
}

// ----------------------------------------------------------------------------

// the calling routine for the mapping for both pCMS and hCMS.

void Event :: MapUnitToPSGeneral(
                 PSMode psmode, 
                 P0_Reference p0refIn,
                 double* unit, 
                 double* jac1, 
                 double* jacrest
              ) {

   switch (frame) {

      case pCMS:
         if (npartons <= 1)
            errf(-1,"Event :: MapUnitToPSGeneral: npartons too small.");
         MapUnitToPSGeneralPCMS(psmode, p0refIn, unit, jac1, jacrest);
         break;

      case hCMS:
         if (npartons > 2)
            errf(-1, "Event::MapUnitToPSGeneral: npartons too large.");
         MapUnitToPSGeneralHCMS(psmode, p0refIn, unit, jac1, jacrest);
         break;

      default:
         errf(-1,"Event :: MapUnitToPSGeneral: no such frame.");
   }
}

// ----------------------------------------------------------------------------

// Phase space for the pCMS:
// maps an array of numbers in [0,1] to a set of phase space variables
// with given CM energy ecm; returns the jacobian factors for
// the first particle and for the rest separately.
// The energy, polar and azimuthal angles can be set explicitly
// for the first particle, if psmode is set to `FIRSTEXPLICIT':
// unit[ 0] = energy (in absolute units), 
// unit[ 1] = polar angle in [0,1],
// unit[2]    not used
// unit[-2] = sin(azimuthal angle)
// unit[-1] = cos(azimuthal angle)
// In this case, the jacobian of the first particle is set to zero
// (because it is undetermined).
// p0ref determines the reference momentum for p0.
// If set to IS_REFERENCE, the first three entries of `unit' are
// as stated above, but now the polar angle and the azimuthal angle
// relative to the z-axis. 
// If set to `ALLUNIT', the mapping is done by a variable transformation
// from [0,1]^(3n-4).
// Polar, normalized and Cartesian coordinates are defined.

// :_CAUTION_: have to do the mapping for v -> 1 for the remaining partons!

void Event :: MapUnitToPSGeneralPCMS(
                 PSMode psmode, 
                 P0_Reference p0refIn,
                 double* unit,
                 double* jac1, 
                 double* jacrest
               ) {

   register int i;
   register int imom;

   double x,yvar,v,z,r[4],s,e1_L,v1,e1max,e2_L,v2,e2max,v12_L,
          v1z,rabs,cosphi,sinphi,coschim,
          ejacobian,vjacobian,lambda_L,energy,emax,
          rdotn,rsquared,rdenom;

   double sinPhi1;
   double cosPhi1;

   double sinPhi2;
   double cosPhi2;

   double t;
   double f; 

   double jacobian;

   nparticles = 0; 

   Particle* k1;
   Particle* k2;
   Particle* k;
   Particle* kfinal;

   p0ref = p0refIn;

   // :_TAWM_:
   e1_L = 0.;
   v1   = 0.;

   r[0] = ecm;
   r[1] = 0.;
   r[2] = 0.;
   r[3] = 0.;

   // the 2*pi-factors in the jacobian:
   *jacrest = pow(TwoPi, (double) 4 - 3 * npartons);

   *jac1 = 1.;

   int ind = 0;

   // first and second momentum

   k1 = AddParticle(0, 'u');

   if (psmode == ALLUNIT) {

      e1max = 0.5 * ecm;

      if (npartons >= 3) {

         x = unit[ind++];
         e1_L =   e1max 
                * variableTransformation(
                     x,
                     1, 
                     vtp -> Get_mapParameter(0),
                     vtp -> Get_mapParameter(6),
                     vtp -> Get_mapParameter(7),
                     vtp -> Get_mapParameter(8),
                     &jacobian
                  );
         ejacobian = e1max * jacobian;
 
      } else {

         ejacobian = 1. / (2. * ecm);
         e1_L = e1max;

      }

      yvar = unit[ind++];
      v = variableTransformation(
             yvar,
             1, 
             vtp -> Get_mapParameter( 1),
             vtp -> Get_mapParameter(13),
             vtp -> Get_mapParameter(14),
             vtp -> Get_mapParameter(15),
             &vjacobian
          );

      z   = unit[ind++];
      v1  = v;
      double phi1 = TwoPi * z;
      sinPhi1 = sin(phi1);
      cosPhi1 = cos(phi1);

      *jac1 *= 0.5 * e1_L * ejacobian * vjacobian * 2. * TwoPi;

   } else if (psmode == FIRSTEXPLICIT) {

      // :_MOD_: should distinguish here between npartons==2 and npartons==3
      // (no energy variable required in the former case)
      e1_L    = unit[ind++];
      v1      = unit[ind++];
      ++ind;
      sinPhi1 = unit[-2];
      cosPhi1 = unit[-1];
      *jac1 = 0.;

   } else {

      // :_TAWM_:
      sinPhi1 = 0;  
      cosPhi1 = 0;

      errf(-1, "Event :: MapUnitToPSGeneral: psmode not defined");

   }

   if (npartons >= 3) {

      yvar = unit[ind++];
      v2 = variableTransformation(
              yvar,
              1, 
              vtp -> Get_mapParameter( 3),
              vtp -> Get_mapParameter(10),
              vtp -> Get_mapParameter(11),
              vtp -> Get_mapParameter(12),
              &vjacobian
           );

      // prevent v2 from being == 1 
      // (this would mean parton[1] || initial in Born term)
      // :_TCUT_:
      if (v2 > 1.-1.e-8)
         v2 = 1.-1.e-8;

      z = unit[ind++];
      double phi2 = TwoPi * z;
      sinPhi2 = sin(phi2);
      cosPhi2 = cos(phi2);

      switch (p0ref) {

         case FS_REFERENCE:
            v12_L = v1;
            break;

         case IS_REFERENCE:
            v12_L = ThreePolarVWCos(v1, v2,  cosPhi1 * cosPhi2
                                           + sinPhi1 * sinPhi2);
            break;

         default:
            errf(-1,"Event :: MapUnitToPSGeneral: p0ref not defined (1)");
            v12_L = 0; // :_TAWM_:
      }
         
      e2max = 0.5 * ecm * (ecm - 2. * e1_L) / (ecm - 2. * e1_L * v12_L);

      if (npartons >= 4) {

         x = unit[ind++];
         t = pow(x, f = vtp -> Get_mapParameter(2));
         e2_L = e2max * t;
         ejacobian = e2max * f * t / x;

      } else {

         e2_L = e2max;
         ejacobian = 1. / (2. * (ecm - 2. * e1_L * v12_L));

      }

      *jacrest *= 0.5 * e2_L * ejacobian * vjacobian * 2. * TwoPi;

      k2 = AddParticle(1, 'u');
      k2 -> SetEnergy(e2_L);
      k2 -> SetMomentum(e2_L);
      k2 -> SetV(v2);
      k2 -> SetSinCosPhi(sinPhi2, cosPhi2);
      k2 -> CalculateCartesianFromPolar();

      // memorize this for later evaluation!
      k1 -> SetSinCosPhiRef(sinPhi1, cosPhi1);

      double sinPhitwiddle;
      double cosPhitwiddle;

      switch (p0ref) {

         case FS_REFERENCE:
            v1z = ThreePolarVWCos(v1, v2, cosPhi1);
            cosPhitwiddle = GetCosPhi(v1, v1z, v2);
            sinPhitwiddle = AbsSinFromCos(cosPhitwiddle);
            if (sinPhi1 > 0.)
               sinPhitwiddle = -sinPhitwiddle;
            sinPhi1 = sinPhi2 * cosPhitwiddle + cosPhi2 * sinPhitwiddle;
            cosPhi1 = cosPhi2 * cosPhitwiddle - sinPhi2 * sinPhitwiddle;
            break;

         case IS_REFERENCE:
            v1z = v1;
            break;

         default:
            errf(-1,"Event::MapUnitToPSGeneral: p0ref not defined (2)");
            v1z = 0; // TAWM
      }

      k1 -> SetEnergy(e1_L);
      k1 -> SetMomentum(e1_L);
      k1 -> SetV(v1z);
      k1 -> SetSinCosPhi(sinPhi1, cosPhi1);
      k1 -> CalculateCartesianFromPolar();

      for (i = 0; i <= 3; ++i)
          r[i] -= k1 -> GetFvect(i) + k2 -> GetFvect(i);

   } else {  // of (npartons >= 3)

      k1 -> SetEnergy(e1_L);
      k1 -> SetMomentum(e1_L);
      k1 -> SetV(v1);
      k1 -> SetSinCosPhi(sinPhi1, cosPhi1);
      k1 -> CalculateCartesianFromPolar();

      for (i = 0; i <= 3; ++i)
          r[i] -= k1 -> GetFvect(i);

   }  // of (npartons >= 3)

   // now determine the remaining momenta (up to the final two)
   for (imom = 3; imom <= npartons - 2; ++imom) {

       k = AddParticle(imom-1, 'u');

       rabs = sqrt(r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);

       x    = unit[ind++];
       yvar = unit[ind++];
       z    = unit[ind++];
   
       t = pow(yvar, f = vtp -> Get_mapParameter(3));
       v = t;
       vjacobian = f * t / yvar;

       double phi = TwoPi*z;
       cosphi = cos(phi);
       sinphi = sin(phi);

       s = 2. * sqrt(v * (1. - v)); // :_CHECK_: check whether well-defined!
       coschim = (  s * (r[1] * cosphi + r[2] * sinphi) 
                  + r[3] *( 1. -2 * v))
                / r[0];

       lambda_L = rabs / r[0];
       emax = 0.5 * r[0] * (1. - lambda_L * lambda_L) / (1. - coschim);

       t = pow(x, f = vtp -> Get_mapParameter(2));
       energy = emax * t;
       ejacobian = emax * f * t / x;

       *jacrest *= 0.5 * energy * ejacobian * vjacobian * 2. * TwoPi;
       
       k -> SetEnergy(energy);
       k -> SetMomentum(energy);
       k -> SetVAndSV(v, s);
       k -> SetSinCosPhi(sinphi, cosphi);
       k -> CalculateCartesianFromPolar();

       for (i = 0; i <= 3; ++i)
           r[i] -= k -> GetFvect(i);
   }

   // the final two momenta, if at least 4 particles

   if (npartons >= 4) {

      k = AddParticle(npartons - 2, 'u');

      rabs = sqrt(r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);

      yvar = unit[ind++];
      t = pow(yvar, f = vtp -> Get_mapParameter(3));
      v = t;
      vjacobian = f * t / yvar;

      z = unit[ind++];
   
      double phi = TwoPi * z;
      cosphi = cos(phi);
      sinphi = sin(phi);

      s = 2. * sqrt(v * (1. - v)); // :_CHECK_: check whether well-defined!
      rdotn = (s * (r[1] * cosphi + r[2] * sinphi) + r[3] *(1. - 2. * v));
      rsquared = r[0] * r[0] -rabs * rabs;
      rdenom = 2. * (r[0] - rdotn);

      energy = rsquared / rdenom;

      *jacrest *= 0.5 * rsquared / (rdenom * rdenom) * vjacobian * 2. * TwoPi;
       
      k -> SetEnergy(energy);
      k -> SetMomentum(energy);
      k -> SetV(v);
      k -> SetSinCosPhi(sinphi, cosphi);
      k -> CalculateCartesianFromPolar();

      kfinal = AddParticle(npartons - 1, 'u');

      for (i = 0; i <= 3; ++i) 
          kfinal -> SetFvect(i, r[i] - k -> GetFvect(i));

   } else {

      kfinal = AddParticle(npartons - 1, 'u');

      for (i = 0; i <= 3; ++i)
          kfinal -> SetFvect(i, r[i]);

      // :_CHECK_: jacrest *= ???
   }   

   kfinal -> CalculatePolarFromCartesian();
}

// ----------------------------------------------------------------------------

// Phase space for the hCMS:
// maps an array of numbers in [0,1] to a set of phase space variables
// with given CM energy ecm; returns the jacobian factors for
// the first particle and for the rest separately.
// The energy, polar and azimuthal angles can be set explicitly
// for the first particle, if psmode is set to `FIRSTEXPLICIT':
// unit[ 0] = energy (in absolute units), 
// unit[ 1] = polar angle in [0,1],
// unit[ 2]   not used
// unit[-2] = sin(azimuthal angle)
// unit[-1] = cos(azimuthal angle)
// In this case, the jacobian of the first particle is set to zero
// (because it is undetermined).
// p0ref determines the reference momentum for p0.
// If set to IS_REFERENCE, the first three entries of `unit' are
// as stated above, but now the polar angle and the azimuthal angle
// relative to the z-axis. 
// If set to `ALLUNIT', the mapping is done by a variable transformation
// from [0,1]^(3n-4).
// Polar, normalized and Cartesian coordinates are defined.

void Event :: MapUnitToPSGeneralHCMS(
                 PSMode psmode, 
                 P0_Reference p0refIn,
                 double* unit, 
                 double* jac1, 
                 double* jacrest
              ) {

   nparticles = 0; 

   register int i;

   double v,x,yvar,z,r[4],e1_L,e1max,
          v1z,ejacobian,vjacobian;

   double t;
   double f;

   Particle* k1;
   Particle* k2;

   p0ref = p0refIn;

   // first momentum
   k1 = AddParticle(0, 'u');

   // the 2*pi-factors in the jacobian:
   *jacrest = pow(TwoPi, (double) 4 - 3 * npartons);

   *jac1 = 1.;

   if (npartons == 1) {

      // 1 outgoing parton: parton model

      if (psmode != ALLUNIT)
         errf(-1, "Event :: MapUnitToPSGeneralHCMS: ALLUNIT not defined");

      *jacrest *= xB / Q2;
      Setxi(xB);

      k1 -> SetEnergy(0.5 * W);
      k1 -> SetMomentum(0.5 * W);
      k1 -> SetVAndSV(1., 0.);
      k1 -> SetSinCosPhi(0., 1.);

   } else {

      // 2 outgoing partons

      int ind = 0;

      double sinPhi;
      double cosPhi;

      if (psmode == ALLUNIT) {

         e1max = 0.5 * W;

         x = unit[ind++];
         t = pow(x, f = vtp -> Get_mapParameter(0));
         e1_L = e1max * t;
         ejacobian = e1max * f * t / x;

         yvar = unit[ind++];
         t = pow(yvar, f = vtp -> Get_mapParameter(1));
         v = t;
         vjacobian = f * t / yvar;

         z = unit[ind++];
         double phi = TwoPi * z;
         sinPhi = sin(phi);
         cosPhi = cos(phi);

         *jac1 *= 0.5 * e1_L * ejacobian * vjacobian * 2. * TwoPi;

      } else if (psmode == FIRSTEXPLICIT) {

         e1_L = unit[ind++];
         v    = unit[ind++];
         ++ind;
         sinPhi = unit[-2];
         cosPhi = unit[-1];
         *jac1 = 0.;
   
      } else {

         errf(-1, "Event :: MapUnitToPSGeneralHCMS: psmode not defined");
         // TAWM
         e1_L   = 0.;
         v      = 0.;
         sinPhi = 0;
         cosPhi = 0;
      } 

      k1 -> SetSinCosPhiRef(sinPhi, cosPhi);

      double eps1 = 2. * e1_L / W;

      // now construct the momentum vectors

      double sinPhi1 = sinPhi;
      double cosPhi1 = cosPhi;

      switch (p0ref) {
         case FS_REFERENCE:
            // for energy close to PS limit evaluate differently!
            if (eps1 <= 0.9)
               v1z = (1. - v) / (1. - eps1 * (2. - eps1) * v);
            else
               v1z =   1. / (1. + (1. - eps1) 
                    * (1. - eps1) * (1. / (1. - v) - 1.));
            *jacrest *= (1. - eps1) / (1. - eps1 * (2. - eps1) * v);
            break;
         case IS_REFERENCE:
            v1z = v;
            break;
         default:
            errf(-1, "Event :: MapUnitToPSGeneralHCMS: p0ref not defined.");
            // :_TAWM_:
            v1z = 0;
      }
      
      *jacrest *= (1. - xB) / (W2 * (1. - eps1 * v));
      Setxi(xB + eps1 * (1. - xB) * (1. - v1z) / (1. - eps1 * v1z));

      k1 -> SetEnergy(e1_L);
      k1 -> SetMomentum(e1_L);
      k1 -> SetV(v1z);
      k1 -> SetSinCosPhi(sinPhi1, cosPhi1);
      k1 -> CalculateCartesianFromPolar();

      // this is actually p1 + p2
      r[0] = proton_energy * (1. - 2. * xB + xi);
      r[1] = 0.;
      r[2] = 0.;
      r[3] = proton_energy * (-1. + xi);

      // second momentum
      k2 = AddParticle(1, 'u');

      for (i = 0; i <= 3; ++i)
          k2 -> SetFvect(i, r[i] - k1 -> GetFvect(i));

      k2 -> CalculatePolarFromCartesian();
   }
}

// ----------------------------------------------------------------------------

// map the unit hypercube to phase space; 
// return the jacobian.

double Event :: MapUnitToPSChoice(
                   P0_Reference p0refIn,
                   double* unit
                ) { 

   double jac1;
   double jacrest;

   MapUnitToPSGeneral(ALLUNIT, p0refIn, unit, &jac1, &jacrest);

   return jac1 * jacrest;
}

// ----------------------------------------------------------------------------

// map the unit hypercube to phase space; 
// return the jacobian.
// use the reference pointer as already specified in the event.

double Event :: MapUnitToPS(double* unit) { 

 return MapUnitToPSChoice(p0ref, unit);
}

// ----------------------------------------------------------------------------

// the calling routine for the inverse mapping for both pCMS and hCMS.

void Event :: MapPSToUnitGeneral(
                 P0_Reference p0refIn,
                 double* unit
              ) {

   switch (frame) {
      case pCMS:
         if (npartons <= 1)
            errf(-1, "Event :: MapPSToUnitGeneral: npartons too small.");
         MapPSToUnitGeneralPCMS(p0refIn, unit);
         break;
      case hCMS:
         if (npartons > 2)
            errf(-1, "Event :: MapPSToUnitGeneral: npartons too large.");
         MapPSToUnitGeneralHCMS(p0refIn, unit);
         break;
      default:
         errf(-1, "Event :: MapPSToUnitGeneral: no such frame.");
   }
}

// ----------------------------------------------------------------------------

// map phase space to the unit hypercube in the pCMS

void Event :: MapPSToUnitGeneralPCMS(
                 P0_Reference p0refIn,
                 double* unit
              ) {

   register int i;
   register int imom;

   double r[4], v,x,rabs,lambda_L,s,coschim,ejacobian,
          v2,v1z;

   double sinPhi;
   double cosPhi;

   Particle* k;
   Particle* k2;

   int ind = 3 * npartons - 5;

   for (i = 0; i <= 3; ++i)
       r[i] = 0.;

   k = particle[npartons-1];
   k -> CalculateCartesianFromPolar();
   for (i = 0; i <= 3; ++i)
       r[i] = k -> GetFvect(i);

   // work backwards; treat the angles of the first momentum differently, 
   // but only if there are not less than three momenta
   // and if the reference momentum is in the final state!
    
   for (imom = npartons - 2; imom >= 0; --imom) {

       k = particle[imom];
       k -> CalculateCartesianFromPolar();

       for (i = 0;i <= 3; ++i)
           r[i] += k -> GetFvect(i);

       // is it the first momentum?
       if (imom == 0 && npartons >= 3 && p0refIn == FS_REFERENCE) {

          k2 = particle[1];
          k2 -> CalculateCartesianFromPolar();
          v2    = k2 -> GetV();
          v1z   = k -> GetV();
          v = V(k, k2);
          cosPhi = GetCosPhi(v1z, v, v2);
          sinPhi = AbsSinFromCos(cosPhi);
          // modify phi?
          if (  k -> GetCosPhi() * k2 -> GetSinPhi()
              - k -> GetSinPhi() * k2 -> GetCosPhi() < 0.)
             sinPhi = - sinPhi;

       } else {

          v = k -> GetV();
          k -> GetSinCosPhi(&sinPhi, &cosPhi);

       }

       unit[ind--] = atan2ZeroTwoPi(sinPhi, cosPhi) / TwoPi;

       if (imom == 0)
          unit[ind--] 
             = variableTransformationInverse(
                  v,
                  1,
                  vtp -> Get_mapParameter(1),
                  vtp -> Get_mapParameter(13),
                  vtp -> Get_mapParameter(14),
                  vtp -> Get_mapParameter(15)
               );
       else
          unit[ind--] 
             = variableTransformationInverse(
                  v,
                  1,
                  vtp -> Get_mapParameter(3),
                  vtp -> Get_mapParameter(10),
                  vtp -> Get_mapParameter(11),
                  vtp -> Get_mapParameter(12)
               );

       // not the final two momenta?
       if (imom <= npartons - 3) {

          rabs = sqrt(r[1] * r[1] + r[2] * r[2] + r[3] * r[3]);
          lambda_L = rabs / r[0];
          s = 2. * sqrt(v * (1. - v)); // check whether well-defined!
          coschim = (  s * (r[1] * cosPhi + r[2] * sinPhi)
                     + r[3] * (1. - 2. * v)) / r[0];
          ejacobian = 0.5 * r[0] * (1. - lambda_L * lambda_L) / (1. - coschim);
          x = k -> GetEnergy() / ejacobian;
          if (imom == 0)
             unit[ind--] 
                = variableTransformationInverse(
                     x,
                     1,
                     vtp -> Get_mapParameter(0),
                     vtp -> Get_mapParameter(6),
                     vtp -> Get_mapParameter(7),
                     vtp -> Get_mapParameter(8)
                  );
          else
             unit[ind--] = pow(x, 1. / vtp -> Get_mapParameter(2));
       }
   }
}

// ----------------------------------------------------------------------------

// map phase space to the unit hypercube in the hCMS

void Event :: MapPSToUnitGeneralHCMS(
                 P0_Reference p0refIn,
                 double* unit) {

   if (npartons == 1) {

      // nothing to determine!

   } else {

      double v;
      Particle* k;

      int ind = 0;

      k = particle[0];

      unit[ind++] = pow(2. * k -> GetEnergy() / W, 
                        1. / vtp -> Get_mapParameter(0)
                    );

      switch (p0refIn) {
         case FS_REFERENCE:
            v = V(k, particle[1]);
            break;
         case IS_REFERENCE:
            v = k -> GetV();
            break;
         default:
            // :_TAWM_:
            v = 0.;
      }

      unit[ind++] = pow(v, 1. / vtp -> Get_mapParameter(1));

      unit[ind++] = k -> CalculatePhi() / TwoPi;
   }
}

// ----------------------------------------------------------------------------

// map phase space to the unit hypercube; 
// use original reference as actual reference

void Event :: MapPSToUnit(double* unit) {

   MapPSToUnitGeneral(p0ref, unit);
}

// ----------------------------------------------------------------------------

// performs soft, collinear and double limits for an event.
// Assumes that the partons occupy the first few entries in the event record
// Normalized, polar and Cartesian coordinates are assumed to be defined.
// Initial state terms are handled.

#define PERFLIMIT 0

void Event :: Limit1(
                 LimitType limittype, 
                 P0_Reference p0refIn,
                 Event* limit,
                 double* jrelative, 
                 double* energy, 
                 double* v
              ) {

   if (disaster -> Get_print_debug()) {
      global -> To_status() << "Event :: Limit1: original event:\n";
      Print(global -> To_status());
   }   

   limit -> Reset();
   CopyPSInformationInto(limit);
   limit -> SetPartonNumber(GetPartonNumber());
   limit -> SetLimitType(limittype);
   switch (p0refIn) {
      case FS_REFERENCE:
         limit -> SetLimit1(particle[0] -> GetLabel());
         limit -> SetLimit2(particle[1] -> GetLabel());
         break;
      case IS_REFERENCE:
         limit -> SetLimit1(particle[0] -> GetLabel());
         limit -> SetLimit2(-1);
         break;
   }

#if PERFLIMIT==1
   Event ratio,collapsed;
   double frac;
#endif

   // :_MOD_: should be improved: limit is possible (eP -> 2 jets)
   if (npartons <= 2 && frame == pCMS)
      errf(-1, "Event :: Limit1: less than 2 final state particles.");

   int ulength = 3 * npartons - 4;
   if (frame == hCMS)
      ++ulength;
   // + 2 to account for sinPhi [-2] and cosPhi [-1] variables
   double* u_L = new double[ulength + 2] + 2;

   Event *edum = tEval
                 -> GetEFree()
                 -> GetDefinedObject(tEval, GetMaxNumber()); 
   CopyPSInformationInto(edum);
   edum -> SetPartonNumber(GetPartonNumber());

   // get parameters in unit cube for the given reference system
   MapPSToUnitGeneral(p0refIn, u_L);

   if (disaster -> Get_print_debug()) {
      global -> To_status() << "Event :: Limit1: u_L:\n";
      for (int i = -2; i < ulength; ++i)
          global -> To_status() << FS("u_L[%d]", i)
                                << FS("%30.22e\n", u_L[i]);
   }   

   // get jacobian and kinematics for the full (``unlimited'') PS
   double jdummy;
   double jbeforelimit;
   edum -> MapUnitToPSGeneral(ALLUNIT, p0refIn, u_L, &jdummy, &jbeforelimit);

   if (disaster -> Get_print_debug()) {
      global -> To_status() << "Event :: Limit1: edum:\n";
      edum -> Print(global -> To_status());
   }   

   double vold;
   double sinPhiOld;
   double cosPhiOld;

//   double eold = edum -> particle[0] -> GetFvect(0);
   double eold = edum -> particle[0] -> GetEnergy();

   double v0, v1, v01_L;
   v0 = edum -> particle[0] -> GetV();
   
   double sinPhi0;
   double cosPhi0;
   edum -> particle[0] -> GetSinCosPhi(&sinPhi0, &cosPhi0);

   double sinPhi0z;
   double cosPhi0z;
   double sinPhi1;
   double cosPhi1;

   if (p0refIn == FS_REFERENCE) {

      vold = v01_L = V(edum -> particle[0], edum -> particle[1]);

      switch (frame) {

         case pCMS:

            v1 = edum -> particle[1] -> GetV();
            cosPhi0z = GetCosPhi(v0, v01_L, v1);
            sinPhi0z = AbsSinFromCos(cosPhi0z);
            edum -> particle[1] -> GetSinCosPhi(&sinPhi1, &cosPhi1);
            if (sinPhi0 * cosPhi1 - sinPhi1 * cosPhi0 >= 0.) 
               sinPhi0z = - sinPhi0z;
            sinPhiOld = sinPhi0z;
            cosPhiOld = cosPhi0z;
            break;
  
         case hCMS:
            sinPhiOld = sinPhi0;
            cosPhiOld = cosPhi0;
            break;

         default:
            // :_TAWM_:
            sinPhiOld = 0.;
            cosPhiOld = 0.;
      }

   } else {

      vold   = v0;
      sinPhiOld = sinPhi0;
      cosPhiOld = cosPhi0;
   }

   // return the actual values of the energy and angle
   *energy = eold;
   *v      = vold;

#if PERFLIMIT==1
   tEval -> To_status() << FS("old event: (%f)\n", jbeforelimit);
   edum -> Print(stdout);
#endif

   u_L[-2] = sinPhiOld;
   u_L[-1] = cosPhiOld;
   u_L[2] = 1e40;

   switch (limittype) {
      case SOFT:
         // soft limit phase space
         u_L[0] = 0.;  // soft limit
         u_L[1] = vold;
         break;
      case COLLINEAR:
         // collinear limit phase space
         u_L[0] = eold;
         u_L[1] = 0.;  // collinear limit
         break;
      case SOFT_AND_COLLINEAR:
         // soft and collinear limit phase space
         u_L[0] = 0.;  // soft limit
         u_L[1] = 0.;  // collinear limit
         break;
      default:
         errf(-1,"Event::Limit1: no limit is impossible!");
   }

   // create the event record in the appropriate limit
   double jafterlimit;
   limit -> MapUnitToPSGeneral(
               FIRSTEXPLICIT,
               p0refIn,
               u_L,
               &jdummy,
               &jafterlimit
            );
   TransferPartonLabelsInto(limit);
   CopyAdditionalParticlesInto(limit);

#if PERFLIMIT==1
   tEval -> To_status() << FS("new event: (%f)\n", jafterlimit);
   limit->Print(stdout);
#endif

   if (fabs(jbeforelimit) < 1.e-50) {
      *jrelative = 0.;
      const char msg1[] = "Event :: Limit1 (1)";
      fpnan -> Set_flag_mth(-1, msg1);
   } else
      *jrelative = jafterlimit / jbeforelimit;

#if PERFLIMIT==1
   // stepwise limit
   collapsed.Define(GetMaxNumber(),global);
   limit->CopyInto(&collapsed);
   // collapse the first two vectors (in the soft case, [0] is zero anyway)
   collapsed.Combine(0,1);
   tEval << "after collapse:\n";
   collapsed.Print(stdout);

   ratio.Define(npartons-1,global);
   
   MapPSToUnit(u);
   edum -> MapUnitToPSGeneral(ALLUNIT, p0ref, u, &jdummy, &jbeforelimit);
   for (frac=1.;frac>=0.0000000001;frac/=10.) {
       if (limittype==SOFT) {
          u[0]=frac*eold;
          u[1]=vold;
          u[2]=phiold;
       } else {
          u[0]=eold;
          u[1]=frac*vold;
          u[2]=phiold;
       }
       edum -> SetNumber(npartons);
       edum -> MapUnitToPSGeneral(FIRSTEXPLICIT, FS_REFERENCE,
                                  u, &jdummy, &jafterlimit);
       // collapse the first two vectors 
       edum -> Combine(0,1);
       ratio.Divide(&collapsed, edum);
       tEval -> To_status()
          << FS("RELATIVE frac=%e", frac)
          << FS(": (%f)\n", jafterlimit);
       ratio.Print(stdout);
   }
#endif

   delete [] (u_L - 2);
  
   tEval -> GetEFree() -> ReturnObject(edum);

   if (disaster -> Get_print_debug()) {
      global -> To_status() << "Event :: Limit1: limit event:\n";
      limit -> Print(global -> To_status());
   }   
}

// ----------------------------------------------------------------------------

// perform the limit of type `limittype' for particles with _labels_ i, j. 
// The event record `in' may be reordered...
// Polar, normalized and Cartesian coordinates are defined.
// If i or j are negative, the limit is performed by means
// of an initial state phase space
// We assume that we know what we do: no check w.r.t. NULL pointers
// We assume that the partons are always the first entries in the
// event record!

void Event :: Limit(
                 LimitType limittype,
                 int i,
                 int j,
                 Event *limit,
                 double* jrelative,
                 double* energy,
                 double* v
              ) {

   int ipos;
   int jpos;

   if (i >= 0 && j >= 0) {

      ipos = mapping[i] -> GetLocation(); 
      jpos = mapping[j] -> GetLocation(); 
      PermuteTwoTo01(ipos, jpos);
      Limit1(limittype, FS_REFERENCE, limit, jrelative, energy, v);

   } else

   if (i < 0) {
      jpos = mapping[j] -> GetLocation(); 
      PermuteOneTo0(jpos);
      Limit1(limittype, IS_REFERENCE, limit, jrelative, energy, v);

   } else {

      ipos = mapping[i] -> GetLocation(); 
      PermuteOneTo0(ipos);
      Limit1(limittype, IS_REFERENCE, limit, jrelative, energy, v);
   }
}

// ----------------------------------------------------------------------------

// map unit variables to the lepton xB and y variables;
// use variable transformation to optimize MC integration
// based on the cuts in *global.
// returns the value of the jacobian

double Event :: MapUnitToxBy(double zx, double zq) {


   double log_xBmin_actual;
   double log_xBmax_actual;
   double Q2min_loc;
   double Q2max_loc;
   double log_Q2min_loc;
   double log_Q2max_loc;
   double Q2loc;
   double trho;
   double tsigma;
   double Amin;
   double Amax;
   double A;
   double Bmin;
   double Bmax;
   double B;

   double jacobian;

   switch (vtp -> Get_mapFunction(0)) {

      case 0:

         // logarithmic transformation

         log_xBmin_actual = log(disaster -> Get_xBmin_actual());
         log_xBmax_actual = log(disaster -> Get_xBmax_actual());

         xB = exp(  zx        * log_xBmax_actual
                  + (1. - zx) * log_xBmin_actual
              );

         Q2min_loc = SH * xB * dmax3(disaster 
                                     -> Get_Q2min_actual() / SH / xB,
                                     disaster
                                     -> Get_W2min_actual() / SH / (1. - xB),
                                     disaster
                                     -> Get_ymin()
                               );

         Q2max_loc = SH * xB * dmin3(disaster
                                     -> Get_Q2max_actual() / SH / xB,
                                     disaster
                                     -> Get_W2max_actual() / SH / (1. - xB),
                                     disaster
                                     -> Get_ymax()
                               );

         if (Q2max_loc < Q2min_loc)
            errf(-1, "Q2max_loc < Q2min_loc");

         log_Q2min_loc = log(Q2min_loc);
         log_Q2max_loc = log(Q2max_loc);
 
         Q2loc = exp(  zq * log_Q2max_loc
                     + (1. - zq) * log_Q2min_loc
                 );

         y = Q2loc / SH / xB;

         if (y > 1.) 
            y = 0.9999;

         jacobian
            =     (log_xBmax_actual - log_xBmin_actual)
                * (log_Q2max_loc - log_Q2min_loc)
                * xB * y;

         break;

      case 1:
   
         // power law transformation

         trho   = vtp -> Get_mapParameter(4);
         tsigma = vtp -> Get_mapParameter(5);

         Amin = pow(disaster -> Get_xBmin_actual(), trho);
         Amax = pow(disaster -> Get_xBmax_actual(), trho);

         A = Amin + zx * (Amax - Amin);

         xB = pow(A, 1. / trho);

         Q2min_loc = SH * xB * dmax3(disaster 
                                     -> Get_Q2min_actual() / SH / xB,
                                     disaster
                                     -> Get_W2min_actual() /SH / (1. - xB),
                                     disaster
                                     -> Get_ymin()
                               );

         Q2max_loc = SH * xB * dmin3(disaster
                                     -> Get_Q2max_actual() / SH / xB,
                                     disaster
                                     -> Get_W2max_actual() / SH / (1. - xB),
                                     disaster
                                     -> Get_ymax()
                               );

         if (Q2max_loc < Q2min_loc)
            errf(-1, "Q2max_loc < Q2min_loc");

         Bmin = pow(Q2min_loc, tsigma - 1.);
         Bmax = pow(Q2max_loc, tsigma - 1.);

         B = Bmin + zq * (Bmax-Bmin);
         Q2loc = pow ( B, 1./(tsigma-1.) );

         y = Q2loc / SH / xB;

         if (y > 1.) 
            y = 0.9999;

         jacobian
            =     1. / SH
                * (Amax - Amin) / trho          / A
                * (Bmax - Bmin) / (tsigma - 1.) / B * Q2loc;

         break;

      default:
         
         errf(-1, "Event :: MapUnitToxBy: parametrization not known");
         jacobian = 0.;  // :_TAWM_:
   }

   return jacobian;
}

// ----------------------------------------------------------------------------

// The lepton phase space jacobian; for dxB dy.

double Event :: LeptonJacobian() {

   Q2 = SH * xB * y;
   W2 = SH * (1. - xB) * y;

   Q = sqrt(Q2);
   W = sqrt(W2);

   return 1. / (16. * Pi * Pi) * SH * y;
}

// ----------------------------------------------------------------------------

// map a unit variables to the momentum fraction xi.
// parametrizes dxi / xi (??) from xB to 1
// use variable transformation to optimize MC integration
// assume xB is defined (lower bound)
// returns the value of the jacobian

double Event :: MapUnitToXi(double zxi) {

   // first do a power law transformation

   double locJacobian;
   double zLocal = 1. - variableTransformationLower(
                           1. - zxi,  
                           1, 
                           vtp -> Get_mapParameter(16),
                           &locJacobian
                        );

   // larger xi_alpha -> zLocal closer to 0

   double log_ximin = log(xB);
   double log_ximax = 0.;

   // zLocal = 0 -> xi=xB: larger alpha -> closer to xB
   Setxi( exp(        zLocal  * log_ximax 
              + (1. - zLocal) * log_ximin) ); // set xi and xi_hard

   return (log_ximax - log_ximin) * xi * locJacobian;
}

// ----------------------------------------------------------------------------

// map unit variables to the splitting variable u
// parametrizes du / u (??) from xi_hard to 1
// use variable transformation to optimize MC integration
// assume xi_hard is defined (lower bound)
// stores the value of the jacobian

void Event :: MapUnitTo_u(double zu) {

   // first do a power law transformation

   double zLocal = 1. - variableTransformationLower(
                           1. - zu,  
                           1, 
                           vtp -> Get_mapParameter(9),
                           &ujacobian
                        );

   // then perform a logarithmic mapping

   double log_min = log(xi_hard);
   double log_max = 0.;

   u = exp(zLocal * log_max + (1. - zLocal) * log_min);

   ujacobian *= (log_max - log_min) * u;
}

// ----------------------------------------------------------------------------

// Add a new particle and return its pointer

Particle* Event :: AddParticle() {

   if (nparticles == maxparticles + addentries) {
      Print(stdout);
      errf(-1, "Event :: AddParticle: event record too small");
   }

   return particle[nparticles++];
}

// ----------------------------------------------------------------------------

// Add a new particle and return it's pointer; 
// set label and type

Particle* Event :: AddParticle(int labelIn, char typeIn) {

   Particle* p = AddParticle();
   p -> SetLabel(labelIn);
   p -> SetType(typeIn);
   mapping[labelIn] = p;
   return p;
}

// ----------------------------------------------------------------------------

// add leptons and initial state to the event record: e+e- annihilation

void Event :: EEAddParticles() {

   double ee = 0.5 * ecm;

   Particle* a;

   // z-axis
   a = AddParticle(-1, 'u');
   a -> SetEnergy(1.);  
   a -> SetMomentum(1.);
   a -> SetVAndSV(0., 0.);
   a -> SetSinCosPhi(0., 1.);

   // incident lepton
   a = AddParticle(-4, 'u');
   a -> SetEnergy(ee);  
   a -> SetMomentum(ee);
   a -> SetVAndSV(0., 0.);
   a -> SetSinCosPhi(0., 1.);

   // incident antilepton
   a = AddParticle(-5, 'u');
   a -> SetEnergy(ee);  
   a -> SetMomentum(ee);
   a -> SetVAndSV(1., 0.);
   a -> SetSinCosPhi(0., 1.);

   // virtual photon
   a = AddParticle(-6, 'u');
   a -> SetEnergy(ecm);  
   a -> SetMomentum(0.);
   a -> SetVAndSV(0., 0.);
   a -> SetSinCosPhi(0., 1.);
}

// ----------------------------------------------------------------------------

// add leptons and initial state to the event record: DIS

void Event :: DISAddParticles() {

   double xi_local;

   switch (frame) {
      case pCMS:
         xi_local = xi;
         break;
      case hCMS:
         xi_local = 1.;
         break;
      default:
         errf(-1, "Event :: DISAddParticles: frame not known");
         xi_local = 0.; // TAWM
   }

   double diff = xi_local - xB;
   if (diff < 0.) 
      errf(-1, "Event::DISAddParticles: diff<0");

   proton_energy = 0.5 * Q / sqrt(xB * diff);
   double el  = proton_energy * (xi_local - xB * y) / y;
   double elp = proton_energy * (xi_local * (1. - y) + xB * y) / y;

   double dl = (xi_local - xB) / (xi_local - xB * y);
   correct_range01(&dl);

   double dlp = (1. - y) * (xi_local - xB) / (xi_local * (1. - y) + xB * y);
   correct_range01(&dlp);

   Particle* a;

   // incident lepton
   a = AddParticle(-4, 'u');
   a -> SetEnergy(el);  
   a -> SetMomentum(el);
   a -> SetV(dl);
   a -> SetSinCosPhi(0., 1.);

   // outgoing lepton
   a = AddParticle(-5, 'u');
   a -> SetEnergy(elp);  
   a -> SetMomentum(elp);
   a -> SetV(dlp);
   a -> SetSinCosPhi(0., 1.);
   
   // incident proton
   a = AddParticle(-7, 'u');
   a -> SetEnergy(proton_energy);  
   a -> SetMomentum(proton_energy);
   a -> SetVAndSV(0., 0.);
   a -> SetSinCosPhi(0., 1.);
}

// ----------------------------------------------------------------------------

// add remnant to the event record
// assumes that DISAddParticles() has been called before
// assumes that proton_energy is already stored before

void Event :: DISAddRemnantAndIncident() {

   Particle* a;

   // incident parton
   a = AddParticle(-1, 'u');
   a -> SetEnergy(xi * proton_energy);  
   a -> SetMomentum(xi * proton_energy);
   a -> SetVAndSV(0., 0.);
   a -> SetSinCosPhi(0., 1.);

   // proton remnant; proton_energy already stored before!
   a = AddParticle(-8, 'u');
   a -> SetEnergy((1.-xi) * proton_energy);  
   a -> SetMomentum((1.-xi) * proton_energy);
   a -> SetVAndSV(0., 0.);
   a -> SetSinCosPhi(0., 1.);
}

// ----------------------------------------------------------------------------

// add remnant and incident parton(s) to the event record

void Event :: EcmAddRemnantAndIncident() {

   DefineEcm(); // now ``xi'' is defined

   switch (n_in_partons) {

      case 1:
         DISAddRemnantAndIncident();
         break;

      default:
         errf(-1, "Event :: AddRemnantAndIncident: illegal!");
   }
}

// ----------------------------------------------------------------------------

// set the reference frame

void Event :: SetFrame(Frame frameI) {

   frame = frameI;
}

// ----------------------------------------------------------------------------

// calculate the phase space variables
// assumes that ``mapping'' is defined

void Event :: CalculatePartonInvariantsLocal() {

   CalculatePartonInvariants(this);
}

// ----------------------------------------------------------------------------

// calculates the list of limits D / sqrt(v)

void Event :: CalculateWList() {

//   tEval -> To_status() << FS("Wlist: calculate %d\n", (int)p0ref);

   register int i, j;

   // prepare list of terms to be included
   iwentries = 0;
   switch (n_in_partons) {
      case 0: // e+e-
         iwlist[iwentries++] = -4; // incident electron
         iwlist[iwentries++] = -5; // incident positron
         break;
      case 1: // DIS
         iwlist[iwentries++] = -4; // incident electron
         iwlist[iwentries++] = -5; // outgoing electron
         if (p0ref == FS_REFERENCE)
            iwlist[iwentries++] = -1; // incident parton
         break;
      case 2: // pp
         iwlist[iwentries++] = -1; // incident parton 1
         iwlist[iwentries++] = -2; // incident parton 2
         break;
   }
   // this determines which additional final state partons are included
   int parton_start = (p0ref == FS_REFERENCE) ? 2 : 1;
   for (i = parton_start; i < npartons; ++i) {
       iwlist[iwentries++] = particle[i] -> GetLabel();
//       tEval -> To_status()
//          << FS("ADDED %d", iwlist[iwentries-1])
//          << FS(" %d\n", iwlist[iwentries-1],parton_start);
   }

   // Now calculate the cosines of the differences in azimuthal
   // angles and the polar angles; always relative to the first parton entry.

   // This is the azimuthal angle of the 0th particle relative to the
   // reference momentum as the z-axis. In the FS case, the x-axis is
   // determined by the momentum of the incident particle. In the IS case, 
   // the x-axis is freely floating.
   double sinPhi0;
   double cosPhi0;
   particle[0] -> GetSinCosPhiRef(&sinPhi0, &cosPhi0);

   for (i = 0; i < iwentries; ++i) {

//       tEval -> To_status() 
//          << FS("i=%d", i)
//          << FS(" iwlist=%d\n", iwlist[i]);

       Particle* p = mapping[iwlist[i]];

       double sinPhiRef;
       double cosPhiRef;

       if (p0ref == FS_REFERENCE) {

          // azimuthal angle with particle 1 as z-axis and incident
          // particle (DIS: parton, e+e-: lepton) as x-axis
          // :_CHECK_: for e+e-, is mapping[-1] the incident lepton??
          switch (n_in_partons) {
             case 0:
             case 1:
                p -> CalculateRelativeAzimuth(
                        particle[1],
                        mapping[-1],
                        &sinPhiRef,
                        &cosPhiRef
                     );
             break;
          }

          // the corresponding polar angle
          p -> vz = V(p, particle[1]);

       } else {

          // azimuthal angle with incident particle as z-axis and 
          // freely floating x-axis
          p -> GetSinCosPhi(&sinPhiRef, &cosPhiRef);
          // the corresponding polar angle
          p -> vz = p -> GetV();
       }

       p -> cdp = cosPhiRef * cosPhi0 + sinPhiRef * sinPhi0;
       double term = p -> vz * (1. - p -> vz);
       correct_range01(&term);
       p -> sq = sqrt(term);
       
       dtlist[iwlist[i]] = p -> GetEnergy() * 2. * p -> sq * p -> cdp;
   }

   register int k, l;
   for (i = 0; i < iwentries - 1; ++i) 
       for (j = i + 1; j < iwentries; ++j) {
           k = iwlist[i];
           l = iwlist[j];
           wlist[k][l]= 2. * (  mapping[k] -> vz 
                                 * mapping[l] -> sq 
                                 * mapping[l] -> cdp
                              - mapping[l] -> vz 
                                 * mapping[k] -> sq 
                                 * mapping[k]->cdp
                             );
           wlist[l][k] = - wlist[k][l];
       }
}

// ----------------------------------------------------------------------------

// set the kinemtical defaults before calling the user routine
// for the incident phase space

void Event :: SetKinematicalDefaults() {

   switch (n_in_partons) {
      case 0:
         break;
      case 1:
         break;
      default:
         errf(-1, "Event::SetKinematicalDefaults: not (yet) defined");
   }
}

// ----------------------------------------------------------------------------

// determines the combined vector h and the ratios lambda or eta, 
// for the final / initial state case

void Event :: DetermineHLambdaEta() {

   Particle* h;

   switch (limitType) {

      case COLLINEAR:
      case SOFT_AND_COLLINEAR:
      case COLLINEAR_ADDED:
      case SOFT_AND_COLLINEAR_ADDED:

         // the cluster:
         // h = pA + pB, if A, B final state
         // h = pI - pA, if I initial state and A final state
         if (mapping[-3] == NULL)
            h = AddParticle(-3, 'u');
         else
            // might exist if a phase space is re-used
            h = mapping[-3];

         switch (p0ref) {

            case FS_REFERENCE:

               mapping[-3] -> Add(mapping[limit1], mapping[limit2]);
               SetLambda(  mapping[limit1] -> GetEnergy()
                         / mapping[-3]     -> GetEnergy()
               );
               break;

            case IS_REFERENCE: 

               mapping[-3] -> Subtract(mapping[limit2], mapping[limit1]);
               SetEta(1. - (  mapping[limit1] -> GetEnergy()
                            / mapping[limit2] -> GetEnergy()
                           )
               );
               break;
         }
         break;

      case SOFT_ADDED:  // :_MOD_: this is too general!

         // put the initial particle into h
         if (mapping[-3] == NULL)
            h = AddParticle(-3, 'u');
         else
            // might exist if a phase space is re-used
            h = mapping[-3];

         mapping[-1] -> CopyInto(h);
         h -> SetLabel(-3);  // sic! (was overwritten)
         eta = 1.;
         
         break;
      
      default:

         // do nothing
         ;

   }  // of switch
}

// ----------------------------------------------------------------------------

// check whether the event satisfies the technical cut condition
// :_MOD_: should make the scale dependent on whether one of the partons 
//       is incident!

int Event :: CheckTechnicalCut(Process* process, char* location) {

   double ref_scale;

   switch ( process
               -> GetGlobal()
               -> disaster
               -> GetScaleChoiceForTechnicalCut() ) {
      case 0:
        ref_scale = 0.5 * Q2;
        break;
      case 1:
        ref_scale = 0.5 * W2;
        break;
      case 2:
        ref_scale = 0.5 * SH;
        break;
      default:
        errf(-1, "Event :: CheckTechnicalCut: scale choice error");
        ref_scale = 0.; // :_TAWM_:
   }

   double cut_scale = process
                         -> GetGlobal()
                         -> disaster
                         -> GetTechnicalCutParameter() 
                    * ref_scale;

   register int i;
   for (i = - n_in_partons; i < npartons; ++i)
       mapping[i] -> CalculateCartesianFromPolar();

   register int j;
   register int ret = TRUE;
   for (i = - n_in_partons; i < npartons - 1; ++i)
       for (j = i + 1; j < npartons; ++j)
           if (DotCartesianMap(i, j) < cut_scale) {
              if (disaster -> Get_print_debug() != 0) {
                 tEval -> To_status() << FS("i = %d ", i)
                                      << FS("j = %d\n", j);
              }
              ret = FALSE;
           }

   if (!ret) {
      disaster -> Set_technicalCut_flag();
      if (disaster -> Get_print_debug() != 0) {
         tEval -> To_error()
            << "*** EVENT DOES NOT FULFILL THE TECHNICAL CUT ***\n";
         tEval -> To_error()
            << FS("error = :%s:\n", location);
         Print(tEval -> To_status());
      }
   }

   return ret;
}

// ----------------------------------------------------------------------------

// SOFT_ADDED:
// create an event record for a soft added subtraction term by renumbering
// a given one, and by adding a fictitious zero-momentum particle
// Assumes that the addedEvent is already defined.
// Assumes that mapping is defined.

// COLLINEAR_ADDED:
// create an event record for a collinear added subtraction term by renumbering
// a given one, and by adding a fictitious zero-momentum particle
// (final state, jangle >= 0) or final-momentum particle 
// (initial state, jangle < 0).
// Assumes that the addedEvent is already defined.
// Assumes that mapping is defined.
// final state: particle #0 is jangle, added particle is iangle
// initial state: added particle is iangle

// SOFT_AND_COLLINEAR_ADDED:
// create an event record for a soft & collinear added subtraction term by 
// renumbering a given one, and by adding a fictitious forward particle
// Assumes that the addedEvent is already defined.
// Assumes that mapping is defined
// At present only for i.s.s. (!) :_CAUTION_:

Event* Event :: CreateAdded(
                   Process* process, 
                   LimitType lt, 
                   int iangle, 
                   int jangle
                ) {
   
   P0_Reference p0r;

   switch (lt) {
      case SOFT_ADDED:
         p0r = FS_REFERENCE;
         break;
      case COLLINEAR_ADDED:
         if (jangle >= 0)
            p0r = FS_REFERENCE;
         else
            p0r = IS_REFERENCE;
         break;
      case SOFT_AND_COLLINEAR_ADDED:
         p0r = IS_REFERENCE;  // :_MOD_: might be redundant
         if (jangle >= 0)
            errf(-1, "Event :: CreateAdded: error; not i.s.s.");
         break;
      default:
         p0r = FS_REFERENCE; // :_TAWM_:
         errf(-1, "Event :: CreateAdded: no such LimitType");
         break;
   }

   Event* limit;
   Particle* p;
   register int i;
   register int pos;

   // try to find the event in the list
   Event* newEvent = tEval 
              -> GetEFree()
              -> GetDefinedObject(tEval, GetMaxNumber() + 1);
   newEvent -> SetFrame(GetFrame());
   newEvent -> SetIn(GetIn());
   newEvent -> SetNpartons(GetPartonNumber()+1); // to be definite
   newEvent -> SetEventType(GetEventType());
   newEvent -> SetLimitType(lt);
   newEvent -> SetP0ref(p0r);
   newEvent -> SetToBeSubtracted(GetToBeSubtracted());

   Event* find = tEval -> process -> GetEventList() -> FindSame(newEvent);

   if (find != NULL) {

      // already present in the list, re-use!

      tEval
         -> GetEFree()
         -> ReturnObject(newEvent);
      limit = find;

   } else {

      // not yet present in the list, construct!

      CopyInto(newEvent);
      newEvent -> SetLimitType(lt);
      newEvent -> SetP0ref(p0r);
      newEvent -> SetLimit1(iangle);
      newEvent -> SetLimit2(jangle);
      process -> GetEventList() -> InsertObject(newEvent);      
      limit = newEvent;

      switch (lt) {

         case SOFT_ADDED:

            limit -> mappingSave[iangle] = NULL;

            // rename partons
            pos = 0;
            for (i = 0; i < limit -> npartons + 1; ++i)
                if (i != iangle) {
                   limit -> mappingSave[i] = limit -> mapping[pos];
                   limit -> mapping[pos++] -> SetLabel(i);
                }
            for (i = 0; i < limit -> npartons + 1; ++i)
                limit -> mapping[i] = limit -> mappingSave[i];

            // now add a fictitious zero-momentum parton
            p = limit -> AddParticle(iangle, 'u');
            p -> SetToZeroMomentum();
            ++ (limit -> npartons);

            limit -> xi_int  = limit -> xi;
            limit -> xi_pden = limit -> xi;
            break;

         case COLLINEAR_ADDED:

               if (jangle >= 0) {

                  // final state collinear term

                  limit -> mappingSave[iangle] = NULL;
                  limit -> mapping[0] -> SetLabel(jangle); 
                  limit -> mappingSave[jangle] = limit -> mapping[0];

                  // rename partons
                  pos = 1;
                  for (i = 0; i < limit -> npartons + 1; ++i)
                      if (i != iangle && i != jangle) {
                          limit -> mappingSave[i] = limit -> mapping[pos];
                          limit -> mapping[pos++] -> SetLabel(i);
                      }
                  for (i = 1; i < limit -> npartons + 1; ++i)
                      limit -> mapping[i] = limit -> mappingSave[i];


                  // now add a fictitious zero-momentum parton
                  // :_MOD_: should better be a splitting (if ff. are appended)
                  p = limit -> AddParticle(iangle, 'u');
                  p -> SetToZeroMomentum();
                  ++ (limit -> npartons);

                  limit -> xi_int  = limit -> xi;
                  limit -> xi_pden = limit -> xi;

               } else {
 
                  // initial state collinear term

                  // rearrange the momentum fraction variables
                  limit -> xi_hard = limit -> xi;
                  limit -> xi      = limit -> xi_hard / limit -> u;
                  limit -> xi_int  = limit -> xi_hard; 
                  limit -> xi_pden = limit -> xi;
      
                  limit -> mappingSave[iangle] = NULL;

                  // rename partons
                  pos = 0;
                  for (i = 0; i < limit -> npartons + 1; ++i)
                      if (i != iangle) {
                         limit -> mappingSave[i] = limit -> mapping[pos];
                         limit -> mapping[pos++] -> SetLabel(i);
                      }
                  for (i = 0; i < limit -> npartons + 1; ++i)
                      limit -> mapping[i] = limit -> mappingSave[i];

                  // now add a fictitious parton
                  p = limit -> AddParticle(iangle, 'u');
                  p -> SetVAndSV(   limit -> mapping[jangle]);
                  p -> SetSinCosPhi(limit -> mapping[jangle]);
                  p -> SetEnergy( (1. / limit -> u - 1.) 
                                 * limit -> mapping[jangle] -> GetEnergy()
                       );
                  p -> SetMomentum(p -> GetEnergy());
                  ++ (limit -> npartons);

                  // rescale the initial-state parton
                  limit -> mapping[jangle] -> DivideEnergyBy(  limit -> u);
                  limit -> mapping[jangle] -> DivideMomentumBy(limit -> u);

                  // rescale the remnant
                  limit 
                  -> mapping[-8] 
                  -> MultiplyEnergyBy(   
                        (1. - limit -> xi) / (1. - limit -> xi_hard)
                     );
                  limit 
                  -> mapping[-8] 
                  -> MultiplyMomentumBy( 
                        (1. - limit -> xi) / (1. - limit -> xi_hard)
                     );
               }

               break;

         case SOFT_AND_COLLINEAR_ADDED:

            // initial state collinear term
            // new particle is soft --> no rescaling of the remnant, 
            //                          mo modification of the xi-variables

            limit -> mappingSave[iangle] = NULL;

            // rename partons
            pos = 0;
            for (i = 0; i < limit -> npartons + 1; ++i)
                if (i != iangle) {
                   limit -> mappingSave[i] = limit -> mapping[pos];
                   limit -> mapping[pos++] -> SetLabel(i);
                }
            for (i = 0; i < limit -> npartons + 1; ++i)
                limit -> mapping[i] = limit -> mappingSave[i];

            // now add a fictitious zero-momentum parton
            p = limit -> AddParticle(iangle, 'u');
            p -> SetToZeroMomentum();
            ++ (limit -> npartons);

            limit -> xi_int  = limit -> xi;
            limit -> xi_pden = limit -> xi;
        
            break;

         default:
            errf(-1, "Event :: CreateAdded: no such LimitType");
            break;
      }

   }

   return limit; 
}

// ----------------------------------------------------------------------------

// Create an event record for an added subtraction term by renumbering
// a given one; it is assumed that the additional particle has
// already been included as mapping[npartons-1] (with already incremented
// npartons).
// final state: particle #0 is jangle, (already) added particle is iangle
// initial state: added particle is iangle
// Exploits that the ``Born term'' fulfilled particle[i] == mapping[i]
// for i < npartons_Born

void Event :: AddedRenumbering(int iangle, int jangle) {

   CHECK_EVENT_RECORD(this)

   register int i, pos;

   Particle* active = mapping[limit1];

   if (jangle >= 0) {

      particle[0] -> SetLabel(jangle);
      mapping[jangle] = particle[0];

      // rename partons
      pos = 1;
      for (i = 0; i < npartons; ++i)
          if (i != iangle && i != jangle) {
             particle[pos] -> SetLabel(i);
             mapping[i] = particle[pos++];
          }
 
      active -> SetLabel(iangle); 
      mapping[iangle] = active;

      SetLimit1(iangle);
      SetLimit2(jangle);

   } else {

      // rename partons
      pos = 0;
      for (i = 0; i < npartons; ++i)
          if (i != iangle) {
             particle[pos] -> SetLabel(i);
             mapping[i] = particle[pos++];
          }

      active -> SetLabel(iangle); 
      mapping[iangle] = active; 

      SetLimit1(iangle);
      SetLimit2(jangle);
   }

   CHECK_EVENT_RECORD(this)
}

// ----------------------------------------------------------------------------

// boost in z-direction such that the vector alpha * P + q == 0

void Event :: ZBoost(double alpha) {

#define CHECKINV 0
#if CHECKINV
   Event *copy = tEval 
                 -> GetEFree()
                 -> GetDefinedObject(tEval, GetMaxNumber());
   CopyInto(copy);
//   tEval -> To_status() << "old record:\n";
//   copy -> Print(stdout);
#endif

   Particle* proton
                = tEval
                  -> GetPFree()
                  -> GetDefinedObject(global -> disaster);
    
   // have to copy proton because it is modified when boosted
   mapping[-7] -> CopyInto(proton);

   Particle* virtualPhoton
                = tEval
                  -> GetPFree()
                  -> GetDefinedObject(global -> disaster);

   virtualPhoton -> LinearCombination(1.0, mapping[-4], -1.0, mapping[-5]);

   Particle* zeropart
                = tEval
                  -> GetPFree()
                  -> GetDefinedObject(global -> disaster);

   zeropart -> LinearCombination(alpha, proton, 1.0, virtualPhoton);

   double c0new = sqrt(Dot(zeropart, zeropart));   
   double oneOverC0new = 1. / c0new;
   double p0new = Dot(proton, virtualPhoton) * oneOverC0new;
   double oneOverP0new = 1. / p0new;

   proton   -> CalculateCartesianFromPolar();
   zeropart -> CalculateCartesianFromPolar();

   register int i;
   for (i = 0; i < nparticles; ++i) {
       particle[i] -> CalculateCartesianFromPolar();
       // have to use intermediate variables...
       double k0new =   DotCartesian(particle[i], zeropart)
                      * oneOverC0new;
       double kznew = k0new -   DotCartesian(particle[i], proton)
                              * oneOverP0new;
       particle[i] -> SetFvect(0, k0new);
       particle[i] -> SetFvect(3, kznew);
       // transverse components [1] and [2] are identical...
       // calculate the other coordinates
       particle[i] -> CalculatePolarFromCartesian();
   }

   tEval -> GetPFree() -> ReturnObject(proton);
   tEval -> GetPFree() -> ReturnObject(zeropart);
   tEval -> GetPFree() -> ReturnObject(virtualPhoton);

#if CHECKINV
//   tEval -> To_sttaus() << "new record:\n";
//   Print(stdout);
   register int k, l;
   double invold, invnew;
//   tEval -> To_sttaus() << "checking...\n";
   double deltaabs = 1.e-8;
   double deltarel = 1.e-8;
   for (k=0; k<nparticles; ++k)
       for (l=0; l<=k; ++l) {
           invold = Dot(copy -> particle[k], copy -> particle[l]);
           invnew = Dot(particle[k], particle[l]);
           if (!(   fabs(invold) < deltaabs && fabs(invnew) < deltaabs
                 || fabs(invold/(invold+invnew)-0.5) < deltarel
                )
              )
           tEval -> To_status() 
              << FS("(%2d", k)
              << FS(", %2d", l)
              << FS(") -> %30.22e", invold)
              << FS(" %30.22e\n", invnew);
       }
   tEval -> GetEFree() -> ReturnObject(copy);
#endif
}

// ----------------------------------------------------------------------------

// boost to the hCMS

void Event :: BoostTo_hCMS() {

   ZBoost(1.0);   
}

// ----------------------------------------------------------------------------

// boost to the Breit system

void Event :: BoostTo_Breit() {

   ZBoost(2. * xB);   
}

// ----------------------------------------------------------------------------

// boost to the laboratory system (DIS)

void Event :: BoostTo_DIS_Lab(double ratio_Ep_Ee) {

   SetUserFrame(lab);

   // determine the boost parameters
   Particle* vE = tEval
                  -> GetPFree() 
                  -> GetDefinedObject(global -> disaster);
   Particle* vZ = tEval  
                  -> GetPFree() 
                  -> GetDefinedObject(global -> disaster);
   Particle* vX = tEval  
                  -> GetPFree() 
                  -> GetDefinedObject(global -> disaster);
      
   vE -> LinearCombination( 
            1.0, 
            GetMapping(-7),      // incident proton
            ratio_Ep_Ee,
            GetMapping(-4)       // incident lepton
         );
   vZ -> LinearCombination(
            1.0, 
            GetMapping(-7), 
            - ratio_Ep_Ee,
            GetMapping(-4) 
         );
   *vX = *(GetMapping(-5));  // the outgoing lepton
   
   BoostEZX(vE, vZ, vX);
   
   tEval -> GetPFree() -> ReturnObject(vE);
   tEval -> GetPFree() -> ReturnObject(vZ);
   tEval -> GetPFree() -> ReturnObject(vX);
}

// ----------------------------------------------------------------------------

// calculate the current thrust

double Event :: CurrentThrust() {

   register int i, k;
   double sum = 0;
   k = 0;
   for (i = 0; i < npartons; ++i) {
       if (mapping[i] -> GetV() > 0.5) {
          sum += mapping[i] -> GetMomentum() * (2 * mapping[i] -> GetV() - 1.);
          ++k;
       }
   }
   sum *= 2. / Q;

   if (k == 0)  // definition a la B. Webber
      sum = 1.;

   return sum;
}

// ----------------------------------------------------------------------------

// Lorentz boost. 
// New coordinates defined by unit vectors.

void Event :: BoostWithCartesianVector(
                 Particle* nE, 
                 Particle* nX, 
                 Particle* nY, 
                 Particle* nZ
              ) {

   for (register int i = 0; i < nparticles; ++i) {
       particle[i] -> CalculateCartesianFromPolar();
       particle[i] -> NewComponentsCartesian(particle[i], nE, nX, nY, nZ);
       particle[i] -> CalculatePolarFromCartesian();
   }
}

// ----------------------------------------------------------------------------

// boost such that vE points into the energy direction and
// vZ points into the z direction.
// (==> required that vE is timelike and vZ is spacelike!)
// Moreover the x-axis is fixed by the component of vX transverse to vZ.
// The y-axis is fixed such that x, y, z form a right-handed system.

void Event :: BoostEZX(
                 Particle* vE, 
                 Particle* vZ, 
                 Particle* vX
              ) {

#undef CHECKINV
#define CHECKINV 0
#if CHECKINV
   Event* copy
             = tEval
               -> GetEFree()
               -> GetDefinedObject(tEval, GetMaxNumber());
   CopyInto(copy);
#endif

   Particle *nE = tEval -> GetPFree() -> GetDefinedObject(tEval);
   Particle *nX = tEval -> GetPFree() -> GetDefinedObject(tEval);
   Particle *nY = tEval -> GetPFree() -> GetDefinedObject(tEval);
   Particle *nZ = tEval -> GetPFree() -> GetDefinedObject(tEval);

   double scale;

   nE -> SetVAndSV(vE);
   nE -> SetSinCosPhi(vE);
   scale = 1 / vE -> CalculateMass();
   nE -> SetEnergy(  scale * vE -> GetEnergy());
   nE -> SetMomentum(scale * vE -> GetMomentum());

   nZ -> SetVAndSV(vZ);
   nZ -> SetSinCosPhi(vZ);
   scale = 1 / vZ -> CalculateImaginaryMass();
   nZ -> SetEnergy(  scale * vZ -> GetEnergy());
   nZ -> SetMomentum(scale * vZ -> GetMomentum());
 
   nE -> CalculateCartesianFromPolar();
   vX -> CalculateCartesianFromPolar();
   nX -> LinearCombinationCartesian(1.0, vX, - DotCartesian(vX, nE), nE);
   
   nZ -> CalculateCartesianFromPolar();
   nX -> LinearCombinationCartesian(1.0, nX, DotCartesian(nX, nZ), nZ);

   nX -> ScaleByCartesian(1 / nX -> CalculateImaginaryMassCartesian());
   
   nY -> EpsilonCartesian(nE, nZ, nX);

   BoostWithCartesianVector(nE, nX, nY, nZ);

#if 0
   Particle *test = tEval -> GetPFree() -> GetDefinedObject(tEval);

   test -> NewComponentsCartesian(nE, nE, nX, nY, nZ);
   test -> Print(stdout);
   tEval -> To_status() << "\n";

   test -> NewComponentsCartesian(nX, nE, nX, nY, nZ);
   test -> Print(stdout);
   tEval -> To_status() << "\n";

   test -> NewComponentsCartesian(nY, nE, nX, nY, nZ);
   test -> Print(stdout);
   tEval -> To_status() << "\n";

   test -> NewComponentsCartesian(nZ, nE, nX, nY, nZ);
   test -> Print(stdout);
   tEval -> To_status() << "\n";

   tEval -> GetPFree() -> ReturnObject(test);
#endif

   tEval -> GetPFree() -> ReturnObject(nE);
   tEval -> GetPFree() -> ReturnObject(nX);
   tEval -> GetPFree() -> ReturnObject(nY);
   tEval -> GetPFree() -> ReturnObject(nZ);

#if CHECKINV
//   tEval -> To_status() << "new record:\n";
//   Print(stdout);
   register int k, l;
   double invold, invnew;
//   tEval -> To_status() << "checking...\n";
   double deltaabs = 1.e-8;
   double deltarel = 1.e-8;
   for (k=0; k<nparticles; ++k)
       for (l=0; l<=k; ++l) {
           invold = Dot(copy -> particle[k], copy -> particle[l]);
           invnew = Dot(particle[k], particle[l]);
           if (!(   fabs(invold) < deltaabs && fabs(invnew) < deltaabs
                 || fabs(invold/(invold+invnew)-0.5) < deltarel
                )
              )
           tEval -> To_status() 
              << FS("(%2d", k)
              << FS(", %2d", l)
              << FS(") -> %30.22e", invold)
              << FS(" %30.22e\n", invnew);
       }
   tEval -> GetEFree() -> ReturnObject(copy);
#endif
}

// ----------------------------------------------------------------------------

int Event :: Cut_Jets_pT_PseudoRapidity(
                double pT_Min,
                double pT_Max,
                double eta_Min,
                double eta_Max
             ) {
   
   register int i;
   register int fulfilled = 1;
   for (i = 0; i < npartons; ++i) {
       if ( ! inRange(mapping[i] -> Calculate_pT(), 
                      pT_Min, pT_Max))
          fulfilled = 0;
       if ( ! inRange(mapping[i] -> Calculate_PseudoRapidity(), 
                      eta_Min, eta_Max))
          fulfilled = 0;
   }
   return fulfilled;
}

// ----------------------------------------------------------------------------

int Event :: Cut_Lepton_pT_PseudoRapidity(
                double pT_Min,
                double pT_Max,
                double eta_Min,
                double eta_Max
             ) {
   
   register int fulfilled = 1;
   if ( ! inRange(mapping[-5] -> Calculate_pT(), 
                  pT_Min, pT_Max))
      fulfilled = 0;
   if ( ! inRange(mapping[-5] -> Calculate_PseudoRapidity(), 
                  eta_Min, eta_Max))
      fulfilled = 0;

   return fulfilled;
}

// ----------------------------------------------------------------------------

int Event :: UseListIsSame(Event* other) {

   return (
         GetNpartons() == other -> GetNpartons()
      && GetN_in_partons() == other -> GetN_in_partons()
      && GetFrame() == other -> GetFrame()
      && GetLimitType() == other -> GetLimitType()
      && GetP0ref() == other -> GetP0ref()
      && GetToBeSubtracted() == other -> GetToBeSubtracted()
      && GetEventType() == other -> GetEventType()
   );
}

// ----------------------------------------------------------------------------

Event* Event :: UserFrameFindEvent(Frame f) {

   HorizontalList <Event> * h;

   h = GetH();
   while (h != NULL && h -> GetObject() -> GetUserFrame() != f) {
      h = h -> GetHorizontalPtr();
   }

   return (h == NULL) ? (Event*) NULL : h -> GetObject();
}

// ----------------------------------------------------------------------------

void Event :: PrestoreGraph(
                 Book* book,
                 double x,
                 int entry
              ) {

   Contribution* c;
   int iev;

   for (c = GetContribution();  
        c != NULL && (iev = c -> GetLocation(), TRUE);
        c = c -> GetSameEvent()) {
       book -> PrestoreGraph(x, iev, entry);
   }
}

// ----------------------------------------------------------------------------

void Event :: PrestoreHistogram(
                 Book* book,
                 double x
              ) {

   Contribution* c;
   int iev;

   for (c = GetContribution();  
        c != NULL && (iev = c -> GetLocation(), TRUE);
        c = c -> GetSameEvent()) {
       book -> PrestoreHistogram(x, iev);
   }
}

// ----------------------------------------------------------------------------

void Event :: PrestoreHistogram(
                 Book* book,
                 double x, 
                 double weight
              ) {

   Contribution* c;
   int iev;

   for (c = GetContribution();  
        c != NULL && (iev = c -> GetLocation(), TRUE);
        c = c -> GetSameEvent()) {
       book -> PrestoreHistogram(x, iev, weight);
   }
}

// ----------------------------------------------------------------------------

double Event :: DotMap(int i, int j) {

   return Dot(mapping[i], mapping[j]);
}

// ----------------------------------------------------------------------------

double Event :: DotCartesianMap(int i, int j) {

   return DotCartesian(mapping[i], mapping[j]);
}

// ----------------------------------------------------------------------------

double Event :: InvariantMassCartesianMap(int i, int j) {

   return InvariantMassCartesian(mapping[i], mapping[j]);
}

// ----------------------------------------------------------------------------

double Event :: VMap(int i, int j) {

   return V(mapping[i], mapping[j]);
}

// ----------------------------------------------------------------------------

double Event :: VNormalMap(int i, int j) {

   return VNormal(mapping[i], mapping[j]);
}

// ----------------------------------------------------------------------------

void Event :: DotAndVMasslessNormalMap(int i, int j, double* d, double* v) {

   DotAndVMasslessNormal(mapping[i], mapping[j], d, v);
}

// ----------------------------------------------------------------------------

void Event :: InvariantAndVMasslessNormalMap(
                 int i, 
                 int j, 
                 double* d, 
                 double* v) {

   InvariantAndVMasslessNormal(mapping[i], mapping[j], d, v);
}

// ----------------------------------------------------------------------------

void Event :: CalculateNormalFromPolarRangeMap(int imin, int imax) {

   for (register int i = imin; i <= imax; ++i)
       mapping[i] -> CalculateNormalFromPolar();
}

// ----------------------------------------------------------------------------

void Event :: CalculateCartesianFromPolarRangeMap(int imin, int imax) {

   for (register int i = imin; i <= imax; ++i)
       mapping[i] -> CalculateCartesianFromPolar();
}

// ----------------------------------------------------------------------------

void Event :: CalculateCartesianAndNormalFromPolarRangeMap(int imin, int imax) {

   for (register int i = imin; i <= imax; ++i)
       mapping[i] -> CalculateCartesianAndNormalFromPolar();
}

// ----------------------------------------------------------------------------

void Event :: CalculatePolarFromCartesianRangeMap(int imin, int imax) {

   for (register int i = imin; i <= imax; ++i)
       mapping[i] -> CalculatePolarFromCartesian();
}

// ----------------------------------------------------------------------------
//
// --> Variable transformation parameters
//
// preliminary: define a size of 10 (should be large enough), 
//              should resize automatically!
//
// ----------------------------------------------------------------------------

VariableTransformationParameters :: VariableTransformationParameters() {
   
   offset = new Array_int(0, 20);
   offsetData = offset -> GetDataPtr();
   n_offset = 0;

   mapFunction   = new Array_int(0, 20);
   mapFunctionData = mapFunction -> GetDataPtr();
   n_mapFunction = 0;
 
   mapParameter = new Array_double(0, 20);
   mapParameterData = mapParameter -> GetDataPtr();
   n_mapParameter = 0;
}

// ----------------------------------------------------------------------------

VariableTransformationParameters :: ~VariableTransformationParameters() {
   
   delete offset;
   delete mapFunction;
   delete mapParameter;
}

// ----------------------------------------------------------------------------

void VariableTransformationParameters 
        :: CopyInto(VariableTransformationParameters* into) {

   offset -> CopyInto(into -> offset);
   mapFunction -> CopyInto(into -> mapFunction);
   mapParameter -> CopyInto(into -> mapParameter);
}

// ----------------------------------------------------------------------------

void VariableTransformationParameters :: Set_offset(int pos, int in) {

   if (pos >= offset -> GetRmin() && pos <= offset -> GetRmax()) {
      offset -> SetData(pos, in);
      n_offset = max(n_offset, pos + 1);
   }
}

// ----------------------------------------------------------------------------

void VariableTransformationParameters :: Set_mapFunction(int pos, int in) {

   if (pos >= mapFunction -> GetRmin() && pos <= mapFunction -> GetRmax()) {
      mapFunction -> SetData(pos, in);
      n_mapFunction = max(n_mapFunction, pos + 1);
   }
}

// ----------------------------------------------------------------------------

void VariableTransformationParameters :: Set_mapParameter(int pos, double in) {

   if (pos >= mapParameter -> GetRmin() && pos <= mapParameter -> GetRmax()) {
      mapParameter -> SetData(pos, in);
      n_mapParameter = max(n_mapParameter, pos + 1);
   }
}

// ----------------------------------------------------------------------------

void VariableTransformationParameters :: PrintParameters(Mstream &s_out) {

   s_out << "variable transformation parameters:\n";

   s_out << "offset:\n";
   offset -> Print(s_out, 0, n_offset - 1);

   s_out << "mapFunction:\n";
   mapFunction -> Print(s_out, 0, n_mapFunction - 1);

   s_out << "mapParameter:\n";
   mapParameter -> Print(s_out, 0, n_mapParameter - 1);
}

// ----------------------------------------------------------------------------
//
// --> Base class for matrix elements
//
// ----------------------------------------------------------------------------

MatrixElement :: MatrixElement() {

   global = NULL;
   name = new String(-1);
   particleAssignment = new String(-1);
   Set_nIntVar(0);
}

// ----------------------------------------------------------------------------

MatrixElement :: ~MatrixElement() {

   delete name;
   delete particleAssignment;
}

// ----------------------------------------------------------------------------

void MatrixElement :: AssignGlobal(Global* globalIn) {

   global = globalIn;
}

// ----------------------------------------------------------------------------

void MatrixElement :: AssignProcess(Process* processIn) {

   process = processIn;
}

// ----------------------------------------------------------------------------

void MatrixElement :: Define() {

   // don't do anything here
}

// ----------------------------------------------------------------------------

Process* MatrixElement :: GetProcess() {

   return process;
}

// ----------------------------------------------------------------------------

String* MatrixElement :: GetName() {

   return name;
}

// ----------------------------------------------------------------------------

char* MatrixElement :: GetParticleAssignment() {

   return particleAssignment -> GetData();
}

// ----------------------------------------------------------------------------
//
// --> Base class for processes
//
// ----------------------------------------------------------------------------

Process :: Process() {

   global = NULL;
   me_list = NULL;

   alpha_s_Order = -1;

   eventList = NULL;
}

// ----------------------------------------------------------------------------

Process :: ~Process() {

   Delete();
}

// ----------------------------------------------------------------------------

void Process :: Delete() {

   if (global != NULL) { // ... if it has been Define()'d
      DeleteMatrixElements();
      global = NULL;
   }

   if (eventList != NULL)   
      delete eventList;
}

// ----------------------------------------------------------------------------

void Process :: DefineProcess(Global* globalIn) {

   global = globalIn;
   disaster = global -> disaster;

   eventList = new EventUseList(global);
   eventList -> DefineName("eventList");
}

// ----------------------------------------------------------------------------

void Process :: DeleteMatrixElements() {

   if (me_list != NULL) { // ... if it has been Define()'d

      // note that this is safe if n_components == 0
      register int ic;
      register int ime;
      for (ic = 0; ic < n_components; ++ic)
          for (ime = 0; ime < n_me[ic]; ++ime)
              delete me_list[ic][ime];  
           
      for (ic = 0; ic < n_components; ++ic) {
          delete [] me_list[ic];
      }
    
      delete [] me_list;
      delete [] to_be_subtracted;
      delete [] n_me;

      delete [] componentFlags;

      delete [] vtp;

      me_list = NULL;
   }
}

// ----------------------------------------------------------------------------

void Process :: CreateNME(int n_componentsI) {

   register int i;

   n_components = n_componentsI;
   n_me = new int [n_components];
   to_be_subtracted = new SubtractionType [n_components];

   vtp = new VariableTransformationParameters [n_components];

   componentFlags = new int [n_components];

   for (i = 0; i < n_components; ++i) {

       // set defaults

       (vtp + i) -> Set_offset(0, 0);  // lepton variable offset
       (vtp + i) -> Set_offset(1, 0);  // xi offset
       (vtp + i) -> Set_offset(2, 0);  // u offset
       (vtp + i) -> Set_offset(3, 0);  // parton phase space offset

       (vtp + i) -> Set_mapFunction(0, global 
                                       -> disaster 
                                       -> Get_xB_y_parametrization()
                    );

       (vtp + i) -> Set_mapParameter(0, global -> disaster -> Get_e_gamma());
       (vtp + i) -> Set_mapParameter(1, global -> disaster -> Get_v_delta());
       (vtp + i) -> Set_mapParameter(2, global -> disaster -> Get_e_gamma1());
       (vtp + i) -> Set_mapParameter(3, global -> disaster -> Get_v_delta1());
       (vtp + i) -> Set_mapParameter(
                       4, 
                       global -> disaster -> Get_xB_y_trho()
                    );
       (vtp + i) -> Set_mapParameter(5, global 
                                        -> disaster 
                                        -> Get_xB_y_tsigma()
                    );
       (vtp + i) -> Set_mapParameter(6, global 
                                        -> disaster 
                                        -> Get_e_gamma_upper()
                    );
       (vtp + i) -> Set_mapParameter( 7, global -> disaster -> Get_e_lambda());
       (vtp + i) -> Set_mapParameter( 8, global -> disaster -> Get_e_mu());
       (vtp + i) -> Set_mapParameter( 9, global -> disaster -> Get_u_alpha());
       (vtp + i) -> Set_mapParameter(
                       10, 
                       global -> disaster -> Get_v_delta1_upper()
                    );
       (vtp + i) -> Set_mapParameter(
                       11, 
                       global -> disaster -> Get_v_delta1_lambda()
                    );
       (vtp + i) -> Set_mapParameter(
                       12, 
                       global -> disaster -> Get_v_delta1_mu()
                    );
       (vtp + i) -> Set_mapParameter(
                       13, 
                       global -> disaster -> Get_v_delta_upper()
                    );
       (vtp + i) -> Set_mapParameter(
                       14, 
                       global -> disaster -> Get_v_delta_lambda()
                    );
       (vtp + i) -> Set_mapParameter(
                       15, 
                       global -> disaster -> Get_v_delta_mu()
                    );
       (vtp + i) -> Set_mapParameter(
                       16, 
                       global -> disaster -> Get_xi_alpha()
                    );

       componentFlags[i] = 1;
   }   
}

// ----------------------------------------------------------------------------

void Process :: CreateOtherArrays() {

   // print variable transformations
   vtp -> PrintParameters(global -> To_status());

   // create the variables holding the pointers to the matrix elements
   me_list = new MatrixElement** [n_components];
   for (register int ic = 0; ic < n_components; ++ic) {
       me_list[ic] = new MatrixElement* [n_me[ic]];
   }
}

// ----------------------------------------------------------------------------

// assign matrix elements;
// check for consistency

void Process :: AssignInformationToMatrixElements() {

   const char* sub_type[] = {"unsubtracted", 
                             "subtracted", 
                             "added subtraction",
                             "collinear initial-state radiation",
                             "error"
                            };

   register int ic;
   register int ime;

   for (ic = 0; ic < n_components; ++ic) {
       int lst;
       switch (to_be_subtracted[ic]) {
          case NoSubtraction: 
             lst = 0;
             break;
          case Subtraction: 
             lst = 1;
             break;
          case AddedSubtraction: 
             lst = 2;
             break;
          case CollinearInitial: 
             lst = 3;
             break;
          default:
             errf(-1, "Process :: AssignInformationToMatrixElements: "
                      "not known");
             lst = 4;  // :_TAWM_:
       }
       global -> To_status() 
          << FS("Component #%3d", ic)
          << FS(": (%s)\n", sub_type[lst]);
       global -> To_status() 
          << FS("   leptonVariableOffset:   %3d\n", 
                (vtp + ic) -> Get_offset(0)
             );
       global -> To_status() 
          << FS("   xiOffset:               %3d\n", 
                (vtp + ic) -> Get_offset(1)
             );
       global -> To_status() 
          << FS("   uOffset:                %3d\n", 
                (vtp + ic) -> Get_offset(2)
             );
       global -> To_status() 
          << FS("   partonPhaseSpaceOffset: %3d\n", 
                (vtp + ic) -> Get_offset(3)
             );
       for (ime = 0; ime < n_me[ic]; ++ime) {
           me_list[ic][ime] -> AssignGlobal(global);
           me_list[ic][ime] -> AssignProcess(this); 
           me_list[ic][ime] -> Define(); 
           global -> To_status() 
              << FS("   %-30s", me_list[ic][ime] -> GetName() -> GetData())
              << FS("  [ %10s", me_list[ic][ime] -> GetParticleAssignment())
              << FS(" ]  (in: %2d", me_list[ic][ime] -> GetIn())
              << FS(", out: %2d", me_list[ic][ime] -> GetOut())
              << FS(", var: %2d)\n", me_list[ic][ime] -> Get_nIntVar()
                  );
           if (me_list[ic][ime] -> GetIn() != me_list[ic][0] -> GetIn())
              errf(-1,"number of incident partons inconsistent");
           if (me_list[ic][ime] -> GetOut() != me_list[ic][0] -> GetOut())
              errf(-1,"number of outgoing partons inconsistent");
           if (   me_list[ic][ime] -> Get_nIntVar() 
               != me_list[ic][0]   -> Get_nIntVar())
              errf(-1, "number of integration variables inconsistent");
       }
       if (me_list[ic][0] -> GetIn() != me_list[0][0] -> GetIn())
          errf(-1, "number of incident partons inconsistent");
   }
   global -> To_status() << "\n";

   SetMinMaxOut();
   Set_Max_nIntVar();
}

// ----------------------------------------------------------------------------

int Process :: GetIn() {

   return me_list[0][0] -> GetIn();
}

// ----------------------------------------------------------------------------

// returns the number of final state partons of a particular component

int Process :: GetOut(int ic) {

   return me_list[ic][0] -> GetOut();
}

// ----------------------------------------------------------------------------

// sets the minimum and maximum number of final state partons

void Process :: SetMinMaxOut() {

   maxout = 0;
   minout = 9999;

   register int ic;
   for (ic = 0; ic < n_components; ++ic) {
       maxout = max(maxout, GetOut(ic));
       minout = min(minout, GetOut(ic));
   }
}

// ----------------------------------------------------------------------------

// returns the maximum number of final state partons

int Process :: GetMaxOut() {

   return maxout;
}

// ----------------------------------------------------------------------------

int Process :: GetNComponents() {

   return n_components;
}

// ----------------------------------------------------------------------------

int Process :: GetNMe(int ic) {

   return n_me[ic];
}

// ----------------------------------------------------------------------------

int Process :: Get_nIntVar(int ic) {

   return me_list[ic][0] -> Get_nIntVar();
}

// ----------------------------------------------------------------------------

// sets the maximum number of integration variables

void Process :: Set_Max_nIntVar() {

   max_nIntVar = -999999;  // nIntVar may be -1 (for q->q, for example)
   register int ic;
   for (ic = 0; ic < n_components; ++ic)
       max_nIntVar = max(max_nIntVar, Get_nIntVar(ic));
}

// ----------------------------------------------------------------------------

// returns the maximum number of integration variables

int Process :: Get_Max_nIntVar() {

   return max_nIntVar;
}

// ----------------------------------------------------------------------------

MatrixElement* Process :: GetPMatrixElement(int ic, int ime) {

   return me_list[ic][ime];
}

// ----------------------------------------------------------------------------

SubtractionType Process :: GetToBeSubtracted(int ic) {

   return to_be_subtracted[ic];
}

// ----------------------------------------------------------------------------

// sum over all matrix elements of a particular component

void Process :: Evaluate(
                   int flag, 
                   int ic, 
                   LimitType limittypeIn, 
                   int ifixIn, 
                   int jfixIn, 
                   Event* eIn,
                   double factor,
                   int lstore
                ) {

   // excluded?
   if (!flag)
      return;

   int ifix, jfix;
   LimitType limittype;

   limittype = limittypeIn;

   Event* newE = NULL;
   Event* e;

   if (disaster -> Get_ps_permute_n()[eIn -> npartons] > 0) {

      newE = disaster
                -> GetEFree()
                -> GetDefinedObject(disaster, eIn -> npartons);
      eIn -> CopyInto(newE);
      newE -> PermuteLabels(newE->npartons, 
                            disaster -> Get_ps_permute()[newE -> npartons]);

      if (limittype != NO_LIMIT) {
         if (ifixIn >= 0) {
            ifix = disaster -> Get_ps_permute()[newE -> npartons][ifixIn];
            newE -> SetLimit1(ifix);
         } else
            ifix = ifixIn;
         if (jfixIn >= 0) {
            jfix = disaster -> Get_ps_permute()[newE -> npartons][jfixIn];
            newE -> SetLimit2(jfix);
         } else
            jfix=jfixIn;
      } else {
         // :_TAWM_:
         ifix = -3;
         jfix = -3;
      }

      e = newE;

   } else {

      e = eIn;
      ifix = ifixIn;
      jfix = jfixIn;

   }

   e -> DetermineHLambdaEta();

   // calculate the invariants, also azimuthal-angle-dependent terms if required
   e -> CalculatePartonInvariantsLocal();

   int ime;
   for (ime = 0; ime < n_me[ic]; ++ime)
       me_list[ic][ime]
       -> Evaluate(
             e, 
             ca -> GetContribution(lstore), 
             factor
          );

   if (newE != NULL)
      global -> disaster -> GetEFree() -> ReturnObject(newE);
}

// ----------------------------------------------------------------------------

// calculate the weight for a given set of random numbers / parameters

void Process :: CalculateWeight(
                   double* unit,
                   ContributionArray* caIn
                ) {

   ca = caIn;

   int tcpass;

   double ustore;
   ustore = 0; // :_TAWM_:

   Event* find;

   // do sum over all contributing components

   for (register int ic = 0; ic < n_components; ++ic) {
   if (Get_componentFlag(ic)) {

       int varpos = 0;
       double xjacobian; 
 
       Event* additional
          = global
            -> disaster
            -> GetEFree()
            -> GetDefinedObject(global -> disaster, 0);

       // set the discrete parameters and check whether this 
       // particular event record is already available
       // choice of frame: hCMS only if DIS or pp, and <= two particles 
       // :CAUTION: if only added terms, this should be hCMS!
       additional -> SetFrame(frame);
       additional -> SetNpartons(0); // to be definite
       additional -> SetIn(GetIn());
       additional -> SetEventType(Additional);
       additional -> SetP0ref(FS_REFERENCE); // to be definite

       switch (to_be_subtracted[ic]) {

          case AddedSubtraction:
          case CollinearInitial:
             additional -> SetToBeSubtracted(Collinear);
             break;

          default:
             additional -> SetToBeSubtracted(to_be_subtracted[ic]);
             break;
      }

      find = GetEventList() -> FindSame(additional);

      if (find != NULL) {

         // already present in the list, re-use!
         global
            -> disaster
            -> GetEFree()
            -> ReturnObject(additional);
         additional = find;
         xjacobian = additional -> GetXJacobian();
         varpos = additional -> GetVarpos();
         ustore = additional -> GetUstore(); // this should be conditional!

      } else {

         // not yet present in the list, construct!
         GetEventList() -> InsertObject(additional);      

         // assign the variable transformation parameters
         additional -> Set_vtp(Get_vtp(ic));

         additional -> SetEcmFull(disaster -> Get_ECMFULL());

         xjacobian = 1.0;

         switch (additional -> n_in_partons) {

            case 0: // *** e+e- annihilation ***
               additional -> DefineEcm();
               break;

            case 1: // *** DIS ep scattering ***

               switch (disaster -> Get_lepton_integration()) {
                  case 0:
                     additional -> SetLeptonVariables(
                                      disaster -> Get_xB_fixed(), 
                                      disaster -> Get_y_fixed()
                                   );
                     break;
                  case 1:
                     xjacobian *= additional 
                                  -> MapUnitToxBy(
                                        unit[  varpos 
                                             + Get_vtp(ic) -> Get_offset(0)
                                        ], 
                                        unit[  1 
                                             + varpos
                                             + Get_vtp(ic) -> Get_offset(0)
                                        ]
                                     );
                     varpos += 2;
                     break;
                  default:
                     errf(-1,"Process :: CalculateWeight:"
                             " no such lepton_integration");
               }

               xjacobian *= additional -> LeptonJacobian();

               // if variables in the pCMS, set xi here (else later)
               if (additional -> GetFrame() == pCMS) {
                  switch (disaster -> Get_xi_integration()) {
                     case 0:
                        additional -> Setxi(disaster -> Get_xi_fixed());
                        break;
                     case 1:
                        xjacobian *= additional 
                                     -> MapUnitToXi(
                                           unit[  varpos++
                                                + Get_vtp(ic) -> Get_offset(1)
                                           ]
                                        );
                        break;
                     default:
                        errf(-1,"Process :: CalculateWeight:"
                                " no such xi_integration");
                  }
                  additional -> DefineEcm();
               }

               if (
                      to_be_subtracted[ic] == AddedSubtraction
                   || to_be_subtracted[ic] == CollinearInitial) {
                  switch(additional -> GetFrame()) {
                     case pCMS:
                        // set u-variable: 
                        //:_CAUTION_: have to supply a jacobian here!
                        additional -> MapUnitTo_u(
                                         unit[  varpos++
                                              + Get_vtp(ic) -> Get_offset(2)
                                         ]
                                      );
                        break;
                     case hCMS:
                        // store value for u-variable 
                        //:_CAUTION_: have to supply a jacobian here!
                        ustore = unit[  varpos++
                                      + Get_vtp(ic) -> Get_offset(2)
                                 ];
                        break;
                     default:
                        errf(-1, "Process :: CalculateWeight: frame not known");
                  }
               }
               break;

            default:
               errf(-1, "Process :: CalculateWeight: "
                        "process not yet implemented");
               varpos = 0; // :_TAWM_:
         }

         switch (additional -> n_in_partons) {
            case 0:
               additional -> EEAddParticles();
               break;
            case 1:
               additional -> DISAddParticles();
               if (additional -> GetFrame() == pCMS)
                  additional -> DISAddRemnantAndIncident();
               break;  
            default:
               errf(-1, "Process :: CalculateWeight: "
                        "process not yet implemented");
         }

         additional -> SetKinematicalDefaults();

         additional -> SetXJacobian(xjacobian);
         additional -> SetVarpos(varpos);
         additional -> SetUstore(ustore);

         additional -> SetIsConstructed(TRUE);
      }

      int flag;
      global -> user -> IncidentPhaseSpace(this, additional, flag);

      if (!flag) {
         global -> To_status() 
            << "Process :: CalculateWeight: "
            << "initial state phase space rejected\n";
         return;
      }
   
      double fs_jacobian, is_jacobian;

      // create outgoing parton variables
      Event *final 
               = global 
                 -> disaster 
                 -> GetEFree()
                 -> GetDefinedObject(
                       global -> disaster, 
                       GetOut(ic)
                    );

      additional -> CopyPSInformationInto(final);
      final -> SetPartonNumber(GetOut(ic));
      final -> SetEventType(Final);

      find = GetEventList() -> FindSame(final);
      if (find != NULL) {
         // already present in the list, re-use!
         global
            -> disaster
            -> GetEFree()
            -> ReturnObject(final);
         final = find;
         fs_jacobian = final -> GetFSJacobian();
      } else {
         // not yet present in the list, construct!
         GetEventList() -> InsertObject(final);      

// :_CAUTION_: should be 1!
#if 1 
         fs_jacobian = final 
                       -> MapUnitToPSChoice(
                             FS_REFERENCE, 
                             unit 
                                + varpos 
                                + Get_vtp(ic) -> Get_offset(3)
                          );
#else
         fs_jacobian = final -> MapUnitToPSChoice(IS_REFERENCE, unit + varpos
                                           + Get_vtp(ic) -> Get_offset(3)
#endif

         additional -> CopyAdditionalParticlesInto(final);

         final -> SetFSJacobian(fs_jacobian);

         final -> SetIsConstructed(TRUE);
      }

      if (disaster -> Get_print_debug()) {
         global -> To_status() << "Process :: CalculateWeight: final\n";
         final -> Print(global -> To_status());
      }
                 
      switch (GetToBeSubtracted(ic)) {

         case Subtraction:
        
            // this term has to be subtracted
         
            switch (global -> subtraction_type) {
   
               case 0:
   
                  final -> DoSumOverSubtractionTerms(
                              this, 
                              ic,
                              xjacobian * fs_jacobian,
                              ca
                           );
                  break;
   
               case 1:
   
                  // sum over final state singular regions
                  final -> DoSumOverDiffSubtractionTerms(
                              this, 
                              ic,
                              xjacobian * fs_jacobian,
                              ca
                           );
                  if (final -> n_in_partons > 0) {
                     // create event for sum over initial state singular regions
                     Event* initial
                          = global
                            -> disaster
                            -> GetEFree()
                            -> GetDefinedObject(global -> disaster,
                                                final -> GetMaxNumber());
                     initial -> SetPartonNumber(final -> GetPartonNumber());
                     additional -> CopyPSInformationInto(initial);
                     is_jacobian = initial  
                                   -> MapUnitToPSChoice(
                                         IS_REFERENCE,
                                         unit 
                                            + varpos 
                                            + Get_vtp(ic) -> Get_offset(3)
                                      );
                     additional -> CopyAdditionalParticlesInto(initial);

                     if (disaster -> Get_print_debug()) {
                        global -> To_status() 
                           << "Process :: CalculateWeight: initial\n";
                        initial -> Print(global -> To_status());
                     }
                 
                     // sum over initial state singular regions
                     initial -> DoSumOverDiffSubtractionTerms(
                                   this,
                                   ic,
                                   xjacobian * is_jacobian,
                                   ca
                                );
                     global
                        -> disaster
                        -> GetEFree()
                        -> ReturnObject(initial);
                  }
                  break;

               default:
                  errf(-1, "Process :: CalculateWeight: "
                           "no such subtraction_type");
            }
            break; 

         case NoSubtraction:

            // this term is unsubtracted

            // add this if necessary
            if (   final -> GetFrame() == hCMS 
                && final -> GetMapping(-1) == NULL) {
               final -> EcmAddRemnantAndIncident();
            }

            // check the technical cut
            tcpass = final -> CheckTechnicalCut(
                                 this, 
                                 "Process :: CalculateWeight (NoSubtraction)"
                              );
            if (tcpass) {
               global -> user -> PhaseSpace(
                                    this, 
                                    final,
                                    flag,
                                    ca -> NextContribution()
                                 );
               Evaluate(flag, ic, NO_LIMIT, -1, -2, final,
                        xjacobian * fs_jacobian,
                        ca -> GetActual());
            }
            break;
  
         case AddedSubtraction:

            // this term is an added subtraction

            if (   final -> GetFrame()==hCMS
                && final -> GetMapping(-1) == NULL) {
               final -> MapUnitTo_u(ustore);
               final -> EcmAddRemnantAndIncident();
            }

            // check the technical cut
            tcpass = final -> CheckTechnicalCut(
                                 this,
                                 "Process :: CalculateWeight (AddedSubtraction)"
                              );
            if (tcpass) {
               final -> DoSumOverAddedSubtractionTerms(
                           this,  
                           ic,
                           xjacobian * fs_jacobian,
                           ca
               );
            }
            break;
  
         case CollinearInitial:

            // this term is a collinear contribution in the initial state

            if (   final -> GetFrame()==hCMS
                && final -> GetMapping(-1) == NULL) {
               final -> MapUnitTo_u(ustore);
               final -> EcmAddRemnantAndIncident();
            }

            // check the technical cut
            tcpass = final -> CheckTechnicalCut(
                                 this,
                                 "Process :: CalculateWeight (CollinearInitial)"
                              );
            if (tcpass) {
               final -> DoSumOverCollinearInitialTerms(
                           this, 
                           ic,
                           xjacobian * fs_jacobian,
                           ca
                        );
            }
            break;
  
         default:

            errf(-1, "Process :: CalculateWeight: unknown subtraction type");

      } // of switch ()

   } // of if
   } // of for (ic)
}

// ----------------------------------------------------------------------------

// copy the offsets from another component

void Process :: Copy_vtp(int to, int from) {

   (vtp + from) -> CopyInto(vtp + to);
}

// ----------------------------------------------------------------------------

void Process :: Set_componentFlag(int i, int val) {

   if (i >= 0 && i < n_components)
      componentFlags[i] = val;
   else {
      global -> To_error() << FS("i=%d", i)
                           << FS(" n_components=%d\n", n_components);
      errf(-1, "Process :: Set_componentFlag: no such component!");
   }
}

// ----------------------------------------------------------------------------

void Process :: Print_componentFlags(Mstream& s_out) {

   s_out << "Component flags:\n";
   for (register int i = 0; i < n_components; ++i)
       s_out << FS("#%3d", i)
             << FS(" --> %3d\n", componentFlags[i]);
}

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
// ============================================================================
//
// --> string class
//
// file:              strng.cc
// created:           29.03.1997
// last modification: 11.12.1997
//
// ============================================================================

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

//#include "container.h"

#include "global.h"  // for errf()
#include "mth.h"

#include "strng.h"

// ============================================================================
//
// -> functions
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> some utility routines to extract parameters from a string
//
// ----------------------------------------------------------------------------

long ReadLong(Global* global, char* pname, char* ins) {
 
   long i;
   int assigned = sscanf(ins, "%ld", &i);
   
   if (assigned < 1) {
      global -> To_error() << FS("string -> :%s:\n",ins);
      errf(-1,"ReadLong.");
   }

   global -> To_status() << FS(":%-42s: --> ", pname)
                         << FS("%30ld\n", i);

   global -> To_log()    << FS("%-42s ", pname)
                         << FS("%30ld\n", i);

   return i;
}

// ----------------------------------------------------------------------------

int ReadInt(Global* global, char* pname, char* ins) {
 
   return (int) ReadLong(global, pname, ins);
}

// ----------------------------------------------------------------------------

int ReadBool(Global* global, char* pname, char* ins)
{
   long in = ReadLong(global, pname, ins);
  
   int ret;
   switch (in) {
      case 0:
         ret = FALSE;
         break;
      case 1:
         ret = TRUE;
         break;
      default:
         ret = FALSE; // :_TAWM_:
         errf(-1,"ReadBool: only 0 and 1 are allowed");
   }

   return ret;
}

// ----------------------------------------------------------------------------

double ReadDouble(Global* global, char* pname,char* ins)
{
   double x;
   int assigned = sscanf(ins, "%lf", &x);
   
   if (assigned < 1) {
      global -> To_error() << FS("string -> :%s:\n",ins);
      errf(-1,"ReadDouble.");
   }

   global -> To_status() << FS(":%-42s: --> ", pname)
                         << FS("%30.22e\n", x);

   global -> To_log()    << FS("%-42s ", pname)
                         << FS("%30.22e\n", x);

   return x;
}

// ----------------------------------------------------------------------------

void ReadString(Global* global, char* pname, char* ins, String* str) {

   // :_TAWM_:
   global = global;
   pname  = pname;

   int i = 0;
   while (ins[i] != '\0' && isspace(ins[i])) 
      ++i;

   int quoted;
   quoted = (ins[i] == '\0') ? FALSE : ins[i]=='\"';
   if (quoted)
      ++i;

   int j = 0;
   while (   ins[i] != '\0'
          && (quoted && ins[i] != '\"' || !quoted && !isspace(ins[i])) 
          && j < str -> maxlength)
      str -> data[j++] = ins[i++];
   str -> data[j++] = '\0';

   if (    quoted && ins[i] != '\"' 
       || !quoted && !(isspace(ins[i]) || ins[i] == '\0') 
      ) {
      global  
      -> To_error() << FS("string (max=%d) -> ", str -> maxlength)
                    << FS(":%s:\n", ins);
      errf(-1,"ReadString: string too long || unterminated");
   }

   global -> To_log() << FS("# ReadString: :%s:\n", str -> data);
}

// ----------------------------------------------------------------------------
//
// --> skips white space and comments on a file
//
// ----------------------------------------------------------------------------

// skips white space and comments (#: skip until end of line)
int SkipWhiteSpace(FILE* in) {

   int ch, ret;

   do {
      while ((ch = fgetc(in)) != EOF && isspace(ch));
      if (ch == '#') {
         // skip until EOL
         while ((ch = fgetc(in)) != EOF && ch != '\n');
      }   
   } while (ch != EOF && ch == '\n');

   if (ch == EOF) {
      ret = TRUE;
   } else {
      ungetc(ch, in);
      ret = FALSE;
   }

   return ret;
}

// ============================================================================
//
// --> classes
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> string class
//
// ----------------------------------------------------------------------------

String :: String(int maxlengthI) {

   data = NULL;
   Define(maxlengthI);
}

// ----------------------------------------------------------------------------

String :: String(const char *in) {

   data = NULL;
   Define(strlen(in));
   CopyFromChar(in);
}

// ----------------------------------------------------------------------------

String :: String(const String *str1, const String *str2) {

   data = NULL;
   Define(str1, str2);
}

// ----------------------------------------------------------------------------

String :: String(const String &str1, const String &str2) {

   data = NULL;
   Define(str1, str2);
}

// ----------------------------------------------------------------------------

String :: ~String() {

   Delete();
}

// ----------------------------------------------------------------------------

void String :: Define(int maxlengthI) {

   if (data != NULL)
      errf(-1,"String :: Define.");
   maxlength = maxlengthI;

   if (maxlength >= 0) {
      data = new char[maxlength+1];
      data[0] = '\0';
   }
   else {
      data = NULL;
   }
}

// ----------------------------------------------------------------------------

void String :: Define(const String *str1, const String *str2) {

   Define(str1 -> DetermineLength() + str2 -> DetermineLength());
   CopyFromString(str1);
   AppendFromString(str2);
}

// ----------------------------------------------------------------------------

void String :: Define(const String &str1, const String &str2) {

   Define(&str1, &str2);
}

// ----------------------------------------------------------------------------

void String :: Delete() {

   if (data != NULL) {
      delete[] data;
      data = NULL;
   }
}

// ----------------------------------------------------------------------------

void String :: Print(FILE* file) {

   fprintf(file, "%s", data);
   errf(-1, "USED String :: Print!");
}

// ----------------------------------------------------------------------------

int String :: DetermineLength() const {
   return data == NULL ? 0 : strlen(data);
}

// ----------------------------------------------------------------------------

void String::CopyFromChar(const char* in) {

   register int i;
   for (i = 0; i < maxlength && in[i] != '\0'; ) {
       data[i] = in[i];
       ++i;
   }
   data[i] = '\0';
}

// ----------------------------------------------------------------------------

void String :: AppendFromChar(const char* in) {

   register int i, j;
   for (i = strlen(data), j=0; i < maxlength && in[j] != '\0'; ) {
       data[i++] = in[j++];
   }
   data[i] = '\0';

   if (in[j] != '\0') {
      gl -> To_error() << FS(":%s: ", in)
                       << FS("%d", maxlength);
      errf(-1,"String :: AppendFromChar: string too long");
   }
}

// ----------------------------------------------------------------------------

void String::DefineAndCopyFromChar(const char* in) {

   Define(strlen(in));
   CopyFromChar(in);
}

// ----------------------------------------------------------------------------

void String :: CopyFromString(const String* in) {

   register int i;
   for (i = 0; i < maxlength && in -> data[i] != '\0'; ) {
       data[i] = in -> data[i];
       ++i;
   }
   data[i] = '\0';

   if (in -> data[i] != '\0') {
      gl -> To_error() << FS(":%s: ", in -> data)
                       << FS("%d", maxlength);
      errf(-1,"String :: CopyFromString: string too long");
   }
}

// ----------------------------------------------------------------------------

void String :: AppendFromString(const String* in) {

   register int i, j;

   for (i = strlen(data), j = 0; i < maxlength && in -> data[j] != '\0'; ) {
       data[i++] = in -> data[j++];
   }
   data[i] = '\0';

   if (in -> data[j] != '\0') {
      gl -> To_error() << FS(":%s: ", in -> data)
                       << FS("%d", maxlength);
      errf(-1,"String :: AppendFromString: string too long");
   }
}

// ----------------------------------------------------------------------------

void String :: DefineAndCopyFromString(const String* in) {

   Define(in -> DetermineLength());
   CopyFromString(in);
}

// ----------------------------------------------------------------------------

void String :: DefineAndCopyFromString(const String &in) {

   Define(in.DetermineLength());
   CopyFromString(&in);
}

// ----------------------------------------------------------------------------

// modify all '\' into '\n'

void String :: BackslashToNewline() {

   char* p = data - 1;
   while (*(++p) != '\0')
      if (*p == '\\')
         *p = '\n';
}

// ----------------------------------------------------------------------------

// reads a string enclosed in '' from a file

void String :: ReadBoundedFromFile(FILE* in) {

   int ch;
   int i = 0;

   while ((ch = fgetc(in)) != EOF && isspace(ch));
   if (ch == EOF)
      errf(-1,"String :: ReadBoundedFromFile: EOF reached.");

   if (ch != '\'')
      errf(-1,"String :: ReadBoundedFromFile: \' expected.");

   while ((ch = fgetc(in)) != EOF && ch!='\'' && i<maxlength) {
      data[i++] = ch;
   } 

   if (ch == EOF)
      errf(-1, "String :: ReadBoundedFromFile: EOF reached.");

   if (ch != '\'')
      errf(-1, "String :: ReadBoundedFromFile: string too long.");

   data[i] = '\0';
}

// ----------------------------------------------------------------------------

// reads a string of non-white characters; terminates if the
// first white character or EOF is encountered;
// also handles quotes

void String :: ReadNonWhiteFromFile(FILE* in) {

   int i = 0;
   int ch;
   int quoted;

   if ((ch = fgetc(in)) != EOF && ch == '\"') {
      quoted = TRUE;
      data[i++] = ch;
   }
   else {
      quoted = FALSE;
      ungetc(ch, in);
   }

   while (   (ch = fgetc(in)) != EOF 
          && (quoted && ch != '\"' || !quoted && !isspace(ch)) 
          && i < maxlength) {
      data[i++] = ch;
   } 

   if (quoted && ch == '\"' && i < maxlength)
      data[i++] = ch;

   data[i] = '\0';

   if (    quoted && ch != '\"'
       || !quoted && !(isspace(ch) || ch == EOF)
      )
      errf(-1, "String :: ReadNonWhiteFromFile: string too long.");
}

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
// ============================================================================
//
// --> user routines
//
// file:              user.cc
// created:           29.03.1997
// last modification: 15.12.1997
//
// ============================================================================

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "global.h"
#include "qcd.h"
#include "disaster.h"
#include "strng.h"
#include "mth.h"
#include "container.h"

#include "user.h"

// ============================================================================
//
// definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

// ============================================================================
//
// --> class-unrelated functions
//
// ============================================================================

// ============================================================================
//
// --> classes
//
// ============================================================================
   
// ----------------------------------------------------------------------------
//
// --> flavour factors and alphaS
//
// ----------------------------------------------------------------------------
   
PDFC :: PDFC() {
}

// ----------------------------------------------------------------------------
   
PDFC :: ~PDFC() {
}

// ----------------------------------------------------------------------------
   
void PDFC :: Define(Disaster* disasterIn) {

   if (IsDefined())
      errf(-1, "PDFC :: Define: already defined");

   disaster = disasterIn;
}

// ----------------------------------------------------------------------------
   
void PDFC :: Reset(Disaster* disasterIn) {

   disaster = disasterIn;
}

// ----------------------------------------------------------------------------
   
void PDFC :: ReturnToFreeList() {
}

// ----------------------------------------------------------------------------
   
void PDFC :: Print(Mstream& s_out) {
   
   s_out = s_out; // :_TAWM_:

   s_out << "PDFC :: Print\n";

   register int i;
   for (i = 0; i < 11; ++i)
       s_out << FS("%16.6e\n", ff[i]);

//   errf(-1, "PDFC :: Print: not implemented");
}

// ----------------------------------------------------------------------------
   
void PDFC :: Print(FILE* out) {

   Mstream s_out;
   s_out.Associate(out);   

   Print(s_out);
}

// ----------------------------------------------------------------------------
   
int PDFC :: UseListIsSame(PDFC* in) {
   
   in = in; // :_TAWM_:

   errf(-1, "PDFC :: UseListIsSame: not implemented");

   return FALSE;
}

// ----------------------------------------------------------------------------
   
long PDFC :: GetIdent() {
   
   errf(-1, "PDFC :: GetIdent: not implemented");

   return -999999;
}

// ----------------------------------------------------------------------------
//
// --> flavour factors
//
// ----------------------------------------------------------------------------
   
FlavourFactors :: FlavourFactors() {
}

// ----------------------------------------------------------------------------
   
FlavourFactors :: ~FlavourFactors() {

   Delete();
}

// ----------------------------------------------------------------------------
   
void FlavourFactors :: Define(Disaster* disasterIn) {

   if (IsDefined())
      errf(-1, "FlavourFactors :: Define: already defined");

   disaster = disasterIn;

   pdfcfl = NULL;
   cache  = NULL;   

   pdf_Server = NULL;
   alphaS_Server = NULL;
   alphaEM_Server = NULL;
}

// ----------------------------------------------------------------------------
   
void FlavourFactors :: Delete() {

   if (cache != NULL)
      DeleteCache();
}

// ----------------------------------------------------------------------------
   
void FlavourFactors :: CalculateFlavourFactors(
                          PDFC* pdfc,
                          double xi,
                          double muF,
                          double muR,
                          int nflavour_pden,
                          int nflavour_out,
                          int pdf_collectionIn,                               
                          int pdf_parametrizationIn,
                          int pdf_setIn,
                          int alphaSVariantIn,
                          int alphaSOrderIn,
                          double alphasLambdaQCD4In,
                          int alphaEMVariantIn,
                          double alphaEMScaleIn
                       ) {

   double* ff = pdfc -> GetFlavourFactorPtr();
   double* alphaSPtr  = pdfc -> GetAlphaSPtr();
   double* alphaEMPtr = pdfc -> GetAlphaEMPtr();

   double* qcharge  = disaster
                      -> GetQuarkCharge()
                      -> GetDataPtr();
   
   double* qcharge2  = disaster
                       -> GetQuarkChargeSquared()
                       -> GetDataPtr();

   double* qcharge2Single = disaster
                            -> GetQuarkChargeSquaredPartialSumSingle()
                            -> GetDataPtr();

   double** qchargeDouble = disaster
                            -> GetQuarkChargePartialSumDouble()
                            -> GetDataPtr();

   double** qcharge2Double = disaster
                             -> GetQuarkChargeSquaredPartialSumDouble()
                             -> GetDataPtr();

   PartonDensity* pdensityX;
   double* pdensity;

   pdensityX 
      = pdf_Server
        -> FindOrCreateAndFillAndCachePartonDensity(
              xi, muF, -1,
              pdf_collectionIn, 
              pdf_parametrizationIn, 
              pdf_setIn
           );

    pdensity = pdensityX -> GetDataPtr();

    register int i;
    for (i = 0; i <= 10; ++i)
        ff[i] = 0.;

    double qfac1, qfac2;

    for (i = nflavour_pden; i >= 1; --i) {
        ff[ 0] += pdensity[ i] * qcharge2[i];
        ff[ 1] += pdensity[-i] * qcharge2[i];
        ff[ 7] += pdensity[ i] 
                     * (qfac2 = qcharge2Double[nflavour_out][i]);
        ff[ 8] += pdensity[-i]  
                     * qfac2;
        ff[ 9] += pdensity[ i] 
                     * qcharge[i]    
                     * (qfac1 = qchargeDouble[nflavour_out][i]);
        ff[10] += pdensity[-i] 
                     * qcharge[i] 
                     * qfac1;
    }
    ff[5] = ff[0] * (nflavour_out - 1);
    ff[6] = ff[1] * (nflavour_out - 1);

    ff[2] = pdensity[0] * qcharge2Single[nflavour_out];

    for (i = min(nflavour_pden, nflavour_out); i >= 1; --i) {
        ff[3] += pdensity[ i] * qcharge2[i];
        ff[4] += pdensity[-i] * qcharge2[i];
    }

    *alphaSPtr = alphaS_Server
                 -> EvaluateAlphaSAutomaticNquarks(
                       muR,
                       alphaSVariantIn,
                       alphaSOrderIn,
                       alphasLambdaQCD4In
                    );

#if 0
    // cross check for PDFLIB alpha_s
    AlphaS asCheck(disaster);
    asCheck.Define();
    asCheck.Set_outputFlag(0);
    double myAlphaS 
              = asCheck 
                .EvaluateAlphaSAutomaticNquarks(
                      muR,
                      1, 
                      2, 
                      alphaS_Server -> GetLambdaQCD_4()
                   );
   disaster -> To_status()
      << FS("%2d", alphaSVariantIn)
      << FS("%2d", alphaSOrderIn)
      << FS("%16.6e", muR)
      << FS("%16.6e", alphaS_Server -> GetLambdaQCD_4())
      << FS("%16.6e", *alphaSPtr)
      << FS("%16.6e\n", myAlphaS);
#endif

    *alphaEMPtr = alphaEM_Server
                  -> EvaluateAlphaEM(alphaEMScaleIn, alphaEMVariantIn);
}

// ----------------------------------------------------------------------------
   
PDFC* FlavourFactors :: FindOrCreateAndFillAndCacheFlavourFactors(
                           double xi,
                           double muF,
                           double muR,
                           int nflavour_pden,
                           int nflavour_out,
                           int pdf_collectionIn,                               
                           int pdf_parametrizationIn, 
                           int pdf_setIn,
                           int alphaSVariantIn,
                           int alphaSOrderIn,
                           double alphasLambdaQCD4In,
                           int alphaEMVariantIn,
                           double alphaEMScaleIn
                        ) {
   
   PDFC* pdfc;

   if (cache == NULL)   
      errf(-1, "FlavourFactors"
               " :: FindOrCreateAndFillAndCacheFlavourFactors: no cache!");
                        
   intKey[0]    = nflavour_pden;
   intKey[1]    = nflavour_out;
   intKey[2]    = pdf_collectionIn;
   intKey[3]    = pdf_parametrizationIn;
   intKey[4]    = pdf_setIn;
   intKey[5]    = alphaSVariantIn;
   intKey[6]    = alphaSOrderIn;
   intKey[7]    = alphaEMVariantIn;
   doubleKey[0] = xi;
   doubleKey[1] = muF;   
   doubleKey[2] = muR;   
   doubleKey[3] = alphasLambdaQCD4In;   
   doubleKey[4] = alphaEMScaleIn;   
               
   // try to find in cache
   pdfc = cache -> FindInCache(intKeyArray, doubleKeyArray);
                        
   if (pdfc == NULL) {
      pdfc = pdfcfl -> GetDefinedObject(disaster);
      CalculateFlavourFactors(
         pdfc,
         xi,
         muF,
         muR,
         nflavour_pden,
         nflavour_out,
         pdf_collectionIn,
         pdf_parametrizationIn,
         pdf_setIn,
         alphaSVariantIn,
         alphaSOrderIn,
         alphasLambdaQCD4In,
         alphaEMVariantIn,
         alphaEMScaleIn
      );
      cache -> StoreInCache(intKeyArray, doubleKeyArray, pdfc);
   }
               
   return pdfc;
}

// ----------------------------------------------------------------------------
      
void FlavourFactors :: CreateCache(int size) {
   
   if (cache != NULL)
      errf(-1, "FlavourFactors :: CreateCache: "
               "cache already exists");
      
   int nInt    = 8;
   int nDouble = 5; 
   cache = new Cache <PDFC> (disaster, size, nInt, nDouble);
   cache -> DefineName("FlavourFactors cache");
   cache -> SetFreeList(pdfcfl);

   intKeyArray    = new Array_int(nInt);
   intKey         = intKeyArray -> GetDataPtr();

   doubleKeyArray = new Array_double(nDouble);
   doubleKey      = doubleKeyArray -> GetDataPtr();
}

// ----------------------------------------------------------------------------

void FlavourFactors :: DeleteCache() {

   delete cache;
   cache = NULL;

   delete intKeyArray;
   intKeyArray = NULL;

   delete doubleKeyArray;
   doubleKeyArray = NULL;
}
 
// ----------------------------------------------------------------------------
//
// --> MC User routines
//
// ----------------------------------------------------------------------------
   
MCUser :: MCUser() {

   global  = NULL;
   gl -> To_status() << "MCUser :: MCUser()\n";
}

// ----------------------------------------------------------------------------
   
MCUser :: ~MCUser() {

   gl -> To_status() << "MCUser :: ~MCUser()\n";
   Delete();
}

// ----------------------------------------------------------------------------
   
void MCUser :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined()) 
      errf(-1, "MCUser :: Define: already defined"); 

   gl -> To_status() << "MCUser :: Define()\n";

   global  = globalIn;
}

// ----------------------------------------------------------------------------
   
void MCUser :: Delete() {

   if (definedObject.Status()) {
      definedObject.ResetDefined();
      gl -> To_status() << "MCUser :: Delete()\n";
   }
}

// ----------------------------------------------------------------------------
//
// --> MC User routine called from Disaster
//
// ----------------------------------------------------------------------------
   
MCUserForDisaster :: MCUserForDisaster() {

   gl -> To_status() << "MCUserForDisaster :: MCUserForDisaster()\n";
}

// ----------------------------------------------------------------------------
   
MCUserForDisaster :: ~MCUserForDisaster() {

   gl -> To_status() << "MCUserForDisaster :: ~MCUserForDisaster()\n";
   Delete();
}

// ----------------------------------------------------------------------------
   
void MCUserForDisaster :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined()) 
      errf(-1, "MCUserForDisaster :: Define: already defined"); 

   gl -> To_status() << "MCUserForDisaster :: Define()\n";

   MCUser :: Define(globalIn);

   flavourFactors = new FlavourFactors;
   flavourFactors -> Define(globalIn);

   pdfcFree = new PDFCFreeList(global);
   pdfcFree -> DefineName("PDFC free list");

   pdFree 
      = global
        -> disaster
        -> Get_qcd()
        -> Get_pdFree();

   alphaS_Server 
      = global
        -> disaster
        -> Get_qcd()
        -> Get_alphaS_Server();

   pdf_Server
      = global
        -> disaster
        -> Get_qcd()
        -> Get_pdf_Server();

   flavourFactors -> SetPDFCFL(pdfcFree);
// DEBUG   flavourFactors -> CreateCache(500);
   flavourFactors -> CreateCache(20);
   flavourFactors -> Set_pdf_Server(pdf_Server);
   flavourFactors -> Set_alphaS_Server(alphaS_Server);

   alphaEM_Server 
      = global
        -> disaster
        -> Get_electroWeak()
        -> Get_alphaEM_Server();

   flavourFactors -> Set_alphaEM_Server(alphaEM_Server);

   uoEnergyPower = 0.5;
   uoAnglePower  = 0.5;

   nOutFlavours  = 5;
   nPdenFlavours = 6;

   variant_alpha_s = 1;          // running coupling
   user_defined_alpha_s = FALSE; // fixed by selected process
   order_alpha_s = 1;            // if not: first order
   lambda_4_alpha_s = 0.2;       // with this lambda_QCD value for 4 flavours
   pdf_collection = 1;           // PDFLIB
   pdf_parametrization = 4;      // | CTEQ leading order
   pdf_set = 29;                 // ^

   universalObservableChoice = 0;

   testCase = FALSE;

   observablesD = 0;

   nFlag = 0;

   mcUserInterface -> Define(global);
}

// ----------------------------------------------------------------------------
   
void MCUserForDisaster :: Delete() {

   if (definedObject.Status()) {

      definedObject.ResetDefined();

      gl -> To_status() << "MCUserForDisaster :: Delete()\n";

      delete flavourFactors;
      delete pdfcFree;

      MCUser :: Delete();
   }
}

// ----------------------------------------------------------------------------
   
void MCUserForDisaster :: Start() {

   mcUserInterface -> Start();
}

// ----------------------------------------------------------------------------
   
void MCUserForDisaster :: End() {

   mcUserInterface -> End();
}

// ----------------------------------------------------------------------------

// called if a parameter name parsed by the MC cannot be resolved
// --> assumed to be a user parameter

int MCUserForDisaster :: SetParameter(char* pname, char* parameter) {

   int ret = 0;

   if (!strcmp(pname,"TEST_CASE")) {
      testCase=ReadBool(global,pname,parameter);
   }
   else

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   // observables

   if (!strcmp(pname, "UNIVERSAL_OBSERVABLE_CHOICE")) {
      universalObservableChoice = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "UNIVERSAL_OBSERVABLE_ENERGY_POWER")) {
      uoEnergyPower = ReadDouble(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "UNIVERSAL_OBSERVABLE_ANGLE_POWER")) {
      uoAnglePower = ReadDouble(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "OBSERVABLES_D")) {
      observablesD = ReadInt(global, pname, parameter);
   }
   else

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   // test: flags for components

   if (!strcmp(pname, "COMPONENT_FLAG_RESET")) {
      nFlag = 0;
   }
   else

   if (!strcmp(pname, "COMPONENT_FLAG_INDEX")) {
      cFlagIndex[nFlag] = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "COMPONENT_FLAG_VALUE")) {
      cFlagValue[nFlag] = ReadInt(global, pname, parameter);      
      ++nFlag;
   }
   else

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   {
      ret = mcUserInterface -> SetParameter(pname, parameter);
   }

   return ret;
}

// ----------------------------------------------------------------------------

void MCUserForDisaster :: StartIntegration1() {

   n_events = 50;

   // modify component flags
   for (int index = 0; index < nFlag; ++index)
       global 
       -> process 
       -> Set_componentFlag(cFlagIndex[index], cFlagValue[index]);
   global
   -> process 
   -> Print_componentFlags(global -> To_status());

   globalLibrary = new Library(global, 20);

   ( adaptationObservable = globalLibrary -> CreateNewBook() )
        -> Define(Graph, Linear, 1, 0.0, 0.0,
                  "adaptationObservable", "adaptationObservable", 
                  n_events);

   if (observablesD) {

      unitCubeHistogram = new Book* [global -> disaster -> Get_dimension()];
      for (int idim = 0; 
           idim < global -> disaster -> Get_dimension(); 
           ++idim) {
         char hname[100];
         char digits[100];
         sprintf(digits, "%02d", idim);
         strcpy(hname, "unitCube_");
         strcat(hname, digits);
         ( unitCubeHistogram[idim]
              = globalLibrary -> CreateNewBook() )
                -> Define(Histogram, Linear, 10, 0.0, 1.0,
                          hname, hname, 
                          n_events);
      }   

   }  // of ``observablesD''
   
   globalLibrary -> ResetStorage();

   mcUserInterface -> StartIntegration1();
}

// ----------------------------------------------------------------------------

void MCUserForDisaster :: StartIntegration2() {

   global -> To_status() << "\n";
   global -> To_status()
      << FS("event         count = %8d\n",
            global -> disaster -> Get_event_counter());
   global -> To_status()
      << FS("fpError       count = %8d\n", 
            global -> disaster -> Get_fpError_counter());
   global -> To_status()
      << FS("technical Cut count = %8d\n", 
            global -> disaster -> Get_technicalCut_counter());
   global -> To_status() << "\n";

   // has already been filled up; prepare for final run
   globalLibrary -> ResetStorage();  
       
   if (testCase) {
      nused = 30;
      int i;
      for (i = 0; i < nused; ++i) {
          sum[i]   = 0.;
          sum_2[i] = 0.;
      }
   }

   mcUserInterface -> StartIntegration2();
}

// ----------------------------------------------------------------------------

void MCUserForDisaster :: EndIntegration(double average, double error) {

   global -> To_status() << "\n";
   global -> To_status()
      << FS("event         count = %8d\n",
            global -> disaster -> Get_event_counter());
   global -> To_status()
      << FS("fpError       count = %8d\n", 
            global -> disaster -> Get_fpError_counter());
   global -> To_status()
      << FS("technical Cut count = %8d\n", 
            global -> disaster -> Get_technicalCut_counter());
   global -> To_status() << "\n";

//   fpnan -> PrintStatus(global -> To_error());

   mcUserInterface -> EndIntegration(average, error);

//   fpnan -> PrintStatus(global -> To_error());      

   if (testCase) {
      int i;
      for (i = 0; i < nused; ++i) {
          double term 
                 = 1. / (global -> VEGAS_pts_in_last_iteration-1)
                    * (  global -> VEGAS_pts_in_last_iteration * sum_2[i]
                       - sum[i] * sum[i]
                      );
          if (term >= 0.)
             err[i] = sqrt(term);
          else {
             global -> To_status()  
                          << "MCUserForDisaster :: EndIntegration: negative\n";
             err[i] = 0.;
          }
          if ( !(fabs(sum[i]) < 1.e-50 || fabs(err[i]) > 1.e50) )
             relative_err[i] = fabs(err[i] / sum[i]);
          else
             relative_err[i] = -9.999e97;

          global -> To_status()
                       << FS("sum[%2d] = ", i) 
                       << FS("%16.6e", sum[i])
                       << FS("+- %16.6e", err[i])
                       << FS(" ( =^= %16.6e %% )\n", 100*relative_err[i]);
      }
   }
 
   // -------------------------------------------------------------------------

   globalLibrary -> CalculateErrors();
   globalLibrary -> Print(global -> To_status());   

   double averageX;
   double errorX;
   double relErrorX;
                            
   adaptationObservable -> GetEntry(-1, &averageX, &errorX, &relErrorX);

   global -> To_status() 
                << "uobs  =  " 
                << FS("%10.3e", averageX)
                << " +-"
                << FS("%10.3e", errorX)
                << "\n";
   
   delete globalLibrary;

   if (observablesD)
      delete [] unitCubeHistogram;

//   fpnan -> PrintStatus(global -> To_error());      
}

// ----------------------------------------------------------------------------

void MCUserForDisaster :: BeginEvent(Process* process) {

   mcUserInterface -> BeginEvent(process);

   globalLibrary -> PrepareForEvent();
}

// ----------------------------------------------------------------------------

void MCUserForDisaster :: EndOfEvent(Process* process) {

   mcUserInterface -> EndOfEvent(process);
}

// ----------------------------------------------------------------------------

double MCUserForDisaster :: EndEvent(
                    Process* process,
                    ContributionArray* ca,
                    int ifFinalRun
                 ) {

   // boost all event records to the Breit frame

   UseListContainer <Event> * ec;
   Event* e;

   // run through the list of all events
   for (ec = process -> GetEventList() -> GetPtrToFirstContainer();
        ec != NULL && (e = ec -> GetDataPtr(), TRUE);
        ec = ec -> GetNext()) {
       if (e -> GetContribution() != NULL)
          BoostToBreitFrame(process, e, ca, ifFinalRun);
   }

   double ret = mcUserInterface -> EndEvent(process, ca, ifFinalRun);


   // -------------------------------------------------------------------------

   if (testCase) {

      if (ifFinalRun) {
    
         register int i;

         for (i = 0; i < nused; ++i)
             contr[i] = 0.;

         contr[3] = event_sum;
         contr[8] = event_sum_q;
         contr[9] = event_sum_g;
         for (i = 0; i <= 10; ++i)
             contr[i+10] = event_ff[i];

         // sum up terms
 
         for (i = 0; i < nused; ++i) {
             sum  [i] += contr[i];
             sum_2[i] += contr[i] * contr[i];
         }          
      }
   }

   // -------------------------------------------------------------------------

   // the weight returned to the integration routine

   if (testCase) {

      // test case: no parton density, no coupling constants
     
      event_sum = 0;
      for (register int l = ca -> GetActual(); l >= 0; --l) {
          Contribution* con = ca -> GetContribution(l);
          double *data = con -> GetData();
          event_sum += data[0];   
      }
      ret = event_sum;
   }

   adaptationObservable -> PrestoreGraph(1.0, 0, 0);
   adaptationObservable -> SPS_Graph(0, 0, ret);

   if (   ifFinalRun
       && observablesD) {

      for (int idim = global -> disaster -> Get_dimension() - 1;
           idim >= 0;
           --idim) {
           unitCubeHistogram[idim]
           -> StoreHistogram(
                 (global -> disaster -> Get_unitCube()) [idim], 
                 ret
              ); 
      }
   }

   globalLibrary -> AddUpEvent();
   
   return ret;
}

// ----------------------------------------------------------------------------

// :_MOD_: the arrays should be safe arrays!

void MCUserForDisaster :: ConvoluteWithDistributions(
                  Process* process,
                  ContributionArray* ca, 
                  double* ren_scale,
                  double* fact_scale, 
                  int* nflavour_out,
                  int* nflavour_pden, 
                  int* nflavour_rs,
                  int* nflavour_fs,
                  int pdf_collectionIn,
                  int pdf_parametrizationIn,
                  int pdf_setIn,
                  int alphaSVariantIn,
                  int alphaSOrderIn, 
                  double alphasLambdaQCD4In,
                  int alphaEMVariantIn,
                  double* alphaEMScaleIn 
               ) {

   // :_TAWM_:
   process = process;

   register int i, l;

   // calculate weight: loop over all contributions

   Contribution* con;
   event_sum_q = 0;
   event_sum_g = 0;

   for (l = ca -> GetActual(); l >= 0; --l) {

       con = ca -> GetContribution(l);

       double xi = con -> Get_xi_U();

       PDFC* pdfc;

#if 0
       global -> To_status()
          << FS("%16.6e ", xi)
          << FS("%16.6e ", fact_scale[l])
          << FS("%16.6e ", ren_scale[l])
          << FS("%16.6e\n", alphaEMScaleIn[l]);
       global -> To_status() << FS("%d\n", l);
       global -> To_status()
          << FS("%d ", nflavour_pden[l])
          << FS("%d ", nflavour_out[l])
          << FS("%d ", pdf_collectionIn)
          << FS("%d ", pdf_parametrizationIn)
          << FS("%d\n", pdf_setIn);
#endif

       pdfc = flavourFactors 
              -> FindOrCreateAndFillAndCacheFlavourFactors(
                    xi,
                    fact_scale[l],
                    ren_scale[l],
                    nflavour_pden[l],
                    nflavour_out[l],
                    pdf_collectionIn,
                    pdf_parametrizationIn,
                    pdf_setIn,
                    alphaSVariantIn,
                    alphaSOrderIn,
                    alphasLambdaQCD4In,
                    alphaEMVariantIn,
                    alphaEMScaleIn[l]
                 );

       double* ff = pdfc -> GetFlavourFactorPtr();
       double alphaS  = *(pdfc -> GetAlphaSPtr());
       double alphaEM = *(pdfc -> GetAlphaEMPtr());

       double afact = 1.;

       // alphaS ^ (con -> GetOrderAlphaS())
       switch(con -> GetOrderAlphaS()) {
          case 2:
             afact *= alphaS;
             // fall through
          case 1:
             afact *= alphaS;
             // fall through
          case 0:
             break;
          default:
             errf(-1, "MCUserForDisaster::ConvoluteWithDistributions: "
                      "order in alphaS too large...");
       }

       // alphaEM ^ (con -> Get_orderAlphaEM())
       switch(con -> Get_orderAlphaEM()) {
          case 2:
             afact *= alphaEM;
             // fall through
          case 1:
             afact *= alphaEM;
             // fall through
          case 0:
             break;
          default:
             errf(-1, "MCUserForDisaster::ConvoluteWithDistributions: "
                      "order in alphaEM too large...");
       }

       double* data = con -> GetData();

       // append factors for the scale dependence
    
       switch (con -> GetSLF()) {
          case RenormalizationScaleLogarithmWithNf:
             afact *= nflavour_rs[l];
             // --- fall-through! ---
          case RenormalizationScaleLogarithm:
             afact *= 2 * log( ren_scale[l] / con -> Get_Q_U() );
             break;
          case FactorizationScaleLogarithmWithNf:
             afact *= nflavour_fs[l];
             // --- fall-through! ---
          case FactorizationScaleLogarithm:
             afact *= 2 * log( fact_scale[l] / con -> Get_Q_U() );
             break;
          default:
             // do nothing
             ;
       }

       double event_ff_c[11];

       for (i = 0; i <= 10; ++i) {
           event_ff[i] += (event_ff_c[i] = ff[i] * data[i] * afact);
       }

       tsumq[l] =   event_ff_c[ 0]
                  + event_ff_c[ 1]
                  + event_ff_c[ 3]
                  + event_ff_c[ 4]
                  + event_ff_c[ 5]
                  + event_ff_c[ 6]
                  + event_ff_c[ 7]
                  + event_ff_c[ 8]
                  + event_ff_c[ 9]
                  + event_ff_c[10]
       ;

       tsumg[l] =   event_ff_c[2];

       tsum[l] = tsumq[l] + tsumg[l];

       event_sum_q += tsumq[l];
       event_sum_g += tsumg[l];
   }

   event_sum = event_sum_q + event_sum_g;
}

// ----------------------------------------------------------------------------

void MCUserForDisaster :: DropEvent(Process* process) {

   mcUserInterface -> DropEvent(process);
}

// ----------------------------------------------------------------------------

void MCUserForDisaster :: AcceptEvent(Process* process) {

   mcUserInterface -> AcceptEvent(process);
}

// ----------------------------------------------------------------------------

// Determine the ``universal observable''
// (to obtain the weight to be returned to the adaptive routine)

double MCUserForDisaster :: UniversalObservable(Event* e) {

   double obs;

   switch (universalObservableChoice) {

      case 0:

         if ( global -> numberOfFinalStatePartonsInBornTerm == 2) {

            Event* univ;

            if (e -> npartons == 3) {
               univ = global
                      -> disaster
                      -> GetEFree()
                      -> GetDefinedObject(
                            global -> disaster, 
                            e -> GetMaxNumber()
                         );
               e -> CopyInto(univ);
               // :_MOD_: this includes the remnant in the 
               // clustering procedure...
               // should be better if we had the incident parton instead...
               univ -> ClusterGeneralOneStep(JADE_algorithm, E_Scheme);
            } else 
               univ = e;

            switch (e -> npartons) {

               case 2:
               case 3:  // :_CAUTION_: for 3 partons??

                  obs = pow(
                              univ -> GetMapping(0) -> GetEnergy()
                            * univ -> GetMapping(1) -> GetEnergy(),
                            uoEnergyPower
                        )
                      * pow(                  
                              univ -> VMap(0,  1)
                            * univ -> VMap(0, -1)
                            * univ -> VMap(1, -1),
                            uoAnglePower
                        );

                  break;

               default: 
                  e -> Print();
                  errf(-1, "MCUserForDisaster :: UniversalObservable: error");
                  obs = 0; // :_TAWM_:
            }

            if (e -> npartons == 3)
               global -> disaster -> GetEFree() -> ReturnObject(univ);

         } else {

            obs = 1.;
         }
         break;

      case 1:

         obs = 1.;         
         break;

      default:

         errf(-1, "MCUserForDisaster :: UniversalObservable: choice not known");
         obs = 0.;  // :_TAWM_:
   }

   return obs;
}

// ----------------------------------------------------------------------------

// the event to determine the ``universal observable''
// (to obtain the weight to be returned to the adaptive routine)

double MCUserForDisaster :: CalculateAdaptationObservable(
                    Process* process,           
                    ContributionArray* ca,           
                    int pdf_collectionIn,             
                    int pdf_parametrizationIn,          
                    int pdf_setIn,          
                    int alphaSVariant_In,          
                    int alphaSOrder_In,          
                    double alphasLambdaQCD4_In,          
                    int alphaEMVariant_In           
                 ) {

   double* adaptObs = new double[n_events];

   // run through the list of all events
   UseListContainer <Event> * ec;
   Event* e;  
   Contribution* c;
   for (ec = process 
             -> GetEventList() 
             -> GetPtrToFirstContainer();
        ec != NULL && (e = ec -> GetDataPtr(), TRUE);
        ec = ec -> GetNext()) {
       for (c = (e -> GetContribution());
            c != NULL;
            c = c -> GetSameEvent()) {
           Event *breitE = e -> UserFrameFindEvent(breit);
           adaptObs[c -> GetLocation()] = UniversalObservable(breitE);
       }
   }
                  
   ConvoluteWithDistributions(
      process, ca,   
      QArray, QArray,
      noutArray, npdenArray,
      nfSwitchArray, nfSwitchArray,
      pdf_collectionIn,
      pdf_parametrizationIn,
      pdf_setIn,
      alphaSVariant_In,
      alphaSOrder_In,
      alphasLambdaQCD4_In,
      alphaEMVariant_In,
      QArray
   ); 

   double sum_L = 0;
   for (register int l = ca -> GetActual(); l >= 0; --l) {
       sum_L += tsum[l] * adaptObs[l];
   }

   delete [] adaptObs;

   return sum_L;
}

// ----------------------------------------------------------------------------

   // create the links between the event and the list of contributions

void MCUserForDisaster :: PhaseSpace(
                  Process* process, 
                  Event* event,
                  int& flag,
                  Contribution* con, 
                  int a
               ) {

   // :_TAWM_:
   process = process;
   a = a;

   flag = TRUE;

   // -------------------------------------------------------------------------

//   gl -> To_status() << "MCUserForDisaster :: PhaseSpace\n";
//   event -> Print(gl -> To_status());

   con -> SetEvent(event);
   if (event -> GetContribution() == NULL) 
      con -> SetSameEvent(NULL);
   else
      con -> SetSameEvent(event -> GetContribution());
   event -> SetContribution(con);
}

// ----------------------------------------------------------------------------

void MCUserForDisaster :: BoostToBreitFrame(
                             Process* process, 
                             Event* event,
                             ContributionArray* ca, 
                             int ifFinalRun
                          ) {

   // :_TAWM_:
   process = process; 
   ca = ca;
   ifFinalRun = ifFinalRun;

   // -------------------------------------------------------------------------

   // boost to the breit frame

   HorizontalList <Event> * firstInH = event -> GetH();

   // event record in the Breit frame
   breit_event
          = global
            -> GetEFree()
            -> GetDefinedObject(global -> disaster, event -> GetMaxNumber());
   event -> CopyInto(breit_event);
   breit_event -> SetUserFrame(breit);
   firstInH -> InsertAfterFirstElement(breit_event -> GetH());
   breit_event -> BoostTo_Breit();

//   gl -> To_status() << " MCUserForDisaster :: BREITFRAME\n";
//   event -> Print(gl -> To_status());

   // -------------------------------------------------------------------------

   // prepare parameters for the flavour factors

   Contribution* c;
   int iev;

   for (c = event -> GetContribution();
       c != NULL && (iev = c -> GetLocation(), TRUE);
       c = c -> GetSameEvent()) {

       // set some quatities to be used later
       c -> Set_nflavour_out_U(nOutFlavours);
       c -> Set_nflavour_pden_U(nPdenFlavours);
       c -> Set_ren_scale_U(event -> Q);
       c -> Set_fact_scale_U(event -> Q);
       c -> Set_xi_U(event -> xi);
       c -> Set_Q_U(event -> Q);
 
       QArray[iev] = event -> Q;
//       QArray[iev] = 2 * event -> Q;  // :_CAUTION_:
       noutArray[iev]  = nOutFlavours;
       npdenArray[iev] = nPdenFlavours;

// :_CAUTION_: fixed to 5 flavours? should not be 5 in general...
#if 0
       nfSwitchArray[iev] = alphaS_Server
                            -> FlavourSwitch(event -> Q);
#else
       nfSwitchArray[iev] = 5;
#endif
   }  // of for

}
   
// ----------------------------------------------------------------------------

// may be used to reject events based on the lepton phase space

void MCUserForDisaster :: IncidentPhaseSpace(
                  Process* process,
                  Event* event,
                  int& flag
               ) {

   // :_TAWM_:
   process = process;
   event = event;

   mcUserInterface -> IncidentPhaseSpace(
                         process,
                         event, 
                         flag 
                      );
}

// ----------------------------------------------------------------------------
//
// --> 
//
// ----------------------------------------------------------------------------
   
MCUserCalls :: MCUserCalls() {

   gl -> To_status() << "MCUserCalls :: MCUserCalls()\n";
}

// ----------------------------------------------------------------------------
   
MCUserCalls :: ~MCUserCalls() {

   gl -> To_status() << "MCUserCalls :: ~MCUserCalls()\n";
   Delete();
}

// ----------------------------------------------------------------------------
   
void MCUserCalls :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined()) 
      errf(-1, "MCUserCalls :: Define: already defined"); 

   gl -> To_status() << "MCUserCalls :: Define()\n";

   glob = globalIn;
}

// ----------------------------------------------------------------------------
   
void MCUserCalls :: Delete() {

   if (definedObject.Status()) {
      definedObject.ResetDefined();
      gl -> To_status() << "MCUserCalls :: Delete()\n";
   }
}

// ----------------------------------------------------------------------------
//
// --> 
//
// ----------------------------------------------------------------------------
   
MCUserCallsOther :: MCUserCallsOther() {

   gl -> To_status() << "MCUserCallsOther :: MCUserCallsOther()\n";
}

// ----------------------------------------------------------------------------
   
MCUserCallsOther :: ~MCUserCallsOther() {

   gl -> To_status() << "MCUserCallsOther :: ~MCUserCallsOther()\n";
   Delete();
}

// ----------------------------------------------------------------------------
   
void MCUserCallsOther :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined()) 
      errf(-1, "MCUserCallsOther :: Define: already defined"); 

   gl -> To_status() << "MCUserCallsOther :: Define()\n";

   glob = globalIn;
}

// ----------------------------------------------------------------------------
   
void MCUserCallsOther :: Delete() {

   if (definedObject.Status()) {
      definedObject.ResetDefined();
      gl -> To_status() << "MCUserCallsOther :: Delete()\n";
   }
}

// ----------------------------------------------------------------------------
//
// --> called from MCUserInterface
//
// ----------------------------------------------------------------------------
   
MCUserCalled :: MCUserCalled() {

   gl -> To_status() << "MCUserCalled :: MCUserCalled()\n";
}

// ----------------------------------------------------------------------------
   
MCUserCalled :: ~MCUserCalled() {

   gl -> To_status() << "MCUserCalled :: ~MCUserCalled()\n";
   Delete();
}

// ----------------------------------------------------------------------------
   
void MCUserCalled :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined()) 
      errf(-1, "MCUserCalled :: Define: already defined"); 

   gl -> To_status() << "MCUserCalled :: Define()\n";

   MCUserForDisaster :: Define(globalIn);
   MCUserCalls       :: Define(globalIn);
}

// ----------------------------------------------------------------------------
   
void MCUserCalled :: Delete() {

   if (definedObject.Status()) {

      definedObject.ResetDefined();

      gl -> To_status() << "MCUserCalled :: Delete()\n";

      MCUserForDisaster :: Delete();
      MCUserCalls       :: Delete();
   }
}

// ----------------------------------------------------------------------------

void MCUserCalled
     :: ConvoluteWithDistributions_User(
           Process* process,
           ContributionArray* ca,  
           double* d1,
           double* d2,
           int* nout,
           int* npden,
           int* nRS,
           int* nFS,
           int pdf_collectionIn,
           int pdf_parametrizationIn,
           int pdf_setIn,
           int alphaSVariant_In,
           int alphaSOrder_In,
           double alphasLambdaQCD4_In,
           int alphaEMVariant_In,
           double* alphaEMScale_In
        ) {

   mcUserInterface
   -> ConvoluteWithDistributions_User(
         process,
         ca,
         d1,
         d2,
         nout, 
         npden,
         nRS,
         nFS,
         pdf_collectionIn,
         pdf_parametrizationIn,
         pdf_setIn,
         alphaSVariant_In,
         alphaSOrder_In,
         alphasLambdaQCD4_In,
         alphaEMVariant_In,
         alphaEMScale_In
      );
}

// ----------------------------------------------------------------------------

double MCUserCalled :: CalculateAdaptationObservable_User(
                          Process* process,           
                          ContributionArray* ca,           
                          int pdf_collectionIn,             
                          int pdf_parametrizationIn,          
                          int pdf_setIn,          
                          int alphaSVariant_In,          
                          int alphaSOrder_In,          
                          double alphasLambdaQCD4_In,          
                          int alphaEMVariant_In           
                       ) {

   return
   mcUserInterface
   -> CalculateAdaptationObservable_User(
         process,
         ca,
         pdf_collectionIn,
         pdf_parametrizationIn,
         pdf_setIn,
         alphaSVariant_In, 
         alphaSOrder_In,
         alphasLambdaQCD4_In,
         alphaEMVariant_In
      );
}

// ----------------------------------------------------------------------------
//
// --> the user routine called to implement pre/post conditions
//
// ----------------------------------------------------------------------------
   
MCUserInterface :: MCUserInterface() {

   gl -> To_status() << "MCUserInterface :: MCUserInterface()\n";
}

// ----------------------------------------------------------------------------
   
MCUserInterface :: ~MCUserInterface() {

   gl -> To_status() << "MCUserInterface :: ~MCUserInterface()\n";
   Delete();
}

// ----------------------------------------------------------------------------
   
void MCUserInterface :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined())
      errf(-1, "MCUserInterface :: Define: already defined");
  
   gl -> To_status() << "MCUserInterface :: Define()\n";
   MCUser :: Define(globalIn);
   MCUserCallsOther :: Define(globalIn);
}

// ----------------------------------------------------------------------------
   
void MCUserInterface :: Delete() {

   if (definedObject.Status()) {

      definedObject.ResetDefined();

      gl -> To_status() << "MCUserInterface :: Delete()\n";

      MCUser :: Delete();
      MCUserCallsOther :: Delete();
   }
}

// ----------------------------------------------------------------------------
   
void MCUserInterface :: Start() {

   mcUserCalled -> Start_User();
}

// ----------------------------------------------------------------------------
   
void MCUserInterface :: End() {

   mcUserCalled -> End_User();
}

// ----------------------------------------------------------------------------

int MCUserInterface :: SetParameter(char* pname, char* parameter) {

   return mcUserCalled -> SetParameter_User(pname, parameter);
}

// ----------------------------------------------------------------------------

void MCUserInterface :: StartIntegration1() {

   mcUserCalled -> StartIntegration1_User();
}

// ----------------------------------------------------------------------------

void MCUserInterface :: StartIntegration2() {

   mcUserCalled -> StartIntegration2_User();
}

// ----------------------------------------------------------------------------

void MCUserInterface :: EndIntegration(double average, double error) {

   mcUserCalled -> EndIntegration_User(average, error);
}

// ----------------------------------------------------------------------------

void MCUserInterface :: BeginEvent(Process* process) {

   mcUserCalled -> BeginEvent_User(process);
}

// ----------------------------------------------------------------------------

void MCUserInterface :: EndOfEvent(Process* process) {

   mcUserCalled -> EndOfEvent_User(process);
}

// ----------------------------------------------------------------------------

double MCUserInterface :: EndEvent(
                         Process* process,
                         ContributionArray* ca,
                         int ifFinalRun
                      ) {

   return mcUserCalled -> EndEvent_User(
                             process,
                             ca,
                             ifFinalRun
                          );
}

// ----------------------------------------------------------------------------

void MCUserInterface :: DropEvent(Process* process) {

   mcUserCalled -> DropEvent_User(process);
}

// ----------------------------------------------------------------------------

void MCUserInterface :: AcceptEvent(Process* process) {

   mcUserCalled -> AcceptEvent_User(process);
}

// ----------------------------------------------------------------------------

// may be used to reject events based on the lepton phase space

void MCUserInterface :: IncidentPhaseSpace(
                       Process* process,
                       Event* event,
                       int& flag
                    ) {

   mcUserCalled -> IncidentPhaseSpace_User(
                      process,
                      event,
                      flag 
                    );
}

// ----------------------------------------------------------------------------

void MCUserInterface 
     :: ConvoluteWithDistributions_User(
           Process* process,
           ContributionArray* ca,  
           double* d1,
           double* d2,
           int* nout,
           int* npden,
           int* nRS,
           int* nFS,
           int pdf_collectionIn,
           int pdf_parametrizationIn,
           int pdf_setIn,
           int alphaSVariant_In,
           int alphaSOrder_In,
           double alphasLambdaQCD4_In,
           int alphaEMVariant_In,
           double* alphaEMScale_In
        ) {

   mcUserForDisaster
   -> ConvoluteWithDistributions(
         process,
         ca,
         d1,
         d2,
         nout, 
         npden,
         nRS,
         nFS,
         pdf_collectionIn,
         pdf_parametrizationIn,
         pdf_setIn,
         alphaSVariant_In,
         alphaSOrder_In,
         alphasLambdaQCD4_In,
         alphaEMVariant_In,
         alphaEMScale_In
      );
}

// ----------------------------------------------------------------------------

double MCUserInterface 
       :: CalculateAdaptationObservable_User(
             Process* process,           
             ContributionArray* ca,           
             int pdf_collectionIn,             
             int pdf_parametrizationIn,          
             int pdf_setIn,          
             int alphaSVariant,          
             int alphaSOrder,          
             double alphasLambdaQCD4,          
             int alphaEMVariant           
          ) {

   return
   mcUserForDisaster
   -> CalculateAdaptationObservable(
         process,
         ca,
         pdf_collectionIn,
         pdf_parametrizationIn,
         pdf_setIn,
         alphaSVariant, 
         alphaSOrder,
         alphasLambdaQCD4,
         alphaEMVariant
      );
}

// ----------------------------------------------------------------------------
//
// --> the user routine called to implement pre/post conditions
//     for the FORTRAN interface
//
// ----------------------------------------------------------------------------
   
MCUserInterface_FORTRAN :: MCUserInterface_FORTRAN() {

   gl -> To_status() 
            << "MCUserInterface_FORTRAN :: MCUserInterface_FORTRAN()\n";
}

// ----------------------------------------------------------------------------
   
MCUserInterface_FORTRAN :: ~MCUserInterface_FORTRAN() {

   gl -> To_status() 
            << "MCUserInterface_FORTRAN :: ~MCUserInterface_FORTRAN()\n";
   Delete();
}

// ----------------------------------------------------------------------------
   
void MCUserInterface_FORTRAN :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined())
      errf(-1, "MCUserInterface_FORTRAN :: Define: already defined");
    
   gl -> To_status() << "MCUserInterface_FORTRAN :: Define()\n";
   MCUserInterface :: Define(globalIn);

   changed_print_counter.Set_maximum(5);
   changed_print_counter.Reset();
}

// ----------------------------------------------------------------------------
   
void MCUserInterface_FORTRAN :: Delete() {

   if (definedObject.Status()) {

      definedObject.ResetDefined();

      gl -> To_status() << "MCUserInterface_FORTRAN :: Delete()\n";

      MCUserInterface :: Delete();
   }
}

// ----------------------------------------------------------------------------
   
void MCUserInterface_FORTRAN :: Start() {

   Control_To_FORTRAN();
   mcUserCalled -> Start_User();
   Control_To_C();
}

// ----------------------------------------------------------------------------
   
void MCUserInterface_FORTRAN :: End() {

   Control_To_FORTRAN();
   mcUserCalled -> End_User();
   Control_To_C();
}

// ----------------------------------------------------------------------------

int MCUserInterface_FORTRAN :: SetParameter(char* pname, char* parameter) {

   Control_To_FORTRAN();
   int ret = mcUserCalled -> SetParameter_User(pname, parameter);
   Control_To_C();
   return ret;
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN :: StartIntegration1() {

   Control_To_FORTRAN();
   mcUserCalled -> StartIntegration1_User();
   Control_To_C();
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN :: StartIntegration2() {

   Control_To_FORTRAN();
   mcUserCalled -> StartIntegration2_User();
   Control_To_C();
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN :: EndIntegration(double average, double error) {

   Control_To_FORTRAN();
   mcUserCalled -> EndIntegration_User(average, error);
   Control_To_C();
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN :: BeginEvent(Process* process) {

   Control_To_FORTRAN();
   mcUserCalled -> BeginEvent_User(process);
   Control_To_C();
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN :: EndOfEvent(Process* process) {

   Control_To_FORTRAN();
   mcUserCalled -> EndOfEvent_User(process);
   Control_To_C();
}

// ----------------------------------------------------------------------------

double MCUserInterface_FORTRAN 
       :: EndEvent(
             Process* process,
             ContributionArray* ca,
             int ifFinalRun
          ) {

   Control_To_FORTRAN();

   double weight 
      = mcUserCalled -> EndEvent_User(
                           process,
                           ca,
                           ifFinalRun
                        );
   Control_To_C();

   return weight;
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN 
     :: DropEvent(Process* process) {

   Control_To_FORTRAN();
   mcUserCalled -> DropEvent_User(process);
   Control_To_C();
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN 
     :: AcceptEvent(Process* process) {

   Control_To_FORTRAN();
   mcUserCalled -> AcceptEvent_User(process);
   Control_To_C();
}

// ----------------------------------------------------------------------------

// may be used to reject events based on the lepton phase space

void MCUserInterface_FORTRAN 
       :: IncidentPhaseSpace(
             Process* process,
             Event* event,
             int& flag
          ) {

   Control_To_FORTRAN();
   mcUserCalled -> IncidentPhaseSpace_User(
                      process,
                      event,
                      flag 
                    );
   Control_To_C();
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN
     :: ConvoluteWithDistributions_User(
           Process* process,
           ContributionArray* ca,  
           double* d1,
           double* d2,
           int* nout,
           int* npden,
           int* nRS,
           int* nFS,
           int pdf_collectionIn,
           int pdf_parametrizationIn,
           int pdf_setIn,
           int alphaSVariant,
           int alphaSOrder,
           double alphasLambdaQCD4,
           int alphaEMVariant,
           double* alphaEMScale
        ) {

   Control_To_C();

   mcUserForDisaster
   -> ConvoluteWithDistributions(
         process,
         ca,
         d1,
         d2,
         nout, 
         npden,
         nRS,
         nFS,
         pdf_collectionIn,
         pdf_parametrizationIn,
         pdf_setIn,
         alphaSVariant,
         alphaSOrder,
         alphasLambdaQCD4,
         alphaEMVariant,
         alphaEMScale
      );

   Control_To_FORTRAN();
}

// ----------------------------------------------------------------------------

double MCUserInterface_FORTRAN
       :: CalculateAdaptationObservable_User(
             Process* process,           
             ContributionArray* ca,           
             int pdf_collectionIn,             
             int pdf_parametrizationIn,          
             int pdf_setIn,          
             int alphaSVariant,          
             int alphaSOrder,          
             double alphasLambdaQCD4,          
             int alphaEMVariant           
          ) {

   Control_To_C();

   double obs =
   mcUserForDisaster
   -> CalculateAdaptationObservable(
         process,
         ca,
         pdf_collectionIn,
         pdf_parametrizationIn,
         pdf_setIn,
         alphaSVariant, 
         alphaSOrder,
         alphasLambdaQCD4,
         alphaEMVariant
      );

   Control_To_FORTRAN();

   return obs;
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN
     :: Control_To_FORTRAN() {

//   global -> To_status() 
//      << "MCUserInterface_FORTRAN :: Control_To_FORTRAN\n";

     fpnan -> Hold();
}

// ----------------------------------------------------------------------------

void MCUserInterface_FORTRAN
     :: Control_To_C() {

//   global -> To_status() 
//      << "MCUserInterface_FORTRAN :: Control_To_C\n";

   Boolean changed; 
//   changed = False;
   changed = fpnan -> Continue();

   if (changed) {
      if (changed_print_counter.IncrementAndCheck()) {   
         global -> To_error()
            << "MCUserInterface_FORTRAN :: Control_To_C: \n"
            << "floating point status has changed in the user routine!\n";
         fpnan -> PrintStatus(global -> To_error());
         fpnan -> PrintSaveStatus(global -> To_error());
      }   
      if (changed_print_counter.BorderlineCase())   
         global -> To_error() 
            << "MCUserInterface_FORTRAN :: Control_To_C: "
            << "no more warning messages\n";   
   }
}

// ----------------------------------------------------------------------------
//
// --> my own user routines
//
// ----------------------------------------------------------------------------
   
MyUserClass :: MyUserClass() {

   gl -> To_status() << "MyUserClass :: MyUserClass()\n";

   mcUserInterface = new MCUserInterface;
   mcUserInterface -> Set_mcUserCalled(this);
   mcUserInterface -> Set_mcUserForDisaster(this);
}

// ----------------------------------------------------------------------------
   
MyUserClass :: ~MyUserClass() {

   gl -> To_status() << "MyUserClass :: ~MyUserClass()\n";
   Delete();

   delete mcUserInterface;
}

// ----------------------------------------------------------------------------
   
void MyUserClass :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined())
      errf(-1, "MyUserClass :: Define: already defined");

   gl -> To_status() << "MyUserClass :: Define()\n";

   MCUserCalled :: Define(globalIn);  

   // initialization of parameters
   cluster_flag=TRUE;
   cluster_algorithm=NONE;
   min_combine=0;
   max_combine=1;
   cluster_scale_choice=0;
   cluster_cut=0.;
   ren_scale_choice=0;
   fact_scale_choice=0;
   decisionMode=FALSE;
   minClusterNumber=-1;
   maxClusterNumber=-1;
   cat = JADE_algorithm;
   rt = E_Scheme;

   ratioProtonElectronInLaboratory = 1.;

   testObservables           = 0;
   observablesA              = 1;
   observablesB              = 1;
   observablesC              = 1;
   observablesE              = 0;
}

// ----------------------------------------------------------------------------
   
void MyUserClass :: Delete() {

   if (definedObject.Status()) {

      definedObject.ResetDefined();

      gl -> To_status() << "MyUserClass :: Delete()\n";

      MCUserCalled :: Delete();
   }
}

// ----------------------------------------------------------------------------
   
void MyUserClass :: Start_User() {
}

// ----------------------------------------------------------------------------
   
void MyUserClass :: End_User() {
}

// ----------------------------------------------------------------------------

// called if a parameter name parsed by the MC cannot be resolved
// --> assumed to be a user parameter

int MyUserClass :: SetParameter_User(char* pname, char* parameter) {

   int ret = 0;

   if (!strcmp(pname,"CLUSTER_FLAG")) {
      cluster_flag=ReadBool(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"CLUSTER_ALGORITHM")) {
      long inclu;
      inclu=ReadInt(global,pname,parameter);
      switch (inclu) {
         case 0:
            cluster_algorithm=NONE;
            break;
         case 1:
            cluster_algorithm=LORENTZ_INVARIANT_1;
            break;
         case 2:
            cluster_algorithm=LORENTZ_INVARIANT_OLD_1;
            break;
         default:
            errf(-1,"MCUser::SetParameter: cluster algorithm not known");
      }
   }
   else

   if (!strcmp(pname,"MIN_COMBINE")) {
      min_combine=ReadInt(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"MAX_COMBINE")) {
      max_combine=ReadInt(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"CLUSTER_SCALE_CHOICE")) {
      cluster_scale_choice=ReadInt(global,pname,parameter);
   } 
   else

   if (!strcmp(pname,"CLUSTER_CUT")) {
      cluster_cut=ReadDouble(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"CLUSTER_ALGORITHM_TYPE")) {
      long inclu;
      inclu=ReadInt(global,pname,parameter);
      switch (inclu) {
         case 1:
            cat=JADE_algorithm;
            break;
         case 2:
            cat=kT_algorithm;
            break;
         default:
            errf(-1,"MCUser::SetParameter: cat not known");
      }
   }
   else

   if (!strcmp(pname,"RECOMBINATION_TYPE")) {
      long inclu;
      inclu=ReadInt(global,pname,parameter);
      switch (inclu) {
         case 1:
            rt=JADE_Scheme;
            break;
         case 2:
            rt=E_Scheme;
            break;
         case 3:
            rt=E0_Scheme;
            break;
         case 4:
            rt=P_Scheme;
            break;
         case 5:
            rt=P0_Scheme;
            break;
         default:
            errf(-1,"MCUser::SetParameter: rt not known");
      }
   }
   else

   if (!strcmp(pname,"REN_SCALE_CHOICE")) {
      ren_scale_choice=ReadInt(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"FACT_SCALE_CHOICE")) {
      fact_scale_choice=ReadInt(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"TEST_CASE")) {
      testCase=ReadBool(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"DECISION_MODE")) {
      decisionMode=ReadInt(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"MIN_CLUSTER_NUMBER")) {
      minClusterNumber=ReadInt(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"MAX_CLUSTER_NUMBER")) {
      maxClusterNumber=ReadInt(global,pname,parameter);
   }
   else

   if (!strcmp(pname,"RATIO_PROTON_ELECTRON_IN_LABORATORY")) {
      ratioProtonElectronInLaboratory = ReadDouble(global, pname, parameter);
   }
   else

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   // parton densities, alpha_s and alpha_em

   if (!strcmp(pname, "PDF_COLLECTION")) {
      pdf_collection = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "PDF_PARAMETRIZATION")) {
      pdf_parametrization = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "PDF_SET")) {
      pdf_set = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "ALPHA_S_VARIANT")) {
      variant_alpha_s = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "ALPHA_S_ORDER")) {
      order_alpha_s = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "USER_DEFINED_ALPHA_S")) {
      user_defined_alpha_s = ReadBool(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "LAMBDA_4_ALPHA_S")) {
      lambda_4_alpha_s = ReadDouble(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "ALPHA_EM_VARIANT")) {
      variant_alpha_em = ReadInt(global, pname, parameter);
   }
   else

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   // observables

   if (!strcmp(pname, "TEST_OBSERVABLES")) {
      testObservables = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "OBSERVABLES_A")) {
      observablesA = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "OBSERVABLES_B")) {
      observablesB = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "OBSERVABLES_C")) {
      observablesC = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "OBSERVABLES_E")) {
      observablesE = ReadInt(global, pname, parameter);
   }
   else

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   // flavours

   if (!strcmp(pname, "N_OUT_FLAVOURS")) {
      nOutFlavours = ReadInt(global, pname, parameter);
   }
   else

   if (!strcmp(pname, "N_PDEN_FLAVOURS")) {
      nPdenFlavours = ReadInt(global, pname, parameter);
   }
   else

   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ret = -1;

   return ret;
}

// ----------------------------------------------------------------------------

void MyUserClass :: StartIntegration1_User() {

   // have to do it here to define n_events
//   MCUser :: StartIntegration1();

   // define set of observables
   switch ( global -> numberOfFinalStatePartonsInBornTerm ) {
      case 1:
         obsSet = 1;
         minScaleJets = 2;
         maxScaleJets = 3;
  //     minScaleJets = 2;
  //     maxScaleJets = 2;
         minJets = 2;
         maxJets = 2;
         minPtJets = 2;
         maxPtJets = 2;
         break;
      case 2:
         obsSet = 2;
         minScaleJets = 3;
         maxScaleJets = 3;
         minJets = 3;
         maxJets = 3;
         minPtJets = 3;
         maxPtJets = 3;
         break;
      case 3:
         obsSet = 2;
         minScaleJets = 4;
         maxScaleJets = 4;
         minJets = 4;
         maxJets = 4;
         minPtJets = 4;
         maxPtJets = 4;
         break;
      default:
         errf(-1, "MCUser::StartIntegration1: obsSet not defined");
   }

   global -> To_status() << "obsSet=" << FS("%d\n", obsSet);

   rsfsscale = TRUE;

   if (   observablesA
       || observablesB
       || observablesC
       || observablesE
       || testObservables
       || universalObservableChoice < 0
      )
      library = new Library(global, 50);
   else
      library = NULL;
   
   if (observablesA) {

      if (obsSet == 1)
         ( total = library -> CreateNewBook() )
             -> Define(Graph, Linear, 1, 0.0, 0.0,
                       "total", "total", 
                       n_events);

   }  // of ``observablesA''
   
   if (   observablesA
       || observablesC) {

      ( t_total = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_total", "t_total", 
                    n_events);
   
      ( t_q = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_q", "t_q", 
                    n_events);
   
      ( t_g = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_g", "t_g", 
                    n_events);
   
   }  // of ``observablesA || observablesC''
   
   if (   observablesA) {

      ( JADE_E_hCMS_cut = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "JADE_E_hCMS_cut", "JADE_E_hCMS_cut", 
                    n_events);
   
#if 0
      ( JADE_E0_hCMS_cut = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "JADE_E0_hCMS_cut", "JADE_E0_hCMS_cut", 
                    n_events);
#endif
   
      ( JADE_P_hCMS_cut = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "JADE_P_hCMS_cut", "JADE_P_hCMS_cut", 
                    n_events);
   
      ( JADE_JADE_hCMS_cut = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "JADE_JADE_hCMS_cut", "JADE_JADE_hCMS_cut", 
                    n_events);
   
      ( kT_E_Breit_cut = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.1, 1.0,
                    "kT_E_Breit_cut", "kT_E_Breit_cut", 
                    n_events);
   
      ( JADE_EBreit_scale = library -> CreateNewBook() )
          -> Define(Graph, Log, 5, 0.5, 2.0,
                    "JADE_EBreit_scale", "JADE_EBreit_scale", 
                    n_events);
   
      ( JADE_Q_scale = library -> CreateNewBook() )
          -> Define(Graph, Log, 5, 0.5, 2.0,
                    "JADE_Q_scale", "JADE_Q_scale", 
                    n_events);
   
      if (rsfsscale) {
         ( JADE_Q_rs_scale = library -> CreateNewBook() )
             -> Define(Graph, Log, 5, 0.5, 2.0,
                       "JADE_Q_rs_scale", "JADE_Q_rs_scale", 
                       n_events);
   
         ( JADE_Q_fs_scale = library -> CreateNewBook() )
             -> Define(Graph, Log, 5, 0.5, 2.0,
                       "JADE_Q_fs_scale", "JADE_Q_fs_scale", 
                       n_events);
      }
   
      ( kT_Q_scale = library -> CreateNewBook() )
          -> Define(Graph, Log, 5, 0.5, 2.0,
                    "kT_Q_scale", "kT_Q_scale", 
                    n_events);
   
      if (rsfsscale) {
         ( kT_Q_rs_scale = library -> CreateNewBook() )
             -> Define(Graph, Log, 5, 0.5, 2.0,
                       "kT_Q_rs_scale", "kT_Q_rs_scale", 
                       n_events);
   
         ( kT_Q_fs_scale = library -> CreateNewBook() )
             -> Define(Graph, Log, 5, 0.5, 2.0,
                       "kT_Q_fs_scale", "kT_Q_fs_scale", 
                       n_events);
      }
   
      ( pT_hCMS = library -> CreateNewBook() )
          -> Define(Histogram, Linear, 60, 0.0, 50,
                    "pT_hCMS", "pT_hCMS", 
                    n_events);
   
      ( kT_Breit = library -> CreateNewBook() )
          -> Define(Histogram, Linear, 60, 0.0, 100,
                    "kT_Breit", "kT_Breit", 
                    n_events);
   
      ( Thrust_Breit = library -> CreateNewBook() )
          -> Define(Histogram, Linear, 60, 0.0, 1.0,
                    "Thrust_Breit", "Thrust_Breit", 
                    n_events);
   
      ( energy_Breit = library -> CreateNewBook() )
          -> Define(Histogram, Linear, 60, 0.0, 100,
                    "energy_Breit", "energy_Breit", 
                    n_events);
   
   }  // of ``observablesA''

   if (   observablesA
       || observablesB
       || universalObservableChoice == -1
      ) {

      ( t_cut_total = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_cut_total", "t_cut_total", 
                    n_events);
   
      ( t_cut_q = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_cut_q", "t_cut_q", 
                    n_events);
   
      ( t_cut_g = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_cut_g", "t_cut_g", 
                    n_events);
   
      ( t_cut_old_total = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_cut_old_total", "t_cut_old_total", 
                    n_events);
   
      ( t_cut_old_q = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_cut_old_q", "t_cut_old_q", 
                    n_events);
   
      ( t_cut_old_g = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "t_cut_old_g", "t_cut_old_g", 
                    n_events);
   
      ( three_q = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "three_q", "three_q", 
                    n_events);
   
      ( three_g = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "three_g", "three_g", 
                    n_events);
   
      ( direct_sum = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "direct_sum", "direct_sum", 
                    n_events);
   
   }  // of ``observablesA || observablesB || ...''
   
   if (testObservables) {

      ( test2History = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test2History", "test2History", 
                    n_events);

      ( test2Other = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test2Other", "test2Other", 
                    n_events);
   
      ( test2Old = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test2Old", "test2Old", 
                    n_events);
   
      ( test3History = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test3History", "test3History", 
                    n_events);
   
      ( test3Other = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test3Other", "test3Other", 
                    n_events);
   
      ( test3Old = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test3Old", "test3Old", 
                    n_events);

      ( test4History = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test4History", "test4History", 
                    n_events);

      ( test4Other = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test4Other", "test4Other", 
                    n_events);

      ( test4Old = library -> CreateNewBook() )
          -> Define(Graph, Log, 10, 0.01, 0.1,
                    "test4Old", "test4Old", 
                    n_events);

   }  // of ``testObservables''

   if (   observablesE
       || universalObservableChoice == -2
      ) {

      ( kt_cmp = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "kt_cmp", "kt_cmp", 
                    n_events);
   
      ( kt_cut_cmp = library -> CreateNewBook() )
          -> Define(Graph, Linear, 1, 0.0, 0.0,
                    "kt_cut_cmp", "kt_cut_cmp", 
                    n_events);
   
   }  // of ``observablesE ...''
   
   if (library != NULL)
      library -> ResetStorage();
}

// ----------------------------------------------------------------------------

void MyUserClass :: StartIntegration2_User() {
}

// ----------------------------------------------------------------------------

void MyUserClass :: EndIntegration_User(double average, double error) {

   // :_TAWM_:
   average = average;
   error   = error;

   // -------------------------------------------------------------------------

   if (library != NULL) {
      library -> CalculateErrors();
      library -> Print(global -> To_status());
//      library -> PrintMathForm(global -> PlotOut);
   }

   double averageQ;
   double errorQ;
   double relErrorQ;
   double averageG;
   double errorG;
   double relErrorG;

   if (   observablesA
       || observablesB
       || universalObservableChoice == -1
      ) {

      t_cut_q -> GetEntry(-1, &averageQ, &errorQ, &relErrorQ);
      t_cut_g -> GetEntry(-1, &averageG, &errorG, &relErrorG);
 
      global -> To_status() 
                   << FS("sigma =  %10.3e",  averageQ)
                   << FS(" +-%10.3e (Q)  ",  errorQ)
                   << FS("%10.3e",           averageG)
                   << FS(" +-%10.3e (G)\n",  errorG);
   }

   if (   observablesA
       || observablesC
      ) {

      t_q -> GetEntry(-1, &averageQ, &errorQ, &relErrorQ);
      t_g -> GetEntry(-1, &averageG, &errorG, &relErrorG);

      global -> To_status() 
                   << FS("nocut =  %10.3e",  averageQ)
                   << FS(" +-%10.3e (Q)  ",  errorQ)
                   << FS("%10.3e",           averageG)
                   << FS(" +-%10.3e (G)\n",  errorG);
   }

   if (   observablesA
       || observablesB
       || universalObservableChoice == -1
      ) {

      three_q -> GetEntry(-1, &averageQ, &errorQ, &relErrorQ);
      three_g -> GetEntry(-1, &averageG, &errorG, &relErrorG);

      global -> To_status() 
                   << FS("three =  %10.3e",  averageQ)
                   << FS(" +-%10.3e (Q)  ",  errorQ)
                   << FS("%10.3e",           averageG)
                   << FS(" +-%10.3e (G)\n",  errorG);
   }

   if (library != NULL)
      delete library;

//   MCUser :: EndIntegration(average, error);
}

// ----------------------------------------------------------------------------

void MyUserClass :: BeginEvent_User(Process* process) {

   // :_TAWM_:
   process = process;

//   MCUser :: BeginEvent(process);

   if (library != NULL)
      if (global -> finalrun || universalObservableChoice < 0)
         library -> PrepareForEvent();
}

// ----------------------------------------------------------------------------

void MyUserClass :: EndOfEvent_User(Process* process) {

   // :_TAWM_:
   process = process;
}

// ----------------------------------------------------------------------------

double MyUserClass :: EndEvent_User(
                         Process* process,
                         ContributionArray* ca,
                         int ifFinalRun
                      ) {

//   gl -> To_status() << "MyUserClass :: EndEvent_User\n";

//   for (int iii = 0; iii < 100; ++iii) 
//       provoke_floating_point_exception();

   double ret = 0; // :_CHECK_:

   // :_TAWM_:
   process = process;
   ifFinalRun = ifFinalRun;

   register int l;

//   ca -> Print(stdout, 2);
//   ca -> Print(global -> To_status(), 2);

   // -------------------------------------------------------------------------

   // calculate observables

   UseListContainer <Event> * ec;
   Event* e;
   // run through the list of all events
   for (ec = process -> GetEventList() -> GetPtrToFirstContainer();
        ec != NULL && (e = ec -> GetDataPtr(), TRUE);
        ec = ec -> GetNext()) {
//       gl -> To_status() << "ec = ec -> GetNext()\n";
       if (e -> GetContribution() != NULL)
          CalculateObservablesForEvent_User(process, e, ca, ifFinalRun);
   }

   // -------------------------------------------------------------------------

   // alpha_s
   if (user_defined_alpha_s) {
      alphaSOrder = order_alpha_s;
      alphasLambdaQCD4 = lambda_4_alpha_s;
   } else {
      alphaSOrder = global
                    -> disaster
                    -> process 
                    -> GetAlphaS_Order();
      alphasLambdaQCD4 = pdf_Server 
                         -> GetPDFLambda_4();
   }

   // -------------------------------------------------------------------------

   register int i;

   Book* book;

   if (       ifFinalRun
           && (
                  testObservables
               || observablesA
               || observablesB
               || observablesC
               || observablesE
              )
       || universalObservableChoice < 0
      ) {

         // ----------------------------------------------------------------------
         //
         // event with fs & rs == Q
         //
         // ----------------------------------------------------------------------
   
         for (l = ca -> GetActual(); l >= 0; --l) {
             nRSArray[0][l] = alphaS_Server
                                   -> FlavourSwitch(QArray[l]);
             nRSArray[0][l] = 5; // :_CAUTION_: just for test purposes!
             nFSArray[0][l] = nRSArray[0][l]; 
         }

         ConvoluteWithDistributions_User(
            process, ca,   
            QArray, QArray, 
            noutArray, npdenArray,
            nRSArray[0], nFSArray[0],
            pdf_collection,
            pdf_parametrization,
            pdf_set,
            variant_alpha_s,           
            alphaSOrder,            
            alphasLambdaQCD4,
            variant_alpha_em,
            QArray
         ); 

   }  // of ``ifFinalRun...''

   if (       ifFinalRun
           && (
                  observablesA
               || observablesB
              )
       || universalObservableChoice == -1
      ) {

         for (l = ca -> GetActual(); l >= 0; --l) {
             t_cut_total -> SPS_Graph( l, 0, tsum  [l] );
             t_cut_q     -> SPS_Graph( l, 0, tsumq [l] );
             t_cut_g     -> SPS_Graph( l, 0, tsumg [l] );
             t_cut_old_total -> SPS_Graph( l, 0, tsum  [l] );
             t_cut_old_q     -> SPS_Graph( l, 0, tsumq [l] );
             t_cut_old_g     -> SPS_Graph( l, 0, tsumg [l] );
             three_q         -> SPS_Graph( l, 0, tsumq [l] );
             three_g         -> SPS_Graph( l, 0, tsumg [l] );
             direct_sum      -> SPS_Graph( l, 0, tsum  [l] );
         }

         if (universalObservableChoice == -1)
            ret = t_cut_total -> GetBeforeSum(0);

   }   

   if (      ifFinalRun
          && observablesE
       || universalObservableChoice == -2
      ) {

         for (l = ca -> GetActual(); l >= 0; --l) {
             kt_cmp     -> SPS_Graph( l, 0, tsum[l] );
             kt_cut_cmp -> SPS_Graph( l, 0, tsum[l] );
         }

         if (universalObservableChoice == -2)
            ret = kt_cmp -> GetBeforeSum(0);
      }

   if (ifFinalRun) {

      if (observablesA) {
   
         if (obsSet == 1) {
            book = total;
            for (l = ca -> GetActual(); l >= 0; --l)
                book -> SPS_Graph(l, 0, tsum[l]);
         }
   
         book = JADE_E_hCMS_cut;
         for (i = 0; i < book -> GetSize(); ++i)
             for (l = ca -> GetActual(); l >= 0; --l)
                 book -> SPS_Graph(l, i, tsum[l]);

      }

      if (   observablesA
          || observablesC) {
   
         for (l = ca -> GetActual(); l >= 0; --l) {
             t_total -> SPS_Graph( l, 0, tsum  [l] );
             t_q     -> SPS_Graph( l, 0, tsumq [l] );
             t_g     -> SPS_Graph( l, 0, tsumg [l] );
         }

      }   

      if (testObservables) {

         book = test2History;
         for (i = 0; i < book -> GetSize(); ++i)
             for (l = ca -> GetActual(); l >= 0; --l) {
                 test2History -> SPS_Graph(l, i, tsum[l]);
                 test2Other   -> SPS_Graph(l, i, tsum[l]);
                 test2Old     -> SPS_Graph(l, i, tsum[l]);
                 test3History -> SPS_Graph(l, i, tsum[l]);
                 test3Other   -> SPS_Graph(l, i, tsum[l]);
                 test3Old     -> SPS_Graph(l, i, tsum[l]);
                 test4History -> SPS_Graph(l, i, tsum[l]);
                 test4Other   -> SPS_Graph(l, i, tsum[l]);
                 test4Old     -> SPS_Graph(l, i, tsum[l]);
             }
      }

      if (observablesA) {

#if 0
         book = JADE_E0_hCMS_cut;
         for (i = 0; i < book -> GetSize(); ++i)
             for (l = ca -> GetActual(); l >= 0; --l)
                 book -> SPS_Graph(l, i, tsum[l]);
#endif
   
         book = JADE_P_hCMS_cut;
         for (i = 0; i < book -> GetSize(); ++i)
             for (l = ca -> GetActual(); l >= 0; --l)
                 book -> SPS_Graph(l, i, tsum[l]);
   
         book = JADE_JADE_hCMS_cut;
         for (i = 0; i < book -> GetSize(); ++i)
             for (l = ca -> GetActual(); l >= 0; --l)
                 book -> SPS_Graph(l, i, tsum[l]);
   
         book = kT_E_Breit_cut;
         for (i = 0; i < book -> GetSize(); ++i)
             for (l = ca -> GetActual(); l >= 0; --l)
                 book -> SPS_Graph(l, i, tsum[l]);
   
         book = pT_hCMS;
         for (l = ca -> GetActual(); l >= 0; --l)
             book -> SPS_Histogram(l, tsum[l]);
   
         book = kT_Breit;
         for (l = ca -> GetActual(); l >= 0; --l)
             book -> SPS_Histogram(l, tsum[l]);
   
         book = Thrust_Breit;
         for (l = ca -> GetActual(); l >= 0; --l)
             book -> SPS_Histogram(l, tsum[l]);
   
         book = energy_Breit;
         for (l = ca -> GetActual(); l >= 0; --l)
             book -> SPS_Histogram(l, tsum[l]);
   
         // -------------------------------------------------------------------
         //
         // scale study
         //
         // -------------------------------------------------------------------
   
         book = JADE_EBreit_scale;
         for (i = 0; i < book -> GetSize(); ++i) {
             for (l = ca -> GetActual(); l >= 0; --l) {
                 rsArray[0][l] = dmax(EBreitArray[l] * book -> GetX(i), 2.5);
                 fsArray[0][l] = rsArray[0][l];
                 nRSArray[0][l] = alphaS_Server
                                       -> FlavourSwitch(rsArray[0][l]);
                 nFSArray[0][l] = alphaS_Server
                                       -> FlavourSwitch(fsArray[0][l]);
             }
             ConvoluteWithDistributions_User(
                process, ca,   
                rsArray[0], fsArray[0], 
                noutArray, npdenArray,
                nRSArray[0], nFSArray[0],
                pdf_collection,
                pdf_parametrization,
                pdf_set,
                variant_alpha_s,           
                alphaSOrder,            
                alphasLambdaQCD4, 
                variant_alpha_em,
                QArray
             ); 
             for (l = ca -> GetActual(); l >= 0; --l)
                 book -> SPS_Graph(l, i, tsum[l]);
         }
   
         // -------------------------------------------------------------------
   
         book = JADE_Q_scale;
         for (i = 0; i < book -> GetSize(); ++i) {
             for (l = ca -> GetActual(); l >= 0; --l) {
                 rsArray[0][l] = QArray[l] * book -> GetX(i);
                 fsArray[0][l] = rsArray[0][l];
                 nRSArray[0][l] = alphaS_Server
                                       -> FlavourSwitch(rsArray[0][l]);
                 nFSArray[0][l] = alphaS_Server
                                       -> FlavourSwitch(fsArray[0][l]);
             }
             ConvoluteWithDistributions_User(
                process, ca,   
                rsArray[0], fsArray[0], 
                noutArray, npdenArray,
                nRSArray[0], nFSArray[0],
                pdf_collection,
                pdf_parametrization,
                pdf_set,
                variant_alpha_s,           
                alphaSOrder,            
                alphasLambdaQCD4, 
                variant_alpha_em,
                QArray
             ); 
             for (l = ca -> GetActual(); l >= 0; --l) {
                 book -> SPS_Graph(l, i, tsum[l]);
                 kT_Q_scale -> SPS_Graph(l, i, tsum[l]);
             }
         }
   
         // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
         if (rsfsscale) {
   
            book = JADE_Q_rs_scale;
            for (i = 0; i < book -> GetSize(); ++i) {
                for (l = ca -> GetActual(); l >= 0; --l) {
                    rsArray[0][l] = QArray[l] * book -> GetX(i);
                    fsArray[0][l] = QArray[l];
                    nRSArray[0][l] = alphaS_Server
                                          -> FlavourSwitch(rsArray[0][l]);
                    nFSArray[0][l] = alphaS_Server
                                          -> FlavourSwitch(fsArray[0][l]);
                }
                ConvoluteWithDistributions_User(
                   process, ca,   
                   rsArray[0], fsArray[0], 
                   noutArray, npdenArray,
                   nRSArray[0], nFSArray[0],
                   pdf_collection,
                   pdf_parametrization,
                   pdf_set,
                   variant_alpha_s,           
                   alphaSOrder,            
                   alphasLambdaQCD4,
                   variant_alpha_em,
                   QArray
                ); 
                for (l = ca -> GetActual(); l >= 0; --l) {
                    book -> SPS_Graph(l, i, tsum[l]);
                    kT_Q_rs_scale -> SPS_Graph(l, i, tsum[l]);
                }
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            book = JADE_Q_fs_scale;
            for (i = 0; i < book -> GetSize(); ++i) {
                for (l = ca -> GetActual(); l >= 0; --l) {
                    rsArray[0][l] = QArray[l];
                    fsArray[0][l] = QArray[l] * book -> GetX(i);
                    nRSArray[0][l] = alphaS_Server
                                        -> FlavourSwitch(rsArray[0][l]);
                    nFSArray[0][l] = alphaS_Server
                                        -> FlavourSwitch(fsArray[0][l]);
                }
                ConvoluteWithDistributions_User(
                   process, ca,   
                   rsArray[0], fsArray[0], 
                   noutArray, npdenArray,
                   nRSArray[0], nFSArray[0],
                   pdf_collection,
                   pdf_parametrization,
                   pdf_set,
                   variant_alpha_s,           
                   alphaSOrder,            
                   alphasLambdaQCD4,
                   variant_alpha_em,
                   QArray
                ); 
                for (l = ca -> GetActual(); l >= 0; --l) {
                    book -> SPS_Graph(l, i, tsum[l]);
                    kT_Q_fs_scale -> SPS_Graph(l, i, tsum[l]);
                }
            }
   
         }
   
      }
   
   } // of ifFinalRun

   // -------------------------------------------------------------------------

   if (universalObservableChoice >= 0) {

      ret = CalculateAdaptationObservable_User(
               process,
               ca,
               pdf_collection,
               pdf_parametrization,
               pdf_set,
               variant_alpha_s,           
               alphaSOrder,            
               alphasLambdaQCD4,
               variant_alpha_em
            );
   }

   // -------------------------------------------------------------------------

   if (   ifFinalRun 
       && library != NULL   
      )
      library -> AddUpEvent();
   
   return ret;
}

// ----------------------------------------------------------------------------

void MyUserClass :: DropEvent_User(Process* process) {

   // :_TAWM_:
   process = process;

//   global -> To_status() << "MyUserClass :: DropEvent_User: event dropped\n";
}

// ----------------------------------------------------------------------------

void MyUserClass :: AcceptEvent_User(Process* process) {

   // :_TAWM_:
   process = process;

//   global -> To_status()
//      << "MyUserClass :: AcceptEvent_User: event accepted\n";
}

// ----------------------------------------------------------------------------

void MyUserClass :: CalculateObservablesForEvent_User(
                       Process* process, 
                       Event* event,
                       ContributionArray* ca, 
                       int ifFinalRun
                    ) {

   // :_TAWM_:
   process = process;
   ca = ca;

   Event* clu = NULL;  
   Event* cop = NULL;

   hCMS_event  = NULL;
   lab_event   = NULL;

   // boost to other frames

   HorizontalList <Event> * firstInH = event -> GetH();

      if (      ifFinalRun
             && (
                    testObservables
                 || observablesA
                 || observablesB
                 || observablesC
                 || observablesE
                )
          || universalObservableChoice < 0
         ) {

            // event record in the hCMS
            hCMS_event 
                   = global
                     -> GetEFree()
                     -> GetDefinedObject(
                           global -> disaster, 
                           event -> GetMaxNumber()
                        );
            event -> CopyInto(hCMS_event);
            hCMS_event -> SetUserFrame(hCMS);
            firstInH -> InsertAfterFirstElement(hCMS_event -> GetH());  
            hCMS_event -> BoostTo_hCMS();
   
            // event record in the Laboratory frame
            lab_event
                   = global
                     -> GetEFree()
                     -> GetDefinedObject(
                           global -> disaster, 
                           event -> GetMaxNumber()
                        );
            event -> CopyInto(lab_event);
            firstInH -> InsertAfterFirstElement(lab_event -> GetH());

            lab_event -> BoostTo_DIS_Lab(ratioProtonElectronInLaboratory);

            clu 
                  = global
                    -> disaster
                    -> GetEFree()
                    -> GetDefinedObject(
                          global -> disaster,  
                          event -> GetMaxNumber()
                       );
      
            cop
                  = global
                    -> disaster
                    -> GetEFree()
                    -> GetDefinedObject(
                          global -> disaster, 
                          event -> GetMaxNumber()
                       );
      }
   
   // -------------------------------------------------------------------------

   Contribution* c;
   int iev;

   // -------------------------------------------------------------------------

      int startCluster;
      int finalCluster;
     
      double cc;
      double clusc2;

      if (      ifFinalRun 
             && (
                    observablesA 
                 || observablesB
                )
          || universalObservableChoice == -1
         ) {
      
         clusc2 = lab_event -> W2;
         cc = 0.02;
         int reduct;
         int etaCut;
         int fcs;
         reduct = lab_event -> ClusterInvariantOld(LORENTZ_INVARIANT_OLD_1,
                                                   sqrt(cc*clusc2), 
                                                   startCluster,
                                                   etaCut, 0);
      
         finalCluster = startCluster - reduct;
         fcs = finalCluster;
      
         if (inRange(finalCluster, minJets, maxJets) && etaCut == 1) {
      
               event -> PrestoreGraph(t_cut_old_total, 1.0, 0);
               event -> PrestoreGraph(t_cut_old_q,     1.0, 0);
               event -> PrestoreGraph(t_cut_old_g,     1.0, 0);
         }
      
         // -------------------------------------------------------------------
      
         clusc2 = lab_event -> W2;
         cc = 0.02;
         lab_event -> ClusterGeneralWithEvent(
                                       clu, 
                                       JADE_algorithm,
                                       E_Scheme,     
                                       cc * clusc2,
                                       startCluster, 
                                       finalCluster);
      
         // -------------------------------------------------------------------

#if 0      
         if (   finalCluster != fcs
             ||    clu -> Cut_Jets_pT_PseudoRapidity(1.0, 300.0, -3.5, 3.5) 
                != etaCut) {
            global -> To_status()
               << FS("other =%d ",
                     clu -> Cut_Jets_pT_PseudoRapidity(1.0, 300.0, -3.5, 3.5))
               << FS("jets=%d ", finalCluster)
               << FS("diff=%d\n\n", 
                      finalCluster - fcs
                     +10*( clu 
                           -> Cut_Jets_pT_PseudoRapidity(1.0, 300.0, -3.5, 3.5)
                          -etaCut));
            
            lab_event -> ClusterInvariantOld(LORENTZ_INVARIANT_OLD_1,
                                             sqrt(cc*clusc2),
                                             startCluster,   
                                             etaCut, 1);
            register int ii;
            for (ii = 0; ii < clu -> npartons; ++ii) {
                global -> To_status()
                   << FS("new: %d ", ii)
                   << FS("%16.6e ", clu -> GetMapping(ii) -> Calculate_pT())
                   << FS("%16.6e\n", 
                         clu -> GetMapping(ii) -> Calculate_PseudoRapidity());
            }
            clu -> Print(stdout);
         }
#endif
       
         if (inRange(finalCluster, minJets, maxJets)) {
      
            if ( clu -> Cut_Jets_pT_PseudoRapidity(1.0, 300.0, -3.5, 3.5) ) {
               event -> PrestoreGraph(t_cut_total, 1.0, 0);
               event -> PrestoreGraph(t_cut_q,     1.0, 0);
               event -> PrestoreGraph(t_cut_g,     1.0, 0);
            }
         }
      
         // -------------------------------------------------------------------

         clusc2 = lab_event -> W2;
         cc = 0.02;
         lab_event -> ClusterGeneralWithEvent(
                                       clu,
                                       JADE_algorithm,
                                       E_Scheme,
                                       cc * clusc2,
                                       startCluster,
                                       finalCluster); 

         int ok1 =    inRange(finalCluster, 3, 3)                
                   && clu -> Cut_Jets_pT_PseudoRapidity(1.0, 300.0, -3.5, 3.5);

         clusc2 = lab_event -> W2;
         cc = 0.01;
         lab_event -> ClusterGeneralWithEvent(
                                       clu,
                                       JADE_algorithm,
                                       E_Scheme,
                                       cc * clusc2,
                                       startCluster,
                                       finalCluster); 

         int ok2 =    inRange(finalCluster, 4, 4)                
                   && clu -> Cut_Jets_pT_PseudoRapidity(1.0, 300.0, -3.5, 3.5);

         if (ok1 && ok2) {
//         if (ok2) {
            event -> PrestoreGraph(three_q, 1.0, 0);      
            event -> PrestoreGraph(three_g, 1.0, 0);      
         }

         event -> PrestoreGraph(direct_sum, 1.0, 0);      

      }  // of ifFinalRun...

      if (      ifFinalRun
             && observablesE
          || universalObservableChoice == -2
         ) {
      
         clusc2 = breit_event -> Q2;
            
         breit_event -> ClusterGeneralWithEvent(
                           clu,
                           kT_algorithm,
                           E_Scheme,
                           1.0 * clusc2,
                           startCluster,
                           finalCluster
                        );

         if (inRange(finalCluster, minJets, maxJets)) {

            event -> PrestoreGraph(kt_cmp, 1.0, 0);      
            clu -> CopyInto(cop);
            cop -> BoostTo_DIS_Lab(ratioProtonElectronInLaboratory); 
            if ( cop -> Cut_Jets_pT_PseudoRapidity(5.0, 300.0, -3.5, 3.5) ) {
               event -> PrestoreGraph(kt_cut_cmp, 1.0, 0);      
            }
         }
      
      }  // of ifFinalRun...


   if (ifFinalRun) {

         // find pointers to the various frames
         // :_MOD_: defined above; not required
         //hCMS_event = event -> UserFrameFindEvent(hCMS);
         //lab_event  = event -> UserFrameFindEvent(lab);
   
      register int i;
   
      double tmom;
      double thrust;
   
      Book* book;
   
      Array_double scale2(10);
      Array_int    nCluster(10);
      int nEntries;

      if (observablesA) {

         // -------------------------------------------------------------------
      
         if (obsSet == 1) {
            book = total;
            event -> PrestoreGraph(book, 1.0, 0);
         }
      
         // -------------------------------------------------------------------
      
         book = JADE_E_hCMS_cut;
         clusc2 = hCMS_event -> W2;
      
         hCMS_event -> ClusterGeneralHistoryWithEvent(
                          clu,
                          JADE_algorithm,
                          E_Scheme,     
                          &scale2, 
                          &nCluster, 
                          nEntries  
                       );
   
         for (i = 0; i < book -> GetSize(); ++i) {
             cc = book -> GetX(i);
             if (inRange(clusterNumber(
                            &scale2, 
                            &nCluster, 
                            nEntries, 
                            cc*clusc2
                         ), 
                 minJets, maxJets)) 
                event -> PrestoreGraph(book, 1.0, i);
         }

      }  // of observablesA
   
      // ----------------------------------------------------------------------
      
      if (testObservables) {

         book = test2History;
         clusc2 = hCMS_event -> W2;
      
         hCMS_event -> ClusterGeneralHistoryWithEvent(
                          clu,
                          JADE_algorithm,
                          E_Scheme,     
                          &scale2, 
                          &nCluster, 
                          nEntries  
                       );
      
         for (i = 0; i < book -> GetSize(); ++i) {
             cc = book -> GetX(i);
             if (inRange(clusterNumber(
                            &scale2, &nCluster, nEntries, cc*clusc2
                         ), 
                 2, 2)) 
                event -> PrestoreGraph(test2History, 1.0, i);
             if (inRange(clusterNumber(
                            &scale2, &nCluster, nEntries, cc*clusc2
                         ), 
                 3, 3)) 
                event -> PrestoreGraph(test3History, 1.0, i);
             if (inRange(clusterNumber(
                            &scale2, &nCluster, nEntries, cc*clusc2
                         ), 
                 4, 4)) 
                event -> PrestoreGraph(test4History, 1.0, i);
         }
      
         // -------------------------------------------------------------------
      
         // this has problems with back-to-back jets!
         book = test2Other;
         clusc2 = hCMS_event -> W2;
         for (i = 0; i < book -> GetSize(); ++i) {
             cc = book -> GetX(i);
             hCMS_event -> ClusterGeneralWithEvent(
                                           clu, 
                                           JADE_algorithm,
                                           E_Scheme,     
                                           cc * clusc2,
                                           startCluster, 
                                           finalCluster);
             if (inRange(finalCluster, 2, 2)) {
                event -> PrestoreGraph(test2Other, 1.0, i);
             }
             if (inRange(finalCluster, 3, 3)) {
                event -> PrestoreGraph(test3Other, 1.0, i);
             }
             if (inRange(finalCluster, 4, 4)) {
                event -> PrestoreGraph(test4Other, 1.0, i);
             }
         }
      
         // -------------------------------------------------------------------
      
         book = test2Old;
         clusc2 = lab_event -> W2;
   
         for (i = 0; i < book -> GetSize(); ++i) {
             cc = book -> GetX(i);
             int reductX;
             int etaCutX;
             int fcsX;
             reductX = lab_event -> ClusterInvariantOld(LORENTZ_INVARIANT_OLD_1,
                                                       sqrt(cc*clusc2), 
                                                       startCluster,
                                                       etaCutX, 0);
      
             finalCluster = startCluster - reductX;
             fcsX = finalCluster;
      
             if (inRange(finalCluster, 2, 2) && etaCutX == 1) {
                event -> PrestoreGraph(test2Old, 1.0, i);
             }
             if (inRange(finalCluster, 3, 3) && etaCutX == 1) {
                event -> PrestoreGraph(test3Old, 1.0, i);
             }
             if (inRange(finalCluster, 4, 4) && etaCutX == 1) {
                event -> PrestoreGraph(test4Old, 1.0, i);
             }
         }
      
         // -------------------------------------------------------------------

      }

      if (observablesA) {

#if 0
         // this has problems with back-to-back jets!
         book = JADE_E0_hCMS_cut;
         clusc2 = hCMS_event -> W2;
         for (i = 0; i < book -> GetSize(); ++i) {
             cc = book -> GetX(i);
             hCMS_event -> ClusterGeneralWithEvent(
                                           clu, 
                                           JADE_algorithm,
                                           E0_Scheme,     
                                           cc * clusc2,
                                           startCluster, 
                                           finalCluster);
             if (finalCluster==3) 
                event -> PrestoreGraph(book, 1.0, i);
         }
#endif
      
         // -------------------------------------------------------------------
      
         book = JADE_P_hCMS_cut;
         clusc2 = hCMS_event -> W2;
      
         hCMS_event -> ClusterGeneralHistoryWithEvent(
                          clu,
                          JADE_algorithm,
                          P_Scheme,     
                          &scale2, 
                          &nCluster, 
                          nEntries  
                       );
      
         for (i = 0; i < book -> GetSize(); ++i) {
             cc = book -> GetX(i);
             if (inRange(clusterNumber(
                            &scale2, &nCluster, nEntries, cc*clusc2
                         ), 
                 minJets, maxJets)) 
                event -> PrestoreGraph(book, 1.0, i);
         }
      
         // -------------------------------------------------------------------
      
         book = JADE_JADE_hCMS_cut;
         clusc2 = hCMS_event -> W2;
      
         hCMS_event -> ClusterGeneralHistoryWithEvent(
                          clu,
                          JADE_algorithm,
                          JADE_Scheme,     
                          &scale2, 
                          &nCluster, 
                          nEntries  
                       );
      
         for (i = 0; i < book -> GetSize(); ++i) {
             cc = book -> GetX(i);
             if (inRange(clusterNumber(
                            &scale2, &nCluster, nEntries, cc*clusc2
                         ), 
                 minJets, maxJets)) 
                event -> PrestoreGraph(book, 1.0, i);
         }
      
         // -------------------------------------------------------------------
      
         book = kT_E_Breit_cut;
         clusc2 = breit_event -> Q2;
      
         breit_event -> ClusterGeneralHistoryWithEvent(
                           clu,
                           kT_algorithm,
                           E_Scheme,     
                           &scale2, 
                           &nCluster, 
                           nEntries  
                        );
      
         for (i = 0; i < book -> GetSize(); ++i) {
             cc = book -> GetX(i);
             if (inRange(clusterNumber(
                            &scale2, &nCluster, nEntries, cc*clusc2
                         ), 
                 minJets, maxJets)) 
                event -> PrestoreGraph(book, 1.0, i);
         }
      
         // -------------------------------------------------------------------
      
         clusc2 = hCMS_event -> W2;
         cc = 0.02;
#if 0
         hCMS_event -> ClusterGeneralWithEvent(
                                       clu, 
                                       JADE_algorithm,
                                       JADE_Scheme,     
                                       cc * clusc2,
                                       startCluster, 
                                       finalCluster);
#else
         hCMS_event -> ClusterGeneralWithEvent(
                                       clu, 
                                       JADE_algorithm,
                                       P_Scheme,     
                                       cc * clusc2,
                                       startCluster, 
                                       finalCluster);
#endif
      
         // -------------------------------------------------------------------
      
         if (inRange(finalCluster, minScaleJets, maxScaleJets)) {
      
            book = JADE_Q_scale;
            for (i = 0; i < book -> GetSize(); ++i)
                event -> PrestoreGraph(book, 1.0, i);
            if (rsfsscale) {
               book = JADE_Q_rs_scale;
               for (i = 0; i < book -> GetSize(); ++i)
                   event -> PrestoreGraph(book, 1.0, i);
               book = JADE_Q_fs_scale;
               for (i = 0; i < book -> GetSize(); ++i)
                   event -> PrestoreGraph(book, 1.0, i);
            }
         }
      
         // boost the clustered event to the Breit frame
         clu -> CopyInto(cop);
         cop -> BoostTo_Breit();
      
         double taverage = 0;
         for (i=0; i<finalCluster-1; ++i)
             taverage += cop -> mapping[i] -> GetEnergy();
         taverage /= (finalCluster-1);
         for (c = event -> GetContribution();
              c != NULL && (iev = c -> GetLocation(), TRUE);
              c = c -> GetSameEvent()) {
             EBreitArray[iev] = taverage;
             c -> Set_EBreit_U(taverage);
         }
      
         // -------------------------------------------------------------------
      
         if (inRange(finalCluster, minScaleJets, maxScaleJets)) {
            book = JADE_EBreit_scale;
            for (i = 0; i < book -> GetSize(); ++i)
                event -> PrestoreGraph(book, 1.0, i);
         }
      
         // -------------------------------------------------------------------
      
         if (inRange(finalCluster, minPtJets, maxPtJets)) {
      
            book = pT_hCMS;
            for (i=0; i<finalCluster-1; ++i) {
               tmom = clu -> mapping[i] -> Get_pT();
               event -> PrestoreHistogram(book, tmom);
            }
         }
      
         // -------------------------------------------------------------------
      
         if (inRange(finalCluster, minScaleJets, maxScaleJets)) {
            book = energy_Breit;
            for (i=0; i<finalCluster-1; ++i) {
                tmom = cop -> mapping[i] -> GetEnergy();
                event -> PrestoreHistogram(book, tmom);
            }
         }
      
         // -------------------------------------------------------------------
      
         clusc2 = breit_event -> Q2;
         cc = 0.5;
         breit_event -> ClusterGeneralWithEvent(
                                    clu, 
                                    kT_algorithm,
                                    E_Scheme,     
                                    cc * clusc2,
                                    startCluster, 
                                    finalCluster);
      
         // -------------------------------------------------------------------
      
         if (inRange(finalCluster, minScaleJets, maxScaleJets)) {
            book = kT_Q_scale;
            for (i = 0; i < book -> GetSize(); ++i)
                event -> PrestoreGraph(book, 1.0, i);
            if (rsfsscale) {
               book = kT_Q_rs_scale;
               for (i = 0; i < book -> GetSize(); ++i)
                   event -> PrestoreGraph(book, 1.0, i);
               book = kT_Q_fs_scale;
               for (i = 0; i < book -> GetSize(); ++i)
                   event -> PrestoreGraph(book, 1.0, i);
            }
         }
      
         // -------------------------------------------------------------------
      
         if (inRange(finalCluster, minPtJets, maxPtJets)) {
      
            book = kT_Breit;
            for (i=0; i<finalCluster-1; ++i) {
                tmom = clu -> mapping[i] -> Get_pT();
                event -> PrestoreHistogram(book, tmom);
            }
         }
      
         // -------------------------------------------------------------------
      
            book = Thrust_Breit;
            thrust = clu -> CurrentThrust();
            event -> PrestoreHistogram(book, thrust);
   
         // -------------------------------------------------------------------
      
      }  // of observablesA

      if (   observablesA
          || observablesC) {

         clusc2 = hCMS_event -> W2;
         cc = 0.02;
         hCMS_event -> ClusterGeneralWithEvent(
                                       clu, 
                                       JADE_algorithm,
                                       E_Scheme,     
                                       cc * clusc2,
                                       startCluster, 
                                       finalCluster);
      
         if (inRange(finalCluster, minJets, maxJets)) {
      
            event -> PrestoreGraph(t_total, 1.0, 0);
            event -> PrestoreGraph(t_q,     1.0, 0);
            event -> PrestoreGraph(t_g,     1.0, 0);
         }

      }  // of observablesA || observablesC

   } // of if (ifFinalRun)

   if (      ifFinalRun
          && (
                 testObservables
              || observablesA
              || observablesB
              || observablesC
              || observablesE
             )
       || universalObservableChoice < 0
      ) {
      global -> disaster -> GetEFree() -> ReturnObject(clu);
      global -> disaster -> GetEFree() -> ReturnObject(cop);
   }
}
   
// ----------------------------------------------------------------------------

// may be used to reject events based on the lepton phase space

void MyUserClass :: IncidentPhaseSpace_User(
                       Process* process,
                       Event* event,
                       int& flag
                    ) {

   // :_TAWM_:
   process = process;
   event = event;

   flag = TRUE;
}

// ----------------------------------------------------------------------------
//
// --> FORTRAN interface class
//
// ----------------------------------------------------------------------------

FORTRAN_User :: FORTRAN_User() {

   gl -> To_status() << "FORTRAN_User :: FORTRAN_User()\n";

   mcUserInterface = new MCUserInterface_FORTRAN;
   mcUserInterface -> Set_mcUserCalled(this);
   mcUserInterface -> Set_mcUserForDisaster(this);
}

// ----------------------------------------------------------------------------
   
FORTRAN_User :: ~FORTRAN_User() {

   gl -> To_status() << "FORTRAN_User :: ~FORTRAN_User()\n";
   Delete();

   delete mcUserInterface;
}

// ----------------------------------------------------------------------------
   
void FORTRAN_User :: Define(Disaster* globalIn) {

   if (definedObject.IsDefined())
      errf(-1, "FORTRAN_User :: Define: already defined");

   gl -> To_status() << "FORTRAN_User :: Define()\n";

   MCUserCalled :: Define(globalIn);
}

// ----------------------------------------------------------------------------
   
void FORTRAN_User :: Delete() {

   if (definedObject.Status()) {

      definedObject.ResetDefined();   
      gl -> To_status() << "FORTRAN_User :: Delete()\n";
 
      MCUserCalled :: Delete();
   }
}

// ----------------------------------------------------------------------------
   
void FORTRAN_User :: Start_User() {

   user1_C(1);
}

// ----------------------------------------------------------------------------
   
void FORTRAN_User :: End_User() {

   user1_C(2);
}

// ----------------------------------------------------------------------------

// called if a parameter name parsed by the MC cannot be resolved
// --> assumed to be a user parameter

int FORTRAN_User :: SetParameter_User(char* pname, char* parameter) {

   // :_TAWM_:
   pname = pname;
   parameter = parameter;

   int ret = -1;
   return ret;
}

// ----------------------------------------------------------------------------

void FORTRAN_User :: StartIntegration1_User() {

   user1_C(3);
}

// ----------------------------------------------------------------------------

void FORTRAN_User :: StartIntegration2_User() {

   user1_C(4);
}

// ----------------------------------------------------------------------------

void FORTRAN_User :: EndIntegration_User(double average, double error) {

   user3_C(average, error);
   user1_C(5);
}

// ----------------------------------------------------------------------------

void FORTRAN_User :: BeginEvent_User(Process* process) {

   // :_TAWM_:
   process = process;

   user1_C(6);
}

// ----------------------------------------------------------------------------

void FORTRAN_User :: EndOfEvent_User(Process* process) {

   // :_TAWM_:
   process = process;

   user1_C(9);
}

// ----------------------------------------------------------------------------

double FORTRAN_User :: EndEvent_User(
                          Process* process,
                          ContributionArray* ca,
                          int ifFinalRun
                      ) {

   // :_TAWM_:
   process = process;
   ifFinalRun = ifFinalRun;

   caSave = ca;
   processSave = process;

   // do this to determine the event record in the Breit frame
   // calculate observables
                    
   UseListContainer <Event> * ec;
   Event* e;
   // run through the list of all events
   for (ec = process -> GetEventList() -> GetPtrToFirstContainer();
        ec != NULL && (e = ec -> GetDataPtr(), TRUE);
        ec = ec -> GetNext()) {
       if (e -> GetContribution() != NULL)
          CalculateObservablesForEvent_User(process, e, ca, ifFinalRun);
   }        
   
   const int cIndex[] = {-8, -7, -5, -4, -1, 0, 1, 2};
   int cSizeStart = 5;

   int nr;
   int nl;
   double fvect[4*21*30];
   int npartons[30];
   double xb[30];
   double q2[30];
   double xi[30];
   double weight[1000];
   int irps[50];
   int ialphas[50];
   int ialphaem[50];
   int lognf[50];

   // copy into new variables

   // parameters of FORTRAN arrays
   const int weightSize    =  11;
   const int particleSize  =  21;
   const int particleStart = -10;

//   UseListContainer <Event> * ec;
//   Event* e;

   Contribution* c;
   int iev;

//   ca -> Print(stdout, 1);

   nl = ca -> GetActual() + 1;

   int eventCo = 0;

   // run through the list of all events 
   for (ec = process
             -> GetEventList()
             -> GetPtrToFirstContainer();
        ec != NULL && (e = ec -> GetDataPtr(), TRUE);
        ec = ec -> GetNext()) {

       if (e -> GetContribution() != NULL) {

          // list of contributions associated with the event
          for (c = e -> GetContribution();
               c != NULL && (iev = c -> GetLocation(), TRUE);
               c = c -> GetSameEvent()) {
              irps[iev] = eventCo + 1;  // FORTRAN offset (array starts at 1)
              ialphas[iev]  = c -> GetOrderAlphaS();
              ialphaem[iev] = c -> Get_orderAlphaEM();
              lognf[iev]    = c -> GetSLF();
              double* base = weight + weightSize * iev;
              double* data = c -> GetData();
              for (register int k = weightSize - 1; k >= 0; --k) {
                  *(base++) = *(data++);
              }
          }

          Event* breitE = e -> UserFrameFindEvent(breit);

          npartons[eventCo] = breitE -> GetNpartons();
          xb[eventCo]       = breitE -> Get_xB();
          q2[eventCo]       = breitE -> Get_Q2();
          xi[eventCo]       = breitE -> Getxi();

          double* base = fvect + 4 * (particleSize * eventCo - particleStart);

          // loop over particles
          for (register int ip = cSizeStart + npartons[eventCo] - 1; 
               ip >= 0; 
               -- ip
              ) {
              breitE 
              -> GetMapping(cIndex[ip]) 
              -> CalculateCartesianFromPolar();
              double* newBase = base + 4 * cIndex[ip];
              double* cartesian = breitE 
                                  -> GetMapping(cIndex[ip]) 
                                  -> GetFvect();
              *(newBase++) = *(cartesian++);
              *(newBase++) = *(cartesian++);
              *(newBase++) = *(cartesian++);
              *(newBase++) = *(cartesian++);
          }
          ++ eventCo;
       }
   }

   nr = eventCo;

   double retp
      = user2_C(  
           nr,
           nl,
           fvect,
           npartons,
           xb,  
           q2,  
           xi,
           weight,
           irps,   
           ialphas,
           ialphaem,
           lognf
     );

   return retp;
}

// ----------------------------------------------------------------------------

void FORTRAN_User :: DropEvent_User(Process* process) {

   // :_TAWM_:
   process = process;

   user1_C(7);
}

// ----------------------------------------------------------------------------

void FORTRAN_User :: AcceptEvent_User(Process* process) {

   // :_TAWM_:
   process = process;  

   user1_C(8);
}

// ----------------------------------------------------------------------------

void FORTRAN_User :: CalculateObservablesForEvent_User(
                       Process* process, 
                       Event* event,
                       ContributionArray* ca, 
                       int ifFinalRun
                    ) {

   // :_TAWM_:
   process = process;
   event = event;
   ca = ca;
   ifFinalRun = ifFinalRun;
}
   
// ----------------------------------------------------------------------------

// may be used to reject events based on the lepton phase space

void FORTRAN_User :: IncidentPhaseSpace_User(
                        Process* process,
                        Event* event,
                        int& flag
                     ) {

   // :_TAWM_:
   process = process;
   event = event;

   flag = TRUE;
}

// ----------------------------------------------------------------------------
// 
// --> individual routines called from FORTRAN
//
// ----------------------------------------------------------------------------

// global variables (...)

   Disaster* disaster_F;
   FORTRAN_User* user_F;
 
// ----------------------------------------------------------------------------

extern "C" {
void disaster_c_from_f_start() {

   disaster_F = new Disaster;  // this also defines ``gl''
   disaster_F -> disaster = disaster_F;

   // output via a FORTRAN routine
   (gl -> To_status()).Associate(print_f_C,  6);
   (gl ->  To_error()).Associate(print_f_C,  6);
   (gl ->    To_log()).Associate(print_f_C, 42);

   disaster_F -> Define();

   gl -> To_status() << "\nDISASTER++ FORTRAN INTERFACE START\n\n";
   
   user_F = new FORTRAN_User;
   user_F -> Define(disaster_F);
   disaster_F -> SetUser(user_F);
   disaster_F -> Start();

   // set a few default parameters
   disaster_F -> PPCommand("SCALE_CHOICE_FOR_TECHNICAL_CUT", "0");
   disaster_F -> PPCommand("TECHNICAL_CUT_PARAMETER", "1.e-10");
   disaster_F -> PPCommand("SUBTRACTION_TYPE", "1");
   disaster_F -> PPCommand("E_GAMMA", "1.0");
   disaster_F -> PPCommand("E_GAMMA1", "1.0");
   disaster_F -> PPCommand("V_DELTA", "2.0");
   disaster_F -> PPCommand("V_DELTA1", "2.0");
   disaster_F -> PPCommand("MC_EPS", "1.e-12");
   disaster_F -> PPCommand("ERROR_DUMP", "0");
   disaster_F -> PPCommand("RATIO_PROTON_ELECTRON_IN_LABORATORY", "29.8845");
   disaster_F -> PPCommand("XI_INTEGRATION", "1");
   disaster_F -> PPCommand("EXPLICIT_FRAME_CHOICE", "-1");
   disaster_F -> PPCommand("PS_PERMUTE_N", "1");
   disaster_F -> PPCommand("PS_PERMUTE_N", "2");
   disaster_F -> PPCommand("PS_PERMUTE_N", "3");
   disaster_F -> PPCommand("PS_PERMUTE_N", "4");
   disaster_F -> PPCommand("PDF_COLLECTION", "1");
   disaster_F -> PPCommand("PDF_PARAMETRIZATION", "5");
   disaster_F -> PPCommand("PDF_SET", "6");
   disaster_F -> PPCommand("USER_DEFINED_ALPHA_S", "1");
   disaster_F -> PPCommand("ALPHA_S_ORDER", "2");
   disaster_F -> PPCommand("ALPHA_S_VARIANT", "1");
   disaster_F -> PPCommand("ALPHA_EM_VARIANT", "1");
   disaster_F -> PPCommand("N_OUT_FLAVOURS", "5");
   disaster_F -> PPCommand("N_PDEN_FLAVOURS", "5");
   disaster_F -> PPCommand("TEST_OBSERVABLES", "0");
   disaster_F -> PPCommand("OBSERVABLES_A", "0");
   disaster_F -> PPCommand("OBSERVABLES_B", "1");
   disaster_F -> PPCommand("OBSERVABLES_C", "1");
   disaster_F -> PPCommand("OBSERVABLES_D", "1");
   disaster_F -> PPCommand("OBSERVABLES_E", "1");
   disaster_F -> PPCommand("UNIVERSAL_OBSERVABLE_CHOICE", "0");
   disaster_F -> PPCommand("PROCESS_INDEX", "2");
   disaster_F -> PPCommand("NUMBER_OF_FINAL_STATE_PARTONS_IN_BORN_TERM", "2");
   disaster_F -> PPCommand("POINTS", "100");
   disaster_F -> PPCommand("ITERATIONS", "10");
   disaster_F -> PPCommand("FACT_PREP", "0.5");
   disaster_F -> PPCommand("FACT_FINAL", "1.5");
   disaster_F -> PPCommand("LEPTON_INTEGRATION", "0");
   disaster_F -> PPCommand("XB_FIXED", "0.01");
   disaster_F -> PPCommand("Y_FIXED", "0.444444");
   disaster_F -> PPCommand("ECM", "300.");

#if 1
   disaster_F -> PPCommand("UNIVERSAL_OBSERVABLE_ENERGY_POWER","1.0");
   disaster_F -> PPCommand("UNIVERSAL_OBSERVABLE_ANGLE_POWER","0.5");
#endif

#if 0
   disaster_F -> PPCommand("SHIFTED_PARTON_VARIABLES", "1");
   disaster_F -> PPCommand("E_GAMMA", "1.0");
   disaster_F -> PPCommand("E_GAMMA_UPPER", "2.0");
   disaster_F -> PPCommand("E_LAMBDA", "0.5");
   disaster_F -> PPCommand("E_MU", "0.5");
   disaster_F -> PPCommand("E_GAMMA1", "1.0");
   disaster_F -> PPCommand("V_DELTA", "2.0");
   disaster_F -> PPCommand("V_DELTA_UPPER", "1.0");
   disaster_F -> PPCommand("V_DELTA_LAMBDA", "0.5");
   disaster_F -> PPCommand("V_DELTA_MU", "0.5");
   disaster_F -> PPCommand("V_DELTA1", "2.5");
   disaster_F -> PPCommand("V_DELTA1_UPPER", "2.5");
   disaster_F -> PPCommand("V_DELTA1_LAMBDA", "0.5");
   disaster_F -> PPCommand("V_DELTA1_MU", "0.5");
   disaster_F -> PPCommand("U_ALPHA", "2.0");
   disaster_F -> PPCommand("XI_ALPHA","0.2");
#endif

//   disaster_F -> PPCommand("");
}
}

// ----------------------------------------------------------------------------

extern "C" {
void disaster_c_from_f_stop() {

   disaster_F -> Terminate();
   delete user_F;
   delete disaster_F;

   Mstream s_out;
   s_out.Associate(print_f_C, 6);   
   s_out << "\nDISASTER++ FORTRAN INTERFACE STOP\n\n";
}
}

// ----------------------------------------------------------------------------

// pName is void* because cfortran.h requires this

extern "C" {
void disaster_c_from_f_setInt(void* pName, int in) {

   int err;
   err = disaster_F -> PPCommandInt((char*) pName, in);

   if (err == -1) {
      gl -> To_status()
               << FS(":%s:\n", (char*) pName);
      errf(-1, "disaster_c_from_f_setInt: parameter not known");  
   }
}
}

// ----------------------------------------------------------------------------

// pName is void* because cfortran.h requires this

extern "C" {
void disaster_c_from_f_setDouble(void* pName, double in) {

   int err;
   err = disaster_F -> PPCommandDouble((char*) pName, in);

   if (err == -1) {
      gl -> To_status()
               << FS(":%s:\n", (char*) pName);
      errf(-1, "disaster_c_from_f_setDouble: parameter not known");  
   }
}
}

// ----------------------------------------------------------------------------

// pName is void* because cfortran.h requires this

extern "C" {
void disaster_c_from_f_executeCommand(void* pName) {

   // set parameters if necessary
  // was trying some diagnostics, but did not get very far
  //cout << (char*) pName << "\n";
  //cout << "strcmp result: " << strcmp((char*) pName, "RUN_MC") << "\n";
  //cout << "strcmp result: " << strcmp((char*) pName, "RUN_M|") << "\n";
  //cout << "strcmp result: " << strcmp("RUN_MC", "RUN_MC") << "\n";
  //cout << "strcmp result: " << strcmp("ABD", "RUN_MC") << "\n";

   if (   ! strcmp((char*) pName, "RUN_MC")
       && (   disaster_F -> Get_process_index() == 1
           || disaster_F -> Get_process_index() == 2
          )
       && disaster_F -> Get_numberOfFinalStatePartonsInBornTerm() == 2
      ) {
      // LO and NLO 2+1 jet production
      gl -> To_status() << "modifying phase space parameters "
                        << "(2+1) jet production\n";
      disaster_F -> PPCommand("SHIFTED_PARTON_VARIABLES", "1");
      disaster_F -> PPCommand("E_GAMMA", "1.0");
      disaster_F -> PPCommand("E_GAMMA_UPPER", "2.0");
      disaster_F -> PPCommand("E_LAMBDA", "0.5");
      disaster_F -> PPCommand("E_MU", "0.5");
      disaster_F -> PPCommand("E_GAMMA1", "1.0");
      // GPS TEMP FIX -- not applied
      disaster_F -> PPCommand("V_DELTA", "2.0");
      disaster_F -> PPCommand("V_DELTA_UPPER", "1.0");
      disaster_F -> PPCommand("V_DELTA_LAMBDA", "0.5");
      disaster_F -> PPCommand("V_DELTA_MU", "0.5");
      // GPS TMP FIX -- not applied
      disaster_F -> PPCommand("V_DELTA1", "2.5");
      disaster_F -> PPCommand("V_DELTA1_UPPER", "2.5");
      disaster_F -> PPCommand("V_DELTA1_LAMBDA", "0.5");
      disaster_F -> PPCommand("V_DELTA1_MU", "0.5");
      disaster_F -> PPCommand("U_ALPHA", "2.0");
      disaster_F -> PPCommand("XI_ALPHA","0.2");
   }


   int err;
   err = disaster_F -> PPCommandInt((char*) pName, 0);

   if (err == -1) {
      gl -> To_status() << FS(":%s:\n", (char*) pName);
      errf(-1, "disaster_c_from_f_executeCommand: parameter not known");  
   }
}
}

// ----------------------------------------------------------------------------

extern "C" {
double disaster_c_from_f_adaptationObservable(
          int pdf_collectionIn,
          int pdf_parametrizationIn,
          int pdf_setIn,
          int alphaSVariant,
          int alphaSOrder,
          double alphasLambdaQCD4,
          int alphaEMVariant
       ) {

   return user_F
          -> CalculateAdaptationObservable_User(
                user_F -> Get_processSave(),
                user_F -> Get_caSave(),
                pdf_collectionIn,
                pdf_parametrizationIn,
                pdf_setIn,
                alphaSVariant,
                alphaSOrder,
                alphasLambdaQCD4,
                alphaEMVariant
             );
}
}

// ----------------------------------------------------------------------------

extern "C" {
void disaster_c_from_f_FlavourFactors(
        int pdf_collectionIn,
        int pdf_parametrizationIn,
        int pdf_setIn,
        int alphaSVariantIn,
        int alphaSOrderIn,
        double alphasLambdaQCD4In,
        int alphaEMVariantIn,
        int nf,
        double* din,
        double* dout
     ) {

   ((MCUserInterface_FORTRAN*) user_F -> Get_mcUserInterface()) 
      -> Control_To_C();

   double xi           = din[0];
   double muF          = din[1];
   double muR          = din[2];
   double alphaEMScale = din[3];

   PDFC* pdfc = user_F 
                -> Get_flavourFactors()
                -> FindOrCreateAndFillAndCacheFlavourFactors(
                      xi,
                      muF,
                      muR,
                      nf,
                      nf,
                      pdf_collectionIn,
                      pdf_parametrizationIn,
                      pdf_setIn,
                      alphaSVariantIn,
                      alphaSOrderIn,
                      alphasLambdaQCD4In,
                      alphaEMVariantIn,
                      alphaEMScale
                   );

   double* ff = pdfc -> GetFlavourFactorPtr();
   dout[ 0] = ff[ 0];
   dout[ 1] = ff[ 1];
   dout[ 2] = ff[ 2];
   dout[ 3] = ff[ 3];
   dout[ 4] = ff[ 4];
   dout[ 5] = ff[ 5];
   dout[ 6] = ff[ 6];
   dout[ 7] = ff[ 7];
   dout[ 8] = ff[ 8];
   dout[ 9] = ff[ 9];
   dout[10] = ff[10];
   dout[11] = *(pdfc -> GetAlphaSPtr());
   dout[12] = *(pdfc -> GetAlphaEMPtr());

   ((MCUserInterface_FORTRAN*) user_F -> Get_mcUserInterface()) 
      -> Control_To_FORTRAN();
}
}

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

#if 0
// ============================ GARBAGE =======================================
// ============================ END OF GARBAGE ================================
#endif
// ============================================================================
//
// --> global classes
//
// file:              global.cc
// created:           29.03.1997
// last modification: 15.12.1997
//
// ============================================================================

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "global.h"
#include "strng.h"
#include "utility.h"
#include "cmb.h"
#include "mth.h"
#include "nan.h"

// ============================================================================
//
// --> definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

Global *gl;

//Event    *gEvent    = NULL;
//Particle *gParticle = NULL;

FPNaN *fpnan;

// test purpose
long coTest;
long coCont;

// ----------------------------------------------------------------------------
//
// --> the error function:
//      0: non-fatal
//     -1: fatal
//
// ----------------------------------------------------------------------------

void errf(int level, const char* text) {
   gl -> To_error() << FS("\n\nDISASTER++ error:\n\n --> %s <--\n\n", text);
   if (level == -1) {
      gl -> To_error() << "Fatal error, program terminates.\n\n";
      provoke_segmentation_fault();
   }
}

// ----------------------------------------------------------------------------
//
// --> provoke a segmentation fault (to enable trace-back)
//
// ----------------------------------------------------------------------------

void provoke_segmentation_fault() {

   gl -> To_error() << "\nPROVOKED SEGMENTATION FAULT\n";
   int* dummy = NULL;
   int i = *dummy;
   i = i;
}

// ----------------------------------------------------------------------------
//
// --> provoke a floating point exception (to enable trace-back)
//
// ----------------------------------------------------------------------------

void provoke_floating_point_exception() {

   gl -> To_error() << "\n\nPROVOKED FLOATING POINT EXCEPTION\n";
   double a = 1.;
   double b = sqrt(0.);  
   double c = a / b;
   c = c;
   gl -> To_status() 
      << FS("provoke_floating_point_exception: number is %16.6e\n\n", c);
}

// ============================================================================
//
// --> FORTRAN interface routines
//
// ============================================================================

// ============================================================================
//
// --> classes
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> the global error class
//
// ----------------------------------------------------------------------------

GlobalError :: GlobalError() {
}

// ----------------------------------------------------------------------------

GlobalError :: ~GlobalError() {
}

// ----------------------------------------------------------------------------
//
// --> the global error class for fatal errors
//
// ----------------------------------------------------------------------------

FatalGlobalError :: FatalGlobalError(const String* messageIn) {

   message = NULL;

   CreateAndSetMessage(messageIn);
}

// ----------------------------------------------------------------------------

FatalGlobalError :: FatalGlobalError(const char* messageIn) {

   message = NULL;

   String copy(messageIn);

   // :_CHECK_: does String() survive sufficiently long?
//   CreateAndSetMessage(&String(messageIn));
   CreateAndSetMessage(&copy);
}

// ----------------------------------------------------------------------------

FatalGlobalError :: ~FatalGlobalError() {
}

// ----------------------------------------------------------------------------

void FatalGlobalError :: CreateAndSetMessage(const String* messageIn) {

   if (message != NULL)
      errf(-1, "FatalGlobalError :: CreateAndSetMessage: "
               "message not empty");

   message = new String(messageIn -> DetermineLength());
   message -> CopyFromString(messageIn);
}

// ----------------------------------------------------------------------------

String* FatalGlobalError :: GetMessageString() {

   return message;
}

// ----------------------------------------------------------------------------

char* FatalGlobalError :: GetMessageChar() {

   return message -> GetData();
}

// ----------------------------------------------------------------------------
//
// --> class for Formatted Strings 
//     (used for printing formatted data via stream io)
//
// ----------------------------------------------------------------------------

// default constructor

FS :: FS() { 

   CreateLargeString();
}

// ----------------------------------------------------------------------------

// destructor

FS :: ~FS() { 

   delete [] string; 
}

// ----------------------------------------------------------------------------

// copy constructor

FS :: FS(const FS& fs) { 

   CreateLargeString();
   strcpy(string, fs.string);
}

// ----------------------------------------------------------------------------

// assignment operator
// :_CAUTION_: string must be defined; check for this!

FS &FS :: operator=(const FS& fs) { 

   strcpy(string, fs.string);
   return *this;
}

// ----------------------------------------------------------------------------

// constructor for double

FS :: FS(char* format, const double d) { 

   CreateLargeString();

//   // first check for NaN
//   int cNaN = intnan8_C(d);   

//   if (cNaN)
//      strcpy(string, "***NAN***");
//   else
      sprintf(string, format, d);
}

// ----------------------------------------------------------------------------

// constructor for int

FS :: FS(char* format, const int d) { 

   CreateLargeString();
   sprintf(string, format, d);
}

// ----------------------------------------------------------------------------

// constructor for long

FS :: FS(char* format, const long d) { 

   CreateLargeString();
   sprintf(string, format, d);
}

// ----------------------------------------------------------------------------

// constructor for unsigned long

FS :: FS(char* format, const unsigned long d) { 

   CreateLargeString();
   sprintf(string, format, d);
}

// ----------------------------------------------------------------------------

// constructor for character array

FS :: FS(char* format, const char* d) { 

   CreateLargeString();
   sprintf(string, format, d);
}

// ----------------------------------------------------------------------------

// constructor for char

FS :: FS(char* format, const char d) { 

   CreateLargeString();
   sprintf(string, format, d);
}

// ----------------------------------------------------------------------------

void FS :: Print(FILE* out) {

   Mstream s_out;
   s_out.Associate(out);

   s_out << FS("Print: %s\n", string);
   s_out.Flush();
}

// ----------------------------------------------------------------------------

void FS :: CreateLargeString() { 

   string = new char[1000]; 
}

// ----------------------------------------------------------------------------

char* FS :: GetString() { 

   return string; 
}

// ----------------------------------------------------------------------------
//
// --> Extend stream io
//
// ----------------------------------------------------------------------------

// operator overloading of the stream io output operator
// specify how to print a formatted quantity

ostream& operator << (ostream& s, FS fs) {

   return s << fs.GetString(); 
}
    
// ----------------------------------------------------------------------------
//
// --> output stream for multiple destinations
//
// ----------------------------------------------------------------------------

Mstream :: Mstream() {

   target = Mstream_None;

   bufferSize = 1000;  // size of the buffer excluding the final '\0'
   bufferPos = 0;
   buffer = new char[bufferSize + 1];
}

// ----------------------------------------------------------------------------

Mstream :: ~Mstream() {

   Flush();
   delete buffer;
}

// ----------------------------------------------------------------------------

void Mstream :: Associate(FILE* out) { 

   Flush();   

   target = Mstream_file;
   fileOut = out; 
}

// ----------------------------------------------------------------------------
   
void Mstream :: Associate(ostream* out) { 
  
   Flush();   

   target = Mstream_ostream;
   ostreamOut = out; 
}

// ----------------------------------------------------------------------------

void Mstream :: Associate(char* out, int maxLength_I) { 

   Flush();   

   target = Mstream_charString;
   charStringOut = out;
   maxLength = maxLength_I;
   position = 0;
}
                               
// ----------------------------------------------------------------------------

void Mstream :: Associate(void (*out)(const char*)) { 

   Flush();   

   target = Mstream_function;
   functionOut = out; 
}

// ----------------------------------------------------------------------------

void Mstream :: Associate(
                   void (*out)(const char*, int, int), 
                   int unitNumber_I
                ) { 

   Flush();   

   target = Mstream_FORTRAN;
   FORTRAN_Out = out; 
   unitNumber = unitNumber_I;
}

// ----------------------------------------------------------------------------

void Mstream :: Associate(Mstream* out) { 

   Flush();   

   target = Mstream_Mstream;
   mstreamOut = out; 
}

// ----------------------------------------------------------------------------

void Mstream :: CharStringOut(const char* chstr) {

   int copied;
   int pos;

   switch (target) {

      case Mstream_file:
         fprintf(fileOut, "%s", chstr);
         break;

      case Mstream_ostream:
         *ostreamOut << chstr;
         break;

      case Mstream_charString:
         copied = copyStringSafe(
                     charStringOut + position, 
                     chstr, 
                     maxLength - position
                  );
         if (copied < 0)
            errf(-1, "Mstream :: CharStringOut: string to small");
         else
            position += copied;
         break;

      case Mstream_function:
         functionOut(chstr);
         break;

      case Mstream_FORTRAN:
         // first copy into the buffer until a (possible) '\n' is found
         pos = 0;         
         while (chstr[pos] != '\0') {
            // output until the next '\n'
            while (   chstr[pos] != '\0' 
                   && chstr[pos] != '\n'
                   && bufferPos < bufferSize
                  ) {
                  buffer[bufferPos++] = chstr[pos++];
            }
            if (bufferPos == bufferSize) {
               // terminate string (examination by means of a debugger...)
               buffer[bufferPos-1] = '\0';
               // :_CAUTION_: if errf(...) uses Mstream: what happens?
               errf(-1, "Mstream :: CharStringOut: buffer too small");
            }
            // output if there is a '\n'
            if (chstr[pos] == '\n') {
               ++pos;
               buffer[bufferPos] = '\0';
               FORTRAN_Out(buffer, bufferPos, unitNumber);
               bufferPos = 0;
            }
         }
         break;

      case Mstream_Mstream:
         *mstreamOut << chstr;
         break;

      default:
         errf(-1, "Mstream :: CharStringOut: unknown target");
   }
}

// ----------------------------------------------------------------------------

void Mstream :: Flush() {

   switch (target) {

      case Mstream_file:
         fflush(fileOut);
         break;

      case Mstream_ostream:
         *ostreamOut << flush;
         break;

      case Mstream_charString:
         break;

      case Mstream_function:
         break;

      case Mstream_FORTRAN:
         if (bufferPos > 0) {
            buffer[bufferPos] = '\0';
            FORTRAN_Out(buffer, bufferPos, unitNumber);
            bufferPos = 0;
         }
         break;

      case Mstream_Mstream:
         mstreamOut -> Flush();
         break;

      default:
         ;  // nothing to be done
   }
}

// ----------------------------------------------------------------------------
//
// --> Extend stream io
//
// ----------------------------------------------------------------------------

// operator overloading of the output operator          

Mstream& operator << (Mstream& m, const char* s) {

//   printf("\n\nlength=%d :%s:\n\n", (int) strlen(s), s);
//   fflush(stdout);

   m.CharStringOut(s); 

   return m; 
}
    
// ----------------------------------------------------------------------------

// specify how to print a formatted quantity

Mstream& operator << (Mstream& s, FS fs) {

   return s << fs.GetString(); 
}
    
// ----------------------------------------------------------------------------
//
// --> base class to check for already defined objects
//
// ----------------------------------------------------------------------------
 
DefinedObject::DefinedObject() {
   
   ResetDefined();
}

// ----------------------------------------------------------------------------
 
DefinedObject::~DefinedObject() {
}

// ----------------------------------------------------------------------------
 
void DefinedObject::ResetDefined() {

   isDefined = FALSE;
}

// ----------------------------------------------------------------------------
 
void DefinedObject::SetDefined() {

   isDefined = TRUE;
}

// ----------------------------------------------------------------------------
 
int DefinedObject::IsDefined() {

   if (!isDefined) {
      SetDefined();
      return FALSE;
   } else {
      return TRUE;
   }  
      
}

// ----------------------------------------------------------------------------
//
// --> named objects
//
// ----------------------------------------------------------------------------
 
NamedObject::NamedObject() {

   name = new String(-1);

   const char noNameTemplate[] = "NamedObject (no name)";
   noName = new char[strlen(noNameTemplate)+1];
   strcpy(noName, noNameTemplate);
}

// ----------------------------------------------------------------------------
  
NamedObject::~NamedObject() {

   delete name;
   delete noName;
}

// ----------------------------------------------------------------------------
  
void NamedObject::DefineName(char *nameIn){ 
   // :_CAUTION_: what happens with the old string contents 
   // if the object was created via ()
   name -> DefineAndCopyFromChar(nameIn);
}

// ----------------------------------------------------------------------------
  
void NamedObject::DefineName(const String *nameIn) { 
   name -> DefineAndCopyFromString(nameIn);
}

// ----------------------------------------------------------------------------
  
void NamedObject::DefineName(const String &nameIn) { 
   name -> DefineAndCopyFromString(&nameIn); 
}

// ----------------------------------------------------------------------------
  
char* NamedObject::GetName() {
   return name -> GetData();
}

// ----------------------------------------------------------------------------
  
char* NamedObject::GetNameSafe() {
   if (name -> GetData() != NULL)
      return name -> GetData();  
   else
      return noName;
}

// ----------------------------------------------------------------------------
//
// --> TrivialList (to trap non-destructed objects...)
//
// ----------------------------------------------------------------------------

TrivialList :: TrivialList(Global* globalIn) 
   : GlobalObject(globalIn)
   {

   maximum = 0;
}

// ----------------------------------------------------------------------------

TrivialList :: ~TrivialList() {

   delete[] index;
   delete[] count;
}

// ----------------------------------------------------------------------------

void TrivialList :: Define(int identifierIn, int maximumIn) {

   if (maximum != 0)
      errf(-1,"TrivialList :: Define: not empty");
   identifier = identifierIn;
   maximum = maximumIn;
   index = new long [maximum];
   count = new long [maximum];
   Clear();
}

// ----------------------------------------------------------------------------

void TrivialList :: Clear() {

   actual = 0;
   unique = 0;
   trap = -1;

//   Print(stdout);
}

// ----------------------------------------------------------------------------

long TrivialList :: Unique() {

   long ret;
   ret = unique++;

   return ret;
}

// ----------------------------------------------------------------------------

long TrivialList::Mark() {

   if (unique == trap) {
      global -> To_status() << FS("identifier=%d\n", identifier);
      global -> To_status() << FS("trap=%ld\n",trap);
      global -> To_status() << FS("*** :%s:\n", GetNameSafe());
      global -> To_status() << "*** TrivialList::Mark: trap found\n\n";
      // and terminate?
      if (is_soft_trap == -1) 
         errf(-1,"TrivialList::Mark: trap");
      else if (is_soft_trap == 1) 
         provoke_floating_point_exception();
   }
   if (actual == maximum) {
      global -> To_status() << FS("identifier=%d\n",identifier);
      global -> To_status() << "*** TrivialList::Mark: overflow\n\n";
      errf(-1, "TrivialList::Mark...");
   }
   count[actual] = 0;
   long ret;
   ret = Unique();
   index[actual++] = ret;

   return ret;
}

// ----------------------------------------------------------------------------

int TrivialList::Check(long id, int situation) {

   int ret = 0;
   register int i;

   if (id == trap) {

      // search for entry      
      for (i = 0; i < actual && ! (index[i] == id); ++i);
      if (i == actual) {
         Print(stdout);
         global -> To_status() << FS("identifier=%d\n",identifier);
         global -> To_status() << FS("id=%ld ", id)
                               << FS("situation=%d\n", situation);
         global -> To_status() << FS("*** :%s:\n", GetNameSafe());
         global -> To_status() << "*** TrivialList::Check: not found\n\n";
         errf(-1, "TrivialList::Check...");
      } else {
         count[i]++;
         if (trap_check_count == -1 || count[i] == trap_check_count) {
            ret = 1;
            global -> To_status() << FS("identifier=%d", identifier)
                                  << FS(" count=%ld", count[i])
                                  << FS(" situation=%d\n", situation);
            global -> To_status() << FS("trap=%ld\n",trap);
            global -> To_status() << FS("*** :%s:\n", GetNameSafe());
            global -> To_status() << "*** TrivialList::Check: trap found\n\n";
            // and terminate?
            if (is_soft_trap == -3) 
               errf(-1,"TrivialList::Check: trap");
            else if (is_soft_trap == 3) 
               provoke_floating_point_exception();
         }
      }
   }

   return ret;
}

// ----------------------------------------------------------------------------

void TrivialList::Remove(long remove) {

   register int i;
   for (i = 0; i < actual && ! (index[i] == remove); ++i);
   if (i == actual) {
      Print(stdout);
      global -> To_status() << FS("identifier=%d\n", identifier);
      global -> To_status() << FS("remove=%ld\n", remove);
      global -> To_status() << FS("*** :%s:\n", GetNameSafe());
      global -> To_status() << "*** TrivialList::Remove: not found\n\n";
      errf(-1, "TrivialList::Remove...");
   } else {
      index[i] = index[--actual];
      count[i] = count[  actual];
      ++ count[i];
   }
}

// ----------------------------------------------------------------------------

void TrivialList::Print(FILE *out) {

   Mstream s_out;
   s_out.Associate(out);

   register int i;
   
   s_out << "-----------------------------------------------------\n";
   s_out << FS("*** :%s:\n", GetNameSafe());
   s_out << FS("*** PRINT: identifier=%d\n", identifier);
   s_out << FS("trap was %ld\n", trap);
   for (i = 0; i < actual; ++i) {
       s_out << FS("%d:", i)
             << FS(" -> %ld\n", index[i]);
   }
   s_out << "End of PRINT.\n";
   s_out << "-----------------------------------------------------\n\n";
}

// ----------------------------------------------------------------------------

void TrivialList::SetTrap(long trapIn, 
                          int is_soft_trapIn, 
                          long trap_check_countIn) {
   trap=trapIn;
   is_soft_trap = is_soft_trapIn;
   trap_check_count = trap_check_countIn;
}

// ----------------------------------------------------------------------------
//
// --> trap objects
//
// ----------------------------------------------------------------------------
 
TrapObject::TrapObject() {

   SetTrapCreate(-2);
   SetTrapGetNoCreate(-2);
   SetTrapReturn(-2); 
}

// ----------------------------------------------------------------------------
 
TrapObject::~TrapObject() {
}

// ----------------------------------------------------------------------------

void TrapObject::GetTrapsFrom(TrapObject* to) {

   SetTrapCreate(to -> GetTrapCreate());
   SetTrapGetNoCreate(to -> GetTrapGetNoCreate());
   SetTrapReturn(to -> GetTrapReturn());
}

// ----------------------------------------------------------------------------

void TrapObject :: Notify(long trap) {

     trap = trap; // :_TAWM_:
}

// ----------------------------------------------------------------------------

Counter :: Counter() {

//   gl -> To_status() << "Counter :: Counter: constructor\n";

   maximum = -1;
   Reset();
}

// ----------------------------------------------------------------------------

Counter :: ~Counter() {

//   gl -> To_status() << "Counter :: ~Counter: destructor\n";
}

// ----------------------------------------------------------------------------

void Counter :: Reset() {

   counter = 0;
}

// ----------------------------------------------------------------------------

Boolean Counter :: Check() {

   if (maximum < 0 || counter <= maximum)
      return True;
   else
      return False;
}

// ----------------------------------------------------------------------------

void Counter :: Increment() {

   ++ counter;
}

// ----------------------------------------------------------------------------

Boolean Counter :: IncrementAndCheck() {

   Increment();
   return Check();
}

// ----------------------------------------------------------------------------

Boolean Counter :: BorderlineCase() {

   if (maximum < 0 || counter != maximum)
      return False;
   else
      return True;
}

// ----------------------------------------------------------------------------
//
// --> base class for objects that have a pointer to Global
//
// ----------------------------------------------------------------------------

GlobalObject :: GlobalObject() {

   global = NULL;
}

// ----------------------------------------------------------------------------

GlobalObject :: ~GlobalObject() {
}

// ----------------------------------------------------------------------------
//
// --> floating point trap signal
//
// ----------------------------------------------------------------------------

FPSignalHandler :: FPSignalHandler(Global* globalIn)
   : GlobalObject(globalIn)
   {

   global -> To_status() << "FPSignalHandler constructor.\n";

   newHandler_function = NULL;

   maxprint = 10; // maximum number of printed messages
   nprint = 0;

   Reset();
}

// ----------------------------------------------------------------------------

FPSignalHandler :: ~FPSignalHandler() {

   global -> To_status() << "FPSignalHandler destructor called.\n";

   // before we finish, we reinstall the old handler
   Reinstall_old();
}

// ----------------------------------------------------------------------------

int FPSignalHandler :: Status() {

   return flag;
}

// ----------------------------------------------------------------------------

void FPSignalHandler :: Set(int sigIN) {

   flag = sigIN;
}

// ----------------------------------------------------------------------------

void FPSignalHandler :: Reset() {

   counter = 0;
   Set(FALSE);
}

// ----------------------------------------------------------------------------

// Does not work if this is a method of the class (even if address 
// assignment is casted!)! Why? Cannot take the absolute address
// of a method of a class! 
// Therefore we have to access the flag via fpnan -> fpsh -> Set(TRUE)
// directly!

void FPSignalHandlerReceive(int sig) {

//   gl -> To_status() << FS("SIGFPE=%d\n", SIGFPE);

   if (sig == SIGFPE) {

      fpnan -> fpsh -> Increment_counter();

      if (
               fpnan -> fpsh -> nprint < fpnan -> fpsh -> maxprint
          && ! fpnan -> fpsh -> Get_flag()
         ) {
         ++ fpnan -> fpsh -> nprint;
         gl -> To_status()
                  << "FPSignalHandlerReceive: "
                  << FS("received floating point error signal = %d\n", sig);
      }
      if (fpnan -> fpsh -> nprint == fpnan -> fpsh -> maxprint) {
         ++ fpnan -> fpsh -> nprint;
         gl -> To_status()
                  << "FPSignalHandlerReceive: "
                  << "no more warning messages printed\n";
      }
      fpnan -> fpsh -> Set(TRUE);

#if 0
      // calling the original handler routine
      fpnan 
      -> fpsh 
      -> Get_p_originalHandler_action() 
      -> sa_handler(sig);
#endif

   } else {

      gl -> To_error()
               << FS("FPSignalHandlerReceive: received signal=%d\n", sig);
      gl -> To_error()
               << "Illegal in this routine. Terminating.\n";
      exit(-1);
   }
}

// ----------------------------------------------------------------------------

void FPSignalHandler :: Install(void (*newHandler_function_In)(int)) {

   global -> To_status() << "FPSignalHandler :: Install called.\n";

   if (newHandler_function != NULL)
      errf(-1, "FPSignalHandler :: Install: Install called twice...");

   newHandler_function = newHandler_function_In;
   
   int error;

   // get the original action
   error = sigaction(
              SIGFPE,
              NULL,
              &originalHandler_action
           );

   if (error != 0)
      errf(-1, "FPSignalHandler :: Install: call (1) of sigaction(...) failed");

//   global -> To_status() 
//      << "originalHandler_action.sa_mask:  "
//      << FS("%016x\n", (unsigned long) originalHandler_action.sa_mask);   
//   global -> To_status() 
//      << "originalHandler_action.sa_flags: "
//      << FS("%016x\n", (unsigned long) originalHandler_action.sa_flags);   

   newHandler_action.sa_handler = newHandler_function;
   newHandler_action.sa_mask    = originalHandler_action.sa_mask;
   newHandler_action.sa_flags   = originalHandler_action.sa_flags;

   struct sigaction dummy_action;

   // set new action
   error = sigaction (
              SIGFPE,
              &newHandler_action,
              &dummy_action
           );

   if (error != 0)
      errf(-1, "FPSignalHandler :: Install: call (2) of sigaction(...) failed");

   Set(FALSE);

//   global -> To_status() 
//      << FS("sigaction=%d\n", dummy);
//   global -> To_status() 
//      << FS("old sa_handler=%ld\n", (long) oldaction.sa_handler);
//   global -> To_status() 
//      << FS("old sa_mask=%ld\n", (unsigned long) oldaction.sa_mask);
//   global -> To_status()
//      << FS("old sa_flags=%d\n", oldaction.sa_flags);

}

// ----------------------------------------------------------------------------

void FPSignalHandler :: Reinstall_old() {

//   global -> To_status() << "FPSignalHandler :: Reinstall_old called.\n";
//   (global -> To_status()).Flush();

   int error;
   struct sigaction dummy_action;

   // set original action
   error = sigaction (
              SIGFPE,
              &originalHandler_action,
              &dummy_action
           );

   if (error != 0)
      errf(-1, "FPSignalHandler :: Reinstall_old: sigaction(...) failed");
}

// ----------------------------------------------------------------------------

void FPSignalHandler :: Reinstall_new() {

//   global -> To_status() << "FPSignalHandler :: Reinstall_new called.\n";
//   (global -> To_status()).Flush();

   if (newHandler_function == NULL)
      errf(-1, "FPSignalHandler :: Reinstall_new: no new handler installed");

   int error;
   struct sigaction dummy_action;

   // set new action
   error = sigaction (
              SIGFPE,
              &newHandler_action,
              &dummy_action
           );

   if (error != 0)
      errf(-1, "FPSignalHandler :: Reinstall_new: sigaction(...) failed");
}

// ----------------------------------------------------------------------------
//
// --> Floating point check
//
// ----------------------------------------------------------------------------

FPNaN_info :: FPNaN_info() {

   gl -> To_status() << "FPNaN_info constructor.\n";

   // have to copy, because direct return of the pointer did not
   // work with DEC cxx!

   const char no_message_errno_T[] 
                 = "FPNaN_info :: Get_msg_errno: no error no_message";
   no_message_errno = new char[strlen(no_message_errno_T)+1];
   strcpy(no_message_errno, no_message_errno_T);   

   const char no_message_mth_T[] 
                 = "--- empty ---";
   no_message_mth = new char[strlen(no_message_mth_T)+1];
   strcpy(no_message_mth, no_message_mth_T);   

   text_mth = NULL;

   Reset();
}

// ----------------------------------------------------------------------------

FPNaN_info :: ~FPNaN_info() {

   gl -> To_status() << "FPNaN_info destructor.\n";

   delete text_mth;   
   delete no_message_errno;
   delete no_message_mth;
}

// ----------------------------------------------------------------------------

void FPNaN_info :: Reset() {

   // reset the local flags

   flag_SIGFPE = FALSE;
   flag_errno  = FALSE;

   flag_NaN    = FALSE;
   counter_NaN = 0;

   flag_mth    = FALSE;
   error_mth   = 0;
   counter_mth = 0;
   delete text_mth;
   text_mth = NULL;
}

// ----------------------------------------------------------------------------

void FPNaN_info :: Copy_text_mth_From(const char* in) {

   delete text_mth;

   if (in == NULL)
      text_mth = NULL;
   else {
      text_mth = new char[strlen(in)+1];
      strcpy(text_mth, in);
   }
}

// ----------------------------------------------------------------------------

const char* FPNaN_info :: Get_msg_errno(int flag, int errno_in) {

   if (flag)
      return strerror(errno_in);
   else 
      return no_message_errno;
}

// ----------------------------------------------------------------------------

const char* FPNaN_info :: Get_msg_mth() {

#if 0
   if (text_mth == NULL) {
      (gl -> To_status()).Flush();
      gl -> To_status() << "\nno_message start\n";
      gl -> To_status() << ":\n:" << no_message_mth << ":\n";
      gl -> To_status() << "\nno_message end\n";
   }
#endif

   if (text_mth != NULL)
      return text_mth;
   else 
      return no_message_mth;
}

// ----------------------------------------------------------------------------

void FPNaN_info :: PrintStatus(Mstream& s_out) {

   s_out <<    "-- FPNaN_info :: PrintStatus:\n"
         << FS("-- flag_SIGFPE: %4d", flag_SIGFPE)
         << FS(" (# = %6d)\n", counter_SIGFPE)
         << FS("-- flag_errno:  %4d", flag_errno)
         << FS(" [%6d]\n", value_errno)
         << FS("-- flag_NaN:    %4d\n", flag_NaN)
         << FS("-- flag_mth:    %4d", flag_mth)
         << FS(" (# = %6d)", counter_mth)
         << FS(" [%6d]\n", error_mth)
         << FS("             %s\n", Get_msg_mth());

   if (flag_errno)
      s_out << FS("-- errno msg: %s\n", 
                  Get_msg_errno(flag_errno, value_errno)
               );
}

// ----------------------------------------------------------------------------

void FPNaN_info :: CopyInto(FPNaN_info* into) {

   into -> flag_SIGFPE    = flag_SIGFPE;
   into -> counter_SIGFPE = counter_SIGFPE;
   into -> flag_errno     = flag_errno;
   into -> value_errno    = value_errno;
   into -> flag_NaN       = flag_NaN;
   into -> counter_NaN    = counter_NaN;
   into -> flag_mth       = flag_mth;
   into -> error_mth      = error_mth;
   into -> counter_mth    = counter_mth;
   into -> Copy_text_mth_From(text_mth);
}

// ----------------------------------------------------------------------------

Boolean operator==(const FPNaN_info& a, const FPNaN_info& b) {

   return 

      (Boolean) (

         a.flag_SIGFPE    == b.flag_SIGFPE
      && a.counter_SIGFPE == b.counter_SIGFPE
      && a.flag_errno     == b.flag_errno
      && a.value_errno    == b.value_errno
      && a.flag_NaN       == b.flag_NaN
      && a.counter_NaN    == b.counter_NaN
      && a.flag_mth       == b.flag_mth
      && a.error_mth      == b.error_mth
      && a.counter_mth    == b.counter_mth

      )
   ;
}

// ----------------------------------------------------------------------------

FPNaN :: FPNaN(Global* globalIn, Boolean if_fpHandler) 
   : GlobalObject(globalIn)
   {

   global -> To_status() << FS("FPNaN constructor. (..., %d)\n", 
                               (int) if_fpHandler);

   fpsh = NULL;

#ifndef DISASTER_FPHANDLER
#define DISASTER_FPHANDLER TRUE
#endif
#if DISASTER_FPHANDLER
   Modify_fpHandler(if_fpHandler);
#endif

   Reset();
}

// ----------------------------------------------------------------------------

FPNaN :: ~FPNaN() {

   global -> To_status() << "FPNaN destructor.\n";

   if (fpsh != NULL)
      delete fpsh;
}

// ----------------------------------------------------------------------------

void FPNaN :: Modify_fpHandler(Boolean if_fpHandler) {

   delete fpsh;
   fpsh = NULL;

   if (if_fpHandler) {
      // install floating point exception handler
      fpsh = new FPSignalHandler(global);
      fpsh -> Install(FPSignalHandlerReceive);
   }
}

// ----------------------------------------------------------------------------

int FPNaN :: Get_flag_SIGFPE() {

   return (fpsh == NULL) ? flag_SIGFPE 
                         : (flag_SIGFPE = flag_SIGFPE || fpsh -> Status());
}

// ----------------------------------------------------------------------------

int FPNaN :: Get_counter_SIGFPE() {

   return (counter_SIGFPE 
           = (fpsh == NULL) ? 0 
                            : fpsh -> Get_counter()
          );
}

// ----------------------------------------------------------------------------

int FPNaN :: Get_flag_errno() {

   return (errno == 0) ? flag_errno : (flag_errno = TRUE);
}

// ----------------------------------------------------------------------------

int FPNaN :: Get_value_errno() {

   return (value_errno = errno);
}

// ----------------------------------------------------------------------------

void FPNaN :: Set_flag_mth(int in, const char* text) {

   ++ counter_mth;

   if ( ! flag_mth ) { 
      error_mth = in;
      flag_mth = TRUE;
      Copy_text_mth_From(text);
#if 0
      global -> To_error() << FS("%d, ", in)
                           << FS("%s\n", text);
      global -> To_error() << "FPNaN :: Set_flag_mth\n";
#endif
   };
}

// ----------------------------------------------------------------------------

void FPNaN :: Reset() {

   onHold = False;

   // reset the SIGFPE handler
   if (fpsh != NULL)
      fpsh -> Reset();

   // this is the errno-facility variable
   errno = 0;

   FPNaN_info :: Reset();
}

// ----------------------------------------------------------------------------

void FPNaN :: LogDoubleStatus(double d) {

   // check the bit pattern
   int newStatus;
   if ( (newStatus = (intnan8_C(d) == 1) ? TRUE : FALSE) )
      ++ counter_NaN;

   flag_NaN = flag_NaN || newStatus;
}

// ----------------------------------------------------------------------------

Boolean FPNaN :: Get_Error() {

   return (Boolean) (
             // in order to inquire the status of the handlers
             Get_flag_SIGFPE() 
          || Get_flag_errno()
          || Get_flag_NaN()  
          || Get_flag_mth()
          )
   ;
}

// ----------------------------------------------------------------------------

void FPNaN :: PrintStatus(Mstream& s_out) {

   int errnoF;
   int errnoV;

   s_out <<    "-- FPNaN :: PrintStatus:\n"
         << FS("-- flag_SIGFPE: %4d", Get_flag_SIGFPE())
         << FS(" (# = %6d)\n", Get_counter_SIGFPE())
         << FS("-- flag_errno:  %4d", (errnoF = Get_flag_errno()))
         << FS(" [%6d]\n", (errnoV = Get_value_errno()))
         << FS("-- flag_NaN:    %4d\n", Get_flag_NaN())
         << FS("-- flag_mth:    %4d", Get_flag_mth())
         << FS(" (# = %6d)", counter_mth)
         << FS(" [%6d]\n", error_mth)
         << FS("             %s\n", Get_msg_mth());

   if (errnoF)
      s_out << FS("-- errno msg: %s\n", 
                  Get_msg_errno(errnoF, errnoV)
               );
}

// ----------------------------------------------------------------------------

// put on hold

void FPNaN :: Hold() {

   if (onHold) 
      errf(-1, "FPNaN :: Hold(): already on hold");

   onHold = True;

   // update these
   (void) Get_flag_SIGFPE();
   (void) Get_counter_SIGFPE();
   (void) Get_flag_errno();
   (void) Get_value_errno();

   CopyInto(&saveStatus);

   // hold the signal handler
   if (fpsh != NULL)
      fpsh -> Reinstall_old();

//   gl -> To_status() << "FPNaN :: Hold()\n";
//   PrintStatus(gl -> To_status());
//   saveStatus.PrintStatus(gl -> To_status());
}

// ----------------------------------------------------------------------------

// continue

Boolean FPNaN :: Continue() {

   if ( ! onHold ) 
      errf(-1, "FPNaN :: Continue(): not on hold");

   // update these
   (void) Get_flag_SIGFPE();
   (void) Get_counter_SIGFPE();
   (void) Get_flag_errno();
   (void) Get_value_errno();

//   gl -> To_status() << "FPNaN :: Continue()\n";

   // compare
   Boolean changed =  (Boolean) (! (saveStatus == *this));

   saveStatus.CopyInto(this);

   // reinstall the signal handler
   if (fpsh != NULL)
      fpsh -> Reinstall_new();

   onHold = False;

   return changed;
}

// ----------------------------------------------------------------------------

int FPNaN :: DoubleStatus(double d) {

   int out;

   if (fpsh != NULL) {
      // check only the status term
      double e;
      e = d; // dummy :_CHECK_: what is this good for???
      e = e; // :_TAWM_:
      out = fpsh -> Status();
   } else {
      // check only the bit pattern
      int cint = intnan8_C(d);
      out = (cint == 1) ? TRUE : FALSE;
   }

   if (out) {
      global -> To_error() << FS("FPE error=%d\n", out);
   }
      
   if (errno != 0) {
      global -> To_error() << FS("errno=%d\n", errno);
      global -> To_error() << FS(":%s:\n", strerror(errno));
      out = TRUE;
//      errf(-1, "errno STOP");
   }
      
   return out;
}

// ----------------------------------------------------------------------------

int FPNaN :: DoubleStatusReset(double d) {

   int out;

   out = DoubleStatus(d);

   Reset();

   return out;
}

// ----------------------------------------------------------------------------
//
// --> collection of global objects
//
// ----------------------------------------------------------------------------

Global :: Global() {

   // for the time being, associate the standard output

   statusStream.Associate(&cout);
   errorStream.Associate(&statusStream);
   logStream.Associate(&statusStream);

   if_fpHandler = True;

//   To_status() << "Global: constructor.\n";

   gl = this; // to be used in routines where ``global''
              // cannot be passed; ``gl'' is a global variable

   pMapSquareLinear=new MapSquareLinear(0,0);
   pPermute=new Permute(0);

   // default values

   rho_lim=0.5;
   alpha_lim=0.4;

   ia=0;
   ib=1;

   i0  =  0;
   i1  =  1;
   i2  = -1;
   i3  =  2;
   i2p = -1;
   i3p =  2;

   process=NULL;
   user=NULL;

   mel_index=0;
   process_index=0;
   numberOfFinalStatePartonsInBornTerm=0;
   subtraction_type=1;

   tlist=new (TrivialList*[10]);

#if CHECKDESTRUCT
#define LENGTH_0  1000
#define LENGTH_1  1000000
#else
#define LENGTH_0  1
#define LENGTH_1  1
#endif

   tEvent = new TrivialList(this);
   tEvent -> Define(0, LENGTH_0);  
   tEvent -> DefineName("TrivialList: Event");
   tlist[0] = tEvent;

   tParticle = new TrivialList(this);
   tParticle -> Define(1, LENGTH_1);  
   tParticle -> DefineName("TrivialList: Particle");
   tlist[1] = tParticle;

   idEvent = new TrivialList(this);
   idEvent -> Define(-1,1);  
}

// ----------------------------------------------------------------------------

Global :: ~Global() {

   To_status() << "Global: destructor called.\n";

   delete pMapSquareLinear;
   delete pPermute;

#if CHECKDESTRUCT
      tEvent.Print(stdout);
      tParticle.Print(stdout);
#endif

   delete[] tlist;

   delete tEvent;
   delete tParticle;
   delete idEvent;

   // flush output (if it is buffered)
   statusStream.Flush();
   errorStream.Flush();
   logStream.Flush();
}

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

// ============================================================================
//
// --> combinatorical classes
//
// file:              cmb.cc
// created:           29.03.1997
// last modification: 14.11.1997
//
// ============================================================================

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "global.h"
#include "cmb.h"
#include "mth.h"

// ============================================================================
//
// --> classes
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> Permutations
//
// ----------------------------------------------------------------------------

Permute :: Permute(int numberI) {

   list = NULL;
   startAt = NULL;
   Define(numberI);
}

// ----------------------------------------------------------------------------

Permute :: ~Permute() {

   register int i, n;

   if (list!=NULL) {
      for (n=0;n<number;++n) {
          if (list[n]!=NULL) {
             for (i=0;i<length[n];++i)  
                 delete [] list[n][i];
             delete list[n];
          }
          if (startAt[n]!=NULL)
             delete startAt[n];
      }
      delete [] list;
      delete [] startAt;
      delete [] length;
   }
}

// ----------------------------------------------------------------------------

void Permute::Define(int numberI) {

   int n;

   number=numberI+1;

   if (list!=NULL)
      errf(-1,"Permute::Define.");
   if (number<=1) {
      length=NULL;
      list=NULL;
   } else {
      length=new int [number];
      list=new (int** [number]);
      startAt=new (int* [number]);
      for (n=0;n<number;++n) {
          list[n]=NULL;
          startAt[n]=NULL;
      }
   }
}

// ----------------------------------------------------------------------------

// creates the list of all permutations

void Permute::CreateList(int n) {

   register int i,j,k,l;
   int br,co;
   int *c,*p,*hits;

   if ( n >= number || n == 0 || list[n] != NULL || startAt[n] != NULL )
      errf(-1,"Permute::CreateList: error.");

   c = new int [n];
   p = new int [n];
   hits = new int [n];

   length[n]=lfactorial(n);
   list[n]=new (int* [length[n]]);
   for (i=0;i<length[n];++i) {
       list[n][i]=new int [n];
   }
   startAt[n]=new (int [length[n]]);

#if 1

   // new version

   for (i=0; i<n; ++i) {
       c[i]=n-i-1;  
       p[i]=i;
   }

   co=0; 
   br = FALSE;

   do {
      // store
      for (i=0;i<n;++i) {
          list[n][co][i]=p[i];
      }

      // determine first position from the right where entry changes
      if (co == 0)
         startAt[n][co]=0;
      else {
         int s = 0;
         while ( s < n && list[n][co][s] == list[n][co-1][s] )
            ++s;
         startAt[n][co]=s;
      }
     
      // create next permutation

      // look for first opportunity to increment
      i=n-2;
      while ( i >= 0 && c[i] == 0 )
         --i;
      
      // terminate?
      if ( i < 0 || i== 0 && c[0] == 0 ) {
         br = TRUE;
      } else {

         // find numbers for the others
         for (j=i; j<n; ++j) {
             if (j==i) {
                c[i]--;
                l=p[i]+1;
             } else {
                c[j]=n-j-1;
                l=0;
             }
             for (k=0; k<n; ++k)
                 hits[k]=FALSE;
             for (k=0; k<j; ++k)
                 hits[p[k]]=TRUE;
//             while (l < n && hits[l] == TRUE)
             while (l < n && hits[l])
                ++l;
             p[j]=l;
         }
      }

      co++;

   } while (!br);

#else

   // old version

   for (i=0;i<n;++i) 
       c[i]=1;  

   br=FALSE;
   co=0; 

   while (!br) {

      // create current permutation
      for (i=0;i<n;++i)
          p[i]=-1;
      for (i=n-1;i>=0;--i) {
          k=-1;
          l=0;
          while (l<c[i]) {               
             ++k;
             if (p[k]==-1)
                ++l;
          }
          p[k]=i;
      }

      for (i=0;i<n;++i) {
          list[n][co][i]=p[i];
          // determine first position from the right where entry changes
          if (co == 0)
             startAt[n][co]=0;
          else {
             int s = 0;
             while ( s < n && list[n][co][s] == list[n][co-1][s] )
                ++s;
             startAt[n][co]=s;
          }
      }
     
      // prepare for next one
      for (i=n-1;i>=0 && (i==n-1 || c[i+1]==1);--i) {
          ++c[i];
          if (c[i]==i+2)
             c[i]=1;
      }
      br=(i==-1);
      ++co;
   }
#endif

   delete[] c;
   delete[] p;
   delete[] hits;
}

// ----------------------------------------------------------------------------

void Permute :: PrintList(FILE* file) {

   Mstream s_out;
   s_out.Associate(file);

   int n,co,i;
   for (n=1;n<number;++n) {
      s_out << FS("------------ Permutation n=%d -------------\n", n);
      for (co=0;co<length[n];++co) {
          s_out << FS("%3d", co)
                << FS(" [ %3d ] --> ", startAt[n][co]);
          for (i=0;i<n;++i)
              s_out << FS("%3d ", list[n][co][i]);
          s_out << "\n";
      }
   }
   s_out << "-------------------------------------------\n";
}

// ----------------------------------------------------------------------------

void Permute::CreateListOfLists(int nmin,int nmax) {
 
   register int n;

   for (n=nmin;n<=nmax;++n)
       CreateList(n);
}

// ----------------------------------------------------------------------------

int **Permute :: GetList(int n) {
 
   if (n>=number || n==0 || list[n]==NULL) {
      gl -> To_error() << FS("number,n=%d", number)
                       << FS(",%d\n", n);
      errf(-1,"Permute::GetList: error.");
   }
   return list[n];
}

// ----------------------------------------------------------------------------

int *Permute::GetStartAt(int n) {
 
   if (n>=number || n==0 || startAt[n]==NULL) {
      gl -> To_error() << FS("number,n=%d", number)
                       << FS(",%d\n", n);
      errf(-1,"Permute::GetStartAt: error.");
   }
   return startAt[n];
}

// ----------------------------------------------------------------------------

void Permute::PrintPartial(FILE* file) {

   Mstream s_out;
   s_out.Associate(file);

   int i;
   s_out << "------------ Partial Fractions -------------\n";
   for (i=0;i<nS;++i) {
       s_out << FS("%3d", i)
             << FS(" --> %14.6e\n", dataS[i]);
   }
   s_out << "--------------------------------------------\n";
}

// ----------------------------------------------------------------------------

// This returns the sum over all permutations of
// x_1 * ... * x_n / ([x_s1,...,x_sn]) with fixed s1 == ifixed.
// Note that this expression has no pole in x_s1, and that the
// sum over all s1 of these expressions yields 1.

double Permute::SumOverPermutations(double *data, int n, int ifixed) {

   register int i;
   register int co;

   // save for possible use by the printing routine
   dataS = data;
   nS = n;

   // get permutation list; also checks whether available
   int** plist    = GetList(n-1);
   int   plength  = length[n-1];
   int*  pStartAt = GetStartAt(n-1);

   // create the pointer list and calculate the numerator
   double* sdata = new double [n];

   sdata[n-1] = data[ifixed];
   double numerator = 1.;
   double* dataPtr = data;
   double* sdataPtr = sdata;
   for (i = 0; i < ifixed; ++i) {
       *(sdataPtr++)   = *dataPtr;
       numerator *= *(dataPtr++);
   }
   dataPtr++;
   for (i = ifixed + 1; i < n; ++i) {
       *(sdataPtr++) = *dataPtr;
       numerator *= *(dataPtr++);
   }

   // sum over all permutations
   double sum = 0;
   double* prodList = (new double [n+1]) + 1;
   double* sumList  = (new double [n+1]) + 1;
   prodList[-1] = 1.;
   sumList [-1] = sdata[n-1];

   double* prodList_n_2 = prodList + n - 2;
   int* prow_ptr;
   register double* sumList_i;
   register double* sumList_i_1;
   register double* prodList_i;
   register double* prodList_i_1;
   int* pStartAtCoPtr = pStartAt;
   for (co = 0; co < plength; ++co) {
       int pStartAtCo = *(pStartAtCoPtr++);
       prow_ptr     = plist[co] + pStartAtCo;
       sumList_i    = sumList  + pStartAtCo;
       sumList_i_1  = sumList_i  - 1;
       prodList_i   = prodList + pStartAtCo;
       prodList_i_1 = prodList_i - 1;
       for (i = pStartAtCo; i < n - 1; ++i) {
           *sumList_i  = *sumList_i_1 + sdata[*(prow_ptr++)];
           *prodList_i = *prodList_i_1 * *sumList_i;
           sumList_i_1 = sumList_i;
           ++sumList_i;
           prodList_i_1 = prodList_i;
           ++prodList_i;
       }
       sum += 1. / *prodList_n_2;
   }

   delete [] sdata;
   delete [] (prodList-1);
   delete [] (sumList-1);

   return numerator * sum;
}

// ----------------------------------------------------------------------------
//
// --> mapping of invariants to a linear list and vice versa
//
// ----------------------------------------------------------------------------

MapSquareLinear :: MapSquareLinear(int numberI, int inI) {

   mapij = NULL;
   Define(numberI, inI);
}

// ----------------------------------------------------------------------------

MapSquareLinear :: ~MapSquareLinear() {

   Delete();
}

// ----------------------------------------------------------------------------

void MapSquareLinear :: Define(int numberI, int inI) {

   register int n;
   number = numberI + 1;  // maximum number of outgoing partons
   in = inI;
   if (mapij != NULL) {
      gl -> To_error() << FS("number=%d\n",number);
      errf(-1,"MapSquareLinear::Define.");
   }
   if (number<=1) {
      mapij=NULL;
   } else {
      mapij=new (int** [number]);
      imap=new (int* [number]);
      jmap=new (int* [number]);
      ijmaplength=new int [number];
      for (n=0;n<number;++n) {
          mapij[n]=NULL;
          imap[n]=NULL;
          jmap[n]=NULL;
          ijmaplength[n]=-1;
      }
   }
}

// ----------------------------------------------------------------------------

void MapSquareLinear :: Delete() {

   register int i,n;
   if (mapij!=NULL) {
      for (n=0;n<number;++n) 
          if (mapij[n]!=NULL) {
             for (i=-in;i<n;++i)
                 delete [] (mapij[n][i]-in);
             delete [] (mapij[n]-in);
             delete [] imap[n];
             delete [] jmap[n];
          }
      delete [] mapij;
      delete [] imap;
      delete [] jmap;
      delete [] ijmaplength;
   }
   mapij=NULL;
}

// ----------------------------------------------------------------------------

void MapSquareLinear::Print(FILE *file) {

   Mstream s_out;
   s_out.Associate(file);

   int i,j,n;
   for (n=1;n<number;++n) {
       s_out << FS("------------ MapSquareLinear n=%3d -------------\n", n);
       s_out << FS("ijmaplength=%d\n", ijmaplength[n]);
       if (mapij[n]!=NULL) 
          for (i=-in;i<n;++i)
              for (j=-in;j<n;++j) {
                  s_out << FS("%3d", i)
                        << FS(",%3d", j)
                        << FS(" --> %3d", mapij[n][i][j]
                         );
                  if (mapij[n][i][j]!=-1)
                     s_out << FS(" -> %3d", imap[n][mapij[n][i][j]])
                           << FS(" %3d",    jmap[n][mapij[n][i][j]]);
                  s_out << "\n";
              }
       else
          s_out << "NULL\n";
       s_out << "------------------------------------------------\n";
   }
}

// ----------------------------------------------------------------------------

void MapSquareLinear::CreateList(int n) {
 
   register int i,j,count;

   if (n>=number || n==0 || mapij[n]!=NULL)
      errf(-1,"MapSquareLinear: error.");

   mapij[n]=new (int* [n+in])+in; // shift to allow for negative indices
   for (i=-in;i<n;++i)
       mapij[n][i]=new int [n+in]+in; // shift to allow for negative indices
   
   ijmaplength[n]=((n+in)*(n+in-1))/2;
   // one less in case of two incident partons
   if (in==2)
      --ijmaplength[n];
   imap[n]=new int [ijmaplength[n]];
   jmap[n]=new int [ijmaplength[n]];

   count=0;
   for (i=-in;i<n-1;++i) {
       for (j=i+1;j<n;++j) {
           if (!(i<0 && j<0)) {
              mapij[n][i][j]=count++;
              mapij[n][j][i]=mapij[n][i][j];
              imap[n][mapij[n][i][j]]=i;
              jmap[n][mapij[n][i][j]]=j;
           }
       }
       mapij[n][i][i]=-1; // no diagonal elements
   }
   mapij[n][n-1][n-1]=-1;
   if (in==2) {
      mapij[n][-1][-2]=-1;
      mapij[n][-2][-1]=-1;
   }
}

// ----------------------------------------------------------------------------

void MapSquareLinear::CreateListOfLists(int nmin,int nmax) {
 
   register int n;

   for (n=nmin;n<=nmax;++n)
       CreateList(n);
}

// ----------------------------------------------------------------------------

int **MapSquareLinear :: GetSquareToLinearSafe(int n) {
 
   if (n >= number || n == 0 || mapij[n] == NULL) {
      gl -> To_error() << FS("number,n=%d", number)
                       << FS(",%d\n", n);
      errf(-1,"MapSquareLinear::GetSquareToLinear: error.");
   }

   return mapij[n];
}

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

#if 0
=========================== GARBAGE ===========================================
=========================== END OF GARBAGE ====================================
#endif
// ============================================================================
//
// --> DISASTER++ matrix elements
//
// file:              me.cc
// created:           05.04.1997
// last modification: 12.12.1997
//
// ============================================================================

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "switch.h"
#include "enum.h"
#include "global.h"
#include "mth.h"
#include "strng.h"
#include "user.h"
#include "mcintg.h"
#include "cmb.h"
#include "qcd.h"
#include "disaster.h"
#include "proc.h"
#include "rest.h"
#include "fold.h"

#include "me.h"

// ============================================================================
//
// definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

// ============================================================================
//
// --> classes
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> DIS Matrix Elements
//
// ----------------------------------------------------------------------------

// static variables

long   DISMatrixElement::lastVirtual;
double DISMatrixElement::virtqSave;
double DISMatrixElement::virtgSave;

// ----------------------------------------------------------------------------

DISMatrixElement::DISMatrixElement(DISSubProcess mel_indexIn) {

   mel_index=mel_indexIn;

   const char *tname[] = {
                           "DISMatrixElement q --> q", 
                           "DISMatrixElement q --> qg", 
                           "DISMatrixElement g --> qa", 
                           "DISMatrixElement q --> qgg", 
                           "DISMatrixElement g --> qag", 
                           "DISMatrixElement q --> qqa", 
                           "DISMatrixElement q --> qg V", 
                           "DISMatrixElement g --> qa V", 
                           "DISMatrixElement q --> q V", 
                           "DISMatrixElement q --> qg A", 
                           "DISMatrixElement g --> qa A", 
                           "DISMatrixElement q --> qgg A", 
                           "DISMatrixElement g --> qag A", 
                           "DISMatrixElement q --> qqa A", 
                           "DISMatrixElement q --> qg RS", 
                           "DISMatrixElement g --> qq RS", 
                           "DISMatrixElement q --> q FS", 
                           "DISMatrixElement q --> qg FS", 
                           "DISMatrixElement g --> qq FS", 
                           "DISMatrixElement q --> qg RSNF", 
                           "DISMatrixElement g --> qq RSNF", 
                           "DISMatrixElement g --> qq FSNF", 
                           "DIS undefined", 
                           "DIS test", 
                         };

   const char *particle_assignment[] = {
                           "qq ?", 
                           "qqg ?", 
                           "gqa ?", 
                           "qqgg ?", 
                           "gqag ?", 
                           "qqqa ?", 
                           "qqg ?!!", 
                           "gqa ?!!", 
                           "qq ?", 
                           "qqg ?", 
                           "gqa ?", 
                           "qqgg ?", 
                           "gqag ?", 
                           "qqqa ?", 
                           "qqg ?", 
                           "gqa ?", 
                           "qq ?", 
                           "qqg ?", 
                           "gqa ?", 
                           "qqg ?", 
                           "gqa ?", 
                           "gqa ?", 
                           "test",
                           "undefined ?"
                         };

   const int nGraphsValue[] = {
                               1, 2, 2, 
                               8, 8, 8, 
                               0, 0, 1, 
                               0, 0, 0, 0, 0, 
                               2, 2, 
                               0, 0, 0, 
                               2, 2, 
                               2, 
                               0, 0
                              };

   const int nSoftLengthValue[]  = {
                                    0, 0, 0, 
                                    0, 0, 0, 
                                    0, 0, 0, 
                                    3, 3, 6, 6, 6, 
                                    0, 0, 
                                    0, 0, 0, 
                                    0, 0, 
                                    0, 
                                    0, 0
                                   };

   int cindex;

   switch (mel_index) {
      case qXq:
         cindex=0;
         break;
      case qXqg:
         cindex=1;
         break;
      case gXqq:
         cindex=2;
         break;
      case qXqgg:
         cindex=3;
         break;
      case gXqqg:
         cindex=4;
         break;
      case qXqqq:
         cindex=5;
         break;
      case qXqgV:
         cindex=6;
         break;
      case gXqqV:
         cindex=7;
         break;
      case qXqV:
         cindex=8;
         break;
      case qXqgA:
         cindex=9;
         break;
      case gXqqA:
         cindex=10;
         break;
      case qXqggA:
         cindex=11;
         break;
      case gXqqgA:
         cindex=12;
         break;
      case qXqqqA:
         cindex=13;
         break;
      case qXqgRS:
         cindex=14;
         break;
      case gXqqRS:
         cindex=15;
         break;
      case qXqFS:
         cindex=16;
         break;
      case qXqgFS:
         cindex=17;
         break;
      case gXqqFS:
         cindex=18;
         break;
      case qXqgRSNF:
         cindex=19;
         break;
      case gXqqRSNF:
         cindex=20;
         break;
      case gXqqFSNF:
         cindex=21;
         break;
      case testDIS:
         cindex=22;
         break;
      case NoDISSubProcess:
      default:
         cindex=23;  // TAWM
         errf(-1,"mel_index not defined (1)");
   }

   // create array for cross section contributions
   nGraphs = nGraphsValue[cindex];
   termArray.Define(nGraphs, nGraphs);
   if (nGraphs > 0)
      term = termArray.GetDataPtrWithNullWarning();
   else
      term = NULL;

   // create array for soft coefficients
   nSoftLength = nSoftLengthValue[cindex];
   if (nSoftLength > 0) {
      softCoefficient = new double [nSoftLength];
      softAngle       = new double [nSoftLength];
      softContributes = new int    [nSoftLength];
      index1          = new int    [nSoftLength];
      index2          = new int    [nSoftLength];
   } else {
      softCoefficient = NULL;
      softAngle       = NULL;
      softContributes = NULL;
      index1          = NULL;
      index2          = NULL;
   }

   name -> DefineAndCopyFromChar(tname[cindex]);

   particleAssignment -> DefineAndCopyFromChar(particle_assignment[cindex]);

   lastVirtual = -1;

   // set number of required integration variables
   switch (mel_index) {
      case qXq: 
      case qXqV: 
         Set_nIntVar(-1);
         break;
      case qXqg: 
      case gXqq: 
      case qXqgV: 
      case gXqqV: 
      case qXqgRS: 
      case gXqqRS: 
      case qXqgRSNF: 
      case gXqqRSNF: 
      case gXqqFSNF: 
         Set_nIntVar(2);
         break;
      case qXqgg: 
      case gXqqg: 
      case qXqqq:
      case testDIS:
         Set_nIntVar(5);
         break;
      case qXqgA: 
      case gXqqA: 
      case qXqFS: 
         Set_nIntVar(0);
         break;
      case qXqggA: 
      case gXqqgA: 
      case qXqqqA:
      case qXqgFS: 
      case gXqqFS: 
         Set_nIntVar(3);
         break;
      case NoDISSubProcess:
      default:    
         errf(-1,"mel_index not defined (nIntVar)");    
   }
}

// ----------------------------------------------------------------------------

DISMatrixElement::~DISMatrixElement() {

   delete [] softCoefficient;
   delete [] softAngle;
   delete [] softContributes;
   delete [] index1;
   delete [] index2;
}

// ----------------------------------------------------------------------------

INLINE void DISMatrixElement::Define() {
}

// ----------------------------------------------------------------------------

INLINE int DISMatrixElement::GetIn() {
   int ret;
   ret=-1; // TAWM
   switch (mel_index) {
      case qXq: 
      case qXqg: 
      case gXqq: 
      case qXqgg: 
      case gXqqg: 
      case qXqqq: 
      case qXqgV: 
      case gXqqV: 
      case qXqV: 
      case qXqgA: 
      case gXqqA: 
      case qXqggA: 
      case gXqqgA: 
      case qXqqqA: 
      case qXqgRS: 
      case gXqqRS: 
      case qXqFS: 
      case qXqgFS: 
      case gXqqFS: 
      case qXqgRSNF: 
      case gXqqRSNF: 
      case gXqqFSNF: 
      case testDIS: 
         ret=1;
         break;
      case NoDISSubProcess:
      default:
         errf(-1,"mel_index not defined (2)");
         ret=1;
   }
   return ret;
}

// ----------------------------------------------------------------------------

INLINE int DISMatrixElement::GetOut() {
   int ret;
   ret=-1; // TAWM
   switch (mel_index) {
      case qXq: 
      case qXqV: 
      case qXqgA: 
      case gXqqA: 
      case qXqFS: 
         ret=1;
         break;
      case qXqg: 
      case gXqq: 
      case qXqgV: 
      case gXqqV: 
      case qXqggA: 
      case gXqqgA: 
      case qXqqqA:
      case qXqgRS: 
      case gXqqRS: 
      case qXqgFS: 
      case gXqqFS: 
      case qXqgRSNF: 
      case gXqqRSNF: 
      case gXqqFSNF: 
         ret=2;
         break;
      case qXqgg: 
      case gXqqg: 
      case qXqqq:
         ret=3;
         break;
      case testDIS:
         ret=3;
         break;
      case NoDISSubProcess:
      default:    
         errf(-1,"mel_index not defined (3)");    
         ret=3;
   }
   return ret;
}

// ----------------------------------------------------------------------------

INLINE int DISMatrixElement::GetSpecialCut() {
   return FALSE;
}

// ----------------------------------------------------------------------------

// assumes that invariants have already been calculated in ``e''

void DISMatrixElement::Evaluate(
                            Event* e,
                            Contribution *con,
                            double factor
                         ) {

   // kinematical variables
   CopyPartonInvariantsFrom(e);

   double sterm[5]; // contents: gluon initiated or quark initiated
   double aterm[5]; // contents: anti-quark initiated

   EvaluateLorentzAndColour(e, sterm);

   double bfactor;

   // order of the fine structure constant
   con -> Set_orderAlphaEM(2);

   // factors of FourPi & order of the strong coupling constant

   switch (mel_index) {

      case qXq: 
         bfactor = 1.;
         con -> SetOrderAlphaS(0);
         break;

      case qXqV: 
      case qXqg: 
      case gXqq: 
      case qXqgA: 
      case gXqqA: 
      case qXqFS: 
         bfactor = FourPi;
         con -> SetOrderAlphaS(1);
         break;

      case qXqgg: 
      case gXqqg: 
      case qXqqq: 
      case qXqggA: 
      case gXqqgA: 
      case qXqqqA: 
      case testDIS:
      case qXqgV: 
      case gXqqV: 
      case qXqgRS: 
      case gXqqRS: 
      case qXqgFS: 
      case gXqqFS: 
      case qXqgRSNF: 
      case gXqqRSNF: 
      case gXqqFSNF: 
         bfactor = FourPi * FourPi;
         con -> SetOrderAlphaS(2);
         break;

      case NoDISSubProcess:
      default:
         errf(-1,"mel_index not defined (5)");
         // :_TAWM_:
         bfactor = 0.;
   }

   // electroweak coupling and cross section formula factors
   
//   double alphaEM = e -> GetEventWithAlphaEM() -> GetAlphaEM();
#if 0
   double alphaEM = global 
                    -> disaster 
                    -> Get_electroWeak() 
                    -> Get_alphaEM_Server()
                    -> EvaluateAlphaEM(e -> Q2);

   bfactor *=   alphaEM * alphaEM * e -> y
              * global -> disaster -> GetGeV_2IntoPb()
              / ( 8. * e -> xi_hard * e -> Q2 * e -> Q2 );
#endif

   bfactor *=   e -> y
              * global -> disaster -> GetGeV_2IntoPb()
              / ( 8. * e -> xi_hard * e -> Q2 * e -> Q2 );

   // compensate for the explicit lepton phase space (already included
   // in the cross section formula
   bfactor *= 16.*Pi*Pi / (e->SH * e->y);

   // -------------------------------------------------------------------------

   // the terms to be summed up individually

   // -------------------------------------------------------------------------

   // type of matrix element: for simplicity here only the index (... wrong)
   
   con -> SetMatrixElement(mel_index);

   // -------------------------------------------------------------------------

   // scale logarithm?

   switch (mel_index) {

      case qXq: 
      case qXqV: 
      case qXqg: 
      case gXqq: 
      case qXqgg: 
      case gXqqg: 
      case qXqqq: 
      case qXqgA: 
      case gXqqA: 
      case qXqggA: 
      case gXqqgA: 
      case qXqqqA: 
      case testDIS:
      case qXqgV: 
      case gXqqV: 
         con -> SetSLF(NoScaleLogarithm);
         break;

      case qXqgRS: 
      case gXqqRS: 
         con -> SetSLF(RenormalizationScaleLogarithm);
         break;

      case qXqgRSNF: 
      case gXqqRSNF: 
         con -> SetSLF(RenormalizationScaleLogarithmWithNf);
         break;

      case qXqFS: 
      case qXqgFS: 
      case gXqqFS: 
         con -> SetSLF(FactorizationScaleLogarithm);
         break;

      case gXqqFSNF: 
         con -> SetSLF(FactorizationScaleLogarithmWithNf);
         break;

      case NoDISSubProcess:
      default:
         errf(-1,"mel_index not defined (a)");
   }

   // -------------------------------------------------------------------------

   // calculation for incident antiquarks; here: copy

   switch (mel_index) {

      case qXq: 
      case qXqV: 
      case qXqg: 
      case qXqgA: 
      case qXqgV: 
      case qXqgRS: 
      case qXqgRSNF: 
         aterm[0] = sterm[0];
         break;

      case qXqgg: 
      case qXqggA: 
      case testDIS:
         aterm[0] = sterm[0];
         aterm[1] = sterm[1];
         break;

      case qXqqq: 
      case qXqqqA: 
         aterm[0] = sterm[0];
         aterm[1] = sterm[1];
         aterm[2] = sterm[2];
         aterm[3] = sterm[3];
         aterm[4] = sterm[4];
         break;

      case gXqq: 
      case gXqqA: 
      case gXqqV: 
      case gXqqg: 
      case gXqqgA: 
      case gXqqRS: 
      case gXqqRSNF: 
      case gXqqFSNF: 
         // do nothing
         break;

      case qXqFS: 
      case qXqgFS: 
      case gXqqFS: 
         aterm[0] = sterm[0];
         break;

      case NoDISSubProcess:
      default:
         errf(-1,"mel_index not defined (a)");
   }

   // -------------------------------------------------------------------------

   // now multiply by kinematics factor and global factor, 
   // and assign to the array

   double *cdata = con -> GetData();

   bfactor *= factor * global -> VEGAS_weight_current;

   switch (mel_index) {

      case qXq: 
      case qXqV: 
      case qXqg: 
      case qXqgA: 
      case qXqgV: 
      case qXqgRS: 
      case qXqgRSNF: 

         sterm[0] *= bfactor;
         aterm[0] *= bfactor;

         cdata[0] += sterm[0];
         cdata[1] += aterm[0];

         break;

      case qXqgg: 
      case qXqggA: 
      case testDIS:

         sterm[0] *= bfactor;
         aterm[0] *= bfactor;
         sterm[1] *= bfactor;
         aterm[1] *= bfactor;

         cdata[0] += sterm[0] + sterm[1];
         cdata[1] += aterm[0] + aterm[1];

         break;

      case qXqqq: 
      case qXqqqA: 

         sterm[0] *= bfactor;
         aterm[0] *= bfactor;
         sterm[1] *= bfactor;
         aterm[1] *= bfactor;
         sterm[2] *= bfactor;
         aterm[2] *= bfactor;
         sterm[3] *= bfactor;
         aterm[3] *= bfactor;
         sterm[4] *= bfactor;
         aterm[4] *= bfactor;

         cdata[ 3] += sterm[0] + sterm[1];
         cdata[ 4] += aterm[0] + aterm[1];
         cdata[ 5] += sterm[2];
         cdata[ 6] += aterm[2];
         cdata[ 7] += sterm[3];
         cdata[ 8] += aterm[3];
         cdata[ 9] += sterm[4];
         cdata[10] += aterm[4];

         break;

      case gXqq: 
      case gXqqA: 
      case gXqqV: 
      case gXqqRS: 
      case gXqqRSNF: 
      case gXqqFSNF: 

         sterm[0] *= bfactor;

         cdata[2] += sterm[0];

         break;

      case gXqqg: 
      case gXqqgA: 

         sterm[0] *= bfactor;
         sterm[1] *= bfactor;

         cdata[2] += sterm[0] + sterm[1];

         break;

      case qXqFS: 

         sterm[0] *= bfactor;
         aterm[0] *= bfactor;
         sterm[1] *= bfactor;

         cdata[0] += sterm[0];  // q
         cdata[1] += aterm[0];  // a
         cdata[2] += sterm[1];  // g

         break;

      case qXqgFS: 

         sterm[0] *= bfactor;
         aterm[0] *= bfactor;
         sterm[1] *= bfactor;

         cdata[0] += sterm[0];  // q
         cdata[1] += aterm[0];  // a
         cdata[2] += sterm[1];  // g

         break;

      case gXqqFS: 

         sterm[0] *= bfactor;
         aterm[0] *= bfactor;
         sterm[1] *= bfactor;

         cdata[3] += sterm[0];  // q
         cdata[4] += aterm[0];  // a
         cdata[7] += sterm[0];  // q
         cdata[8] += aterm[0];  // a
         cdata[2] += sterm[1];  // g

         break;

      case NoDISSubProcess:
      default:
         errf(-1,"mel_index not defined (b)");
   }

   // -------------------------------------------------------------------------

#if 0
   if (global -> subtraction_type == 0) {
      if (   global->tstore == global->tcatch 
          && global->lstore != NO_LIMIT
         ) {
         switch (global->tstore) {
            case 1:
               if (!(   global->istore == global->icatch
                     && global->jstore == global->jcatch))
                  sum=0.;
               break;
            case 2:
               if (!(global->istore == global->icatch))
                  sum=0.;
               break;
            case 3:
               if (!(   global->istore == global->icatch
                     && global->jstore == global->jcatch))
                  sum=0.;
              break;
         }
      } else
      if (   global->tcatch != 0
          && global->lstore != NO_LIMIT
         )
         sum=0.;
   } else {
      if (global->tstore == global->tcatch) {
         switch (global->tstore) {
            case 1:
               if (!(   global->istore == global->icatch
                     && global->jstore == global->jcatch))
                  sum=0.;
               break;
            case 2:
               if (!(global->istore == global->icatch))
                  sum=0.;
               break;
            case 3:
               if (!(   global->istore == global->icatch
                     && global->jstore == global->jcatch))
                  sum=0.;
              break;
         }
      } else
      if (global->tcatch != 0)
         sum=0.;
   }
#endif
}

// ----------------------------------------------------------------------------

void DISMatrixElement::EvaluateLorentzAndColour(Event* e, double* sterm) {

   int limit1Loc = -1;
   int limit2Loc = -1;
   LimitType limitTypeLoc = e -> GetLimitType();

   if (limitTypeLoc != NO_LIMIT) {

      limit1Loc    = e->GetLimit1();
      limit2Loc    = e->GetLimit2();

      if (   limit1Loc==limit2Loc
          || limit1Loc < 0
         ) {
         global -> To_error() 
            << FS("limittype, limit1, limit2 %d", (int) limitTypeLoc)
            << FS(" %d", limit1Loc)
            << FS(" %d\n", limit2Loc);
         errf(-1,"DISMatrixElement::Evaluate: limit1, limit2");
      }
   }

   // -------------------------------------------------------------------------

   // evaluate matrix elements

   double lorentz;

   switch (limitTypeLoc) {
      
      case NO_LIMIT:

         // matrix element without limit

         // fill Lorentz structure and prepare for sum
         switch (mel_index) {
            case qXq: 
            case qXqV: 
               FillLorentzLqXq();
               break;
            case qXqg: 
            case qXqgRS: 
            case qXqgRSNF: 
               FillLorentzLqXqg();
               break;
            case gXqq: 
            case gXqqRS: 
            case gXqqRSNF: 
            case gXqqFSNF: 
               FillLorentzLgXqq();
               break;
            case qXqgg: 
            case testDIS:
               FillLorentzLqXqgg();
               break;
            case gXqqg: 
               FillLorentzLgXqqg();
               break;
            case qXqqq: 
               FillLorentzLqXqqq();
               break;
            case qXqgV: 
            case gXqqV: 
               // if possible, re-use calculated value...
               if ( lastVirtual != e->GetId() ) {
                  FillLorentzLpXppV(e, virtq, virtg);
                  virtqSave = virtq;
                  virtgSave = virtg;
                  lastVirtual = e -> GetId();
               } else {
                  virtq = virtqSave;
                  virtg = virtgSave;
               }
               break;
            case NoDISSubProcess:
            default:
               errf(-1,"DISMatrixElement::Evaluate: matrix element not known");
         }
      
         switch (mel_index) {
   
            case qXq: 
               sterm[0] = term[0][0];
               break;
            case qXqV: 
               sterm[0] = term[0][0];
               sterm[0] *=   global -> disaster -> qcd -> cf
                           * (-8-PiSquared/3)/(TwoPi*FourPi);
               break;
            case qXqg: 
            case gXqq: 
               sterm[0] = term[0][0] + 2*term[1][0] + term[1][1];
               break;
            case qXqgRS: 
            case gXqqRS: 
               sterm[0] = term[0][0] + 2*term[1][0] + term[1][1];
               sterm[0] *=   11./6. * global -> disaster -> qcd -> ncolour
                           / (TwoPi * FourPi);
               break;
            case qXqgRSNF: 
            case gXqqRSNF: 
               sterm[0] = term[0][0] + 2*term[1][0] + term[1][1];
               sterm[0] *=   -1./3.
                           / (TwoPi * FourPi);
               break;      
            case gXqqFSNF: 
               sterm[0] = term[0][0] + 2*term[1][0] + term[1][1];
               sterm[0] *= (-1) * global 
                                     -> disaster 
                                     -> qcd 
                                     -> PSplitCounterTermDeltaNf(Q_G_from_G)
                                / (TwoPi * FourPi);
               //sterm[0] /= (-1./3.) / TwoPi;
               break;      
            case qXqgg: 
            case gXqqg: 
               sterm[0] =
                    term[0][0] +   term[2][2] +   term[4][4]      // I*I
                  + term[1][1] +   term[3][3] +   term[5][5]      // II*II
                  + 2 * (  term[2][0] + term[4][0] + term[4][2]   // I*I
                         + term[3][1] + term[5][1] + term[5][3]   // II*II
                         + term[1][0] + term[3][0] + term[5][0]   // I*II
                         + term[2][1] + term[3][2] + term[5][2]   // I*II
                         + term[4][1] + term[4][3] + term[5][4]   // I*II
                        )
               ;
               sterm[1] =
                    2 * (
                           term[1][0] + term[3][0] + term[5][0]   // I*II
                         + term[2][1] + term[3][2] + term[5][2]
                         + term[4][1] + term[4][3] + term[5][4]
                         + term[6][0] + term[7][0]               // I*III
                         + term[6][2] + term[7][2] 
                         + term[6][4] + term[7][4]
                         - (      term[6][6]                // -2 III*III
                            + 2 * term[7][6]  
                            +     term[7][7]
                               )
                         - (  term[6][1] + term[7][1]         // - II*III
                            + term[6][3] + term[7][3] 
                            + term[6][5] + term[7][5]
                           )
                        )
               ;
              break;
            case testDIS: 
               sterm[0] =
                      term[0][0]
               ;
               break;
            case qXqqq: 
               sterm[0] =
                      term[0][0] + term[1][1]
                  +   term[2][2] + term[3][3]
                  +   term[4][4] + term[5][5]
                  +   term[6][6] + term[7][7]
                  + 2 * (  term[1][0]
                         + term[3][2]
                         + term[5][4]
                         + term[7][6]
                         + term[2][0] + term[2][1] + term[3][0] + term[3][1]
                         + term[6][4] + term[6][5] + term[7][4] + term[7][5] 
                        )
               ;
               sterm[1] =
                  2 * (  term[4][0] + term[4][1] + term[5][0] + term[5][1]
                       + term[4][2] + term[4][3] + term[5][2] + term[5][3]
                       + term[6][0] + term[6][1] + term[7][0] + term[7][1]
                       + term[6][2] + term[6][3] + term[7][2] + term[7][3]
                      )
               ;
               sterm[2] =
                    term[0][0] + 2*term[1][0] + term[1][1]
               ;
               sterm[3] =
                    term[2][2] + 2*term[3][2] + term[3][3]
               ;
               sterm[4] =
                    2 * (term[2][0] + term[2][1] + term[3][0] + term[3][1])
               ;
               break;
            case qXqgV: 
               sterm[0] = virtq;
               break;
            case gXqqV: 
               sterm[0] = virtg;
               break;
            case NoDISSubProcess:
            default:
               errf(-1,"mel_index not defined (7)");
         }
         break; 
   
      case SOFT:
      case COLLINEAR:
      case SOFT_AND_COLLINEAR:

         // matrix element with limit
   
         double prefactor;
   
         switch (mel_index) {
   
            case qXq: 
            case qXqg: 
            case gXqq: 
            case qXqgg: 
            case gXqqg: 
            case qXqqq: 
            case testDIS:
   
               // additional factor, to be multipled by the limit cross sections
               // calculated via Mathematica
               switch (limitTypeLoc) {
                  case SOFT:
                     prefactor=1.;
                     break;
                  case COLLINEAR:
                     switch (e->GetP0ref()) {
                        case FS_REFERENCE:
                           prefactor=0.25/(eh*eh*lambda*(1.-lambda));
                           break;
                        case IS_REFERENCE:
                           prefactor=0.25*eta*eta/(eh*eh*(1.-eta));
                           break;
                        default:
                           prefactor=0.; // :_TAWM_:
                     }
                     break;
                  case SOFT_AND_COLLINEAR:
                     prefactor=0.25;
                     break;
                  default:
                     errf(-1, "DISMatrixElement::EvaluateLorentzAndColour:"
                              " impossible");
                     prefactor=0.; // :_TAWM_:
               }
               break;
   
            default:
            errf(-1, "cannot construct subtraction term");
            prefactor=0.;  // :_TAWM_:
         }
   
         switch (mel_index) {
   
            case qXq: 
               errf(-1, "cannot get subtraction for q -> q term!!");
               break;
            case qXqg: 
               sterm[0] = GetLorentzSubtraction(
                            11, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[0] *= prefactor;
               break;
            case gXqq: 
               sterm[0] = GetLorentzSubtraction(
                            12, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[0] *= prefactor;
               break;
            case qXqgg: 
               sterm[0] = GetLorentzSubtraction(
                            1, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[1] = GetLorentzSubtraction(
                            2, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[0] *= prefactor;
               sterm[1] *= prefactor;
               break;
            case testDIS:
               sterm[0] = GetLorentzSubtraction(
                            100, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[1] = GetLorentzSubtraction(
                            101, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[0] *= prefactor;
               sterm[1] *= prefactor;
               break;
            case gXqqg: 
               sterm[0] = GetLorentzSubtraction(
                            3, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[1] = GetLorentzSubtraction(
                            4, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[0] *= prefactor;
               sterm[1] *= prefactor;
              break;
            case qXqqq: 
               sterm[0] = GetLorentzSubtraction(
                            8, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[1] = GetLorentzSubtraction(
                            9, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[2] = GetLorentzSubtraction(
                            5, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[3] = GetLorentzSubtraction(
                            6, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[4] = GetLorentzSubtraction(
                            7, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[0] *= prefactor;
               sterm[1] *= prefactor;
               sterm[2] *= prefactor;
               sterm[3] *= prefactor;
               sterm[4] *= prefactor;
               break;
            case NoDISSubProcess:
            default:
               errf(-1,"mel_index not defined (8)");
         }
         break;
   
      case SOFT_ADDED:
      case COLLINEAR_ADDED:
      case SOFT_AND_COLLINEAR_ADDED:

         switch (mel_index) {
            case qXqgA: 
               sterm[0] = EvaluateLorentzAddedSubtraction(
                            e, 11, limitTypeLoc, limit1Loc, limit2Loc
                         );
               break;
            case gXqqA: 
               sterm[0] = EvaluateLorentzAddedSubtraction(
                            e, 12, limitTypeLoc, limit1Loc, limit2Loc
                         );
               break;
            case qXqggA: 
               sterm[0] = EvaluateLorentzAddedSubtraction(
                            e, 1, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[1] = EvaluateLorentzAddedSubtraction(
                            e, 2, limitTypeLoc, limit1Loc, limit2Loc
                         );
               break;
            case gXqqgA: 
               sterm[0] = EvaluateLorentzAddedSubtraction(
                            e, 3, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[1] = EvaluateLorentzAddedSubtraction(
                            e, 4, limitTypeLoc, limit1Loc, limit2Loc
                         );
              break;
            case qXqqqA: 
               sterm[0] = EvaluateLorentzAddedSubtraction(
                            e, 8, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[1] = EvaluateLorentzAddedSubtraction(
                            e, 9, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[2] = EvaluateLorentzAddedSubtraction(
                            e, 5, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[3] = EvaluateLorentzAddedSubtraction(
                            e, 6, limitTypeLoc, limit1Loc, limit2Loc
                         );
               sterm[4] = EvaluateLorentzAddedSubtraction(
                            e, 7, limitTypeLoc, limit1Loc, limit2Loc
                         );
               break;

            case qXqFS: 
               lorentz = GetBornTerm(qXq);
               sterm[0] = EvaluateLorentzCollinearInitial(
                             e, Q_q_from_q, limitTypeLoc, lorentz
                         );
               sterm[1] = EvaluateLorentzCollinearInitial(
                             e, Q_q_from_G, limitTypeLoc, lorentz
                         );
               break;

            case qXqgFS: 
               lorentz = GetBornTerm(qXqg);
               sterm[0] = EvaluateLorentzCollinearInitial(
                             e, Q_q_from_q, limitTypeLoc, lorentz
                         );
               sterm[1] = EvaluateLorentzCollinearInitial(
                             e, Q_q_from_G, limitTypeLoc, lorentz
                         );
               break;

            case gXqqFS: 
               lorentz = GetBornTerm(gXqq);
               sterm[0] = EvaluateLorentzCollinearInitial(
                             e, Q_G_from_q, limitTypeLoc, lorentz
                         );
               sterm[1] = EvaluateLorentzCollinearInitial(
                             e, Q_G_from_G, limitTypeLoc, lorentz
                         );
               break;

            case NoDISSubProcess:
            default:
               errf(-1,"mel_index not defined (8)");
         }
         break;
   
      default:
         errf(-1, "subtraction not defined...");

   } // of switch (limitTypeLoc)

   // -------------------------------------------------------------------------

   // add colour and symmetry factors

   switch (mel_index) {

      case qXq: 
      case qXqV: 

         sterm[0] *= global -> disaster -> ncolour;
         sterm[0] *= global -> disaster -> quarkfactor;

         break;

      case qXqFS: 

         sterm[0] *= global -> disaster -> ncolour;
         sterm[0] *= global -> disaster -> quarkfactor;

         sterm[1] *= global -> disaster -> ncolour;
         sterm[1] *= global -> disaster -> quarkfactor;

         break;

      case qXqg: 
      case qXqgA: 
      case qXqgRS: 
      case qXqgRSNF: 

         sterm[0] *= global -> disaster -> ncolour * global -> disaster -> cf;
         sterm[0] *= global -> disaster -> quarkfactor;

         break;

      case qXqgFS: 

         sterm[0] *= global -> disaster -> ncolour * global -> disaster -> cf;
         sterm[0] *= global -> disaster -> quarkfactor;

         sterm[1] *= global -> disaster -> ncolour * global -> disaster -> cf;
         sterm[1] *= global -> disaster -> quarkfactor;

         break;

      case gXqq: 
      case gXqqA: 
      case gXqqRS: 
      case gXqqRSNF: 
      case gXqqFSNF: 

         sterm[0] *= global -> disaster -> ncolour * global -> disaster -> cf;
         sterm[0] *= global -> disaster -> gluonfactor;

         break;

      case gXqqFS: 

         sterm[0] *= global -> disaster -> ncolour * global -> disaster -> cf;
         sterm[0] *= global -> disaster -> gluonfactor;

         sterm[1] *= global -> disaster -> ncolour * global -> disaster -> cf;
         sterm[1] *= global -> disaster -> gluonfactor;

         break;

      case qXqgg: 
      case gXqqg: 
      case qXqggA: 
      case gXqqgA: 
      case testDIS:

         sterm[0] *=   global -> disaster -> ncolour 
                     * global -> disaster -> cf
                     * global -> disaster -> cf;
         sterm[1] *=   - 0.5 
                     * global -> disaster -> ncolour 
                     * global -> disaster -> ncolour 
                     * global -> disaster -> cf;

         if (   mel_index == qXqgg 
             || mel_index == qXqggA
             || mel_index == testDIS) {

            sterm[0] *= global -> disaster -> quarkfactor;
            sterm[1] *= global -> disaster -> quarkfactor;

            sterm[0] *= 0.5;
            sterm[1] *= 0.5;

         } else {

            sterm[0] *= global -> disaster -> gluonfactor;
            sterm[1] *= global -> disaster -> gluonfactor;

         }

         break;

      case qXqqq: 
      case qXqqqA: 
            
         sterm[0] *= 0.5 * global -> disaster -> ncolour 
                         * global -> disaster -> cf;
         sterm[1] *=  global -> disaster -> ncolour 
                    * global -> disaster -> cf
                    *(        global -> disaster -> cf 
                      - 0.5 * global -> disaster -> ncolour);
         sterm[2] *= 0.5 * global -> disaster -> ncolour 
                         * global -> disaster -> cf;
         sterm[3] *= 0.5 * global -> disaster -> ncolour 
                         * global -> disaster -> cf;
         sterm[4] *= 0.5 * global -> disaster -> ncolour 
                         * global -> disaster -> cf;
          
         sterm[0] *= global -> disaster -> quarkfactor;
         sterm[1] *= global -> disaster -> quarkfactor;
         sterm[2] *= global -> disaster -> quarkfactor;
         sterm[3] *= global -> disaster -> quarkfactor;
         sterm[4] *= global -> disaster -> quarkfactor;

         sterm[0] *= 0.5;
         sterm[1] *= 0.5;

         break;

      case qXqgV: 
      case gXqqV: 

         sterm[0] *= 1.;

         break;

      case NoDISSubProcess:

      default:
         errf(-1,"mel_index not defined (9)");
   }

}

// ----------------------------------------------------------------------------

void DISMatrixElement::FillLorentzLpXppV(Event* e, double& xq, double& xg) {

   double par[11];

   par[0]=si0;
   par[1]=si1;
   par[2]=s01;
   par[3]=s0k;
   par[4]=e->xB;
   par[5]=e->y;
   par[6]=e->Q2;
   par[7]=e->xi;
//   par[8]=e->GetEventWithAlphaS()->GetAlphaS(); // a_s == 1 in virtual.f!
   par[8]=1234567890;
   par[9]=e->Q2;
   par[10]=5.;

   evalvirtC(
      par,
//      0, 1, 3,
//      0, 1, 1,
      0, 1, 15,
      xq, xg
   );

   xq = xq/FourPi;
   xg = xg/FourPi;
}

// ----------------------------------------------------------------------------

// Calculate the matrix elements for the added subtraction terms.
// The factors from the splitting functions, and from the other 
// contributions are appended.
// Assumes that always ifix > jfix for the collinear limit

double DISMatrixElement::EvaluateLorentzAddedSubtraction(
                             Event* e, int list, 
                             LimitType lt, int ifix, int jfix) {

   // first get the coefficients and evaluate the Born term
   GetLorentzAddedSubtraction(
      list, 
      (lt == SOFT_AND_COLLINEAR_ADDED) ? COLLINEAR_ADDED : lt,  
      ifix, jfix);

   register int i;

   double xi_hard_loc;
   double u_loc;
   double kappa;

   double xB_loc = e -> xB;

   double sum = 0.;

   double term_u_independent, term_u_dependent;
  
   double tw, th;
   double safeSpence, zw;

   QCD* qcd = global -> disaster -> qcd;

   switch (lt) {

      case SOFT_ADDED:

         switch (e -> GetFrame()) {
            case pCMS:
               tw = log( e -> xi / xB_loc - 1.);
               break;
            case hCMS:
               tw = log( e -> W2 / e -> Q2 );
               break;
            default:
               errf(-1, "DISMatrixElement::EvaluateLorentzAddedSubtraction "
                        "frame not known");
               tw = 0; // :_TAWM_:
         }

         // :_MOD_: introduce a safe spence function globally!
         for (i = 0; i < nCoefficient; ++i) {
             if (softContributes[i]) {
                if (fabs(1. - softAngle[i]) < 1.e-10) {
                   zw = - log(softAngle[i]);
                   // first terms in power series
                   safeSpence = zw * (1. - 0.25 * zw); 
                } else {
                   safeSpence = RealSpence(1 - softAngle[i]);
                }
           
                sum +=   softCoefficient[i] 
                       * (  (2. * tw + log(softAngle[i])) * log(softAngle[i])
                          +  2. * safeSpence
                         );
             }
         }
         break;

      case COLLINEAR_ADDED:

         if (collinearContributes) {

            if (jfix >= 0) {

               // final state collinear singularity
               th = 2. * log(2. * e -> GetMapping(-3) -> GetEnergy() / e -> Q);
               sum +=   collinearCoefficient
                      * 0.25 
                      *(  (0.5*th*th-2./3.*PiSquared)
                           *(  qcd->QSplitAt0(splittingFunction)
                             + qcd->QSplitAt1(splittingFunction))
                        + 2 * th * qcd->QSplitRegularIntEps0(splittingFunction)
                        - 2 *      qcd->QSplitRegularIntEps1(splittingFunction)
                       );

            } else {

               // initial state collinear singularity
               xi_hard_loc = e -> xi_hard;
               u_loc       = e -> u;

               switch (e -> GetFrame()) {

                  case pCMS:

                     kappa = xB_loc / xi_hard_loc;

                     term_u_dependent = -0.5 * 1/u_loc * (
                       log(  kappa * u_loc * (1-kappa*u_loc)
                           / ((1-u_loc)*(1-u_loc)))
                       * (  qcd->QSplitHalfRegularEps0(splittingFunction, u_loc)
                          + 1/(1-u_loc) * qcd->QSplitAt1(splittingFunction))
                      + qcd->QSplitHalfRegularEps1(splittingFunction, u_loc)
                     );

                     term_u_independent = 0;

                     break;

                  case hCMS:

                     term_u_dependent = -0.5 * 1/u_loc *(
                       log(  u_loc * u_loc * (1-xB_loc)
                           / ((1-u_loc) * (1-u_loc) * xB_loc))
                       * (  qcd->QSplitHalfRegularEps0(splittingFunction, u_loc)
                          + 1/(1-u_loc) * qcd->QSplitAt1(splittingFunction))
                      + qcd->QSplitHalfRegularEps1(splittingFunction, u_loc)
                     );

                     term_u_independent = 0;

                     break;

                  default:

                     errf(-1, "DISMatrixElement::EvaluateLorentzAdded"
                              "Subtraction: frame not known");
                     term_u_dependent   = 0; // :_TAWM_:
                     term_u_independent = 0;
               }

               sum +=   collinearCoefficient
                      * (  e->ujacobian * term_u_dependent 
                         +                term_u_independent
                        );
            }
         }
         break;

      case SOFT_AND_COLLINEAR_ADDED:

         if (collinearContributes) {
            if (jfix >= 0) {
            
               // final state soft & collinear singularity
               errf(-1,"DISMatrixElement::EvaluateLorentzAddedSubtraction:"
                       " SCA illegal");
            } else {

               // initial state soft & collinear singularity

               xi_hard_loc = e -> xi_hard;
               u_loc       = e -> u;
               switch (e -> GetFrame()) {
                  case pCMS:
                     kappa = xB_loc / xi_hard_loc;
                     term_u_dependent = - (-0.5)
                            *  log(  kappa * u_loc * (1-kappa*u_loc)
                                   / ((1-u_loc)*(1-u_loc)))              
                            / (1-u_loc) 
                            * qcd->QSplitAt1(splittingFunction);
                     term_u_independent = 0.5 * (
                              RealSpence(-kappa/(1-kappa)*(1-xi_hard_loc))
                            + RealSpence(1-xi_hard_loc)
                            - log(1-xi_hard_loc)*log(  kappa*(1-kappa)
                                                     / (1-xi_hard_loc) )
                            + 0.25 * pow(log(kappa*(1-kappa)), 2)
                          ) * qcd->QSplitAt1(splittingFunction);
                     break;

                  case hCMS:
                     term_u_dependent = - (-0.5)
                            * log(  u_loc * u_loc * (1-xB_loc)
                                  / ((1-u_loc) * (1-u_loc) * xB_loc))
                            / (1-u_loc) 
                            * qcd->QSplitAt1(splittingFunction);
                     term_u_independent = 0.5 * (
                              2 * RealSpence(1-xB_loc)
                            + 0.25 * pow(log(xB_loc*(1-xB_loc)), 2)
                          ) * qcd->QSplitAt1(splittingFunction);
                     break;
                  default:
                     errf(-1, "DISMatrixElement::EvaluateLorentzAdded"
                              "Subtraction: frame not known");
                     term_u_dependent   = 0; // :_TAWM_:
                     term_u_independent = 0;
               }
               sum +=   collinearCoefficient
                      * (  e->ujacobian * term_u_dependent 
                         +                term_u_independent
                        );
            }
         }
         break;
      default:
         errf(-1, "DISMatrixElement::EvaluateLorentzAddedSubtraction:"
                  " illegal limit type");
   }

   return sum * bornTerm / (TwoPi * FourPi);
}

// ----------------------------------------------------------------------------

// Calculate the matrix elements for the initial state collinear terms.
// The factors from the splitting functions, and from the other 
// contributions are appended.
// :_MOD_: modify as the AddedSubtraction-term! (function before)


double DISMatrixElement::EvaluateLorentzCollinearInitial(
                             Event* e, SplittingFunction sf, 
                             LimitType lt, double lorentz) {

   double xi_hard_loc;
   double u_loc;

   double sum = 0.;

   double term_u_independent, term_u_dependent;
   double term_L; 
  
   QCD *qcd = global -> disaster -> qcd;

   xi_hard_loc = e -> xi_hard;

   switch (lt) {

      case SOFT_ADDED:

         sum += (-1) * qcd -> PSplitCounterTermDelta(sf, xi_hard_loc);
         break;

      case COLLINEAR_ADDED:

         // initial state collinear singularity
         u_loc = e -> u;
         term_u_dependent
            = (-1) * 1/u_loc
                   * (  qcd -> PSplitCounterTermSubtracted(sf, u_loc)
                      + qcd -> PSplitCounterTermRegular(sf, u_loc));
         term_u_independent = 0;
         sum +=   e->ujacobian * term_u_dependent 
                +                term_u_independent;
         break;

      case SOFT_AND_COLLINEAR_ADDED:

         // initial state soft & collinear singularity
         u_loc = e -> u;
         term_L = (-1)* (- qcd -> PSplitCounterTermSubtracted(sf, u_loc));
         term_L *= e -> ujacobian;
         sum += term_L;
         break;

      default:
         errf(-1, "DISMatrixElement::EvaluateLorentzCollinearInitial:"
                  " illegal limit type");
   }

   return sum * lorentz / (TwoPi * FourPi);
}

// ----------------------------------------------------------------------------
//
// --> DIS process: NoSubtraction Test
//
// ----------------------------------------------------------------------------

DISProcessNoSubtraction_Test::DISProcessNoSubtraction_Test() {
}

// ----------------------------------------------------------------------------

DISProcessNoSubtraction_Test::~DISProcessNoSubtraction_Test() {
}

// ----------------------------------------------------------------------------

void DISProcessNoSubtraction_Test::Define(Global *globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessNoSubtraction_Test::AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 1;

   // *** fill in number of components here
   CreateNME(1);

   // *** fill in number of matrix elements per component here
   n_me             [0] = 1;
   to_be_subtracted [0] = NoSubtraction;

   CreateOtherArrays();

   // *** fill in matrix elements here
   me_list [0][0] = new DISMatrixElement((DISSubProcess)global->mel_index); 

   AssignInformationToMatrixElements();

   // frame to be used for phase space
   if (GetIn() >= 1 && GetMaxOut() <= 2)
      frame=hCMS;
   else
      frame=pCMS;

   global -> numberOfFinalStatePartonsInBornTerm = GetMinOut();
}

// ----------------------------------------------------------------------------
//
// --> DIS process: Subtraction Test
//
// ----------------------------------------------------------------------------

DISProcessSubtraction_Test::DISProcessSubtraction_Test() {
}

// ----------------------------------------------------------------------------

DISProcessSubtraction_Test::~DISProcessSubtraction_Test() {
}

// ----------------------------------------------------------------------------

void DISProcessSubtraction_Test::Define(Global *globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessSubtraction_Test::AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 2;

   // *** fill in number of components here
   CreateNME(1);

   // *** fill in number of matrix elements per component here
   n_me             [0] = 1;
   to_be_subtracted [0] = Subtraction;

   CreateOtherArrays();

   // *** fill in matrix elements here
   me_list [0][0] = new DISMatrixElement((DISSubProcess)global->mel_index); 

   AssignInformationToMatrixElements();

   // frame to be used for phase space
   switch (global->mel_index) {
      case qXqg:
      case gXqq:
         frame = hCMS;
         break;
      case qXqgg:
      case gXqqg:
      case qXqqq:
         frame = pCMS;
         break;
      case qXq:
      case qXqgV:
      case gXqqV:
      case qXqV:
      case qXqgA:
      case gXqqA:
      case qXqggA:
      case gXqqgA:
      case qXqqqA:
      default:
         errf(-1, "(frame) no such mel_index");
   }

   global -> numberOfFinalStatePartonsInBornTerm = GetMinOut();
}

// ----------------------------------------------------------------------------
//
// --> DIS process: AddedSubtraction Test
//
// ----------------------------------------------------------------------------

DISProcessAddedSubtraction_Test::DISProcessAddedSubtraction_Test() {
}

// ----------------------------------------------------------------------------

DISProcessAddedSubtraction_Test::~DISProcessAddedSubtraction_Test() {
}

// ----------------------------------------------------------------------------

void DISProcessAddedSubtraction_Test::Define(Global *globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessAddedSubtraction_Test::AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 2;

   // *** fill in number of components here
   CreateNME(1);

   // *** fill in number of matrix elements per component here
   n_me             [0] = 1;
   to_be_subtracted [0] = AddedSubtraction;

   CreateOtherArrays();

   // *** fill in matrix elements here
   me_list [0][0] = new DISMatrixElement((DISSubProcess)global->mel_index); 

   AssignInformationToMatrixElements();

   // frame to be used for phase space
   switch (global->mel_index) {
      case qXqgA:
      case gXqqA:
         frame = hCMS;
         break;
      case qXqggA:
      case gXqqgA:
      case qXqqqA:
         frame = pCMS;
         break;
      case qXq:
      case qXqgV:
      case gXqqV:
      case qXqV:
      case qXqg:
      case gXqq:
      case qXqgg:
      case gXqqg:
      case qXqqq:
      default:
         errf(-1, "(frame) no such mel_index");
   }

   global -> numberOfFinalStatePartonsInBornTerm = GetMinOut();
}

// ----------------------------------------------------------------------------
//
// --> DIS process: CollinearInitial Test
//
// ----------------------------------------------------------------------------

DISProcessCollinearInitial_Test::DISProcessCollinearInitial_Test() {
}

// ----------------------------------------------------------------------------

DISProcessCollinearInitial_Test::~DISProcessCollinearInitial_Test() {
}

// ----------------------------------------------------------------------------

void DISProcessCollinearInitial_Test::Define(Global *globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessCollinearInitial_Test::AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 2;

   // *** fill in number of components here
   CreateNME(1);

   // *** fill in number of matrix elements per component here
   n_me             [0] = 1;
   to_be_subtracted [0] = CollinearInitial;

   CreateOtherArrays();

   // *** fill in matrix elements here
   me_list [0][0] = new DISMatrixElement((DISSubProcess)global->mel_index); 

   AssignInformationToMatrixElements();

   // frame to be used for phase space
   switch (global->mel_index) {
      case qXqFS:
         frame = hCMS;
         break;
      case qXqgFS:
      case gXqqFS:
         frame = pCMS;
         break;
      default:
         errf(-1, "(frame) no such mel_index");
   }

   global -> numberOfFinalStatePartonsInBornTerm = GetMinOut();
}

// ----------------------------------------------------------------------------
//
// --> DIS process: LO
//
// ----------------------------------------------------------------------------

DISProcessLO::DISProcessLO() {
}

// ----------------------------------------------------------------------------

DISProcessLO::~DISProcessLO() {
}

// ----------------------------------------------------------------------------

void DISProcessLO::Define(Global *globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessLO::AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 1;

   // *** fill in number of components here
   CreateNME(1);

   // *** fill in number of matrix elements per component here
   switch (global->numberOfFinalStatePartonsInBornTerm) {
      case 1:
         n_me             [0] = 1;
         to_be_subtracted [0] = NoSubtraction;
         break;
      case 2:
         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;
         Get_vtp(0) -> Set_mapParameter(
                          1, 
                          Get_vtp(0) -> Get_mapParameter(3)
                       );
         Get_vtp(0) -> Set_mapParameter(
                          13, 
                          Get_vtp(0) -> Get_mapParameter(10)
                       );
         Get_vtp(0) -> Set_mapParameter(
                          14, 
                          Get_vtp(0) -> Get_mapParameter(11)
                       );
         Get_vtp(0) -> Set_mapParameter(
                          15, 
                          Get_vtp(0) -> Get_mapParameter(12)
                       );
         break;
      case 3:
         n_me             [0] = 3;
         to_be_subtracted [0] = NoSubtraction;
         break;
      default:
         errf(-1, "DISProcessLO::AssignMatrixElements: process not known");
   }

   CreateOtherArrays();

   // *** fill in matrix elements and frames here
   switch (global->numberOfFinalStatePartonsInBornTerm) {
      case 1:
         me_list [0][0] 
            = new DISMatrixElement(qXq); 
         frame = hCMS;
         break;
      case 2:
         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 
         frame = pCMS;
         break;
      case 3:
         me_list [0][0] 
            = new DISMatrixElement(qXqgg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqqg); 
         me_list [0][2] 
            = new DISMatrixElement(qXqqq); 
         frame = pCMS;
         break;
   }

   AssignInformationToMatrixElements();
}

// ----------------------------------------------------------------------------
//
// --> DIS process: NLO
//
// ----------------------------------------------------------------------------

DISProcessNLO::DISProcessNLO() {
}

// ----------------------------------------------------------------------------

DISProcessNLO::~DISProcessNLO() {
}

// ----------------------------------------------------------------------------

void DISProcessNLO::Define(Global *globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessNLO::AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 2;

   // *** fill in number of components in CreateNME()

   // *** fill in number of matrix elements per component here
   switch (global->numberOfFinalStatePartonsInBornTerm) {

      case 1:

         CreateNME(5);

         n_me             [0] = 1;
         to_be_subtracted [0] = NoSubtraction;

         n_me             [1] = 1;
         to_be_subtracted [1] = NoSubtraction;
         Copy_vtp(1, 0);

         n_me             [2] = 2;
         to_be_subtracted [2] = AddedSubtraction;

         n_me             [3] = 2;
         to_be_subtracted [3] = Subtraction;

         n_me             [4] = 1;
         to_be_subtracted [4] = CollinearInitial;
         Copy_vtp(4, 2);

         break;

      case 2:

         CreateNME(8);

         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;
         if (global -> disaster -> Get_shifted_parton_variables() == 1)
            Get_vtp(0) -> Set_offset(3, 3);
         Get_vtp(0) -> Set_mapParameter(
                          1, 
                          Get_vtp(3) -> Get_mapParameter(3)
                       );
         Get_vtp(0) -> Set_mapParameter(
                          13, 
                          Get_vtp(3) -> Get_mapParameter(10)
                       );
         Get_vtp(0) -> Set_mapParameter(
                          14, 
                          Get_vtp(3) -> Get_mapParameter(11)
                       );
         Get_vtp(0) -> Set_mapParameter(
                          15, 
                          Get_vtp(3) -> Get_mapParameter(12)
                       );

         n_me             [1] = 2;
         to_be_subtracted [1] = NoSubtraction;
         Copy_vtp(1, 0);

         n_me             [2] = 3;
         to_be_subtracted [2] = AddedSubtraction;
         if (global -> disaster -> Get_shifted_parton_variables() == 1)
            Get_vtp(2) -> Set_offset(3, 2);
         Get_vtp(2) -> Set_mapParameter(
                          1, 
                          Get_vtp(3) -> Get_mapParameter(3)
                       );
         Get_vtp(2) -> Set_mapParameter(
                          13, 
                          Get_vtp(3) -> Get_mapParameter(10)
                       );
         Get_vtp(2) -> Set_mapParameter(
                          14, 
                          Get_vtp(3) -> Get_mapParameter(11)
                       );
         Get_vtp(2) -> Set_mapParameter(
                          15, 
                          Get_vtp(3) -> Get_mapParameter(12)
                       );

         n_me             [3] = 3;
         to_be_subtracted [3] = Subtraction;

         n_me             [4] = 2;
         to_be_subtracted [4] = CollinearInitial;
         Copy_vtp(4, 2);

         n_me             [5] = 1;
         to_be_subtracted [5] = NoSubtraction;
         Copy_vtp(5, 0);

         n_me             [6] = 2;
         to_be_subtracted [6] = NoSubtraction;
         Copy_vtp(6, 0);

         n_me             [7] = 2;
         to_be_subtracted [7] = NoSubtraction;
         Copy_vtp(7, 0);

         break;

      default:
         errf(-1, "DISProcessNLO::AssignMatrixElements: process not known");
   }

   CreateOtherArrays();

   // *** fill in matrix elements and frames here
   switch (global->numberOfFinalStatePartonsInBornTerm) {

      case 1:

         me_list [0][0] 
            = new DISMatrixElement(qXq); 

         me_list [1][0] 
            = new DISMatrixElement(qXqV); 

         me_list [2][0] 
            = new DISMatrixElement(qXqgA); 
         me_list [2][1] 
            = new DISMatrixElement(gXqqA); 

         me_list [3][0] 
            = new DISMatrixElement(qXqg); 
         me_list [3][1] 
            = new DISMatrixElement(gXqq); 

         me_list [4][0]
            = new DISMatrixElement(qXqFS); 

         frame = hCMS;

         break;

      case 2:

         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 

         me_list [1][0] 
            = new DISMatrixElement(qXqgV); 
         me_list [1][1] 
            = new DISMatrixElement(gXqqV); 

         me_list [2][0] 
            = new DISMatrixElement(qXqggA); 
         me_list [2][1] 
            = new DISMatrixElement(gXqqgA); 
         me_list [2][2] 
            = new DISMatrixElement(qXqqqA); 

         me_list [3][0] 
            = new DISMatrixElement(qXqgg); 
         me_list [3][1] 
            = new DISMatrixElement(gXqqg); 
         me_list [3][2] 
            = new DISMatrixElement(qXqqq); 

         me_list [4][0] 
            = new DISMatrixElement(qXqgFS); 
         me_list [4][1] 
            = new DISMatrixElement(gXqqFS); 

         me_list [5][0] 
            = new DISMatrixElement(gXqqFSNF); 

         me_list [6][0] 
            = new DISMatrixElement(qXqgRS); 
         me_list [6][1] 
            = new DISMatrixElement(gXqqRS); 

         me_list [7][0] 
            = new DISMatrixElement(qXqgRSNF); 
         me_list [7][1] 
            = new DISMatrixElement(gXqqRSNF); 

         frame = pCMS;

         break;
   }

   AssignInformationToMatrixElements();
}

// ----------------------------------------------------------------------------
//
// --> DIS process: NLO light...
//
// ----------------------------------------------------------------------------

DISProcessNLO_Light::DISProcessNLO_Light() {
}

// ----------------------------------------------------------------------------

DISProcessNLO_Light::~DISProcessNLO_Light() {
}

// ----------------------------------------------------------------------------

void DISProcessNLO_Light::Define(Global *globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessNLO_Light::AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 2;

   // *** fill in number of components in CreateNME()

   // *** fill in number of matrix elements per component here
   switch (global->numberOfFinalStatePartonsInBornTerm) {

      case 1:

         CreateNME(4);

         n_me             [0] = 1;
         to_be_subtracted [0] = NoSubtraction;

         n_me             [1] = 1;
         to_be_subtracted [1] = NoSubtraction;

         n_me             [2] = 2;
         to_be_subtracted [2] = AddedSubtraction;

         n_me             [3] = 2;
         to_be_subtracted [3] = Subtraction;

         break;

      case 2:

         CreateNME(4);

         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;

         n_me             [1] = 2;
         to_be_subtracted [1] = NoSubtraction;

         n_me             [2] = 3;
         to_be_subtracted [2] = AddedSubtraction;

         n_me             [3] = 3;
         to_be_subtracted [3] = Subtraction;

         break;

      default:
         errf(-1, "DISProcessNLO_Light::AssignMatrixElements: "
                  "process not known");
   }

   CreateOtherArrays();

   // *** fill in matrix elements and frames here
   switch (global->numberOfFinalStatePartonsInBornTerm) {

      case 1:

         me_list [0][0] 
            = new DISMatrixElement(qXq); 

         me_list [1][0] 
            = new DISMatrixElement(qXqV); 

         me_list [2][0] 
            = new DISMatrixElement(qXqgA); 
         me_list [2][1] 
            = new DISMatrixElement(gXqqA); 

         me_list [3][0] 
            = new DISMatrixElement(qXqg); 
         me_list [3][1] 
            = new DISMatrixElement(gXqq); 

         frame = hCMS;

         break;

      case 2:

         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 

         me_list [1][0] 
            = new DISMatrixElement(qXqgV); 
         me_list [1][1] 
            = new DISMatrixElement(gXqqV); 

         me_list [2][0] 
            = new DISMatrixElement(qXqggA); 
         me_list [2][1] 
            = new DISMatrixElement(gXqqgA); 
         me_list [2][2] 
            = new DISMatrixElement(qXqqqA); 

         me_list [3][0] 
            = new DISMatrixElement(qXqgg); 
         me_list [3][1] 
            = new DISMatrixElement(gXqqg); 
         me_list [3][2] 
            = new DISMatrixElement(qXqqq); 

         frame = pCMS;

         break;
   }

   AssignInformationToMatrixElements();
}

// ----------------------------------------------------------------------------
//
// --> DIS process: LO test
//
// ----------------------------------------------------------------------------

DISProcessLO_Test :: DISProcessLO_Test() {
}

// ----------------------------------------------------------------------------

DISProcessLO_Test :: ~DISProcessLO_Test() {
}

// ----------------------------------------------------------------------------

void DISProcessLO_Test :: Define(Global* globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessLO_Test :: AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 1;

   // *** fill in number of components here
   CreateNME(1);

   // *** fill in number of matrix elements per component here
   switch (global -> mel_index) {
      case -2:
         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case 1:
         n_me             [0] = 1;
         to_be_subtracted [0] = NoSubtraction;
         global -> numberOfFinalStatePartonsInBornTerm = 1;
         break;
      case 2:
         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case 3:
         n_me             [0] = 3;
         to_be_subtracted [0] = NoSubtraction;
         global -> numberOfFinalStatePartonsInBornTerm = 3;
         break;
      default:
         errf(-1, "DISProcessLO::AssignMatrixElements: process not known");
   }

   CreateOtherArrays();

   // *** fill in matrix elements and frames here
   switch (global -> mel_index) {
      case -2:
         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 
         frame = pCMS;
         break;
      case 1:
         me_list [0][0] 
            = new DISMatrixElement(qXq); 
         frame = hCMS;
         break;
      case 2:
         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 
         frame = hCMS;
         break;
      case 3:
         me_list [0][0] 
            = new DISMatrixElement(qXqgg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqqg); 
         me_list [0][2] 
            = new DISMatrixElement(qXqqq); 
         frame = pCMS;
         break;
   }

   AssignInformationToMatrixElements();
}

// ----------------------------------------------------------------------------
//
// --> DIS process: NLO Test
//
// ----------------------------------------------------------------------------

DISProcessNLO_Test::DISProcessNLO_Test() {
}

// ----------------------------------------------------------------------------

DISProcessNLO_Test::~DISProcessNLO_Test() {
}

// ----------------------------------------------------------------------------

void DISProcessNLO_Test::Define(Global *globalIn) {

   DefineProcess(globalIn);
}

// ----------------------------------------------------------------------------

void DISProcessNLO_Test::AssignMatrixElements() {

   // delete old definition, if any
   DeleteMatrixElements();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 2;

   // *** fill in number of components here
   switch (global->mel_index) {
      case -10:
         CreateNME(1);
         break;
      case -9:
         CreateNME(3);
         break;
      case -8:
         CreateNME(5);
         break;
      case -7:
         CreateNME(2);
         break;
      case -6:
         CreateNME(3);
         break;
      case -5:
         CreateNME(3);
         break;
      default:
         CreateNME(1);
   }

   // *** fill in number of matrix elements per component here
   switch (global->mel_index) {
      case -10:
         n_me             [0] = 2;
         to_be_subtracted [0] = CollinearInitial;
         break;
      case -9:
         n_me             [0] = 1;
         to_be_subtracted [0] = NoSubtraction;
         n_me             [1] = 1;
         to_be_subtracted [1] = AddedSubtraction;
         n_me             [2] = 1;
         to_be_subtracted [2] = Subtraction;
         break;
      case -8:
         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;
         n_me             [1] = 2;
         to_be_subtracted [1] = CollinearInitial;
         n_me             [2] = 1;
         to_be_subtracted [2] = NoSubtraction;
         n_me             [3] = 2;
         to_be_subtracted [3] = NoSubtraction;
         n_me             [4] = 2;
         to_be_subtracted [4] = NoSubtraction;
         break;
      case -7:

         n_me             [0] = 1;
         to_be_subtracted [0] = NoSubtraction;

         n_me             [1] = 1;
         to_be_subtracted [1] = CollinearInitial;

         break;

      case -6:
         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;
         n_me             [1] = 2;
         to_be_subtracted [1] = CollinearInitial;
         n_me             [2] = 1;
         to_be_subtracted [2] = NoSubtraction;
         break;
      case -5:
         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;
         n_me             [1] = 2;
         to_be_subtracted [1] = NoSubtraction;
         n_me             [2] = 2;
         to_be_subtracted [2] = NoSubtraction;
         break;
      case -4:
         n_me             [0] = 3;
         to_be_subtracted [0] = AddedSubtraction;
         break;
      case -2:
         n_me             [0] = 2;
         to_be_subtracted [0] = NoSubtraction;
         break;
      case 1:
         n_me             [0] = 1;
         to_be_subtracted [0] = Subtraction;
         break;
      case 2:
         n_me             [0] = 2;
         to_be_subtracted [0] = Subtraction;
         break;
      case 3:
         n_me             [0] = 3;
         to_be_subtracted [0] = Subtraction;
         break;
      default:
         global -> To_error() 
            << FS("*** %d ***\n", 
                  global -> numberOfFinalStatePartonsInBornTerm
                 );
         errf(-1, "DISProcessNLO_Test::AssignMatrixElements: "
                  "process not known");
   }

   CreateOtherArrays();

   // *** fill in matrix elements and frames here
   switch (global->mel_index) {
      case -10:
         me_list [0][0] 
            = new DISMatrixElement(qXqgFS); 
         me_list [0][1] 
            = new DISMatrixElement(gXqqFS); 
         frame = pCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case -9:
         me_list [0][0] 
            = new DISMatrixElement(qXq); 
         me_list [1][0] 
            = new DISMatrixElement(qXqgA); 
         me_list [2][0] 
            = new DISMatrixElement(qXqg); 
         frame = hCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case -8:
         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 
         me_list [1][0] 
            = new DISMatrixElement(qXqgFS); 
         me_list [1][1] 
            = new DISMatrixElement(gXqqFS); 
         me_list [2][0] 
            = new DISMatrixElement(gXqqFSNF); 
         me_list [3][0] 
            = new DISMatrixElement(qXqgRS); 
         me_list [3][1] 
            = new DISMatrixElement(gXqqRS); 
         me_list [4][0] 
            = new DISMatrixElement(qXqgRSNF); 
         me_list [4][1] 
            = new DISMatrixElement(gXqqRSNF); 
         frame = pCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case -7:

         me_list [0][0] 
            = new DISMatrixElement(qXq); 

         me_list [1][0] 
            = new DISMatrixElement(qXqFS); 

         frame = hCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 1;

         break;
      case -6:
         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 
         me_list [1][0] 
            = new DISMatrixElement(qXqgFS); 
         me_list [1][1] 
            = new DISMatrixElement(gXqqFS); 
         me_list [2][0] 
            = new DISMatrixElement(gXqqFSNF); 
         frame = pCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case -5:
         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 
         me_list [1][0] 
            = new DISMatrixElement(qXqgRS); 
         me_list [1][1] 
            = new DISMatrixElement(gXqqRS); 
         me_list [2][0] 
            = new DISMatrixElement(qXqgRSNF); 
         me_list [2][1] 
            = new DISMatrixElement(gXqqRSNF); 
         frame = pCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case -4:
         me_list [0][0] 
            = new DISMatrixElement(qXqggA); 
         me_list [0][1] 
            = new DISMatrixElement(gXqqgA); 
         me_list [0][2] 
            = new DISMatrixElement(qXqqqA); 
         frame = pCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case -2:
         me_list [0][0] 
            = new DISMatrixElement(qXqgV); 
         me_list [0][1] 
            = new DISMatrixElement(gXqqV); 
         frame = pCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
      case 1:
         errf(-1, "cannot subtract a 1+1 parton matrix element!");
//         me_list [0][0] 
//            = new DISMatrixElement(qXq); 
         break;
      case 2:
         me_list [0][0] 
            = new DISMatrixElement(qXqg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqq); 
         frame = hCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 1;
         break;
      case 3:
         me_list [0][0] 
            = new DISMatrixElement(qXqgg); 
         me_list [0][1] 
            = new DISMatrixElement(gXqqg); 
         me_list [0][2] 
            = new DISMatrixElement(qXqqq); 
         frame = pCMS;
         global -> numberOfFinalStatePartonsInBornTerm = 2;
         break;
   }

   AssignInformationToMatrixElements();
}

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

#if 0
// ++++++++++++++++++++++++++ GARBAGE +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++ END OF GARBAGE +++++++++++++++++++++++++++++++++++++++++
#endif
// ============================================================================
//
// --> DISASTER++ QCD classes
//
// file:              qcd.cc
// created:           05.04.1997
// last modification: 15.12.1997
//
// ============================================================================

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "container.h"

#include "switch.h"
#include "global.h"
#include "mth.h"
#include "strng.h"
#include "user.h"
#include "mcintg.h"
#include "cmb.h"
#include "disaster.h"
#include "proc.h"
#include "rest.h"
#include "pdflib.h"

#include "qcd.h"

// ============================================================================
//
// definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

const char *QName[] = {
               "NoSplittingFunction",
               "Q_q_from_q",
               "Q_G_from_q",
               "Q_q_from_G",
               "Q_G_from_G"
            };

// ============================================================================
//
// --> classes
//
// ============================================================================

int clusterNumber(Array_double *scale2,                                
                  Array_int *nCluster,              
                  int nEntries,                              
                  double testscale2) {

   register int i = 0;
   while ( i < nEntries-1 && scale2 -> GetData(i) <= testscale2 )
      ++i;

   return nCluster -> GetData(i);
}

// ----------------------------------------------------------------------------
//
// --> booking device
//
// assignment: -3 is under, -2 is over, -1 is total
//
// ----------------------------------------------------------------------------

//Book :: Book() {
//
//   defined = FALSE;
//}

// ----------------------------------------------------------------------------

Book :: Book(Global* global_I) {

   global = global_I;
   defined = FALSE;
}

// ----------------------------------------------------------------------------

Book :: ~Book() {

   Delete();
}

// ----------------------------------------------------------------------------

void Book :: Define(
                BookMode bookModeIn, 
                ScaleType scaleTypeIn, 
                int sizeIn,
                double smin, 
                double smax, 
                char* nameIn, 
                char* idIn, 
                int prestoreSizeIn
           ) {

   if (defined)

      errf(-1, "Book::Define: already defined");

   else {

      defined = TRUE;

      bookMode = bookModeIn;
      scaleType = scaleTypeIn;
      size = sizeIn;
      name = new char [strlen(nameIn)+1];
      strcpy(name, nameIn);
      id = new char [strlen(idIn)+1];
      strcpy(id, idIn);
      prestoreSize = prestoreSizeIn;

      switch (bookMode) {
         case Histogram:
            xarray = new Array_double(size+1);
            prestore = new Array_double(0, prestoreSize-1);
            psData = prestore -> GetDataPtr();
            prestoreFlag = new Array_int(0, prestoreSize-1);
            psFData = prestoreFlag -> GetDataPtr();
            prestoreFactor = new Array_double(0, prestoreSize-1);
            psFactData = prestoreFactor -> GetDataPtr();
            break;
         case Graph:
            xarray = new Array_double(size);
            prestore2 = new TwoD_Array_double (0, prestoreSize-1, 
                                               0, size-1);
            ps2Data = prestore2 -> GetDataPtr();
            break;
      }

      yarray     = new Array_double(-3, size-1);
      y2array    = new Array_double(-3, size-1);
      errors     = new Array_double(-3, size-1);
      relerrors  = new Array_double(-3, size-1);
      zarray     = new Array_double(-3, size-1);

      xdata          = xarray    -> GetDataPtr();
      ydata          = yarray    -> GetDataPtr();
      y2data         = y2array   -> GetDataPtr();
      errorsData     = errors    -> GetDataPtr();
      relerrorsData  = relerrors -> GetDataPtr();
      zdata          = zarray    -> GetDataPtr();

      switch (scaleType) {
         case Arbitrary:
            break;
         case Linear:
            xarray -> LinearScale(smin, smax);
            break;
         case Log:
            xarray -> LogScale(smin, smax);
            break;
      }

      as_counter=0;

//      toBeAdded = TRUE;
   }
}

// ----------------------------------------------------------------------------

void Book::Delete() {

   if (!defined)
      errf(-1, "Book::Delete: not defined");
   else {
      delete xarray;
      delete yarray;
      delete y2array;
      delete errors;
      delete relerrors;
      delete zarray;
      delete name;
      delete id;
      switch (bookMode) {
         case Histogram:
            delete prestore;
            delete prestoreFlag;
            delete prestoreFactor;
            break;
         case Graph:
            delete prestore2;
            break;
      }
   }
}

// ----------------------------------------------------------------------------

void Book::ResetStorage() {


   if ( scaleType == Arbitrary && as_counter != as_counter_max ) {
      global -> To_status()
                   << FS("%d ", as_counter)
                   << FS("%d ", as_counter_max)
                   << FS(":%s: ", name)
                   << FS(":%s:\n", id);
      errf(-1, "Book::ResetStorage: not enough entries for array");
   }

   yarray     -> SetToZero();
   y2array    -> SetToZero();
   errors     -> SetToZero();
   relerrors  -> SetToZero();

   counter = 0;
}

// ----------------------------------------------------------------------------

void Book :: CalculateErrors() {

   int ll;
   switch (bookMode) {
      case Histogram:
         ll = -3;
         break;
      case Graph:
         ll = -1;
         break;
      default:
         errf(-1, "Book::CalculateErrors: no such bookMode");
         ll = 0; // :_TAWM_:
   }

   double sq;
   for (register int i = ll; i < size; ++i) {
       sq = counter * y2data[i] - ydata[i] * ydata[i];
       if (counter < 2 || sq < 0.)
          errorsData[i] = 0;
       else
          errorsData[i] = sqrt(1./(counter-1.) * sq);
       if ( !(fabs(ydata[i]) < 1.e-50 || fabs(errorsData[i]) > 1.e50) )
          relerrorsData[i] = fabs(errorsData[i]/ydata[i]);
       else
          relerrorsData[i] = 0;
   }

}

// ----------------------------------------------------------------------------

void Book::New_X_Entry(double x) {

   if (scaleType != Arbitrary)
      errf(-1, "Book::New_X_Entry: illegal");   

   switch (bookMode) {
      case Histogram:
         as_counter_max = size+1;
         break;
      case Graph:
         as_counter_max = size;
         break;
   }        
   if (as_counter >= as_counter_max)
      errf(-1, "Book::New_X_Entry: too many entries in as_counter");

   xarray->SetData(as_counter++, x);
}

// ----------------------------------------------------------------------------

void Book::PrepareForEvent() {

   zarray -> SetToZero();

   switch (bookMode) {
      case Histogram:
//         prestore  -> SetToZero();
         prestoreFlag -> SetToZero();
         break;
      case Graph:
         prestore2 -> SetToZero();
         break;
   }
}

// ----------------------------------------------------------------------------

void Book::AddUpEvent() {

   ++counter; 

   register int i;

   for (i = -3; i < size; ++i) {
       ydata[i]  += zdata[i];
       y2data[i] += zdata[i] * zdata[i];
   }
}

// ----------------------------------------------------------------------------

void Book::StoreHistogram(double x, double y) {

   if (bookMode != Histogram)
      errf(-1, "Book::StoreHistogram: wrong book mode");

   int found, loc;
   found = xarray -> FindBinLocation(x, loc);

   zdata[-1] += y;  // total

   if (found) {
//      global -> To_status() << FS("(%16.6e ", x)
//                            << FS("%16.6e) ", y)
//                            << FS("storing in %16.6e ", xdata[loc])
//                            << FS("%16.6e\n", xdata[loc+1]);
      zdata[loc] += y;
//      global -> To_status() << FS("%s ", name)
//                            << FS("%16.6e\n", zdata[loc]);
//      if (loc == 0)
//         global -> To_status() << FS("storing 0 = %16.6e ", x)
//                               << FS("in %s\n", name);
   }
   else 
      if (loc==1)
         zdata[-2] += y;  // over
      else 
         zdata[-3] += y;  // under
}

// ----------------------------------------------------------------------------

void Book::StoreGraph(int loc, double y) {

   if (bookMode != Graph)
      errf(-1, "Book::StoreGraph: wrong book mode");

   if (loc < 0 || loc >= size)
      errf(-1, "Book::StoreGraph: unknown location");

   zdata[loc] += y;

   zdata[-1]  += y;  // total
}

// ----------------------------------------------------------------------------

void Book::PrestoreHistogram(double x, int k) {

   if (bookMode != Histogram)
      errf(-1, "wrong Book mode");

   if (k < 0 || k >= prestoreSize)
      errf(-1, "Book :: PrestoreHistogram (2): wrong index");

   psData[k]     = x;
   psFData[k]    = TRUE;
   psFactData[k] = 1.;
}

// ----------------------------------------------------------------------------

void Book::PrestoreHistogram(double x, int k, double weight) {

   if (bookMode != Histogram)
      errf(-1, "wrong Book mode");

   if (k < 0 || k >= prestoreSize)
      errf(-1, "Book :: PrestoreHistogram (3): wrong index");

   psData[k]     = x;
   psFData[k]    = TRUE;
   psFactData[k] = weight;
}

// ----------------------------------------------------------------------------

void Book::PrestoreGraph(double x, int event, int entry) {

   if (bookMode != Graph)
      errf(-1, "wrong Book mode");

   if (   event < 0 || event >= prestoreSize
       || entry < 0 || entry >= size
      )
      errf(-1, "Book::PrestoreGraph: wrong index");

   ps2Data[event][entry] = x;
}

// ----------------------------------------------------------------------------

void Book::SPS_Histogram(int k, double y) {

   if (k < 0 || k >= prestoreSize)
      errf(-1, "Book::SPS_Histogram: wrong index");

   if (psFData[k])
      StoreHistogram(psData[k], y * psFactData[k]);
}

// ----------------------------------------------------------------------------

void Book::SPS_Graph(int event, int entry, double y) {

   if (   event < 0 || event >= prestoreSize
       || entry < 0 || entry >= size
      )
      errf(-1, "Book::SPS_Graph: wrong index");

   StoreGraph(entry, y * ps2Data[event][entry]);
}

// ----------------------------------------------------------------------------

void Book :: Print(Mstream& out) {

   register int i;

   global -> To_status()
                << FS(":%s:\n", name);
   global -> To_status()
                << FS("# of events included: %d\n", counter);

   switch (bookMode) {

      case Histogram:
         
         for (i = -3; i < size; ++i) {

             if (i >= 0)
                out
                    << FS("%3d ",        i)
                    << FS("--> [%11.4e", xdata[i])
                    << FS(",%11.4e] ",   xdata[i+1])
                    << "-->";
             else
                out
                    << FS("%3d --> ", i);

             out 
                 << FS("%11.4e ",      ydata[i])
                 << FS("+- %11.4e ",   errorsData[i])
                 << FS("[%11.4e%%]\n", 100 * relerrorsData[i]);
         }

         break;

      case Graph:
         
         for (i = -1; i < size; ++i) {

             if (i >= 0)
                out
                    << FS("%3d --> ",   i)
                    << FS("%14.4e -->", xdata[i]);
             else
                out  
                    << FS("%3d --> ", i);

             out 
                 << FS("%11.4e ",      ydata[i])
                 << FS("+- %11.4e ",   errorsData[i])
                 << FS("[%11.4e%%]\n", 100 * relerrorsData[i]);
         }

         break;
   }
}

// ----------------------------------------------------------------------------

void Book::PrintMathForm(FILE* out, int numid) {

   Mstream s_out;
   s_out.Associate(out);

   const char *name1[] = {
                          "Arbitrary",
                          "linearPlot",
                          "logPlot"
                         };

   const char *name2[] = {
                          "histogramPlot",
                          "graphPlot"
                         };

   char str[1000];

   int stp;
   switch (scaleType) {
      case Arbitrary:
         stp = 0;
         break;
      case Linear:
         stp = 1;
         break;
      case Log:
         stp = 2;
         break;
      default:
         errf(-1, "Book::PrintMathForm: no such scaleType");
         stp = 0; // :_TAWM_: 
   }

   register int i;

//   s_out << "++bookNumber;\n";
//   s_out << FS("finalBookNumber = %d;\n\n", numid);

   s_out << FS("identId [jobName,%3d]", numid)
         << FS(" = \"%s\";\n", id);

   s_out << FS("name    [jobName,%3d]", numid)
         << FS(" = \"%s\";\n", name);

   s_out << FS("modeId  [jobName,%3d]", numid)
         << FS(" = %s;\n", bookMode == Histogram ? name2[0] : name2[1]);

   s_out << FS("scaleId [jobName,%3d]", numid)
         << FS(" = %s;\n", name1[stp]);

   s_out << FS("size    [jobName,%3d]", numid)
         << FS(" = %d;\n", size);

   s_out << FS("nEvents [jobName,%3d]", numid)
         << FS(" = %d;\n\n", counter);

   switch (bookMode) {

      case Histogram:
         
         for (i=-3; i<size; ++i) {
             if (i>=0) {
                mfffp(xdata[i], str);
                s_out << "x0  [jobName,"
                      << FS("%3d", numid)
                      << ","
                      << FS("%3d", i)
                      << "] = (*"
                      << FS("%14.6e", xdata[i])
                      << " *)"
                      << FS("%s", str)
                      << ";\n";
                mfffp(xdata[i+1], str);
                s_out << "x1  [jobName,"
                      << FS("%3d", numid)
                      << ","
                      << FS("%3d", i)
                      << "] = (*"
                      << FS("%14.6e", xdata[i+1])
                      << " *)"
                      << FS("%s", str)
                      << ";\n";
             }
             mfffp(ydata[i], str);
             s_out << "y   [jobName," 
                   << FS("%3d", numid)
                   << ","
                   << FS("%3d", i)
                   << "] = (*"
                   << FS("%14.6e", ydata[i])
                   << " *)"
                   << FS("%s", str)
                   << ";\n";
             mfffp(errorsData[i], str);
             s_out << "err [jobName,"
                   << FS("%3d", numid)
                   << ","
                   << FS("%3d", i)
                   << "] = (*"
                   << FS("%14.6e", errorsData[i])
                   << " *)"
                   << FS("%s", str)
                   << ";\n\n";
         }
         mfffp(xdata[size], str);
         s_out << "x0  [jobName,"
               << FS("%3d", numid)
               << ","
               << FS("%3d", size)
               << "] = (*"
               << FS("%14.6e", xdata[size])
               << " *)"
               << FS("%s", str)
               << ";\n";
         mfffp(xdata[0], str);
         s_out << "x1  [jobName,"
               << FS("%3d", numid)
               << ","
               << FS("%3d", -1)
               << "] = (*"
               << FS("%14.6e", xdata[0])
               << " *)"
               << FS("%s", str)
               << ";\n";
         break;

      case Graph:
         
         for (i=-1; i<size; ++i) {
             if (i>=0) {
                mfffp(xdata[i], str);
                s_out << "x   [jobName,"
                      << FS("%3d", numid)
                      << ","
                      << FS("%3d", i)
                      << "] = (*"
                      << FS("%14.6e", xdata[i])
                      << " *)"
                      << FS("%s", str)
                      << ";\n";
             }
             mfffp(ydata[i], str);
             s_out << "y   [jobName,"
                   << FS("%3d", numid)
                   << ","
                   << FS("%3d", i)
                   << "] = (*"
                   << FS("%14.6e", ydata[i])
                   << " *)"
                   << FS("%s", str)
                   << ";\n";
             mfffp(errorsData[i], str);
             s_out << "err [jobName,"
                   << FS("%3d", numid)
                   << ","
                   << FS("%3d", i)
                   << "] = (*"
                   << FS("%14.6e", errorsData[i])
                   << " *)"
                   << FS("%s", str)
                   << ";\n\n";
         }

         break;
   }
}

// ----------------------------------------------------------------------------

void Book :: GetEntry(  
                int i,
                double* average,
                double* error,
                double* relError
             ) {
   
   if (i >= size) 
      errf(-1, "Book :: GetEntry: index too large");

   int ll;
   switch (bookMode) {
      case Histogram:
         ll = -3;
         break;
      case Graph:
         ll = -1;
         break;
      default:
         global -> To_status() << FS("bookMode: %d\n", (int) bookMode);
         global -> To_status() << FS("%s\n", name);
         errf(-1, "Book :: GetEntry: no such bookMode");
         ll = 0; // :_TAWM_:
   }

   if (i < ll)  {
      global -> To_status() 
                   << "i, ll: "
                   << FS("%d, ", i)
                   << FS("%d\n", ll);
      global -> To_status() << FS("%s\n", name);
      errf(-1, "Book :: GetEntry: index too small");
   }   

   *average  = ydata[i];
   *error    = errorsData[i];
   *relError = relerrorsData[i];
}

// ----------------------------------------------------------------------------
//
// --> manager for the booking device
//
// ----------------------------------------------------------------------------

Library::Library(Global* global_I, int maxBooksIn) {

   defined = FALSE;
   global = global_I;
   Define(maxBooksIn);   
}

// ----------------------------------------------------------------------------

Library::~Library() {
   Delete();
}

// ----------------------------------------------------------------------------

void Library::Define(int maxBooksIn) {

   if (defined)
      errf(-1, "Library::Define: already defined");   
   else {
      maxBooks = maxBooksIn;
      nBooks = 0; 
      books = new ArrayOfPointers <Book*> (0, maxBooks-1); 
      bookData = books -> GetDataPtr();
      books -> SetToNull();
      defined = TRUE;
   }
}

// ----------------------------------------------------------------------------

void Library::Delete() {

   if (!defined)
      errf(-1, "Library::Delete: not defined");   
   else {
      books -> DeleteAllEntries();
      delete books;
      defined = FALSE;
   }
}

// ----------------------------------------------------------------------------

Book* Library :: CreateNewBook() {

   Book *newBook;
   if (nBooks == maxBooks) {
      errf(-1, "Library::CreateNewBook: too many books requested");   
      newBook = NULL; // :_TAWM_:
   } else {
      newBook = new Book(global);
      books -> SetData(nBooks, newBook);
      ++nBooks;
   }
   return newBook;
}

// ----------------------------------------------------------------------------

void Library::Print(Mstream& out) {

   for (register int i = 0; i<nBooks; ++i) {
       out << "\n\n";
       bookData[i] -> Print(out);
       out << "\n\n";
   }
}

// ----------------------------------------------------------------------------

void Library::PrintMathForm(FILE *out) {

   global -> To_status() << FS("bookNumber[jobName] = %d;\n\n", nBooks);

   for (register int i = 0; i<nBooks; ++i)
       bookData[i] -> PrintMathForm(out, i);
}

// ----------------------------------------------------------------------------

void Library::ResetStorage() {

   for (register int i = 0; i<nBooks; ++i)
       bookData[i] -> ResetStorage();
}

// ----------------------------------------------------------------------------

void Library::CalculateErrors() {

   for (register int i = 0; i<nBooks; ++i)
       bookData[i] -> CalculateErrors();
}

// ----------------------------------------------------------------------------

void Library::PrepareForEvent() {

   for (register int i = 0; i<nBooks; ++i)
       bookData[i] -> PrepareForEvent();
}

// ----------------------------------------------------------------------------

void Library::AddUpEvent() {

   for (register int i = 0; i<nBooks; ++i)
       bookData[i] -> AddUpEvent();
}

// ----------------------------------------------------------------------------
//
// --> Cross section contributions
//
// ----------------------------------------------------------------------------

Contribution::Contribution() {

   defined=FALSE;
}

// ----------------------------------------------------------------------------

Contribution::Contribution(Global *globalIn, int lengthIn) {

   defined=FALSE;
   Define(globalIn, lengthIn);
}

// ----------------------------------------------------------------------------

Contribution::~Contribution() {

   Delete();
}

// ----------------------------------------------------------------------------

void Contribution::Define(Global *globalIn, int lengthIn) {

   if (defined)
      errf(-1, "Contribution::Define: already defined");

   global=globalIn;
   length=lengthIn;

   clist = new Array <double> (length);
   data = clist -> GetDataPtr();

   defined=TRUE;
}

// ----------------------------------------------------------------------------

void Contribution::Delete() {

   if (!defined)
      errf(-1, "Contribution::Delete: not defined");

   delete clist;

   defined=FALSE;
}

// ----------------------------------------------------------------------------

void Contribution::Reset(Global *globalIn) {
   global=globalIn;
}

// ----------------------------------------------------------------------------

void Contribution::ClearEntry() {

   for (register int i=0; i<length; ++i)
       data[i] = 0.;
   orderAlphaS = 0;
   slf = NoScaleLogarithm;
   matrixElement = 0;
}

// ----------------------------------------------------------------------------

void Contribution::Print(FILE* out, int mode) {

   Mstream s_out;
   s_out.Associate(out);

   s_out << FS("M(%d) ", matrixElement)
         << FS("A(%d) ", orderAlphaS)
         << FS("S(%d) ", slf)
         << FS("L(%d) ", length)
         << FS("NO(%d) ", nflavour_out_U)
         << FS("NP(%d) ", nflavour_pden_U)
         << FS("R(%10.3e) ", ren_scale_U)
         << FS("F(%10.3e) ", fact_scale_U)
         << FS("X(%10.3e)\n", xi_U);

   if (mode == 1) {
      for (register int i=0; i<length; ++i)
          s_out << FS("        %4d", i)
                << FS(" --> %16.6e\n", data[i]);
   }
   if (mode == 2)
      s_out << FS("id(%ld) ", GetEvent() -> GetIdent())
            << FS("nxt(%d)\n",
                  (GetSameEvent() == NULL) 
                     ? -1 
                     : GetSameEvent() -> GetLocation()
               );
}

// ----------------------------------------------------------------------------

void Contribution::Print(FILE* out) {

   // :_TAWM_: 
   out = out;

   errf(-1, "Contribution::Print: not yet implemented");
}

// ----------------------------------------------------------------------------
   
int Contribution::UseListIsSame(Contribution* other) {
   // TAWM
   other = other;

   errf(-1, "Contribution::UseListIsSame: not yet implemented"); 
   return 0;
}

// ----------------------------------------------------------------------------
//
// --> Array of contributions
//
// ----------------------------------------------------------------------------

ContributionArray :: ContributionArray() {

   defined = FALSE;
}

// ----------------------------------------------------------------------------
  
ContributionArray::ContributionArray(
                      Global *globalIn, 
                      int maxContributionsIn,
                      int lengthIn) {

   defined=FALSE;
   Define(globalIn, maxContributionsIn, lengthIn);
}

// ----------------------------------------------------------------------------

ContributionArray::~ContributionArray() {

   Delete();
}

// ----------------------------------------------------------------------------

void ContributionArray::Define(
                           Global *globalIn, 
                           int maxContributionsIn,
                           int lengthIn) {

   if (defined)
      errf(-1, "ContributionArray::Define: already defined"); 

   global = globalIn;
   maxContributions=maxContributionsIn;
   length = lengthIn;

   entry = new Array <Contribution> (maxContributions);
   cptr  = entry -> GetDataPtr();

   register int i;
   for (i=0; i<maxContributions; ++i) {
       cptr[i].Define(global, length);
       cptr[i].SetLocation(i);
   }
   ResetActual();

   defined=TRUE;
}

// ----------------------------------------------------------------------------

void ContributionArray::Delete() {

   if (!defined)
      errf(-1, "ContributionArray::Delete: not defined");

   delete entry;

   defined=FALSE;
}

// ----------------------------------------------------------------------------

void ContributionArray::Reset(Global *globalIn) {

   global = globalIn;

   ResetActual();
}

// ----------------------------------------------------------------------------

void ContributionArray::Print(FILE* out, int mode) {

   Mstream s_out;
   s_out.Associate(out);

   s_out << "\n======================================================\n";
   s_out << FS("actual(%d) ", actual)
         << FS("max(%d)\n", maxContributions);
   for (register int i=0; i<=actual; ++i) {
       s_out << FS("%d: ", i);
       cptr[i].Print(out, mode);
   }
   s_out << "======================================================\n";
}

// ----------------------------------------------------------------------------

void ContributionArray::Print(FILE* out) {

   // :_TAWM_: 
   out = out;

   errf(-1, "ContributionArray::Print: not yet implemented");
}

// ----------------------------------------------------------------------------

int ContributionArray::GetActual_SetToZero_Increment() {

   cptr[actual].ClearEntry();
   actual++;
   return actual-1;   
   // :_CAUTION_: in ``return actual++'', is ``actual'' actually incremented?
}

// ----------------------------------------------------------------------------

int ContributionArray::Increment_SetToZero_GetActual() {

   actual++;
   cptr[actual].ClearEntry();
   return actual;   
}

// ----------------------------------------------------------------------------

Contribution *ContributionArray::NextContribution() {

   int i = Increment_SetToZero_GetActual();
   return cptr+i;   
}

// ----------------------------------------------------------------------------

int ContributionArray::UseListIsSame(ContributionArray* other) {

   // TAWM
   other = other;

   errf(-1, "ContributionArray::UseListIsSame: not yet implemented");   
   return 0;
}

// ----------------------------------------------------------------------------
//
// --> Parton lists
//
// ----------------------------------------------------------------------------

PartonList::PartonList() {

   defined=FALSE;
}

// ----------------------------------------------------------------------------
  
PartonList::PartonList(Global *globalIn, int nflavoursIn) {

   defined=FALSE;
   Define(globalIn, nflavoursIn);
}

// ----------------------------------------------------------------------------

PartonList::~PartonList() {

   Delete();
}

// ----------------------------------------------------------------------------

void PartonList::Define(Global *globalIn, int nflavoursIn) {

   if (!defined) {
      defined=TRUE;

      global = globalIn;
      nflavours = nflavoursIn;
      data = new Array_double;
      data -> Define(-nflavours,nflavours);
   } else 
      errf(-1, "PartonList::Define: already defined");
}

// ----------------------------------------------------------------------------

void PartonList::Delete() {

   if (defined) {

      defined=FALSE;

      delete data;
   }
}

// ----------------------------------------------------------------------------

void PartonList::Reset(Global *globalIn) {

   global = globalIn;
}

// ----------------------------------------------------------------------------
   
void PartonList :: Print(Mstream& s_out) {

   s_out << "PartonList:\n";
   data -> Print(s_out);
}

// ----------------------------------------------------------------------------
   
void PartonList::Print(FILE* out) {

   Mstream s_out;
   s_out.Associate(out);

   Print(s_out);
}

// ----------------------------------------------------------------------------
   
int PartonList::UseListIsSame(PartonList* other) {

   // :_TAWM_:
   other = other;

   errf(-1, "PartonList::UseListIsSame: not yet implemented"); 
   return 0;
}

// ----------------------------------------------------------------------------
   
void PartonList::ReturnToFreeList() {
}

// ----------------------------------------------------------------------------
//
// --> Array of parton lists
//
// :MOD: no longer required
//
// ----------------------------------------------------------------------------

PartonListArray::PartonListArray() {

   defined=FALSE;
}

// ----------------------------------------------------------------------------
  
PartonListArray::PartonListArray(Global *globalIn, int lengthIn) {

   defined=FALSE;
   Define(globalIn, lengthIn);
}

// ----------------------------------------------------------------------------

PartonListArray::~PartonListArray() {

   Delete();
}

// ----------------------------------------------------------------------------

void PartonListArray::Define(Global *globalIn, int lengthIn) {

   if (!defined) {
      defined=TRUE;

      global = globalIn;
      length = lengthIn;
      data = new Array <PartonList>;
      data -> Define(length);
      ptr = GetDataPtr();

      register int i;
      for (i=0; i<length; ++i)
          ptr[i].Define(global, 6); // :MOD: This (6) should be more general!

   } else 
      errf(-1, "PartonListArray::Define: already defined");
}

// ----------------------------------------------------------------------------

void PartonListArray::Delete() {

   if (defined) {

      defined=FALSE;

      delete data;
   }
}

// ----------------------------------------------------------------------------

void PartonListArray::Reset(Global *globalIn) {

   global = globalIn;
}

// ----------------------------------------------------------------------------
//
// --> parton densities
//
// ----------------------------------------------------------------------------

PartonDensity::PartonDensity() {

   defined=FALSE;
}

// ----------------------------------------------------------------------------

PartonDensity::~PartonDensity() {

   Delete();
}

// ----------------------------------------------------------------------------

void PartonDensity::Define(Global *globalIn, int nflavoursIn) {

   if (!defined) {
      defined=TRUE;

      nflavours = nflavoursIn;
      PartonList::Define(globalIn, nflavours);

   } else

      errf(-1, "PartonDensity::Define: already defined");
}

// ----------------------------------------------------------------------------

void PartonDensity::Delete() {

   if (defined) {
      defined=FALSE;
   } else
     errf(-1, "PartonDensity::Delete: not defined");
}

// ----------------------------------------------------------------------------

void PartonDensity::Reset(Global *globalIn) {

   global = globalIn;
   PartonList::Reset(global);  
}

// ----------------------------------------------------------------------------
   
void PartonDensity :: Print(Mstream& s_out) {

   s_out << "PartonDensity:\n";
   PartonList :: Print(s_out);
}

// ----------------------------------------------------------------------------
   
void PartonDensity::Print(FILE* out) {

   Mstream s_out;
   s_out.Associate(out);

   Print(s_out);
}

// ----------------------------------------------------------------------------
   
int PartonDensity::UseListIsSame(PartonDensity* other) {
   // TAWM
   other = other;

   errf(-1, "PartonDensity::UseListIsSame: not yet implemented"); 
   return 0;
}

// ----------------------------------------------------------------------------
   
void PartonDensity::ReturnToFreeList() {
}

// ----------------------------------------------------------------------------
//
// --> parton density parametrizations
//
// ----------------------------------------------------------------------------

PartonDensityParametrization :: PartonDensityParametrization(
                                   Global* globalIn
                                ) {

   Set_global(globalIn);

   defined = FALSE;
   collection = -1;
   parametrization =-1;
   set = -1;
}

// ----------------------------------------------------------------------------

PartonDensityParametrization :: ~PartonDensityParametrization() {

   Delete();
}

// ----------------------------------------------------------------------------

void PartonDensityParametrization :: Define() {

   if (!defined) {

//      global -> To_status() << "PartonDensityParametrization :: Define\n";

      defined = TRUE;
      max_flavours = 6;
      pdfl  = NULL;
      cache = NULL;
   }
}

// ----------------------------------------------------------------------------

void PartonDensityParametrization :: Delete() {

   global -> To_status() << "PartonDensityParametrization :: Delete\n";

   if (defined) {
      defined = FALSE;

      if (cache != NULL) {
         global -> To_status() 
            << "PartonDensityParametrization :: Delete (cache)\n";
         DeleteCache();
      }

   } else

      errf(-1, "PartonDensityParametrization :: Delete: error");
}

// ----------------------------------------------------------------------------

void PartonDensityParametrization :: FillPartonDensity(
        PartonDensity* pdensity, 
        double x, 
        double scale, 
        int nquarks
     ) {

   switch (collection) { 

      case 0:
         FillPartonDensity_0(pdensity, x, scale, nquarks);
         break;

      case 1:
         FillPartonDensityPDFLIB(pdensity, x, scale, nquarks);
         break;

      default:
         errf(-1, "PartonDensityParametrization :: Fill...: unknown coll.");
   }
}

// ----------------------------------------------------------------------------

void PartonDensityParametrization :: FillPartonDensity(
                                        PartonDensity* pdensity, 
                                        double x,    
                                        double scale, 
                                        int nquarks,
                                        int collectionIn,
                                        int parametrizationIn, 
                                        int setIn
                                     ) { 

   InitializeParametrization(collectionIn, parametrizationIn, setIn);
   FillPartonDensity(pdensity, x, scale, nquarks);   
}

// ----------------------------------------------------------------------------

PartonDensity* PartonDensityParametrization
                  :: CreateAndFillPartonDensity(
                        double x, 
                        double scale, 
                        int nquarks,
                        int collectionIn,
                        int parametrizationIn, 
                        int setIn
                     ) {

   PartonDensity* pdensity;

   if (pdfl != NULL)
      pdensity = pdfl -> GetDefinedObject(global, max_flavours);
   else {
      pdensity = new PartonDensity;
      pdensity -> Define(global, max_flavours);
   }

   FillPartonDensity(
      pdensity, x, scale, nquarks,
      collectionIn, parametrizationIn, setIn
   );   

   return pdensity;
}

// ----------------------------------------------------------------------------

PartonDensity* PartonDensityParametrization
                  :: FindOrCreateAndFillAndCachePartonDensity(
                        double x, 
                        double scale, 
                        int nquarks,
                        int collectionIn,
                        int parametrizationIn, 
                        int setIn
                     ) {

   PartonDensity* pdensity;

   if (cache == NULL)
      errf(-1, "PartonDensityParametrization"
               " :: FindOrCreateAndFillAndCachePartonDensity: no cache!");

   intKey[0]    = nquarks;
   intKey[1]    = collectionIn;
   intKey[2]    = parametrizationIn;
   intKey[3]    = setIn;
   doubleKey[0] = x;
   doubleKey[1] = scale;

   // try to find in cache
   pdensity = cache -> FindInCache(intKeyArray, doubleKeyArray);

   if (pdensity == NULL) {
      pdensity
         = CreateAndFillPartonDensity(
              x, scale, nquarks,
              collectionIn, parametrizationIn, setIn
           );   
      cache -> StoreInCache(intKeyArray, doubleKeyArray, pdensity);
   }

//   pdensity -> Print();
   
   return pdensity;
}

// ----------------------------------------------------------------------------

// test parton densities

void PartonDensityParametrization::FillPartonDensity_0(
        PartonDensity *pdensity, 
        double x, double scale, int nquarks
     ) {

   // TAWM
   x=x;
   scale=scale;
   nquarks=nquarks;

   double *data=pdensity->GetDataPtr();

   register int i;

   for (i = -max_flavours; i <= max_flavours; ++i)
       data[i]=0.;

   switch (set) {
      case 0: // set all densities to 1
         for (i = -max_flavours; i <= max_flavours; ++i)
             data[i]=1.;
         break;
      case 1: // set all densities to 1-x (i.e. valence and sea)
         for (i = -max_flavours; i <= max_flavours; ++i)
             data[i]=1.-x;
         data[1] *= 2;
         data[2] *= 2;
         break;
      case 2: // set all densities to 1 (i.e. valence and sea)
         for (i = -max_flavours; i <= max_flavours; ++i)
             data[i]=1.;
         data[1] *= 2;
         data[2] *= 2;
         break;
      case 3: // set all densities to 1-x (i.e. valence and sea)
         for (i = -max_flavours; i <= max_flavours; ++i)
             data[i]=pow(1.-x, 2);
         data[1] *= 2;
         data[2] *= 2;
         break;
      case 4: // set all densities to 1-x (i.e. valence and sea)
         for (i = -max_flavours; i <= max_flavours; ++i)
             data[i]=pow(1.-x, 5);
         data[1] *= 2;
         data[2] *= 2;
         break;
      case -2: // set all quarks to 1
         for (i = -max_flavours; i <= max_flavours; ++i)
             data[i]=1.;
         data[0] = 0.;
         break;
      case -3: // set the gluon to 1
         for (i = -max_flavours; i <= max_flavours; ++i)
             data[i]=0.;
         data[0] = 1.;
         break;
      default:
         errf(-1, "PartonDensityParametrization::Fill..._0: unknown set");
   }
}

// ----------------------------------------------------------------------------

// PDFLIB

void PartonDensityParametrization::FillPartonDensityPDFLIB(
        PartonDensity *pdensity, 
        double x, double scale, int nquarks
     ) {

   // TAWM
   nquarks=nquarks;

   double *data=pdensity->GetDataPtr();

   pdfgetC(x, scale, data);
}

// ----------------------------------------------------------------------------

void PartonDensityParametrization :: InitializeCollection() {

   switch (collection) { 

      case 0:   // test distributions
         nflavours = 6;
         order = 1;
         lambda_4 = 0.2;
         x_min = 0.0;
         x_max = 1.0; 
         scale_min = 0.0; 
         scale_max = 1.e15;
         break;

      case 1:   // PDFLIB
         pdfsetparC(
            GetPDFParametrization(), GetPDFSet(),
            nflavours, order,
            lambda_4,
            x_min, x_max, scale_min, scale_max
         );
         break;

      default:
         errf(-1, "PartonDensityParametrization::Fill...: unknown coll.");
   }

   global -> To_status() 
      << "      Parton density parametrization:\n";
   global -> To_status() 
      << FS("      Collection:      %4d\n", collection);
   global -> To_status() 
      << FS("      Parametrization: %4d\n", parametrization);
   global -> To_status() 
      << FS("      Set:             %4d\n", set);
   global -> To_status() 
      << FS("      # of flavours:   %4d\n", nflavours);
   global -> To_status() 
      << FS("      order:           %4d\n", order);
   global -> To_status() 
      << FS(       "      %16.6e",   x_min)
      << FS(" <=   x   <= %16.6e\n", x_max);
   global -> To_status() 
      << FS(       "      %16.6e",         scale_min)
      << FS(" <=  mu_f <= %16.6e [GeV]\n", scale_max);
   global -> To_status() << "\n";
}

// ----------------------------------------------------------------------------

void PartonDensityParametrization :: InitializeParametrization(
                                        int collectionIn,
                                        int parametrizationIn,
                                        int setIn
                                     ) {

   if (   collectionIn      != collection
       || parametrizationIn != parametrization
       || setIn             != set 
      ) {
      collection = collectionIn;
      parametrization = parametrizationIn;
      set = setIn;
      InitializeCollection();
   }
}

// ----------------------------------------------------------------------------

void PartonDensityParametrization::CreateCache(int size) {

   if (cache != NULL)
      errf(-1, "PartonDensityParametrization::CreateCache: "
               "cache already exists");

   int nInt = 4;
   int nDouble = 2;
   cache = new Cache <PartonDensity> (global, size, nInt, nDouble);
   cache -> DefineName("PartonDensityParametrization cache");
   cache -> SetFreeList(pdfl);
   global -> To_status() << "PartonDensityParametrization::CreateCache\n";
#if 0
   delete cache;
   cache = new Cache <PartonDensity> (global, size, nInt, nDouble);
   cache -> DefineName("PartonDensityParametrization cache");
   cache -> SetFreeList(pdfl);
   global -> To_status() << "PartonDensityParametrization::CreateCache (2)\n";
#endif
   intKeyArray = new Array_int(nInt);
   intKey = intKeyArray -> GetDataPtr();
   doubleKeyArray = new Array_double(nDouble);
   doubleKey = doubleKeyArray -> GetDataPtr();
}

// ----------------------------------------------------------------------------

void PartonDensityParametrization :: DeleteCache() {

   delete cache;
   cache = NULL;

   delete intKeyArray;
   intKeyArray = NULL;

   delete doubleKeyArray;
   doubleKeyArray = NULL;
}

// ----------------------------------------------------------------------------
//
// --> coupling constant evaluation
//
// ----------------------------------------------------------------------------

CouplingConstantEvaluation::CouplingConstantEvaluation(Global* globalIn) {

   Set_global(globalIn);

   defined=FALSE;

   outputFlag = 1; // default: output of all parameter changes
}

// ----------------------------------------------------------------------------

CouplingConstantEvaluation::~CouplingConstantEvaluation() {

   Delete();
}

// ----------------------------------------------------------------------------

void CouplingConstantEvaluation::Define() {

   if (!defined) {
      defined=TRUE;
      SetVariant(1);
      SetOrder(-1);
   }
}

// ----------------------------------------------------------------------------

void CouplingConstantEvaluation::Delete() {

   if (defined) {
      defined=FALSE;
   }
}

// ----------------------------------------------------------------------------

double CouplingConstantEvaluation::EvaluateCoupling(
        double scale, int nquarks
     ) {

   // :_TAWM_:
   scale=scale;
   nquarks=nquarks;

   double coupling;

   switch(variant) {
      case 0:
         coupling = 1.;
         break;
      default:
         errf(-1, "CouplingConstantEvaluation::EvaluateCoupling: variant?");
         coupling = -1.; // :_TAWM_:
   }

   return coupling;
}

// ----------------------------------------------------------------------------
//
// --> the strong coupling constant
//
// ----------------------------------------------------------------------------

AlphaS :: AlphaS(Global* globalIn)
   : CouplingConstantEvaluation(globalIn)
   {
}

// ----------------------------------------------------------------------------

AlphaS :: ~AlphaS() {

   Delete();
}

// ----------------------------------------------------------------------------

void AlphaS :: Define() {

   if (!GetDefined()) {

      CouplingConstantEvaluation :: Define();

      SetMcharm  (  1.35);
      SetMbottom (  5.00);
      SetMtop    (170.00);

      SetLambdaQCD_4(-1.0);
   }
}

// ----------------------------------------------------------------------------

void AlphaS :: Delete() {

   if (GetDefined()) {
   }
}

// ----------------------------------------------------------------------------

double AlphaS :: EvaluateAlphaS(double scale, int nquarks) {

   double coupling;

   switch(variant) {

      case -1:
         coupling = 0.1;
         break;

      case 0:
         coupling 
            = CouplingConstantEvaluation :: EvaluateCoupling(scale, nquarks);
         break;

      case 1:
         coupling = EvaluateRunningCoupling(scale, nquarks);
         break;

      case 2:
         coupling = pdfasC(scale);
         break;

      default:
         errf(-1, "AlphaS :: EvaluateAlphaS: variant not known");
         coupling = -1.; // :_TAWM_:
   }

   return coupling;
}

// ----------------------------------------------------------------------------

double AlphaS :: EvaluateAlphaSAutomaticNquarks(double scale) {

   return EvaluateAlphaS(scale, FlavourSwitch(scale));
}

// ----------------------------------------------------------------------------

double AlphaS :: EvaluateAlphaSAutomaticNquarks(
                    double scale,
                    int variantIn,
                    int orderIn,
                    double lambdaQCD4In
                 ) {

   InitializeAlphaS(variantIn, orderIn, lambdaQCD4In);

   return EvaluateAlphaSAutomaticNquarks(scale);
}

// ----------------------------------------------------------------------------

void AlphaS :: CalculateLambdaValues() {

#if 0
         global -> To_status() << FS("%e\n", mCharm);
         global -> To_status() << FS("%e\n", mBottom);
         global -> To_status() << FS("%e\n", mTop); 
         global -> To_status() << "\n";
#endif

   switch (order) {

      case 1:
         lambdaQCD[3] =   lambdaQCD[4] 
                        * pow( mCharm      / lambdaQCD[4], 2./27. );
         lambdaQCD[5] =   lambdaQCD[4]  
                        * pow( lambdaQCD[4] / mBottom,     2./23. );
         lambdaQCD[6] =   lambdaQCD[5] 
                        * pow( lambdaQCD[5] / mTop,        2./21. );
      break;

      case -2:
      case  2:
         lambdaQCD[3] = lambdaQCD[4] * pow( mCharm / lambdaQCD[4], 2./27. )
                      *pow( 2. * log( mCharm / lambdaQCD[4] ), 107./2025. );
         lambdaQCD[5] = lambdaQCD[4] * pow( lambdaQCD[4] / mBottom, 2./23. )
                      *pow( 2. * log( mBottom / lambdaQCD[4] ), -963./13225. );
         lambdaQCD[6] = lambdaQCD[5] * pow( lambdaQCD[5] / mTop, 2./21. )
                      *pow( 2. * log( mTop / lambdaQCD[5] ), -321./3381. );
         break;

      default:
         errf(-1, "AlphaS :: CalculateLambdaValues: order not known");
   }
}

// ----------------------------------------------------------------------------

double AlphaS :: EvaluateRunningCoupling(double scale, int nquarks) {

   if (nquarks < 3 || nquarks > 6)
      errf(-1, "AlphaS :: EvaluateRunningCoupling: nquarks out of range");

   double scalelog = 2. * log(scale / lambdaQCD[nquarks]);

   if (scalelog <= 0.)
      errf(-1, "AlphaS :: EvaluateRunningCoupling: scale too small");

   double alpha_s;

   double one_over_scalelog = 1. / scalelog;

   switch (order) {

      case 1:
         alpha_s = 12. * Pi / (33. - 2. * nquarks) * one_over_scalelog;
         break;

      case 2:
         alpha_s =  12. * Pi / (33. - 2. * nquarks) * one_over_scalelog
                  * (1. - 6. * (153. - 19. * nquarks)
                       /((33. - 2. * nquarks) * (33. - 2. * nquarks))
                          * log(scalelog) * one_over_scalelog
                    );
         break;

      case 3:
         alpha_s =  12. * Pi / (33. - 2. * nquarks) * one_over_scalelog
                  / (1. + 6. * (153. - 19. * nquarks)
                       /((33. - 2. * nquarks) * (33. - 2. * nquarks))
                          * log(scalelog) * one_over_scalelog
                    );
         break;

      default:
         errf(-1, "AlphaS :: EvaluateRunningCoupling: order not known"); 
         alpha_s=0.;  // :_TAWM_:
   }

   return alpha_s;
}

// ----------------------------------------------------------------------------

// calculate the number of flavours depending on the scale

int AlphaS :: FlavourSwitch(double scale) {

   int nflavours;

   if (scale < mCharm)
      nflavours = 3;
   else
   if (scale < mBottom)
      nflavours = 4;
   else
   if (scale < mTop)
      nflavours = 5;
   else
      nflavours = 6;

   return nflavours;
}

// ----------------------------------------------------------------------------

double AlphaS :: EvaluateRunningCouplingAutomaticNquarks(double scale) {

   return EvaluateRunningCoupling(scale, FlavourSwitch(scale));
}

// ----------------------------------------------------------------------------

void AlphaS :: InitializeAlphaS(
                  int variantIn,
                  int orderIn,   
                  double lambdaQCD4In  
               ) {

   if (   variantIn != variant
       || orderIn   != order
       ||    variantIn == 1 
          && ! compareDouble(lambdaQCD4In, lambdaQCD[4], 1.e-8, 0)
      ) {

      variant      = variantIn;
      order        = orderIn;
      lambdaQCD[4] = lambdaQCD4In;

      if (outputFlag == 1)
         global -> To_status() << FS("      alpha_s: Variant=%d", variant)
                               << FS(", Order=%d\n", order);
   
      switch (variant) {

         case -1:
            if (outputFlag == 1)
               global -> To_status() << "alpha_s = 0.1\n";
            break;

         case 0:
            if (outputFlag == 1)
               global -> To_status() <<  "alpha_s = 1.0\n";
            break;

         case 1:
            CalculateLambdaValues();
            if (outputFlag == 1) {
               for (int i = 3; i <= 6; ++i)
                   global -> To_status() 
                      << FS("         Lambda[%2d] = ", i)
                      << FS("%16.6e\n", lambdaQCD[i]);
               global -> To_status() << "\n";
            }
            break;

         case 2:
            pdfassetparC(order, &(lambdaQCD[4]));
            if (outputFlag == 1) {
               global -> To_status() 
                  << "initialization of alpha_s from PDFLIB\n";
               global -> To_status()  
                  << FS("         Lambda[4] = %16.6e\n\n", lambdaQCD[4]);
            }
            break;

         default:
            errf(-1, "AlphaS :: InitializeAlphaS: no such variant");
      }
   }
}

// ----------------------------------------------------------------------------
//
// --> the fine structure constant
//
// ----------------------------------------------------------------------------

AlphaEM :: AlphaEM(Global* globalIn)
   : CouplingConstantEvaluation(globalIn)
   {
}

// ----------------------------------------------------------------------------

AlphaEM :: ~AlphaEM() {

   Delete();
}

// ----------------------------------------------------------------------------

void AlphaEM :: Define() {

   if (!GetDefined()) {
      CouplingConstantEvaluation :: Define();
      fine_structure_constant = 1./137.0359895;
   }
}

// ----------------------------------------------------------------------------

void AlphaEM :: Delete() {

   if (GetDefined()) {
//      CouplingConstantEvaluation::Delete();  // Delete automatically called!
   }
}

// ----------------------------------------------------------------------------

double AlphaEM::EvaluateAlphaEM(double scale, int variantIn) {

   double coupling; 

   switch(variantIn) {
      case 0:
         coupling = CouplingConstantEvaluation :: EvaluateCoupling(scale, -1);
         break;
      case 1:
         coupling = fine_structure_constant;
         break;
      case 2:
         coupling = 1. / 137.;
         break;
      default:
         errf(-1, "AlphaEM :: EvaluateAlphaEM (2): variantIn not known");
         coupling = -1.; // :_TAWM_:
   }

   return coupling;
}

// ----------------------------------------------------------------------------

double AlphaEM :: EvaluateAlphaEM(double scale) {

   return EvaluateAlphaEM(scale, variant);
}

// ----------------------------------------------------------------------------
//
// --> QCD
//
// ----------------------------------------------------------------------------

QCD :: QCD() {

   defined = FALSE;
}

// ----------------------------------------------------------------------------

QCD :: ~QCD() {

   Delete();
}

// ----------------------------------------------------------------------------

void QCD :: Define(Global* globalIn) {

   if (!defined) { 

      defined = TRUE;

      global = globalIn;

//      global -> To_status() << "QCD :: Define\n";

      pdFree = new PartonDensityFreeList(global, 0, 10);
      pdFree -> DefineName("PartonDensity free list");

      alphaS_Server = new AlphaS(global);
      alphaS_Server -> Define();

      pdf_Server = new PartonDensityParametrization(global);
      pdf_Server -> Define();
      pdf_Server -> SetPDFL(pdFree);
      pdf_Server -> CreateCache(20);
// DEBUG_2      pdf_Server -> CreateCache(500);

      ncolour = 3.0;
      cf = (ncolour * ncolour - 1.) / (2. * ncolour);

      quarkfactor = 1. / ncolour;
      gluonfactor = 1. / (ncolour * ncolour - 1.);

      nflavours = 6;
      quark_charge = new PartonList;
      quark_charge -> Define(this, nflavours);
      double* qc = quark_charge -> GetDataPtr();
      qc[-6] = -2./3.;
      qc[-5] =  1./3.;
      qc[-4] = -2./3.;
      qc[-3] =  1./3.;
      qc[-2] = -2./3.;
      qc[-1] =  1./3.;
      qc[ 0] =  0.;
      qc[ 1] = -1./3.;
      qc[ 2] =  2./3.;
      qc[ 3] = -1./3.;
      qc[ 4] =  2./3.;
      qc[ 5] = -1./3.;
      qc[ 6] =  2./3.;

      register int i, j, n;

      quarkChargeSquaredArray = new Array_double(-nflavours, nflavours); 
      quarkChargeSquared = quarkChargeSquaredArray -> GetDataPtr();
      
      quarkChargeProductArray 
         = new TwoD_Array_double(-nflavours, nflavours,
                                 -nflavours, nflavours); 
      quarkChargeProduct = quarkChargeProductArray -> GetDataPtr();

      for (i = -nflavours; i <= nflavours; ++i) {
          quarkChargeSquared[i] = qc[i] * qc[i];
          for (j = -nflavours; j <= nflavours; ++j)
              quarkChargeProduct[i][j] = qc[i] * qc[j];
      }

      quarkChargeSquaredPartialSumSingleArray 
         = new Array_double(0, nflavours); 
      quarkChargeSquaredPartialSumSingle 
         = quarkChargeSquaredPartialSumSingleArray -> GetDataPtr();

      double sum;      

      for (n = 0; n <= nflavours; ++n) {
          sum = 0;
          for (i = n; i >= 1; --i)
              sum += quarkChargeSquared[i];
          quarkChargeSquaredPartialSumSingle[n] = sum;
      }

      quarkChargePartialSumDoubleArray 
         = new TwoD_Array_double(0, nflavours,
                                 0, nflavours); 
      quarkChargePartialSumDouble 
         = quarkChargePartialSumDoubleArray -> GetDataPtr();

      quarkChargeSquaredPartialSumDoubleArray 
         = new TwoD_Array_double(0, nflavours,
                                 0, nflavours); 
      quarkChargeSquaredPartialSumDouble 
         = quarkChargeSquaredPartialSumDoubleArray -> GetDataPtr();

      double sum1;      
      double sum2;      

      for (n = 0; n <= nflavours; ++n)
          for (i = 0; i <= nflavours; ++i) {
             sum1 = 0;
             sum2 = 0;
             for (j = n; j >= 1; --j)
                 if (j != i) {
                    sum1 += qc[j];
                    sum2 += quarkChargeSquared[j];
                 }
             quarkChargePartialSumDouble[n][i]        = sum1;
             quarkChargeSquaredPartialSumDouble[n][i] = sum2;
          }
   }
}

// ----------------------------------------------------------------------------

void QCD :: Delete() {

   if (defined) {

      global -> To_status() << "QCD :: Delete\n";

      defined = FALSE;

      delete quark_charge;
      delete quarkChargeSquaredArray;
      delete quarkChargeProductArray;
      delete quarkChargeSquaredPartialSumSingleArray;
      delete quarkChargePartialSumDoubleArray;
      delete quarkChargeSquaredPartialSumDoubleArray;

      delete pdf_Server;
 
      delete alphaS_Server;

      delete pdFree;
   }
}

// ----------------------------------------------------------------------------

double QCD :: QSplitRegularEps0(SplittingFunction sf, double u) {

   double value; 

   switch (sf) {
      case Q_q_from_q:
         value = - 1. - u;
         break;
      case Q_G_from_q:
         value = - 2. + u;
         break;
      case Q_q_from_G:
         value = 1. - 2. * u * (1. - u);
         break;
      case Q_G_from_G:
         value = u * (1. - u) - 2.;
         break;
      default: 
         errf(-1, "QSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: QSplitRegularEps1(SplittingFunction sf, double u) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = - 1. + u;
         break;
      case Q_G_from_q:
         value = - u;
         break;
      case Q_q_from_G:
         value = - 2. * u * (1. - u);
         break;
      case Q_G_from_G:
         value = 0.;
         break;
      default: 
         errf(-1, "QSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: QSplitHalfRegularEps0(SplittingFunction sf, double u) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = - 1. - u;
         break;
      case Q_G_from_q:
         value = 2. / u - 2. + u;
         break;
      case Q_q_from_G:
         value = 1. - 2. * u * (1. - u);
         break;
      case Q_G_from_G:
         value = 1. / u + u * (1. - u) - 2.;
         break;
      default: 
         errf(-1, "QSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD::QSplitHalfRegularEps1(SplittingFunction sf, double u) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = - 1. + u;
         break;
      case Q_G_from_q:
         value = - u;
         break;
      case Q_q_from_G:
         value = - 2. * u *(1. - u);
         break;
      case Q_G_from_G:
         value = 0.;
         break;
      default: 
         errf(-1, "QSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: QSplitAt0(SplittingFunction sf) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = 0.;
         break;
      case Q_G_from_q:
         value = 2.;
         break;
      case Q_q_from_G:
         value = 0.;
         break;
      case Q_G_from_G:
         value = 1.;
         break;
      default: 
         errf(-1, "QSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: QSplitAt1(SplittingFunction sf) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = 2.;
         break;
      case Q_G_from_q:
         value = 0.;
         break;
      case Q_q_from_G:
         value = 0.;
         break;
      case Q_G_from_G:
         value = 1.;
         break;
      default: 
         errf(-1, "QSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: QSplitRegularIntEps0(SplittingFunction sf) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = - 1.5;
         break;
      case Q_G_from_q:
         value = - 1.5;
         break;
      case Q_q_from_G:
         value = 2. / 3.;
         break;
      case Q_G_from_G:
         value = - 11. / 6.;
         break;
      default: 
         errf(-1, "QSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: QSplitRegularIntEps1(SplittingFunction sf) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = - 6.5; //+2./3.*PiSquared; // regular part only!
         break;
      case Q_G_from_q:
         value = - 6.5; //+2./3.*PiSquared;
         break;
      case Q_q_from_G:
         value = 23. / 9.;
         break;
      case Q_G_from_G:
         value = -67. / 9.; //+2./3.*PiSquared;
         break;
      default: 
         errf(-1, "QSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: PSplitCounterTermDelta(SplittingFunction sf, double lower) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = cf * (1.5 + 2. * log(1. - lower));
         break;
      case Q_G_from_q:
         value = 0.;
         break;
      case Q_q_from_G:
         value = 0.;
         break;
      case Q_G_from_G:
         value = 2. * ncolour * (11. / 12. + log(1. - lower));
         break;
      default: 
         errf(-1, "PSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: PSplitCounterTermSubtracted(SplittingFunction sf, double u) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = cf * 2./ (1. - u);
         break;
      case Q_G_from_q:
         value = 0;
         break;
      case Q_q_from_G:
         value = 0;
         break;
      case Q_G_from_G:
         value = 2. * ncolour / (1. - u);
         break;
      default: 
         errf(-1, "PSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: PSplitCounterTermRegular(SplittingFunction sf, double u) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = cf * (- 1. - u);
         break;
      case Q_G_from_q:
         value = cf * (2. / u - 2. + u);
         break;
      case Q_q_from_G:
         value = 0.5 * (1. - 2. * u * (1. - u));
         break;
      case Q_G_from_G:
         value = 2 * ncolour * (1. / u + u * (1. - u) - 2.);
         break;
      default: 
         errf(-1, "PSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------

double QCD :: PSplitCounterTermDeltaNf(SplittingFunction sf) {

   double value;

   switch (sf) {
      case Q_q_from_q:
         value = 0.;
         break;
      case Q_G_from_q:
         value = 0.;
         break;
      case Q_q_from_G:
         value = 0.;
         break;
      case Q_G_from_G:
         value = - 1. / 3.;
         break;
      default: 
         errf(-1, "PSplit: tried to evaluate for unknown splitting function");
         value = 0.; // :_TAWM_:
   }

   return value;
}

// ----------------------------------------------------------------------------
//
// --> electroweak couplings, etc.
//
// ----------------------------------------------------------------------------

ElectroWeak::ElectroWeak() {
   defined=FALSE;
}

// ----------------------------------------------------------------------------

ElectroWeak::~ElectroWeak() {
   Delete();
}

// ----------------------------------------------------------------------------

void ElectroWeak::Define(Global *globalIn) {

   if (!defined) { 

      defined=TRUE;

      global=globalIn;

      alphaEM_Server = new AlphaEM(global);
      alphaEM_Server -> Define();
   }
}

// ----------------------------------------------------------------------------

void ElectroWeak::Delete() {

   if (defined) {

      defined = FALSE;
      delete alphaEM_Server;
   }
}

// ----------------------------------------------------------------------------
//
// --> Theory Evaluation
//
// ----------------------------------------------------------------------------

TheoryEvaluation::TheoryEvaluation() {

   defined = FALSE;
}

// ----------------------------------------------------------------------------

TheoryEvaluation::~TheoryEvaluation() {

   Delete();
}

// ----------------------------------------------------------------------------

void TheoryEvaluation::Define() {

   if (!defined) {

      defined = TRUE;

      QCD :: Define(this);   

      ElectroWeak :: Define(this);   

      GeV_2IntoPb = 3.8937966e8;

      cFree  = new ContributionFreeList(this, 21);
      cFree -> DefineName("Contribution free list");

      caFree = new ContributionArrayFreeList(this, 101, 21);
      caFree -> DefineName("ContributionArray free list");

      pFree = new ParticleFreeList(this);
      pFree -> DefineName("Particle free list");

      eFree = new EventFreeListArray(this, 0, 10);
      eFree -> DefineName("Event free list");

//      pFree -> SetTrapCreate(51);
//      pFree -> SetTrapGetNoCreate(51);
//      pFree -> SetTrapReturn(51);
   }
}

// ----------------------------------------------------------------------------

void TheoryEvaluation::Delete() {

   if (defined) {

      defined = FALSE;

      delete cFree;
      delete caFree;

      delete eFree;  // have to destruct event free list first, because
                     // events contain particles...
      delete pFree;
   }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

#if 0
// ++++++++++++++++++++++++++ GARBAGE +++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++ END OF GARBAGE +++++++++++++++++++++++++++++++++++++++++
#endif
// ============================================================================
//
// --> container classes
//
// file:              container.cc
// created:           10.06.1997
// last modification: 16.11.1997
//
// ============================================================================

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "global.h"  // for errf()

#include "container.h"

// ============================================================================
//
// ---> global variables
//
// ============================================================================

long uCount = 0;

// ============================================================================
//
// ---> instantiations
//
// ============================================================================

instantiate_Template_1(Array, int)
instantiate_Template_1(Array, double)

instantiate_Template_1(TwoD_Array, double)

// ============================================================================
//
// ---> classes
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> Array of doubles with some additional methods
//
// ----------------------------------------------------------------------------

Array_double :: Array_double() {
}

// ----------------------------------------------------------------------------

Array_double :: Array_double(int rminIn, int rmaxIn)
   : Array <double> (rminIn, rmaxIn) { 
}

// ----------------------------------------------------------------------------

Array_double :: Array_double(int rsize) 
   : Array <double> (rsize) {
}

// ----------------------------------------------------------------------------

Array_double :: ~Array_double() {
}

// ----------------------------------------------------------------------------

void Array_double :: SetToZero() {
   for (register int i = rmin; i <= rmax; ++i)
       data[i] = 0.;
}

// ----------------------------------------------------------------------------

void Array_double :: LinearScale(double smin, double smax) {
   if (rmin == rmax)
      data[rmin] = smin;
   else
      for (register int i = rmin; i <= rmax; ++i)
          data [i] = smin + ((double) (i-rmin)) / (rmax-rmin) * (smax-smin);
}

// ----------------------------------------------------------------------------

void Array_double :: LogScale(double smin, double smax) {
   double lsmin = log(smin);
   double lsmax = log(smax);
   if (rmin == rmax)
      data[rmin] = smin;
   else
      for (register int i = rmin; i <= rmax; ++i)
          data [i] = exp(  lsmin 
                         + ((double) (i-rmin)) / (rmax-rmin) * (lsmax-lsmin)
                        );
}

// ----------------------------------------------------------------------------

// returns TRUE if ok, else FALSE.
// if FALSE, the location is -1 if under and 1 if over,
// and 0 if only one point is available
// if TRUE, the location is such that the lower bound is assumed to 
// be included in the bin

int Array_double :: FindBinLocation(double x, int &loc) {

   if (data[rmin] <= data[rmax]) {
      if (x < data[rmin]) {
         loc = -1;
         return FALSE;
      } 
      if (x >= data[rmax]) {
         loc = 1;
         return FALSE;
      }
   } else {
      if (x < data[rmax]) {
         loc = -1;
         return FALSE;
      } 
      if (x >= data[rmin]) {
         loc = 1;
         return FALSE;
      }
   }

   if (rmin == rmax) {
      loc = 0;
      return FALSE;
   }
   
   register int l;

   if (data[rmin] <= data[rmax]) {
      l = rmin;
      while (l <= rmax-2 && x >= data[l+1])
         ++l;
      loc = l;
   } else {
      l = rmax;
      while (l >= rmin+2 && x > data[l-1])
         --l;
      loc = l-1;
   }

   return TRUE;
}

// ----------------------------------------------------------------------------

void Array_double :: Print(Mstream& s_out, int i0, int i1) {

   for (register int i = i0; i <= i1; ++i)
       s_out << FS("%3d -> ", i)
             << FS("%30.22e\n", data[i]);
}

// ----------------------------------------------------------------------------

void Array_double :: Print(FILE* out, int i0, int i1) {

   Mstream s_out;
   s_out.Associate(out); 

   Print(s_out, i0, i1);
}

// ----------------------------------------------------------------------------
//
// --> Array of ints with some additional methods
//
// ----------------------------------------------------------------------------

Array_int :: Array_int() {
}

// ----------------------------------------------------------------------------

Array_int :: Array_int(int rminIn, int rmaxIn)
   : Array <int> (rminIn, rmaxIn) { 
}

// ----------------------------------------------------------------------------

Array_int :: Array_int(int rsize) 
   : Array <int> (rsize) {
}

// ----------------------------------------------------------------------------

Array_int :: ~Array_int() {
}

// ----------------------------------------------------------------------------

void Array_int :: SetToZero() {
   for (register int i = rmin; i <= rmax; ++i)
       data[i] = 0;
}

// ----------------------------------------------------------------------------

void Array_int :: Print(Mstream& s_out, int i0, int i1) {

   for (register int i = i0; i <= i1; ++i)
       s_out << FS("%3d -> ", i)
             << FS("%8d\n", data[i]);
}

// ----------------------------------------------------------------------------

void Array_int :: Print(FILE* out, int i0, int i1) {

   Mstream s_out;
   s_out.Associate(out);

   Print(s_out, i0, i1);
}

// ----------------------------------------------------------------------------
//
// --> two-dimensional array of doubles with some additional methods
//
// ----------------------------------------------------------------------------

TwoD_Array_double :: TwoD_Array_double() {
}

// ----------------------------------------------------------------------------

TwoD_Array_double :: TwoD_Array_double(int r0minIn, int r0maxIn, 
                                       int r1minIn, int r1maxIn
                                      )
   : TwoD_Array <double> (r0minIn, r0maxIn,
                          r1minIn, r1maxIn
                         ) { 
}

// ----------------------------------------------------------------------------

TwoD_Array_double :: ~TwoD_Array_double() {
}

// ----------------------------------------------------------------------------

void TwoD_Array_double :: SetToZero() {
   for (register int i = r0min; i <= r0max; ++i)
       for (register int j = r1min; j <= r1max; ++j)
           data[i][j] = 0.;
}

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
// ============================================================================
//
// --> some tiny utility routines
//
// file:              utility.cc
// created:           16.10.1997
// last modification: 16.10.1997
//
// ============================================================================

#include <string.h>

// ============================================================================
//
// definitions
//
// ============================================================================

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================
    
// ============================================================================
//
// --> class-unrelated functions
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// copy a string, at most `n' non-'\0' characters
//
// returns +- the number of copied non-'\0' characters if ok, <0 on error
//
// ----------------------------------------------------------------------------

int copyStringSafe(char* destination, const char* source, int n) {

   register int i = 0;

   while (*source != '\0' && i < n)
         destination[i++] = *source++;
   destination[i] = '\0';

   return (*source == '\0') ? i : -i;
}

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
// ============================================================================
//
// --> Interface of FORTRAN routines to C++ for NaN
//
// file:              nan_cc.cc
// created:           16.10.1997
// last modification: 11.12.1997
//
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fortran.h"

#include "nan.h"

// ----------------------------------------------------------------------------

// wrapper for FORTRAN function for the detection of `NaN's

   PROTOCCALLSFFUN1(INT, INTNAN8, intnan8,
                    DOUBLE
   )

#define intnan8_Macro(arg1) \
        CCALLSFFUN1(INTNAN8, intnan8, \
                    DOUBLE, \
                    arg1 \
        )

// the following is the function (with C linkage!) to be called.
        
int intnan8_C(double x) { 
   return intnan8_Macro(x);
}
 
// ----------------------------------------------------------------------------

// wrapper for FORTRAN function for the detection of `NaN's

   PROTOCCALLSFFUN1(INT, INTNAN4, intnan4,
                    FLOAT
   )

#define intnan4_Macro(arg1) \
        CCALLSFFUN1(INTNAN4, intnan4, \
                    FLOAT, \
                    arg1 \
        )

// the following is the function (with C linkage!) to be called.
        
int intnan4_C(float x) { 
   return intnan4_Macro(x);
}
 
// ----------------------------------------------------------------------------

// dummy to avoid warning messages (static functions...)

void dummyfunctionNaN() {
   static char c[]="a";
   c2fstrv(c,c,0,0);
   f2cstrv(c,c,0,0);
   kill_trailing(c,'a');
   vkill_trailing(c,0,0,'a');
   num_elem(c,0,0,0);
}

// ============================================================================
//        
// --> End of file.
//      
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

// DIS tree matrix elements
// Created on {1997, 12, 16, 13, 18, 51}
#include <math.h>
#include <stdio.h>
#include "mathdefs.h"
#include "enum.h"
#include "global.h"
#include "proc.h"
#include "me.h"
void DISMatrixElement::FillLorentzLqXq(
                        ) {
term[0][0] = (
8*Q2inv*si0 + 16*s0k*sik
);
}
void DISMatrixElement::FillLorentzLqXqg(
                        ) {
term[0][0] = (
(16*Q2inv*si1 + 32*s1k*sik)/s01
);
term[1][0] = (
(Q2inv*si0*(-16*s01 + 16*si0 + 16*si1) + 
     s0k*(-16*s1k*si0 + 16*s0k*si1) + 
     sik*(16*s1k*si0 + s0k*(-16*s01 + 32*si0 + 16*si1) - 
        16*s01*sik))/(s01*si1)
);
term[1][1] = (
(16*Q2inv*s01 + 32*s0k*s1k)/si1
);
}
void DISMatrixElement::FillLorentzLgXqq(
                        ) {
term[0][0] = (
(16*Q2inv*si1 + 32*s1k*sik)/si0
);
term[1][0] = (
(s0k*(s1k*(32*s01 - 16*si0 - 16*si1) + 16*s0k*si1) + 
     (-16*s01*s0k - 16*s01*s1k)*sik + 
     Q2inv*(-16*s01*si0 - 16*s01*si1 + 16*POWER2(s01)) + 
     16*si0*POWER2(s1k))/(si0*si1)
);
term[1][1] = (
(16*Q2inv*si0 + 32*s0k*sik)/si1
);
}
void DISMatrixElement::FillLorentzLqXqgg(
                        ) {
term[0][0] = (
(Q2inv*((32*s01 + 32*s12)*si0 - 32*s02*si1 + 
        (32*s01 + 32*s12)*si2) + 
     (s0k*(64*s01 + 64*s12) - 64*s02*s1k + 
        (64*s01 + 64*s12)*s2k)*sik)/
   (s01*POWER2(s01 + s02 + s12))
);
term[1][0] = (
(Q2inv*(32*s12*si0 + (-32*s01 - 32*s02)*si1 + 
        (-32*s01 - 32*s02)*si2) + 
     (64*s0k*s12 + (-64*s01 - 64*s02)*s1k + 
        (-64*s01 - 64*s02)*s2k)*sik)/
   (s01*s02*(s01 + s02 + s12))
);
term[2][0] = (
(32*s0k*s1k*si2 + s1k*(s2k*(-32*si0 - 32*si1) + 
        32*s1k*si2) + Q2inv*
      (32*si0*si1 + si1*
         (-32*s02 - 32*s12 + 32*si1 + 32*si2)) + 
     sik*(s2k*(32*s01 + 32*si1) + s0k*(-32*s12 + 32*si1) + 
        s1k*(-32*s12 + 32*si0 + 64*si1 + 32*si2) + 
        (-32*s01 - 32*s12)*sik))/(s01*(s01 + s02 + s12)*si2)
);
term[3][0] = (
(s1k*(-32*s02*s1k*si0 - 32*s02*s2k*si0) + 
     s0k*(-32*s0k*s12*si1 + (32*s01 + 32*s02)*s2k*si1 + 
        s1k*(32*s12*si0 + 32*s02*si1 - 32*s01*si2)) + 
     sik*(s1k*(-32*s01*si0 + 32*s02*si1 - 32*s01*si2) + 
        s2k*(s01*(-32*s01 - 32*s02) + 
           (32*s01 + 64*s02 + 32*s12)*si0 + 32*s02*si1 - 
           32*s01*si2) + 
        s0k*(32*s01*s12 + (-32*s01 - 32*s12)*si1 - 
           32*s01*si2) + s01*(32*s01 + 32*s02 + 32*s12)*sik)
       + Q2inv*(si2*(-16*s01*s02 - 16*s01*si2) + 
        si0*(-16*s02*s12 + (-32*s01 - 16*s12)*si1 + 
           (32*s02 + 16*s12)*si2) + 
        si1*(16*s02*si1 + (-16*s01 + 16*s02)*si2 + 
           16*POWER2(s02))))/(s01*s02*(s01 + s02 + s12)*si1)
);
term[4][0] = (
(s1k*(32*s02*s1k*si2 + s2k*
         (32*s12*si0 - 32*s02*si1 - 32*s01*si2)) + 
     s0k*(s0k*(si1*(32*s02 + 32*s12 - 32*si1 - 64*si2) + 
           si0*(-32*si1 - 32*si2) + (32*s02 - 32*si2)*si2)\
         + s1k*(si0*(-32*s02 - 32*s12 + 32*si0 + 32*si1) - 
           32*si1*si2 + (32*s01 + 64*s02 - 32*si2)*si2) + 
        s2k*(si1*(-32*s01 - 64*s02 - 32*s12 + 32*si1 + 
              32*si2) + si0*
            (-32*s02 + 32*si0 + 64*si1 + 32*si2))) + 
     sik*(s2k*(s01*(32*s01 + 32*s02) + 
           si0*(32*s02 + 32*s12 - 32*si0 - 32*si1) - 
           64*s01*si1 - 32*s01*si2) + 
        s1k*(-32*s01*s02 + s02*(-32*s02 - 32*s12) + 
           si0*(32*s01 + 64*s02 - 32*si0 - 32*si1) + 
           64*s02*si1 + 32*s02*si2) + 
        s0k*(s01*(-32*s02 - 32*s12) + 
           s02*(-32*s02 - 32*s12) + 
           si0*(32*s01 + 96*s02 - 64*si0 - 96*si1 - 
              32*si2) + si1*
            (64*s01 + 128*s02 + 32*s12 - 32*si1 - 32*si2) + 
           (64*s01 + 64*s02 + 32*s12)*si2) + 
        (s01*(-32*s01 - 64*s02) + s02*(-32*s02 - 32*s12) + 
           (32*s01 + 32*s02)*si0 + (32*s01 + 32*s02)*si1)*
         sik) + Q2inv*(s01*(16*s02 + 16*s12)*si2 + 
        si1*(s02*(-16*s02 - 16*s12) + 16*s02*si1 - 
           16*s01*si2) + 
        si0*(-32*s01*s02 + s02*(-32*s02 - 48*s12) + 
           si0*(32*s01 + 64*s02 + 64*s12 - 32*si0 - 
              64*si1 - 32*si2) + 
           si1*(64*s02 + 48*s12 - 32*si1 - 32*si2) - 
           32*s01*si2 - 16*POWER2(s12))) + 
     32*s01*si1*POWER2(s2k))/
   (s01*(s01 + s02 + s12)*(s12 - si1 - si2)*si2)
);
term[5][0] = (
(Q2inv*si0*(s01*(-32*s01 - 32*s02 - 64*s12) + 
        si0*(64*s01 + 32*s02 + 64*s12 - 32*si0 - 64*si1 - 
           32*si2) + si1*
         (64*s01 + 32*s02 + 64*s12 - 32*si1 - 32*si2) + 
        32*s01*si2) + s0k*
      (s1k*si0*(-32*s01 + 32*si0 + 32*si1) + 
        s2k*si0*(-32*s01 + 32*si0 + 32*si1) + 
        s0k*(si0*(-32*si1 - 32*si2) + 
           si1*(32*s01 - 32*si1 - 32*si2) + 32*s01*si2)) + 
     sik*(s1k*si0*(32*s01 - 32*si0 - 32*si1) + 
        s2k*si0*(32*s01 - 32*si0 - 32*si1) + 
        s0k*(s01*(-32*s01 - 32*s02) + 
           si0*(96*s01 + 32*s02 - 64*si0 - 96*si1 - 
              32*si2) + si1*
            (64*s01 + 32*s02 - 32*si1 - 32*si2) + 32*s01*si2
           ) + (s01*(-32*s01 - 32*s02) + 
           (32*s01 + 32*s02)*si0 + (32*s01 + 32*s02)*si1)*
         sik))/(s01*(s01 + s02 + s12)*si1*(s12 - si1 - si2))
);
term[6][0] = (
(Q2inv*((s02*(-16*s02 - 16*s12) + s01*(-16*s01 + 16*s12))*
         si1 + (s01*(-64*s01 - 64*s02 - 32*s12) + 
           s02*(-32*s02 - 32*s12))*si2 + 
        si0*(16*s02*s12 + s01*(-16*s01 + 16*s02 + 32*s12) + 
           16*POWER2(s12))) + 
     sik*((s02*(-32*s02 - 32*s12) + s01*(-32*s01 + 32*s12))*
         s1k + (s01*(-128*s01 - 128*s02 - 64*s12) + 
           s02*(-64*s02 - 64*s12))*s2k + 
        s0k*(32*s02*s12 + s01*(-32*s01 + 32*s02 + 64*s12) + 
           32*POWER2(s12))))/
   (s01*s12*POWER2(s01 + s02 + s12))
);
term[7][0] = (
(s1k*(32*s02*s1k*si2 + s2k*
         (32*s12*si0 - 32*s02*si1 - 32*s01*si2)) + 
     s0k*(s2k*((16*s01 - 16*s02 + 48*s12)*si0 + 
           (-16*s01 - 48*s02 + 16*s12)*si1) + 
        s0k*((-16*s01 + 16*s02 + 16*s12)*si1 + 
           (-16*s01 + 16*s02 - 48*s12)*si2) + 
        s1k*((16*s01 - 16*s02 - 16*s12)*si0 + 
           (16*s01 + 48*s02 - 48*s12)*si2)) + 
     Q2inv*(si2*(s01*(16*s02 + 24*s12) - 8*s01*si2) + 
        si1*(s02*(-16*s02 - 24*s12) + 24*s02*si1 + 
           (-24*s01 + 8*s02)*si2) + 
        si0*(s02*(-16*s02 - 48*s12) + 
           s01*(16*s01 - 16*s12) + 
           (-16*s01 + 16*s02 - 16*s12)*si0 + 
           (48*s02 - 8*s12)*si1 + 
           (16*s01 + 32*s02 + 24*s12)*si2 - 24*POWER2(s12)))
       + sik*(s2k*(s01*(16*s01 + 16*s02 - 32*s12) + 
           (32*s02 + 32*s12)*si0 - 48*s01*si1 - 16*s01*si2)\
         + (s02*(-32*s02 - 32*s12) + 
           s01*(-16*s01 - 48*s02 + 16*s12))*sik + 
        s1k*(-16*s01*s02 + (16*s01 + 48*s02 - 16*s12)*si0 + 
           48*s02*si1 + 16*s02*si2 - 16*POWER2(s02)) + 
        s0k*(s01*(16*s01 + 16*s12) + 
           s02*(-16*s02 + 16*s12) + 
           (-32*s01 + 32*s02 - 32*s12)*si0 + 
           (-16*s01 + 48*s02)*si1 + 
           (32*s01 + 32*s02 + 16*s12)*si2 + 16*POWER2(s12)))
       + 32*s01*si1*POWER2(s2k))/
   (s01*s12*(s01 + s02 + s12)*(s12 - si1 - si2))
);
term[1][1] = (
(Q2inv*((32*s02 + 32*s12)*si0 + (32*s02 + 32*s12)*si1 - 
        32*s01*si2) + (s0k*(64*s02 + 64*s12) + 
        (64*s02 + 64*s12)*s1k - 64*s01*s2k)*sik)/
   (s02*POWER2(s01 + s02 + s12))
);
term[2][1] = (
(-32*s01*s1k*s2k*si0 + s0k*
      (-32*s0k*s12*si2 + (32*s01 + 32*s02)*s1k*si2 + 
        s2k*(32*s12*si0 - 32*s02*si1 + 32*s01*si2)) + 
     Q2inv*(si1*(-16*s01*s02 - 16*s02*si1 + 
           (16*s01 - 16*s02)*si2) + 
        si0*(-16*s01*s12 + (32*s01 + 16*s12)*si1 + 
           (-32*s02 - 16*s12)*si2) + 
        si2*(16*s01*si2 + 16*POWER2(s01))) + 
     sik*(s2k*(-32*s02*si0 - 32*s02*si1 + 32*s01*si2) + 
        s0k*(32*s02*s12 - 32*s02*si1 + 
           (-32*s02 - 32*s12)*si2) + 
        (32*s01*s02 + s02*(32*s02 + 32*s12))*sik + 
        s1k*(-32*s01*s02 + (64*s01 + 32*s02 + 32*s12)*si0 - 
           32*s02*si1 + 32*s01*si2 - 32*POWER2(s02))) - 
     32*s01*si0*POWER2(s2k))/(s01*s02*(s01 + s02 + s12)*si2)
);
term[3][1] = (
(32*s0k*s2k*si1 + s1k*s2k*(-32*si0 - 32*si2) + 
     Q2inv*(32*si0*si2 + 32*si1*si2 + 
        si2*(-32*s01 - 32*s12 + 32*si2)) + 
     sik*(s1k*(32*s02 + 32*si2) + s0k*(-32*s12 + 32*si2) + 
        s2k*(-32*s12 + 32*si0 + 32*si1 + 64*si2) + 
        (-32*s02 - 32*s12)*sik) + 32*si1*POWER2(s2k))/
   (s02*(s01 + s02 + s12)*si1)
);
term[4][1] = (
(Q2inv*si0*(-32*s01*s02 + s02*(-32*s02 - 64*s12) + 
        si0*(32*s01 + 64*s02 + 64*s12 - 32*si0 - 32*si1 - 
           64*si2) + si1*(32*s02 - 32*si2) + 
        (32*s01 + 64*s02 + 64*s12 - 32*si2)*si2) + 
     s0k*(s1k*si0*(-32*s02 + 32*si0 + 32*si2) + 
        s2k*si0*(-32*s02 + 32*si0 + 32*si2) + 
        s0k*(si1*(32*s02 - 32*si2) + 
           si0*(-32*si1 - 32*si2) + (32*s02 - 32*si2)*si2))\
      + sik*(s1k*si0*(32*s02 - 32*si0 - 32*si2) + 
        s2k*si0*(32*s02 - 32*si0 - 32*si2) + 
        sik*(-32*s01*s02 + (32*s01 + 32*s02)*si0 + 
           (32*s01 + 32*s02)*si2 - 32*POWER2(s02)) + 
        s0k*(-32*s01*s02 + 
           si0*(32*s01 + 96*s02 - 64*si0 - 32*si1 - 
              96*si2) + si1*(32*s02 - 32*si2) + 
           (32*s01 + 64*s02 - 32*si2)*si2 - 32*POWER2(s02)))
     )/(s02*(s01 + s02 + s12)*(s12 - si1 - si2)*si2)
);
term[5][1] = (
(s1k*(32*s02*s1k*si2 + s2k*
         (32*s12*si0 - 32*s02*si1 - 32*s01*si2)) + 
     s0k*(s0k*(si1*(32*s01 - 32*si1 - 64*si2) + 
           si0*(-32*si1 - 32*si2) + 
           (32*s01 + 32*s12 - 32*si2)*si2) + 
        s2k*(si1*(64*s01 + 32*s02 - 32*si1 - 32*si2) + 
           si0*(-32*s01 - 32*s12 + 32*si0 + 32*si2)) + 
        s1k*(32*si1*si2 + 
           si2*(-64*s01 - 32*s02 - 32*s12 + 32*si2) + 
           si0*(-32*s01 + 32*si0 + 32*si1 + 64*si2))) + 
     sik*(s2k*(s01*(-32*s01 - 32*s02 - 32*s12) + 
           32*s01*si1 + si0*
            (64*s01 + 32*s02 - 32*si0 - 32*si2) + 64*s01*si2
           ) + s0k*(s01*(-32*s01 - 32*s02 - 32*s12) - 
           32*s02*s12 + si0*
            (96*s01 + 32*s02 - 64*si0 - 32*si1 - 96*si2) + 
           si1*(64*s01 + 64*s02 + 32*s12 - 32*si2) + 
           (128*s01 + 64*s02 + 32*s12 - 32*si2)*si2) + 
        sik*(s01*(-32*s01 - 64*s02 - 32*s12) + 
           (32*s01 + 32*s02)*si0 + (32*s01 + 32*s02)*si2 - 
           32*POWER2(s02)) + 
        s1k*(32*s01*s02 - 32*s02*si1 + 
           si0*(32*s01 + 32*s12 - 32*si0 - 32*si2) - 
           64*s02*si2 + 32*POWER2(s02))) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s12) + 16*s01*si2) + 
        si1*(16*s01*s02 + 16*s02*s12 - 16*s02*si2) + 
        si0*(s01*(-32*s01 - 32*s02 - 48*s12) + 
           si0*(64*s01 + 32*s02 + 64*s12 - 32*si0 - 
              32*si1 - 64*si2) + si1*(-32*s02 - 32*si2) + 
           (64*s01 + 48*s12 - 32*si2)*si2 - 16*POWER2(s12)))
       + 32*s01*si1*POWER2(s2k))/
   (s02*(s01 + s02 + s12)*si1*(s12 - si1 - si2))
);
term[6][1] = (
(sik*((s02*(128*s02 + 64*s12) + 
           s01*(64*s01 + 128*s02 + 64*s12))*s1k + 
        (s02*(32*s02 - 32*s12) + s01*(32*s01 + 32*s12))*
         s2k + s0k*(s02*(32*s02 - 64*s12) + 
           s01*(-32*s02 - 32*s12) - 32*POWER2(s12))) + 
     Q2inv*((s02*(64*s02 + 32*s12) + 
           s01*(32*s01 + 64*s02 + 32*s12))*si1 + 
        (s02*(16*s02 - 16*s12) + s01*(16*s01 + 16*s12))*
         si2 + si0*(s02*(16*s02 - 32*s12) + 
           s01*(-16*s02 - 16*s12) - 16*POWER2(s12))))/
   (s02*s12*POWER2(s01 + s02 + s12))
);
term[7][1] = (
(s1k*(-32*s02*s1k*si2 + s2k*
         (-32*s12*si0 + 32*s02*si1 + 32*s01*si2)) + 
     s0k*(s2k*((16*s01 - 16*s02 + 16*s12)*si0 + 
           (-48*s01 - 16*s02 + 48*s12)*si1) + 
        s0k*((-16*s01 + 16*s02 + 48*s12)*si1 + 
           (-16*s01 + 16*s02 - 16*s12)*si2) + 
        s1k*((16*s01 - 16*s02 - 48*s12)*si0 + 
           (48*s01 + 16*s02 - 16*s12)*si2)) + 
     sik*(s2k*(s01*(16*s01 + 16*s02) + 
           (-48*s01 - 16*s02 + 16*s12)*si0 - 16*s01*si1 - 
           48*s01*si2) + 
        s1k*(-16*s01*s02 + s02*(-16*s02 + 32*s12) + 
           (-32*s01 - 32*s12)*si0 + 16*s02*si1 + 48*s02*si2)
          + (s02*(16*s02 - 16*s12) + 
           s01*(32*s01 + 48*s02 + 32*s12))*sik + 
        s0k*(s01*(16*s01 - 16*s12) + 
           s02*(-16*s02 - 16*s12) + 
           (-32*s01 + 32*s02 + 32*s12)*si0 + 
           (-32*s01 - 32*s02 - 16*s12)*si1 + 
           (-48*s01 + 16*s02)*si2 - 16*POWER2(s12))) + 
     Q2inv*(si2*(s01*(16*s01 + 24*s12) - 24*s01*si2) + 
        si1*(-16*s01*s02 - 24*s02*s12 + 8*s02*si1 + 
           (-8*s01 + 24*s02)*si2) + 
        si0*(s02*(-16*s02 + 16*s12) + 
           s01*(16*s01 + 48*s12) + 
           (-16*s01 + 16*s02 + 16*s12)*si0 + 
           (-32*s01 - 16*s02 - 24*s12)*si1 + 
           (-48*s01 + 8*s12)*si2 + 24*POWER2(s12))) - 
     32*s01*si1*POWER2(s2k))/
   (s02*s12*(s01 + s02 + s12)*(s12 - si1 - si2))
);
term[2][2] = (
(32*Q2inv*s12 + 64*s1k*s2k)/(s01*si2)
);
term[3][2] = (
(s0k*(s2k*si0*(32*s01 - 32*si0 + 32*si1) + 
        s1k*si0*(32*s02 - 32*si0 + 32*si2) + 
        s0k*si0*(-32*s12 + 32*si1 + 32*si2)) + 
     Q2inv*si0*(32*s01*s02 - 32*s01*si2 + 
        si1*(-32*s02 + 32*si2) + 
        si0*(-32*s01 - 32*s02 + 32*si0 + 32*si1 + 32*si2))\
      + sik*(s2k*si0*(32*s01 + 32*si0 + 32*si1) + 
        s1k*si0*(32*s02 + 32*si0 + 32*si2) + 
        s0k*si0*(-32*s01 - 32*s02 - 64*s12 + 64*si0 + 
           32*si1 + 32*si2) + 
        (-32*s01 - 32*s02 - 32*s12)*si0*sik) - 
     64*s1k*s2k*POWER2(si0))/(s01*s02*si1*si2)
);
term[4][2] = (
(s1k*s2k*(32*s02 - 32*si0) + 
     Q2inv*(32*s01*s02 + s02*(32*s02 + 32*s12) - 
        32*s02*si0 - 32*s02*si1) + 
     s0k*(s2k*(32*s01 + 64*s02 + 32*s12 - 32*si0) + 
        s1k*(32*s02 + 32*si2) + s0k*(-32*s12 + 32*si2)) + 
     (s0k*(-32*s02 - 32*s12) + 32*s01*s2k)*sik - 
     32*s01*POWER2(s2k))/(s01*(s12 - si1 - si2)*si2)
);
term[5][2] = (
(s1k*(32*s1k*si0*si2 + 32*s2k*si0*si2) + 
     s0k*(s0k*si1*(-32*s12 + 32*si1 + 32*si2) + 
        s1k*(32*s02*si1 - 32*si0*si1 - 32*s01*si2) + 
        s2k*(-32*s01*si2 + si1*(32*s02 + 32*si1 + 32*si2) + 
           si0*(-32*s12 + 32*si1 + 64*si2))) + 
     sik*(s0k*(32*s01*s12 + 
           (-32*s01 - 32*s02 - 32*s12)*si1) + 
        s2k*(-32*s01*si1 - 32*s01*si2) + 
        s1k*(-32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        32*s01*s12*sik) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s02) - 16*s01*si2) + 
        si0*(16*s01*s12 - 16*s02*s12 - 32*s01*si1 + 
           (32*s02 + 16*s12)*si2) + 
        si1*(16*s01*s02 + 16*s02*si2 + 16*POWER2(s02))))/
   (s01*si1*si2*(-s12 + si1 + si2))
);
term[6][2] = (
(s1k*(s2k*((16*s01 + 16*s02 + 16*s12)*si0 + 32*s01*si1) - 
        32*s01*s1k*si2) + 
     s0k*(s0k*(-32*s01 - 16*s02 - 48*s12)*si2 + 
        (-32*s01 + 16*s02 - 16*s12)*s1k*si2 + 
        s2k*((32*s01 + 16*s02 + 48*s12)*si0 + 
           (16*s01 - 32*s02)*si1 + 32*s01*si2)) + 
     sik*(s2k*(s01*(-32*s01 + 16*s02 + 16*s12) + 
           (-16*s01 - 32*s02)*si0 + 
           (-32*s01 - 32*s02)*si1 + 16*s01*si2) + 
        s1k*(32*s02*s12 + s01*(16*s02 + 32*s12) + 
           (-32*s01 - 16*s02 - 16*s12)*si0 + 
           (-64*s01 - 32*s02)*si1 + (-32*s01 - 16*s02)*si2)\
         + (s01*(32*s01 + 32*s02 + 16*s12) + 
           s02*(32*s02 + 32*s12))*sik + 
        s0k*(s01*(32*s02 + 16*s12) + 
           s02*(16*s02 + 32*s12) + 
           (-64*s01 - 32*s02 - 32*s12)*si0 + 
           (-32*s01 - 16*s02 + 16*s12)*si1 + 
           (-16*s01 - 32*s02 - 16*s12)*si2 - 16*POWER2(s12))
        ) + Q2inv*(si2*(s01*(-16*s02 - 16*s12) + 
           8*s01*si2) + si1*
         (s01*(16*s02 + 32*s12) + (-32*s01 - 16*s02)*si1 + 
           (-32*s01 - 24*s02)*si2) + 
        si0*(s01*(32*s02 + 16*s12) + 
           s02*(16*s02 + 32*s12) + 
           (-32*s01 - 16*s02 - 16*s12)*si0 + 
           (-32*s01 - 16*s02)*si1 + 
           (-16*s01 - 32*s02 - 8*s12)*si2 + 16*POWER2(s12)))
       - 32*s01*si0*POWER2(s2k))/
   (s01*s12*(s01 + s02 + s12)*si2)
);
term[7][2] = (
(s1k*(32*s1k*si0*si2 + s2k*
         (si0*(16*s12 - 16*si1 - 16*si2) + 32*s02*si2)) + 
     s0k*(s2k*(si0*(16*s12 - 16*si1 - 32*si2) + 
           si1*(16*s01 + 32*s02 + 32*s12 - 16*si2) + 
           (32*s01 + 64*s02 + 32*s12)*si2) + 
        s1k*(si1*(32*s02 - 16*si2) + 
           si0*(-32*si1 - 16*si2) + 
           si2*(-16*s01 + 32*s02 + 16*s12 + 32*si2)) + 
        s0k*(si2*(-16*s12 + 32*si2) + 
           si1*(-32*s12 + 32*si1 + 32*si2))) + 
     sik*(s2k*(-16*s01*s12 - 16*s01*si1 + 32*s01*si2) + 
        s1k*(32*s02*si1 + (-32*s01 - 16*s02)*si2 + 
           si0*(-48*s12 + 16*si1 + 32*si2)) + 
        (48*s01*s12 - 16*s01*si1 - 32*s01*si2)*sik + 
        s0k*(16*s01*s12 - 16*s02*s12 + 
           (-16*s01 - 32*s02 - 16*s12)*si2 + 
           si1*(-32*s01 - 16*s02 - 32*s12 + 16*si1 + 
              32*si2) + si0*(-32*s12 + 32*si1 + 64*si2) - 
           16*POWER2(s12))) + 
     Q2inv*((s01*(-8*s01 + 32*s02 - 16*s12) + 
           s02*(32*s02 + 32*s12))*si2 + 
        si1*(24*s01*s02 + (16*s01 - 16*s02)*si2 + 
           16*POWER2(s02)) + 
        si0*(8*s01*s12 + (-16*s01 - 32*s02 - 16*s12)*si2 + 
           si1*(-32*s01 - 16*s02 - 32*s12 + 16*si1 + 
              32*si2) + si0*(-16*s12 + 16*si1 + 32*si2) + 
           16*POWER2(s12))) - 32*s01*si2*POWER2(s2k))/
   (s01*s12*(s12 - si1 - si2)*si2)
);
term[3][3] = (
(32*Q2inv*s12 + 64*s1k*s2k)/(s02*si1)
);
term[4][3] = (
(32*s1k*s2k*si0*si1 + s0k*
      (s2k*(-32*s02*si1 + 32*s01*si2 - 32*si0*si2) + 
        s0k*(32*si1*si2 + si2*(-32*s12 + 32*si2)) + 
        s1k*(si2*(32*s01 + 32*si2) + 
           si1*(-32*s02 + 32*si2) + 
           si0*(-32*s12 + 64*si1 + 32*si2))) + 
     sik*(s2k*(-32*s12*si0 - 32*s02*si1 + 32*s01*si2) + 
        s1k*(-32*s02*si1 - 32*s02*si2) + 
        s0k*(32*s02*s12 + 
           (-32*s01 - 32*s02 - 32*s12)*si2) + 32*s02*s12*sik
        ) + Q2inv*(s01*(16*s01 + 16*s02)*si2 + 
        si0*(-16*s01*s12 + 16*s02*s12 + 
           (32*s01 + 16*s12)*si1 - 32*s02*si2) + 
        si1*(-16*s01*s02 - 16*s02*si1 + 16*s01*si2 - 
           16*POWER2(s02))) + 32*si0*si1*POWER2(s2k))/
   (s02*si1*si2*(-s12 + si1 + si2))
);
term[5][3] = (
(s1k*(32*s02*s1k + s2k*(-32*s01 + 32*si0)) + 
     s0k*(s1k*(-64*s01 - 32*s02 - 32*s12 + 32*si0) + 
        s2k*(-32*s01 - 32*si1) + s0k*(32*s12 - 32*si1)) + 
     Q2inv*(s01*(-32*s01 - 32*s02 - 32*s12) + 32*s01*si0 + 
        32*s01*si2) + (s0k*(32*s01 + 32*s12) - 32*s02*s1k)*
      sik)/(s02*si1*(-s12 + si1 + si2))
);
term[6][3] = (
(s0k*(s0k*(16*s01 + 32*s02 + 48*s12)*si1 + 
        (-16*s01 + 32*s02 + 16*s12)*s2k*si1 + 
        s1k*((-16*s01 - 32*s02 - 48*s12)*si0 - 32*s02*si1 + 
           (32*s01 - 16*s02)*si2)) + 
     s1k*(32*s02*s1k*si0 + 
        s2k*((-16*s01 - 16*s02 - 16*s12)*si0 - 32*s02*si2))\
      + Q2inv*(si2*(-16*s01*s02 - 32*s02*s12 + 
           (16*s01 + 32*s02)*si2) + 
        si1*(16*s01*s02 + 16*s02*s12 - 8*s02*si1 + 
           (24*s01 + 32*s02)*si2) + 
        si0*(s01*(-16*s01 - 32*s02 - 32*s12) - 16*s02*s12 + 
           (16*s01 + 32*s02 + 16*s12)*si0 + 
           (32*s01 + 16*s02 + 8*s12)*si1 + 
           (16*s01 + 32*s02)*si2 - 16*POWER2(s12))) + 
     sik*(s1k*(-16*s01*s02 + s02*(32*s02 - 16*s12) + 
           (32*s01 + 16*s02)*si0 - 16*s02*si1 + 
           (32*s01 + 32*s02)*si2) + 
        s2k*(s01*(-16*s02 - 32*s12) - 32*s02*s12 + 
           (16*s01 + 32*s02 + 16*s12)*si0 + 
           (16*s01 + 32*s02)*si1 + (32*s01 + 64*s02)*si2) + 
        (s01*(-32*s01 - 32*s02 - 32*s12) + 
           s02*(-32*s02 - 16*s12))*sik + 
        s0k*(s01*(-16*s01 - 32*s02 - 32*s12) - 16*s02*s12 + 
           (32*s01 + 64*s02 + 32*s12)*si0 + 
           (32*s01 + 16*s02 + 16*s12)*si1 + 
           (16*s01 + 32*s02 - 16*s12)*si2 + 16*POWER2(s12)))
       + 32*s02*si1*POWER2(s2k))/
   (s02*s12*(s01 + s02 + s12)*si1)
);
term[7][3] = (
(s1k*(32*s02*s1k*si1 + s2k*
         (-32*s01*si1 + si0*(-16*s12 + 16*si1 + 16*si2))) + 
     s0k*(s0k*(si1*(16*s12 - 32*si1 - 32*si2) + 
           (32*s12 - 32*si2)*si2) + 
        s1k*((-32*s01 - 16*s02 - 32*s12)*si2 + 
           si1*(-64*s01 - 32*s02 - 32*s12 + 16*si2) + 
           si0*(-16*s12 + 32*si1 + 16*si2)) + 
        s2k*(-32*s01*si2 + 
           si1*(-32*s01 + 16*s02 - 16*s12 - 32*si1 + 
              16*si2) + si0*(16*si1 + 32*si2))) + 
     Q2inv*(s01*(-16*s01 - 24*s02)*si2 + 
        si1*(s01*(-32*s01 - 32*s02 - 32*s12) + 
           s02*(8*s02 + 16*s12) + (16*s01 - 16*s02)*si2) + 
        si0*(-8*s02*s12 + 
           si1*(32*s01 + 16*s02 + 16*s12 - 32*si2) + 
           si0*(16*s12 - 32*si1 - 16*si2) + 
           (16*s01 + 32*s02 + 32*s12 - 16*si2)*si2 - 
           16*POWER2(s12))) + 
     sik*(s2k*((16*s01 + 32*s02)*si1 + 
           si0*(48*s12 - 32*si1 - 16*si2) - 32*s01*si2) + 
        s1k*(16*s02*s12 - 32*s02*si1 + 16*s02*si2) + 
        (-48*s02*s12 + 32*s02*si1 + 16*s02*si2)*sik + 
        s0k*(16*s01*s12 - 16*s02*s12 + 
           si1*(32*s01 + 16*s02 + 16*s12 - 32*si2) + 
           si0*(32*s12 - 64*si1 - 32*si2) + 
           (16*s01 + 32*s02 + 32*s12 - 16*si2)*si2 + 
           16*POWER2(s12))) - 32*si0*si1*POWER2(s2k))/
   (s02*s12*si1*(s12 - si1 - si2))
);
term[4][4] = (
(s0k*(64*s2k*si1 + s1k*(64*s12 - 64*si2)) + 
     Q2inv*(32*s01*s12 + 32*s02*si1 - 32*s01*si2 + 
        si0*(-32*s12 + 32*si2)) + s0k*(-64*s12 + 64*si2)*sik
     )/(si2*POWER2(s12 - si1 - si2))
);
term[5][4] = (
(Q2inv*(-32*s12*si0 + (32*s01 + 32*s02)*si1 + 
        (32*s01 + 32*s02)*si2) + 
     s0k*(s1k*(64*si1 + 64*si2) + s2k*(64*si1 + 64*si2)) - 
     64*s0k*s12*sik)/(si1*si2*(-s12 + si1 + si2))
);
term[6][4] = (
(s1k*(-32*s02*s1k*si2 + s2k*
         (-32*s12*si0 + 32*s02*si1 + 32*s01*si2)) + 
     s0k*(s2k*(si0*(-16*s12 - 48*si1 - 16*si2) + 
           si1*(16*s01 + 48*s02 - 16*si1 - 16*si2)) + 
        s1k*(si0*(32*s12 - 32*si1) + 16*si1*si2 + 
           si2*(-16*s01 - 48*s02 + 32*s12 + 16*si2)) + 
        s0k*(si2*(16*s12 + 16*si2) + 
           si1*(-32*s12 + 32*si1 + 48*si2))) + 
     sik*(s2k*(48*s01*s12 + 48*s01*si1 + 
           si0*(-16*s12 + 16*si1 - 16*si2) + 16*s01*si2) + 
        s1k*(-16*s02*s12 - 48*s02*si1 + 
           si0*(48*s12 + 16*si1 - 16*si2) - 16*s02*si2) + 
        (-48*s01*s12 + 16*s02*s12 + 
           (-16*s01 - 16*s02)*si1 + (16*s01 + 16*s02)*si2)*
         sik + s0k*(16*s01*s12 + 
           si1*(-32*s01 - 48*s02 + 16*s12 + 16*si1) + 
           si0*(32*s12 + 32*si1 - 32*si2) + 
           (-32*s01 + 16*s02 + 16*s12 - 16*si2)*si2 - 
           16*POWER2(s12))) + 
     Q2inv*(s01*(-8*s01 - 24*s02 - 24*s12)*si2 + 
        si1*(8*s01*s02 + s02*(24*s02 + 24*s12) - 
           16*s02*si1 + 16*s01*si2) + 
        si0*(24*s01*s12 - 8*s02*s12 + 
           si1*(-32*s01 - 48*s02 - 48*s12 + 16*si1) + 
           si0*(16*s12 + 16*si1 - 16*si2) + 
           (-16*s01 - 16*s12 - 16*si2)*si2 + 24*POWER2(s12))
        ) - 32*s01*si1*POWER2(s2k))/
   (s12*(s01 + s02 + s12)*(s12 - si1 - si2)*si2)
);
term[7][4] = (
(s0k*(s1k*(si1*(64*s12 - 64*si1 - 128*si2) + 
           (64*s12 - 128*si2)*si2) + 
        s2k*((32*s12 - 32*si1)*si1 + (-32*s12 - 32*si2)*si2)
        ) + Q2inv*(si1*(32*s01*s12 + 16*s02*s12 + 
           (-32*s01 - 16*s02)*si1 - 64*s01*si2) + 
        si2*(32*s01*s12 - 16*s02*s12 + 
           (-64*s01 - 16*s02)*si2) + 
        si0*(si1*(16*s12 - 16*si2) + 
           si2*(32*s12 + 16*si2) - 16*POWER2(s12))) + 
     s0k*sik*(si1*(32*s12 - 32*si2) + 
        si2*(64*s12 + 32*si2) - 32*POWER2(s12)))/
   (s12*si2*POWER2(s12 - si1 - si2))
);
term[5][5] = (
(Q2inv*(32*s02*s12 - 32*s02*si1 + si0*(-32*s12 + 32*si1) + 
        32*s01*si2) + s0k*
      (s2k*(64*s12 - 64*si1) + 64*s1k*si2) + 
     s0k*(-64*s12 + 64*si1)*sik)/
   (si1*POWER2(-s12 + si1 + si2))
);
term[6][5] = (
(s1k*(32*s02*s1k*si2 + s2k*
         (32*s12*si0 - 32*s02*si1 - 32*s01*si2)) + 
     s0k*(s0k*(si1*(-16*s12 - 16*si1 - 48*si2) + 
           (32*s12 - 32*si2)*si2) + 
        s2k*(si1*(48*s01 + 16*s02 - 32*s12 - 16*si1 - 
              16*si2) + si0*(-32*s12 + 32*si2)) + 
        s1k*(16*si1*si2 + si2*(-48*s01 - 16*s02 + 16*si2) + 
           si0*(16*s12 + 16*si1 + 48*si2))) + 
     Q2inv*(si2*(s01*(-24*s01 - 8*s02 - 24*s12) + 
           16*s01*si2) + 
        si1*(24*s01*s02 + s02*(8*s02 + 24*s12) - 
           16*s02*si2) + 
        si0*(8*s01*s12 - 24*s02*s12 + 
           si1*(16*s02 + 16*s12 + 16*si1) + 
           si0*(-16*s12 + 16*si1 - 16*si2) + 
           (48*s01 + 32*s02 + 48*s12 - 16*si2)*si2 - 
           24*POWER2(s12))) + 
     sik*(s2k*(16*s01*s12 + 16*s01*si1 + 
           si0*(-48*s12 + 16*si1 - 16*si2) + 48*s01*si2) + 
        s1k*(-48*s02*s12 - 16*s02*si1 + 
           si0*(16*s12 + 16*si1 - 16*si2) - 48*s02*si2) + 
        (-16*s01*s12 + 48*s02*s12 + 
           (-16*s01 - 16*s02)*si1 + (16*s01 + 16*s02)*si2)*
         sik + s0k*(-16*s02*s12 + 
           si1*(-16*s01 + 32*s02 - 16*s12 + 16*si1) + 
           si0*(-32*s12 + 32*si1 - 32*si2) + 
           (48*s01 + 32*s02 - 16*s12 - 16*si2)*si2 + 
           16*POWER2(s12))) + 32*s01*si1*POWER2(s2k))/
   (s12*(s01 + s02 + s12)*si1*(s12 - si1 - si2))
);
term[7][5] = (
(s0k*(s1k*(si1*(32*s12 + 32*si1) + 
           si2*(-32*s12 + 32*si2)) + 
        s2k*(si2*(-64*s12 + 64*si2) + 
           si1*(-64*s12 + 128*si1 + 128*si2))) + 
     s0k*sik*(-32*s12*si2 + 
        si1*(-64*s12 - 32*si1 + 32*si2) + 32*POWER2(s12)) + 
     Q2inv*(si1*(16*s01*s12 - 32*s02*s12 + 
           (16*s01 + 64*s02)*si1 + 64*s02*si2) + 
        si2*(-16*s01*s12 - 32*s02*s12 + 
           (16*s01 + 32*s02)*si2) + 
        si0*(-16*s12*si2 + 
           si1*(-32*s12 - 16*si1 + 16*si2) + 16*POWER2(s12))
        ))/(s12*si1*POWER2(s12 - si1 - si2))
);
term[6][6] = (
(Q2inv*((s01*(-32*s02 - 64*s12) + s02*(32*s02 - 32*s12))*
         si1 + (s01*(32*s01 - 32*s02 - 32*s12) - 
           64*s02*s12)*si2 + 
        si0*(s01*(-20*s01 - 104*s02 - 72*s12) + 
           s02*(-20*s02 - 72*s12) + 12*POWER2(s12))) + 
     sik*((s01*(-64*s02 - 128*s12) + s02*(64*s02 - 64*s12))*
         s1k + (s01*(64*s01 - 64*s02 - 64*s12) - 
           128*s02*s12)*s2k + 
        s0k*(s01*(-40*s01 - 208*s02 - 144*s12) + 
           s02*(-40*s02 - 144*s12) + 24*POWER2(s12))))/
   (POWER2(s12)*POWER2(s01 + s02 + s12))
);
term[7][6] = (
(s1k*(-64*s02*s1k*si2 + s2k*
         (-64*s12*si0 + 64*s02*si1 + 64*s01*si2)) + 
     s0k*(s2k*(-64*s12*si0 + 
           (-32*s01 + 32*s02 + 32*s12)*si1) + 
        s0k*(64*s12*si1 + 64*s12*si2) + 
        s1k*(-64*s12*si0 + (32*s01 - 32*s02 + 32*s12)*si2))\
      + Q2inv*(si2*(s01*(16*s01 - 16*s02) - 16*s01*si2) + 
        si1*(-16*s01*s02 - 16*s02*si1 + 
           (16*s01 + 16*s02)*si2 + 16*POWER2(s02)) + 
        si0*(-28*s01*s12 - 28*s02*s12 + 64*s12*si0 + 
           (-20*s01 - 52*s02 + 28*s12)*si1 + 
           (-52*s01 - 20*s02 + 28*s12)*si2 - 76*POWER2(s12))
        ) + sik*(s2k*(32*s01*s12 + 64*s12*si0 + 
           32*s01*si1 - 32*s01*si2) + 
        s1k*(32*s02*s12 + 64*s12*si0 - 32*s02*si1 + 
           32*s02*si2) + (-64*s01*s12 - 64*s02*s12)*sik + 
        s0k*(8*s01*s12 + 8*s02*s12 + 128*s12*si0 + 
           (-40*s01 - 104*s02 - 8*s12)*si1 + 
           (-104*s01 - 40*s02 - 8*s12)*si2 + 40*POWER2(s12))
        ) - 64*s01*si1*POWER2(s2k))/
   ((s01 + s02 + s12)*(s12 - si1 - si2)*POWER2(s12))
);
term[7][7] = (
(s0k*(s1k*((-64*s12 - 64*si2)*si2 + 
           si1*(-128*s12 + 64*si2)) + 
        s2k*(-128*s12*si2 + si1*(-64*s12 - 64*si1 + 64*si2))
        ) + s0k*sik*(si1*(144*s12 - 40*si1 - 208*si2) + 
        (144*s12 - 40*si2)*si2 + 24*POWER2(s12)) + 
     Q2inv*(si2*(-32*s01*s12 - 64*s02*s12 - 32*s01*si2) + 
        si1*(-64*s01*s12 - 32*s02*s12 - 32*s02*si1 + 
           (32*s01 + 32*s02)*si2) + 
        si0*(si1*(72*s12 - 20*si1 - 104*si2) + 
           (72*s12 - 20*si2)*si2 + 12*POWER2(s12))))/
   (POWER2(s12)*POWER2(s12 - si1 - si2))
);
}
void DISMatrixElement::FillLorentzLgXqqg(
                        ) {
term[0][0] = (
(s0k*s1k*(64*si0 + 64*si2) + s1k*s2k*(64*si0 + 64*si2) + 
     Q2inv*((32*s01 + 32*s12)*si0 - 32*s02*si1 + 
        (32*s01 + 32*s12)*si2) - 64*s02*s1k*sik)/
   (si0*POWER2(-s02 + si0 + si2))
);
term[1][0] = (
(s1k*s2k*(64*s02 - 64*si0) + 64*s0k*s1k*si2 + 
     Q2inv*(32*s02*s12 - 32*s02*si1 + 
        si0*(-32*s12 + 32*si1) + 32*s01*si2) + 
     s1k*(-64*s02 + 64*si0)*sik)/(s02*si0*(s02 - si0 - si2))
);
term[2][0] = (
(s0k*s1k*(-32*si1 - 32*si2) + 
     Q2inv*si1*(-32*s01 - 32*s02 - 32*s12 + 32*si1 + 
        32*si2) + s1k*(s2k*(32*si0 - 32*si1) + 
        s1k*(32*si0 + 32*si2)) + 
     sik*(32*s0k*s12 + s2k*(-32*s01 + 32*si1) + 
        s1k*(-32*s01 - 32*s12 + 64*si1 + 32*si2) - 
        32*s12*sik))/(s12*si0*(-s02 + si0 + si2))
);
term[3][0] = (
(s1k*(s1k*si0*(-32*s02 + 32*si0 + 32*si2) + 
        s2k*(64*s01*s02 + 
           si0*(-32*s01 - 32*s02 + 32*s12 + 32*si0) - 
           32*s02*si1 - 32*s01*si2)) + 
     s0k*(s2k*(32*s02*si1 - 32*si0*si1) + 32*s0k*si1*si2 + 
        s1k*(si0*(32*s12 - 32*si1 - 32*si2) - 32*si1*si2))\
      + sik*(-32*s01*s02*s2k + 
        s1k*((-32*s01 - 32*s12)*si0 + 32*s02*si1) + 
        s0k*(32*s12*si0 - 32*s02*si1 - 32*s01*si2) + 
        32*s01*s02*sik) + 
     Q2inv*(32*s01*s02*s12 + s01*(-16*s02 - 16*s12)*si2 + 
        si1*(s02*(16*s02 - 16*s12) + 16*s02*si1 - 
           16*s01*si2) + 
        si0*(-16*s02*s12 + (-32*s01 - 16*s12)*si1 + 
           16*POWER2(s12))))/(s02*si0*si1*(s02 - si0 - si2))
);
term[4][0] = (
(s1k*(s2k*(s01*(32*s01 + 32*s02) + 
           si0*(32*s02 + 32*s12 - 32*si0 - 64*si1) - 
           32*s01*si1 - 32*s01*si2) + 
        s1k*(-32*s01*s02 + 
           si0*(32*s01 + 64*s02 - 32*si0 - 32*si1) + 
           32*s02*si1 + 32*s02*si2 - 32*POWER2(s02))) + 
     Q2inv*(si2*(s01*(-64*s01 - 48*s02) + 16*s01*si2) + 
        si0*(16*s02*s12 + s01*(-32*s01 - 32*s02 + 32*s12) - 
           16*s12*si1 - 16*s12*si2) + 
        si1*(s01*(-64*s01 - 64*s02 - 32*s12) + 
           (32*s01 + 16*s02)*si1 + (48*s01 + 16*s02)*si2 - 
           16*POWER2(s02)) + 
        s01*(s01*(32*s01 + 64*s02 + 32*s12) + 
           32*POWER2(s02))) + 
     s0k*(s2k*(s01*(32*s01 + 32*s02 + 32*s12) + 
           32*si0*si1 + si1*
            (-64*s01 - 64*s02 - 32*s12 + 32*si1 + 32*si2))\
         + s1k*(s01*(64*s01 + 96*s02 + 32*s12) + 
           s02*(32*s02 + 64*s12) + (-32*s02 - 32*s12)*si2 + 
           si1*(-96*s01 - 128*s02 - 32*s12 + 32*si1 + 
              32*si2) + si0*
            (-32*s01 - 32*s02 - 64*s12 + 64*si1 + 32*si2))\
         + s0k*(-32*s01*s12 - 32*s02*s12 + 
           si1*(32*s01 + 32*s02 + 64*s12 - 32*si1 - 
              32*si2) - 32*POWER2(s12))) + 
     sik*(s2k*(32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        s1k*(s01*(-32*s01 - 64*s02) + 
           s02*(-32*s02 - 32*s12) + (32*s01 + 32*s02)*si0 + 
           (32*s01 + 64*s02)*si1 + 32*s02*si2) - 
        32*s02*s12*sik + 
        s0k*(s01*(-32*s01 - 32*s02) + 64*s02*s12 - 
           32*s12*si0 + (32*s01 - 32*s12)*si1 + 
           32*s01*si2 + 32*POWER2(s12))) - 
     32*si0*si1*POWER2(s2k))/
   (s12*si0*(s12 - si1 - si2)*(-s02 + si0 + si2))
);
term[5][0] = (
(sik*(s0k*(32*s01*si0 + 32*s01*si1 - 32*POWER2(s01)) + 
        s1k*(32*s01*si0 + 32*s01*si1 - 32*POWER2(s01))) + 
     Q2inv*(si1*(s01*(-64*s01 - 32*s02 - 32*s12) + 
           32*s01*si1 + 64*s01*si2) + 
        si0*(s01*(-64*s01 - 32*s02 - 32*s12) + 32*s01*si0 + 
           64*s01*si1 + 64*s01*si2) + 
        (32*s01 + 32*s02 + 32*s12)*POWER2(s01) - 
        64*si2*POWER2(s01)) + 
     s1k*(s1k*(-32*s01*s02 + 
           si0*(32*s01 + 32*s02 - 32*si0 - 32*si1) + 
           32*s02*si1) + 
        s2k*(-32*s01*si0 - 32*s01*si1 + 32*POWER2(s01))) + 
     s0k*(s0k*(-32*s01*s12 + si0*(32*s12 - 32*si1) + 
           (32*s01 + 32*s12 - 32*si1)*si1) + 
        s1k*(s01*(64*s01 + 32*s02 + 32*s12) + 
           si1*(-96*s01 - 32*s02 - 32*s12 + 32*si1) + 
           si0*(-96*s01 - 32*s02 - 32*s12 + 32*si0 + 64*si1)
           ) + s2k*(-32*s01*si0 - 32*s01*si1 + 
           32*POWER2(s01))))/
   (si0*si1*(-s02 + si0 + si2)*(-s12 + si1 + si2))
);
term[6][0] = (
(s0k*s1k*(si2*(-32*s02 + 32*si2) + 
        si0*(-32*s02 - 32*si0 + 64*si2)) + 
     s1k*s2k*(si0*(128*s02 - 128*si0 - 64*si2) + 
        64*s02*si2 - 64*POWER2(s02)) + 
     s1k*sik*(si0*(32*si0 - 32*si2) - 32*s02*si2 + 
        32*POWER2(s02)) + 
     Q2inv*(si2*(-16*s01*s02 + 32*s02*s12 + 16*s01*si2) + 
        si0*(-16*s01*s02 + 64*s02*s12 + 
           si0*(-16*s01 - 64*s12 + 16*si1) + 
           (32*s01 - 32*s12)*si2 - 16*si1*si2) - 
        32*s12*POWER2(s02) + 
        si1*(-16*s02*si2 + 16*POWER2(s02))))/
   (si0*si2*POWER2(-s02 + si0 + si2))
);
term[7][0] = (
(s0k*(s0k*(16*s02*s12 + si0*(16*s12 - 16*si1) + 
           48*s12*si2 + si1*(-16*s02 + 16*si2)) + 
        s2k*(-16*s01*s02 + si0*(-16*s01 - 16*si1) - 
           48*s01*si2 + si1*(48*s02 + 16*si2)) + 
        s1k*(-32*s01*s02 + s02*(-16*s02 - 32*s12) + 
           48*s02*si1 + si2*
            (-32*s01 - 16*s02 + 16*s12 + 16*si2) + 
           si0*(-32*s01 + 32*s12 + 16*si0 + 16*si1 + 16*si2)
           )) + sik*(s1k*
         (48*s01*s02 + s02*(16*s02 + 16*s12) + 
           (-16*s01 - 16*s02)*si0 - 48*s02*si1 + 16*s01*si2)
          + s2k*(-32*s12*si0 - 32*s02*si1 + 32*s01*si2) + 
        s0k*(16*s01*s02 - 48*s02*s12 + 
           (16*s01 + 16*s12)*si0 + (-16*s01 - 48*s12)*si2)\
         + 32*s02*s12*sik) + 
     s1k*(s2k*(-32*s01*s02 + 
           si0*(-16*s02 - 16*s12 + 16*si0 + 48*si1 - 
              32*si2) + 32*s01*si2) + 
        s1k*(si0*(-48*s02 + 16*si0 - 16*si2) - 32*s02*si2 + 
           32*POWER2(s02))) + 
     Q2inv*(s01*(-16*s01*s02 + s02*(-16*s02 - 32*s12)) + 
        si2*(s01*(-16*s01 + 48*s02 + 24*s12) - 
           24*s01*si2) + 
        si1*(48*s01*s02 + s02*(16*s02 + 8*s12) - 
           24*s02*si1 + (8*s01 - 24*s02)*si2) + 
        si0*(-16*s02*s12 + s01*(-16*s01 + 16*s12) + 
           16*s01*si0 + 24*s12*si1 + 
           (-16*s01 + 24*s12)*si2 - 8*POWER2(s12))) + 
     32*si0*si1*POWER2(s2k))/
   (si0*si2*(-s02 + si0 + si2)*(-s12 + si1 + si2))
);
term[1][1] = (
(64*s1k*s2k*si0 + s0k*s1k*(64*s02 - 64*si2) + 
     Q2inv*(32*s01*s02 + 32*s12*si0 - 32*s01*si2 + 
        si1*(-32*s02 + 32*si2)) + s1k*(-64*s02 + 64*si2)*sik
     )/(s02*POWER2(s02 - si0 - si2))
);
term[2][1] = (
(s0k*(32*s0k*s12*si2 + s2k*
         (-32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        s1k*(32*s02*s12 - 32*s02*si1 + 
           (-32*s02 - 32*s12)*si2)) + 
     s1k*(s2k*(32*s01*s02 + 32*s12*si0 - 32*s02*si1) + 
        s1k*(32*s02*si0 + 32*s02*si2 - 32*POWER2(s02))) + 
     sik*(-32*s01*s2k*si0 + 
        s0k*(-32*s02*s12 + 32*s12*si0) + 
        s1k*(32*s01*s02 + (-64*s01 - 32*s02 - 32*s12)*si0 + 
           32*s02*si1 - 32*s01*si2 + 32*POWER2(s02))) + 
     Q2inv*(32*s01*s02*s12 - 16*s01*s12*si2 + 
        si1*(-16*s02*s12 + 16*s02*si1 - 16*s01*si2) + 
        si0*(16*s12*si0 + (-32*s01 - 16*s02 - 16*s12)*si1 - 
           16*s01*si2 + 16*POWER2(s12))) + 
     32*s01*si0*POWER2(s2k))/(s02*s12*si0*(s02 - si0 - si2))
);
term[3][1] = (
(s0k*(-32*s2k*si1 + s1k*(-32*s12 + 32*si2)) + 
     s1k*(s1k*(32*s02 - 32*si2) + 
        s2k*(-32*s01 - 64*s12 + 32*si1 + 32*si2)) + 
     ((-32*s02 + 32*s12)*s1k + (32*s01 + 32*s12)*s2k)*sik + 
     Q2inv*(-32*s01*s12 + 32*s12*si0 + 32*s12*si1 + 
        32*s12*si2 - 32*POWER2(s12)) - 32*si1*POWER2(s2k))/
   (s02*si1*(s02 - si0 - si2))
);
term[4][1] = (
(s1k*(s01*(32*s01 + 32*s02 + 32*s12)*s2k + 
        s1k*(-32*s01*s02 + s02*(-32*s02 - 32*s12) + 
           (32*s01 + 32*s02 + 32*s12)*si0)) + 
     (s01*s0k*(-32*s01 - 32*s02 - 32*s12) + 
        s01*(-32*s01 - 32*s02 - 32*s12)*s1k)*sik + 
     Q2inv*(s01*(-32*s01 - 32*s02 - 32*s12)*si0 + 
        s01*(-32*s01 - 32*s02 - 32*s12)*si1 + 
        s01*(-64*s01 - 64*s02 - 64*s12)*si2 + 
        s01*(s02*(32*s02 + 64*s12) + 
           s01*(32*s01 + 64*s02 + 64*s12) + 32*POWER2(s12)))
       + s0k*(s01*(32*s01 + 32*s02 + 32*s12)*s2k + 
        s0k*(-32*s01*s12 - 32*s02*s12 + 
           (32*s01 + 32*s02 + 32*s12)*si1 - 32*POWER2(s12))\
         + s1k*(s02*(32*s02 + 64*s12) + 
           s01*(64*s01 + 96*s02 + 96*s12) + 
           (-32*s01 - 32*s02 - 32*s12)*si0 + 
           (-32*s01 - 32*s02 - 32*s12)*si1 + 32*POWER2(s12))
        ))/(s02*s12*(s02 - si0 - si2)*(s12 - si1 - si2))
);
term[5][1] = (
(s1k*(s1k*(-32*s01*s02 + s02*(-32*s02 - 32*s12) + 
           si0*(32*s01 + 64*s02 + 32*s12 - 32*si0 - 32*si2))
          + s2k*(s01*(32*s01 + 32*s02 + 32*s12) + 
           si0*(-64*s01 - 32*s02 - 64*s12 + 32*si0 + 
              32*si1 + 32*si2))) + 
     sik*(s2k*(32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        s1k*(s01*(-32*s01 - 32*s12) + 
           s02*(32*s02 + 64*s12) + (32*s01 - 32*s02)*si0 - 
           32*s02*si1 + 32*s01*si2) - 32*s02*s12*sik + 
        s0k*(s01*(-32*s01 - 64*s12) - 32*s02*s12 + 
           (32*s01 + 64*s12)*si0 + (32*s01 + 32*s12)*si1 + 
           32*s12*si2 - 32*POWER2(s12))) + 
     Q2inv*(si2*(s01*(-64*s01 - 48*s12) + 16*s01*si2) + 
        si1*(s01*(-32*s01 + 32*s02 - 32*s12) + 16*s02*s12 - 
           16*s02*si2) + 
        si0*(s01*(-64*s01 - 32*s02 - 64*s12) + 
           (32*s01 + 16*s12)*si0 - 16*s02*si1 + 
           (48*s01 + 16*s12)*si2 - 16*POWER2(s12)) + 
        s01*(s01*(32*s01 + 32*s02 + 64*s12) + 
           32*POWER2(s12))) + 
     s0k*(s2k*(s01*(32*s01 + 32*s12) + 
           si0*(-32*s01 - 64*si1) + 
           (32*s02 + 32*s12 - 32*si1)*si1 - 32*s01*si2) + 
        s0k*(-32*s01*s12 + si0*(32*s12 - 32*si1) + 
           (32*s01 + 64*s12 - 32*si1)*si1 + 32*s12*si2 - 
           32*POWER2(s12)) + 
        s1k*(64*s02*s12 + s01*(64*s01 + 32*s02 + 96*s12) + 
           (-32*s02 - 32*s12)*si2 + 
           si1*(-32*s01 - 64*s02 - 32*s12 + 32*si2) + 
           si0*(-96*s01 - 32*s02 - 128*s12 + 32*si0 + 
              64*si1 + 32*si2) + 32*POWER2(s12))) - 
     32*si0*si1*POWER2(s2k))/
   (s02*si1*(s02 - si0 - si2)*(-s12 + si1 + si2))
);
term[6][1] = (
(s1k*s2k*(si0*(-32*si0 - 32*si2) - 32*s02*si2 - 
        32*POWER2(s02)) + 
     s0k*s1k*(si2*(-64*s02 + 32*si2) + 
        si0*(-32*s02 + 32*si2) - 32*POWER2(s02)) + 
     s1k*sik*(-64*s02*si2 + 
        si0*(-128*s02 + 64*si0 + 64*si2) + 128*POWER2(s02))\
      + Q2inv*(si2*(-32*s01*s02 - 16*s02*s12 + 
           16*s01*si2) + 
        si0*(-16*s01*s02 + si0*(-16*s12 + 32*si1) + 
           (16*s01 - 16*s12)*si2 + si1*(-64*s02 + 32*si2))\
         - 16*s01*POWER2(s02) - 16*s12*POWER2(s02) + 
        si1*(-32*s02*si2 + 64*POWER2(s02))))/
   (s02*si2*POWER2(s02 - si0 - si2))
);
term[7][1] = (
(s0k*(s1k*(-32*s01*s02 + s02*(-16*s02 - 16*s12) + 
           si0*(-32*s01 - 48*s12 + 16*si0 + 32*si1 - 
              16*si2) + (32*s01 + 16*s02 - 16*si2)*si2 + 
           si1*(-32*s02 + 16*si2)) + 
        s0k*(16*s02*s12 + si0*(16*s12 - 16*si1) + 
           16*s12*si2 + si1*(-16*s02 + 48*si2)) + 
        s2k*(-16*s01*s02 + si0*(-16*s01 - 48*si1) - 
           16*s01*si2 + si1*(16*s02 + 48*si2))) + 
     sik*(s2k*(32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        s1k*(s02*(16*s02 + 48*s12) + 
           (32*s01 - 16*s02)*si0 - 16*s02*si1 + 
           (32*s01 + 32*s02)*si2) + 
        s0k*(16*s01*s02 - 16*s02*s12 + 
           (16*s01 + 48*s12)*si0 + (-48*s01 - 16*s12)*si2)\
         - 32*s02*s12*sik) + 
     s1k*(s2k*(16*s01*s02 + 
           si0*(-48*s01 - 16*s02 - 48*s12 + 16*si0 + 
              16*si1) + 16*s01*si2) + 
        s1k*(si0*(48*s02 - 32*si0 - 32*si2) - 16*s02*si2 - 
           16*POWER2(s02))) + 
     Q2inv*(si2*(s01*(16*s01 - 16*s02 + 8*s12) + 
           24*s01*si2) + 
        si1*(-16*s01*s02 + 24*s02*s12 - 8*s02*si1 + 
           (24*s01 - 24*s02)*si2) + 
        s01*(-16*s01*s02 - 16*POWER2(s02)) + 
        si0*(s01*(-16*s01 - 48*s12) + 
           (16*s01 + 16*s12)*si0 + 
           (32*s01 - 16*s02 + 8*s12)*si1 + 
           (48*s01 + 24*s12)*si2 - 24*POWER2(s12))) - 
     32*si0*si1*POWER2(s2k))/
   (s02*(s02 - si0 - si2)*si2*(-s12 + si1 + si2))
);
term[2][2] = (
(32*Q2inv*si2 + 64*s2k*sik)/(s12*si0)
);
term[3][2] = (
(Q2inv*(s01*(-32*s01 - 32*s02 - 32*s12)*si0 + 
        s01*(-32*s01 - 32*s02 - 32*s12)*si1 + 
        (32*s01 + 32*s02 + 32*s12)*POWER2(s01)) + 
     sik*(s01*s0k*(-32*s01 - 32*s02 + 32*s12) + 
        s01*(-32*s01 + 32*s02 - 32*s12)*s1k - 
        64*s2k*POWER2(s01)) + 
     s1k*(s1k*(-32*s01*s02 + 32*s01*si0 + 32*s01*si2) + 
        s2k*(32*s01*si0 - 32*s01*si1 + 32*POWER2(s01))) + 
     s0k*(s1k*(s01*(64*s01 + 32*s02 + 32*s12) - 
           32*s01*si0 - 32*s01*si1 - 64*s01*si2) + 
        s0k*(-32*s01*s12 + 32*s01*si1 + 32*s01*si2) + 
        s2k*(-32*s01*si0 + 32*s01*si1 + 32*POWER2(s01))))/
   (s02*s12*si0*si1)
);
term[4][2] = (
(-32*s1k*s2k*si0 + s0k*(s0k*(32*s12 - 32*si2) + 
        s1k*(-32*s02 + 32*si2) + 
        s2k*(-32*s01 - 64*s02 + 32*si0 + 32*si2)) + 
     (s0k*(32*s02 - 32*s12) + (32*s01 + 32*s02)*s2k)*sik + 
     Q2inv*(-32*s01*s02 + 32*s02*si0 + 32*s02*si1 + 
        32*s02*si2 - 32*POWER2(s02)) - 32*si0*POWER2(s2k))/
   (s12*si0*(s12 - si1 - si2))
);
term[5][2] = (
(s1k*(s2k*si0*(32*s12 - 32*si1) + 32*s1k*si0*si2) + 
     s0k*(s1k*(si1*(32*s02 - 32*si2) + 
           si0*(-32*si1 - 32*si2)) + 
        s0k*si1*(-32*s12 + 32*si1 + 32*si2) + 
        s2k*(64*s01*s12 - 32*s12*si0 + 
           si1*(-32*s01 + 32*s02 - 32*s12 + 32*si1) - 
           32*s01*si2)) + 
     sik*(-32*s01*s12*s2k + 
        s0k*(32*s12*si0 + (-32*s01 - 32*s02)*si1) + 
        s1k*(-32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        32*s01*s12*sik) + 
     Q2inv*(32*s01*s02*s12 + s02*(16*s02 - 16*s12)*si1 + 
        s01*(-16*s02 - 16*s12)*si2 + 
        si0*(-16*s02*s12 + 16*s12*si0 + 
           (-32*s01 - 16*s02)*si1 - 16*s01*si2 + 
           16*POWER2(s12))))/(s12*si0*si1*(s12 - si1 - si2))
);
term[6][2] = (
(s0k*(s2k*(16*s01*s02 + 32*s02*si1 + 
           si0*(-32*s01 - 32*s12 + 16*si1) - 48*s01*si2) + 
        s0k*(-16*s02*s12 + 32*s12*si0 + 48*s12*si2) + 
        s1k*(32*s01*s02 + s02*(16*s02 + 32*s12) + 
           si1*(-16*s02 - 16*si2) + 
           (-32*s01 - 32*s02 - 16*s12 - 16*si2)*si2 + 
           si0*(-64*s01 - 32*s02 - 16*s12 + 32*si1 + 16*si2)
           )) + sik*(s2k*
         (-16*s01*s02 + si0*(16*s01 - 32*si1) + 16*s01*si2)\
         + s1k*(-16*s01*s02 - 16*s02*s12 + 32*s02*si1 + 
           si0*(32*s01 + 16*s02 + 32*s12 - 64*si1 - 
              32*si2) + (16*s01 + 32*s02)*si2) + 
        s0k*(-16*s02*s12 - 32*s12*si0 - 16*s12*si2) + 
        32*s12*si0*sik) + 
     s1k*(s2k*(32*s01*s02 - 32*s02*si1 + 
           si0*(-16*s01 - 16*s02 + 16*s12 - 32*si0 + 
              32*si1 + 16*si2)) + 
        s1k*(si0*(32*s02 - 32*si0 - 16*si2) + 32*s02*si2 - 
           32*POWER2(s02))) + 
     Q2inv*(s01*(16*s01*s02 + s02*(16*s02 + 32*s12)) + 
        si1*(-16*s01*s02 - 24*s02*s12 + 16*s02*si1) + 
        si2*(s01*(-16*s01 - 32*s02 - 8*s12) + 16*s01*si2) + 
        si0*(s01*(-32*s01 - 32*s02 - 16*s12) + 16*s02*s12 + 
           si1*(32*s01 + 16*s02 + 32*s12 - 32*si1 - 
              32*si2) + (16*s01 - 16*s12)*si2 + 
           8*POWER2(s12))) + 32*s01*si0*POWER2(s2k))/
   (s12*si0*si2*(-s02 + si0 + si2))
);
term[7][2] = (
(s1k*(s2k*si0*(-32*s12 - 16*si1 + 16*si2) + 
        s1k*si0*(-32*s12 + 16*si1 + 48*si2)) + 
     Q2inv*(s01*(-32*s01*s12 - 32*s02*s12) + 
        si1*(s02*(16*s02 + 16*s12) + 
           s01*(16*s01 + 16*s02 + 32*s12) - 16*s01*si1 - 
           32*s01*si2) + 
        si2*(32*s02*s12 + s01*(16*s01 + 16*s12) - 
           16*s01*si2) + 
        si0*(16*s01*s12 + 32*s02*s12 + 8*s12*si0 + 
           (-32*s01 - 24*s02 + 16*s12)*si1 + 
           (-8*s01 + 16*s12)*si2) - 32*s12*POWER2(s02)) + 
     sik*(s1k*(32*s01*s12 + 16*s02*s12 - 32*s12*si0 + 
           (-16*s01 + 32*s02)*si1 - 48*s01*si2) + 
        s2k*(16*s01*s12 + 32*s02*s12 - 16*s01*si1 - 
           16*s01*si2) + 32*s01*s12*sik + 
        s0k*(16*s01*s12 + 32*s02*s12 + 16*s12*si0 + 
           (-32*s01 - 32*s02 - 16*s12)*si1 - 16*s12*si2 - 
           32*POWER2(s12))) + 
     s0k*(s2k*(-32*s01*s12 - 64*s02*s12 + 
           si0*(32*s12 - 16*si1) + 
           si1*(16*s01 + 32*s02 + 16*s12 - 32*si2) + 
           (16*s01 + 32*s12)*si2) + 
        s1k*(-64*s01*s12 - 32*s02*s12 + 
           si1*(32*s01 + 16*s02 + 32*s12 - 16*si1 - 
              32*si2) + si0*(16*s12 - 32*si1 - 16*si2) + 
           si2*(32*s01 - 16*s02 + 16*s12 + 16*si2)) + 
        s0k*(-16*s12*si2 + 
           si1*(-32*s12 + 32*si1 + 32*si2) + 32*POWER2(s12))
        ) - 32*s12*si0*POWER2(s2k))/
   (s12*si0*(s12 - si1 - si2)*si2)
);
term[3][3] = (
(32*Q2inv*si2 + 64*s2k*sik)/(s02*si1)
);
term[4][3] = (
(s1k*(32*s02*s1k*si2 + s2k*
         (32*s12*si0 - 32*s02*si1 - 32*s01*si2)) + 
     Q2inv*(32*s01*s02*s12 - 16*s01*s02*si2 + 
        si0*(-16*s02*s12 + 16*s12*si0 + 
           (-32*s01 - 16*s02 - 16*s12)*si1 - 16*s01*si2) + 
        si1*(16*s02*si1 - 16*s01*si2 + 16*POWER2(s02))) + 
     s0k*(s2k*(32*s01*s12 - 32*s12*si0 + 32*s02*si1) + 
        s1k*(32*s02*s12 - 32*s12*si0 + 
           (-32*s02 - 32*s12)*si2) + 
        s0k*(32*s12*si1 + 32*s12*si2 - 32*POWER2(s12))) + 
     sik*(-32*s01*s2k*si1 + 
        s1k*(-32*s02*s12 + 32*s02*si1) + 
        s0k*(32*s01*s12 + 32*s12*si0 + 
           (-64*s01 - 32*s02 - 32*s12)*si1 - 32*s01*si2 + 
           32*POWER2(s12))) + 32*s01*si1*POWER2(s2k))/
   (s02*s12*si1*(s12 - si1 - si2))
);
term[5][3] = (
(Q2inv*si0*(-32*s01 - 32*s02 - 32*s12 + 32*si0 + 32*si2) + 
     s0k*(s2k*(-32*si0 + 32*si1) + s1k*(-32*si0 - 32*si2) + 
        s0k*(32*si1 + 32*si2)) + 
     sik*(32*s02*s1k + s2k*(-32*s01 + 32*si0) + 
        s0k*(-32*s01 - 32*s02 + 64*si0 + 32*si2) - 
        32*s02*sik))/(s02*si1*(-s12 + si1 + si2))
);
term[6][3] = (
(s0k*(s0k*(-16*si0*si1 + si1*(32*s02 - 48*si2)) + 
        s2k*(16*si0*si1 + si1*(32*s02 - 16*si2)) + 
        s1k*(64*s01*s02 + 32*s02*s12 + 
           (-32*s01 - 16*s02 + 16*s12 - 16*si2)*si2 + 
           si1*(-16*s02 + 16*si2) + 
           si0*(-32*s01 - 32*s02 - 16*s12 + 16*si0 + 
              32*si1 + 32*si2))) + 
     sik*(s2k*(-16*s01*s02 - 32*s02*s12 + 16*s01*si0 + 
           16*s01*si2) + 
        s0k*(-32*s01*s02 - 16*s02*s12 + 
           (16*s01 - 32*s12)*si0 + 32*s02*si1 + 48*s01*si2)\
         + s1k*(-16*s01*s02 + s02*(32*s02 - 32*s12) + 
           (32*s01 + 16*s02 + 32*s12)*si0 - 16*s02*si1 + 
           16*s02*si2) - 32*s01*s02*sik) + 
     s1k*(s2k*(32*s01*s02 + 64*s02*s12 - 32*s02*si1 + 
           (-16*s01 - 32*s02)*si2 + 
           si0*(-16*s01 - 16*s02 - 32*s12 + 16*si1 + 32*si2)
           ) + s1k*(si0*(32*s02 - 32*si0 - 32*si2) + 
           16*s02*si2 - 32*POWER2(s02))) + 
     Q2inv*(s01*(32*s01*s02 + 32*s02*s12) + 
        si2*(s01*(-16*s01 - 16*s02) - 32*s02*s12 + 
           16*s01*si2) + 
        si1*(-16*s01*s02 - 32*s02*s12 - 8*s02*si1 + 
           (8*s01 - 16*s02)*si2) + 
        si0*(s01*(-16*s01 - 32*s02 - 16*s12) - 16*s02*s12 + 
           16*s01*si0 + (32*s01 - 16*s02 + 24*s12)*si1 + 
           32*s01*si2 - 16*POWER2(s12)) + 32*s02*POWER2(s12)
        ) + 32*s02*si1*POWER2(s2k))/
   (s02*si1*(s02 - si0 - si2)*si2)
);
term[7][3] = (
(s1k*(s2k*(-16*s01*s12 + si0*(-32*s12 - 16*si1) + 
           (32*s01 + 32*s02)*si1 + 48*s01*si2) + 
        s1k*(16*s02*s12 - 32*s02*si1 - 48*s02*si2)) + 
     sik*(s2k*(16*s01*s12 - 16*s01*si1 + 32*si0*si1 - 
           16*s01*si2) + 
        s1k*(16*s02*s12 + 32*s02*si1 + 16*s02*si2) + 
        s0k*(16*s01*s12 + 16*s02*s12 + 
           si0*(-32*s12 + 64*si1) + 
           (-16*s01 - 32*s12)*si2 + 
           si1*(-32*s01 - 32*s02 - 16*s12 + 32*si2)) - 
        32*s02*si1*sik) + 
     Q2inv*(si2*(s01*(16*s01 + 8*s02 + 32*s12) - 
           16*s01*si2) + 
        si1*(s02*(-8*s02 - 16*s12) + 
           s01*(32*s01 + 16*s02 + 32*s12) + 
           (-16*s01 + 16*s02)*si2) + 
        si0*(16*s01*s12 + 24*s02*s12 + 
           si0*(-16*s12 + 32*si1) + 
           si1*(-32*s01 - 32*s02 - 16*s12 + 32*si2)) + 
        s01*(-16*s01*s12 - 32*s02*s12 - 16*POWER2(s12))) + 
     s0k*(s2k*(-32*s01*s12 + si0*(32*s12 - 32*si1) + 
           si1*(16*s01 - 16*s02 + 16*s12 + 32*si1 - 16*si2))
          + s1k*(-32*s01*s12 - 32*s02*s12 + 
           si1*(64*s01 + 16*s02 + 32*s12 - 16*si2) + 
           si2*(32*s01 + 16*s02 + 32*s12 + 16*si2) + 
           si0*(16*s12 - 32*si1 + 16*si2) - 16*POWER2(s12))\
         + s0k*(-32*s12*si2 + 
           si1*(-32*s12 + 32*si1 + 16*si2) + 32*POWER2(s12))
        ) - 32*s01*si1*POWER2(s2k))/
   (s02*si1*si2*(-s12 + si1 + si2))
);
term[4][4] = (
(s0k*(64*s2k*si1 + s1k*(64*s12 - 64*si2)) + 
     Q2inv*(32*s01*s12 + 32*s02*si1 - 32*s01*si2 + 
        si0*(-32*s12 + 32*si2)) + s0k*(-64*s12 + 64*si2)*sik
     )/(s12*POWER2(s12 - si1 - si2))
);
term[5][4] = (
(Q2inv*(32*s02*s12 - 32*s02*si1 + si0*(-32*s12 + 32*si1) + 
        32*s01*si2) + s0k*
      (s2k*(64*s12 - 64*si1) + 64*s1k*si2) + 
     s0k*(-64*s12 + 64*si1)*sik)/(s12*si1*(s12 - si1 - si2))
);
term[6][4] = (
(s1k*(s2k*(-16*s01*s12 - 16*s01*si1 - 16*s01*si2 + 
           si0*(16*s12 - 48*si1 + 48*si2)) + 
        s1k*(16*s02*s12 + 16*s02*si1 + 16*s02*si2 + 
           si0*(-16*s12 - 16*si1 + 48*si2))) + 
     Q2inv*(si2*(s01*(16*s01 + 8*s02 - 16*s12) + 
           24*s01*si2) + 
        si0*(-16*s01*s12 + 24*s02*s12 - 8*s12*si0 + 
           (32*s01 + 8*s02 - 16*s12)*si1 + 
           (24*s01 - 24*s12)*si2) + 
        si1*(s01*(-16*s01 - 48*s02) + 
           (16*s01 + 16*s02)*si1 + (48*s01 + 24*s02)*si2 - 
           24*POWER2(s02)) + 
        s01*(-16*s01*s12 - 16*POWER2(s12))) + 
     s0k*(s2k*(16*s01*s12 + 16*si0*si1 + 
           si1*(-48*s01 - 48*s02 - 16*s12 + 16*si1) + 
           16*s01*si2) + 
        s0k*(si1*(48*s12 - 32*si1 - 32*si2) - 16*s12*si2 - 
           16*POWER2(s12)) + 
        s1k*(-32*s01*s12 - 16*s02*s12 + 
           si1*(-32*s01 - 48*s02 + 16*si1 - 16*si2) + 
           (32*s01 + 16*s12 - 16*si2)*si2 + 
           si0*(-32*s12 + 32*si1 + 16*si2) - 16*POWER2(s12))
        ) + sik*(s2k*(32*s12*si0 + 32*s02*si1 - 
           32*s01*si2) + 
        s1k*(16*s01*s12 - 16*s02*s12 + 
           (16*s01 + 48*s02)*si1 + (-48*s01 - 16*s02)*si2)\
         - 32*s02*s12*sik + 
        s0k*(48*s02*s12 - 16*s12*si0 + 
           (32*s01 - 16*s12)*si1 + (32*s01 + 32*s12)*si2 + 
           16*POWER2(s12))) - 32*si0*si1*POWER2(s2k))/
   (s12*(s02 - si0 - si2)*(s12 - si1 - si2)*si2)
);
term[7][4] = (
(s0k*sik*(si1*(128*s12 - 64*si1 - 64*si2) + 64*s12*si2 - 
        128*POWER2(s12)) + 
     Q2inv*(si2*(32*s01*s12 + 16*s02*s12 - 16*s01*si2) + 
        si1*(16*s01*s12 + 16*s02*si1 + 
           (-16*s01 + 16*s02)*si2) + 
        si0*(si1*(64*s12 - 32*si1 - 32*si2) + 32*s12*si2 - 
           64*POWER2(s12)) + 16*s01*POWER2(s12) + 
        16*s02*POWER2(s12)) + 
     s0k*(s1k*(si1*(32*s12 - 32*si2) + 
           (64*s12 - 32*si2)*si2 + 32*POWER2(s12)) + 
        s2k*(32*s12*si2 + si1*(32*si1 + 32*si2) + 
           32*POWER2(s12))))/
   (s12*si2*POWER2(s12 - si1 - si2))
);
term[5][5] = (
(Q2inv*(-32*s12*si0 + (32*s01 + 32*s02)*si1 + 
        (32*s01 + 32*s02)*si2) + 
     s0k*(s1k*(64*si1 + 64*si2) + s2k*(64*si1 + 64*si2)) - 
     64*s0k*s12*sik)/(si1*POWER2(-s12 + si1 + si2))
);
term[6][5] = (
(s1k*(s1k*(16*s02*s12 + 16*s02*si1 + 48*s02*si2 + 
           si0*(-16*s12 - 16*si1 + 16*si2)) + 
        s2k*(-16*s01*s12 - 16*s01*si1 - 48*s01*si2 + 
           si0*(48*s12 - 16*si1 + 16*si2))) + 
     sik*(s2k*(-32*s12*si0 - 32*s02*si1 + 32*s01*si2) + 
        s1k*(16*s01*s12 - 48*s02*s12 + 
           (16*s01 + 16*s02)*si1 + (-16*s01 - 48*s02)*si2)\
         + 32*s02*s12*sik + 
        s0k*(48*s01*s12 + 16*s02*s12 - 48*s12*si0 + 
           (-16*s01 - 16*s12)*si1 + 16*s01*si2 + 
           16*POWER2(s12))) + 
     Q2inv*(si2*(s01*(-16*s01 + 24*s02 + 48*s12) - 
           24*s01*si2) + 
        si1*(s01*(-16*s01 + 16*s02) + 
           s02*(-8*s02 - 16*s12) + 16*s01*si1 + 
           (-16*s01 + 24*s02)*si2) + 
        s01*(-16*s01*s12 - 32*s02*s12 - 16*POWER2(s12)) + 
        si0*(48*s01*s12 + 8*s02*s12 - 24*s12*si0 + 
           24*s02*si1 + (8*s01 - 24*s12)*si2 + 
           16*POWER2(s12))) + 
     s0k*(s2k*(-32*s01*s12 + 48*si0*si1 + 
           si1*(-16*s02 - 16*s12 + 16*si1 - 32*si2) + 
           32*s01*si2) + 
        s1k*(-32*s01*s12 - 32*s02*s12 + 
           si0*(48*s12 + 16*si1) + 
           si2*(-32*s01 + 16*s02 - 16*s12 + 16*si2) + 
           si1*(-32*s01 + 32*s02 + 16*si1 + 16*si2) - 
           16*POWER2(s12)) + 
        s0k*(si1*(-48*s12 + 16*si1 - 16*si2) - 32*s12*si2 + 
           32*POWER2(s12))) + 32*si0*si1*POWER2(s2k))/
   (si1*(s02 - si0 - si2)*si2*(-s12 + si1 + si2))
);
term[7][5] = (
(s0k*sik*(32*s12*si2 + si1*(-32*si1 + 32*si2) - 
        32*POWER2(s12)) + 
     Q2inv*(si2*(16*s01*s12 - 32*s02*s12 - 16*s01*si2) + 
        si1*(16*s01*s12 - 64*s02*s12 + 
           (16*s01 + 64*s02)*si1 + (-32*s01 + 32*s02)*si2)\
         + si0*(16*s12*si2 + si1*(-16*si1 + 16*si2) - 
           16*POWER2(s12)) + 32*s02*POWER2(s12)) + 
     s0k*(s1k*(si1*(32*s12 + 32*si1 - 64*si2) + 
           (32*s12 - 32*si2)*si2) + 
        s2k*(-64*s12*si2 + 
           si1*(-128*s12 + 128*si1 + 64*si2) + 
           64*POWER2(s12))))/
   (si1*si2*POWER2(-s12 + si1 + si2))
);
term[6][6] = (
(s1k*s2k*(si0*(64*s02 + 64*si0 - 64*si2) + 128*s02*si2) + 
     s1k*sik*(-64*s02*si2 + si0*(-64*s02 + 128*si2) - 
        64*POWER2(s02)) + 
     s0k*s1k*(si0*(208*s02 - 40*si0 - 144*si2) + 
        si2*(144*s02 + 24*si2) - 40*POWER2(s02)) + 
     Q2inv*(si2*(72*s01*s02 + 64*s02*s12 + 12*s01*si2) + 
        si0*(104*s01*s02 + 32*s02*s12 + 
           (-20*s01 + 32*s12)*si0 + 
           (-72*s01 - 32*s12)*si2 + si1*(-32*s02 + 64*si2))\
         + si1*(-32*s02*si2 - 32*POWER2(s02)) - 
        20*s01*POWER2(s02)))/
   (POWER2(s02 - si0 - si2)*POWER2(si2))
);
term[7][6] = (
(s1k*(s2k*(si0*(32*s12 + 32*si1 - 32*si2) - 64*s01*si2) + 
        s1k*(64*s02*si2 - 64*si0*si2)) + 
     s0k*(s2k*(32*si0*si1 + si1*(32*s02 - 32*si2) - 
           64*s01*si2) + s0k*(64*s12*si2 - 64*si1*si2) + 
        s1k*(-40*s02*s12 + si1*(104*s02 - 8*si2) + 
           si0*(104*s12 - 40*si1 - 8*si2) + 
           (-128*s01 + 8*s02 + 8*s12 - 40*si2)*si2)) + 
     sik*(s2k*(-64*s12*si0 - 64*s02*si1 + 64*s01*si2) + 
        s1k*(-32*s02*s12 - 32*s02*si1 + 
           (64*s01 - 32*s02)*si2) + 
        s0k*(-32*s02*s12 - 32*s12*si0 + 
           (64*s01 - 32*s12)*si2) + 64*s02*s12*sik) + 
     Q2inv*(-20*s01*s02*s12 + 
        si1*(52*s01*s02 + s02*(16*s02 - 16*s12) - 
           16*s02*si1 + 28*s01*si2) + 
        si2*(s01*(-64*s01 - 28*s02 - 28*s12) + 
           76*s01*si2) + 
        si0*(52*s01*s12 - 16*s02*s12 - 16*s12*si0 + 
           (-20*s01 + 16*s02 + 16*s12)*si1 + 28*s01*si2 + 
           16*POWER2(s12))) + 64*si0*si1*POWER2(s2k))/
   ((s02 - si0 - si2)*(-s12 + si1 + si2)*POWER2(si2))
);
term[7][7] = (
(s0k*(s2k*(si1*(64*s12 + 64*si1 - 64*si2) + 128*s12*si2) + 
        s1k*(si1*(208*s12 - 40*si1 - 144*si2) + 
           si2*(144*s12 + 24*si2) - 40*POWER2(s12))) + 
     s0k*sik*(-64*s12*si2 + si1*(-64*s12 + 128*si2) - 
        64*POWER2(s12)) + 
     Q2inv*(si2*(72*s01*s12 + 64*s02*s12 + 12*s01*si2) + 
        si1*(104*s01*s12 + 32*s02*s12 + 
           (-20*s01 + 32*s02)*si1 + (-72*s01 - 32*s02)*si2)\
         + si0*(-32*s12*si2 + si1*(-32*s12 + 64*si2) - 
           32*POWER2(s12)) - 20*s01*POWER2(s12)))/
   (POWER2(s12 - si1 - si2)*POWER2(si2))
);
}
void DISMatrixElement::FillLorentzLqXqqq(
                        ) {
term[0][0] = (
(sik*(s0k*(64*s02*s12 + s01*(128*s02 + 64*s12)) + 
        (s01*(-64*s01 + 64*s02) + 64*s02*s12)*s2k + 
        s1k*(s01*(64*s02 + 64*s12) - 64*POWER2(s02))) + 
     Q2inv*((32*s02*s12 + s01*(64*s02 + 32*s12))*si0 + 
        (s01*(-32*s01 + 32*s02) + 32*s02*s12)*si2 + 
        si1*(s01*(32*s02 + 32*s12) - 32*POWER2(s02))))/
   (POWER2(s12)*POWER2(s01 + s02 + s12))
);
term[1][0] = (
(s1k*(64*s02*s1k*si2 + s2k*
         (64*s12*si0 - 64*s02*si1 - 64*s01*si2)) + 
     s0k*(s2k*(32*s12*si0 + 
           (32*s01 - 32*s02 - 32*s12)*si1) + 
        s1k*(32*s12*si0 + 
           (-32*s01 + 32*s02 - 32*s12)*si2) + 
        s0k*(-32*s12*si1 - 32*s12*si2)) + 
     sik*(s2k*(-32*s01*s12 - 32*s12*si0 - 32*s01*si1 + 
           32*s01*si2) + 
        s0k*(-64*s12*si0 + 64*s02*si1 + 64*s01*si2) + 
        s1k*(-32*s02*s12 - 32*s12*si0 + 32*s02*si1 - 
           32*s02*si2) + (32*s01*s12 + 32*s02*s12)*sik) + 
     Q2inv*(si2*(s01*(-16*s01 + 16*s02) + 16*s01*si2) + 
        si1*(16*s01*s02 + 16*s02*si1 + 
           (-16*s01 - 16*s02)*si2 - 16*POWER2(s02)) + 
        si0*(16*s01*s12 + 16*s02*s12 - 32*s12*si0 + 
           (32*s02 - 16*s12)*si1 + (32*s01 - 16*s12)*si2 + 
           32*POWER2(s12))) + 64*s01*si1*POWER2(s2k))/
   ((s01 + s02 + s12)*(s12 - si1 - si2)*POWER2(s12))
);
term[2][0] = (
(s1k*(-32*s02*s1k*si2 + s2k*
         ((32*s02 - 32*s12)*si0 + 32*s02*si1 + 32*s02*si2))\
      + s0k*(s2k*(-32*s12*si0 - 32*s02*si1 - 32*s01*si2) + 
        s1k*(32*s02*si0 + 32*s12*si2) + 
        s0k*(-32*s02*si1 + 32*s12*si2)) + 
     sik*(s2k*(32*s02*s12 - 32*s01*si0) + 
        s0k*(32*s01*s02 + 32*s02*s12 + 32*s12*si0 - 
           32*s02*si1) + 32*s01*s02*sik + 
        s1k*(32*s12*si0 - 32*s02*si1 + 32*s01*si2 - 
           64*POWER2(s02))) + 
     Q2inv*((32*s02*s12 + s01*(48*s02 + 16*s12))*si2 + 
        si1*(s02*(16*s02 + 16*s12) - 16*s02*si1 + 
           16*s01*si2) + 
        si0*(32*s01*s02 + 16*s02*s12 + 16*s12*si0 + 
           (-16*s02 + 16*s12)*si1 - 16*s01*si2 - 
           16*POWER2(s12))) + 
     (32*s01*si0 - 32*s02*si1)*POWER2(s2k))/
   (s12*(s01 + s02 + s12)*si0*(s02 - si0 - si2))
);
term[3][0] = (
(s1k*(s2k*((-32*s01 + 32*s12)*si0 - 32*s01*si1 - 
           32*s01*si2) + s1k*(-32*s02*si0 + 32*s01*si2)) + 
     s0k*(s2k*(-32*s01*si0 - 32*s12*si1) + 
        s1k*(32*s12*si0 + 32*s02*si1 + 32*s01*si2) + 
        s0k*(-32*s12*si1 + 32*s01*si2)) + 
     sik*(s1k*(-32*s01*s12 + 32*s02*si0) + 
        s0k*(s01*(-32*s02 - 32*s12) - 32*s12*si0 + 
           32*s01*si2) - 32*s01*s02*sik + 
        s2k*(-32*s12*si0 - 32*s02*si1 + 32*s01*si2 + 
           64*POWER2(s01))) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s12) + 16*s01*si2) + 
        si1*(s01*(-48*s02 - 32*s12) - 16*s02*s12 - 
           16*s02*si2) + 
        si0*(s01*(-32*s02 - 16*s12) - 16*s12*si0 + 
           16*s02*si1 + (16*s01 - 16*s12)*si2 + 
           16*POWER2(s12))) + 32*s01*si1*POWER2(s2k))/
   (s12*(s01 + s02 + s12)*si0*(s01 - si0 - si1))
);
term[4][0] = (
(s1k*s2k*(32*s02*si0 + 32*s02*si2) + 
     Q2inv*((32*s01*s02 + 32*s02*s12)*si0 + 
        (32*s01*s02 + 32*s02*s12)*si2) + 
     s0k*(-32*s02*s0k*si1 - 64*s02*s2k*si1 + 
        s1k*(32*s02*si0 + 32*s02*si2)) + 
     sik*(s0k*(32*s01*s02 + 32*s02*s12) + 
        (32*s01*s02 + 32*s02*s12)*s2k - 64*s1k*POWER2(s02))\
      - 32*s02*si1*POWER2(s2k))/
   (s12*(s01 + s02 + s12)*(s02 - si0 - si2)*si2)
);
term[5][0] = (
(s0k*(32*s0k*s12*si1 + (-32*s01 + 32*s12)*s2k*si1 + 
        s1k*(-32*s12*si0 + 32*s01*si2)) + 
     s1k*(-32*s02*s1k*si2 + 
        s2k*(-32*s12*si0 + 32*s02*si1 + 32*s01*si2)) + 
     sik*(s2k*(32*s12*si0 - 32*s01*si2) + 
        s0k*(32*s12*si0 - 32*s02*si1 - 32*s01*si2) + 
        s1k*(32*s02*s12 + 32*s02*si2) - 32*s02*s12*sik) + 
     Q2inv*(si2*(s01*(16*s01 + 16*s12) - 16*s01*si2) + 
        si1*(-16*s01*s02 - 16*s02*s12 + 16*s02*si2) + 
        si0*(-16*s01*s12 + 16*s12*si0 - 16*s02*si1 + 
           (-16*s01 + 16*s12)*si2 - 16*POWER2(s12))) - 
     32*s01*si1*POWER2(s2k))/
   (s12*(s01 + s02 + s12)*(s12 - si1 - si2)*si2)
);
term[6][0] = (
(sik*(s0k*(64*s01*s02 + 64*s02*s12) + 
        (64*s01*s02 + 64*s02*s12)*s2k - 64*s1k*POWER2(s02))\
      + Q2inv*((32*s01*s02 + 32*s02*s12)*si0 + 
        (32*s01*s02 + 32*s02*s12)*si2 - 32*si1*POWER2(s02)))
    /(s01*s12*POWER2(s01 + s02 + s12))
);
term[7][0] = (
(s1k*(32*s02*s1k*si0 + s2k*(-32*s12*si0 + 32*s01*si2)) + 
     s0k*(32*s0k*s12*si1 + (-32*s01 + 32*s12)*s2k*si1 + 
        s1k*(-32*s12*si0 - 32*s02*si1 + 32*s01*si2)) + 
     sik*(s1k*(-32*s01*s02 - 32*s02*si0) + 
        s0k*(32*s12*si0 - 32*s01*si2) + 
        s2k*(32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        32*s01*s02*sik) + 
     Q2inv*(si2*(s01*(16*s01 + 16*s12) - 16*s01*si2) + 
        si1*(16*s01*s02 + 16*s02*s12 + 16*s02*si2) + 
        si0*(-16*s01*s12 + 16*s12*si0 - 16*s02*si1 + 
           (-16*s01 + 16*s12)*si2 - 16*POWER2(s12))) - 
     32*s01*si1*POWER2(s2k))/
   (s01*s12*(s01 + s02 + s12)*(s01 - si0 - si1))
);
term[1][1] = (
(Q2inv*(si2*(32*s02*s12 + 32*s01*si2) + 
        si1*(32*s01*s12 + 32*s02*si1 + 
           (-32*s01 - 32*s02)*si2) + 
        si0*(-32*s12*si2 + si1*(-32*s12 + 64*si2))) + 
     s0k*(-64*s12*si2 + si1*(-64*s12 + 128*si2))*sik + 
     s0k*(s2k*(si1*(64*si1 - 64*si2) + 64*s12*si2) + 
        s1k*(si1*(64*s12 - 64*si2) + 64*POWER2(si2))))/
   (POWER2(s12)*POWER2(s12 - si1 - si2))
);
term[2][1] = (
(s1k*(-32*s02*s1k*si2 + s2k*
         (si0*(-32*s12 - 32*si2) + (32*s01 + 32*s02)*si2))\
      + s0k*(32*s0k*si1*si2 + 
        s2k*(-32*si0*si1 + 32*s12*si2) + 
        s1k*(-32*s12*si0 - 32*s02*si1 + 
           si2*(32*s01 + 64*si2))) + 
     sik*(s2k*(32*s12*si0 + 32*s02*si1 + 32*s01*si2) + 
        s1k*(32*s02*s12 + 32*si0*si2) + 
        s0k*(32*s12*si0 + (-32*s01 - 32*s12)*si2 + 
           32*si1*si2) + (-32*s02*s12 - 32*s01*si2)*sik) + 
     Q2inv*(si2*(32*s02*s12 + s01*(16*s01 + 16*s12) - 
           16*s01*si2) + 
        si1*(-16*s01*s02 + 16*s02*s12 - 48*s02*si2) + 
        si0*(-16*s01*s12 + 16*s12*si0 + 
           (-16*s01 - 16*s12)*si2 + 
           si1*(-16*s02 + 32*si2) - 16*POWER2(s12))) + 
     (-32*si0*si1 - 32*s01*si2)*POWER2(s2k))/
   (s12*si0*(s12 - si1 - si2)*(-s02 + si0 + si2))
);
term[3][1] = (
(s1k*(s2k*((-32*s01 - 32*s02)*si1 + 
           si0*(32*s12 + 32*si1)) + 
        s1k*(32*s02*si1 + 32*si0*si2)) + 
     s0k*(-32*s0k*si1*si2 + 
        s2k*(32*s12*si0 + (-32*s02 - 64*si1)*si1 + 
           32*s01*si2) + s1k*(-32*s12*si1 + 32*si0*si2)) + 
     sik*(s2k*(-32*s01*s12 - 32*si0*si1) + 
        s0k*(-32*s12*si0 + 
           si1*(32*s02 + 32*s12 - 32*si2)) + 
        s1k*(-32*s12*si0 - 32*s02*si1 - 32*s01*si2) + 
        (32*s01*s12 + 32*s02*si1)*sik) + 
     Q2inv*(s01*(16*s02 - 16*s12)*si2 + 
        si1*(s02*(-16*s02 - 16*s12) - 32*s01*s12 + 
           16*s02*si1 + 48*s01*si2) + 
        si0*(16*s02*s12 - 16*s12*si0 + 
           si1*(16*s02 + 16*s12 - 32*si2) + 16*s01*si2 + 
           16*POWER2(s12))) + 32*s01*si1*POWER2(s2k))/
   (s12*si0*(-s01 + si0 + si1)*(s12 - si1 - si2))
);
term[4][1] = (
(s1k*(32*s02*s1k*si2 + s2k*(32*s12*si0 - 32*s01*si2)) + 
     s0k*(-32*s0k*si1*si2 + 
        s1k*(32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        s2k*(32*si0*si1 - 32*si1*si2)) + 
     sik*(s0k*(-32*s12*si0 + 32*s01*si2) + 
        s2k*(-32*s12*si0 - 32*s02*si1 + 32*s01*si2) + 
        s1k*(-32*s02*s12 - 32*s02*si2) + 32*s02*s12*sik) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s12) + 16*s01*si2) + 
        si1*(16*s01*s02 - 16*s02*s12 + 16*s02*si2) + 
        si0*(16*s01*s12 - 16*s12*si0 + 16*s02*si1 + 
           (16*s01 - 16*s12)*si2 + 16*POWER2(s12))) + 
     32*si0*si1*POWER2(s2k))/
   (s12*(s02 - si0 - si2)*(s12 - si1 - si2)*si2)
);
term[5][1] = (
(Q2inv*(si0*si1*(32*s12 - 32*si2) + 
        si1*(-32*s01*s12 - 32*s02*si1 + 32*s01*si2)) + 
     s0k*si1*(64*s12 - 64*si2)*sik + 
     s0k*(s1k*si1*(-64*s12 + 64*si2) - 64*s2k*POWER2(si1)))/
   (s12*si2*POWER2(s12 - si1 - si2))
);
term[6][1] = (
(s0k*(-32*s0k*s12*si1 + (32*s01 - 32*s12)*s2k*si1 + 
        s1k*(32*s12*si0 - 32*s01*si2)) + 
     s1k*(32*s02*s1k*si2 + 
        s2k*(32*s12*si0 - 32*s02*si1 - 32*s01*si2)) + 
     sik*(s2k*(-32*s12*si0 + 32*s01*si2) + 
        s0k*(-32*s12*si0 + 32*s02*si1 + 32*s01*si2) + 
        s1k*(-32*s02*s12 - 32*s02*si2) + 32*s02*s12*sik) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s12) + 16*s01*si2) + 
        si1*(16*s01*s02 + 16*s02*s12 - 16*s02*si2) + 
        si0*(16*s01*s12 - 16*s12*si0 + 16*s02*si1 + 
           (16*s01 - 16*s12)*si2 + 16*POWER2(s12))) + 
     32*s01*si1*POWER2(s2k))/
   (s01*s12*(s01 + s02 + s12)*(s12 - si1 - si2))
);
term[7][1] = (
(s1k*(32*s02*s1k*si1 + s2k*(-32*s01*si1 + 32*si0*si1)) + 
     Q2inv*(si0*si1*(32*s12 - 32*si2) + 
        si1*(-32*s01*s12 + 32*s01*si2)) + 
     sik*(-64*s02*s1k*si1 + s2k*(32*s01*si1 - 32*si0*si1) + 
        s0k*si1*(32*s12 - 32*si2) + 32*s02*si1*sik) + 
     s0k*(s1k*si1*(-32*s12 + 32*si2) - 64*s2k*POWER2(si1)))/
   (s01*s12*(s01 - si0 - si1)*(s12 - si1 - si2))
);
term[2][2] = (
(s1k*s2k*(si0*(64*s02 - 64*si2) + 128*s02*si2) + 
     s0k*s1k*(64*s02*si0 + si2*(64*s02 + 64*si2)) + 
     Q2inv*(si2*(32*s01*s02 + 64*s02*s12 + 32*s01*si2) + 
        si0*(32*s01*s02 + 32*s02*s12 - 32*s12*si2 + 
           32*si1*si2) + si1*(-32*s02*si2 - 32*POWER2(s02)))
       + s1k*sik*(-64*s02*si2 + 64*si0*si2 - 64*POWER2(s02))
     )/(POWER2(si0)*POWER2(-s02 + si0 + si2))
);
term[3][2] = (
(s1k*(s1k*si0*(-32*s02 + 32*si2) + 
        s2k*(64*s12*si0 - 64*s02*si1 - 64*s01*si2)) + 
     s0k*(s2k*(si0*(32*s12 + 32*si1) + 
           si1*(-32*s02 - 32*si2)) - 64*s0k*si1*si2 + 
        s1k*(-32*s01*si2 - 32*si1*si2 + 
           si0*(32*s12 + 32*si2))) + 
     sik*(s1k*(32*s01*s02 + (32*s02 - 32*s12)*si0 + 
           32*s02*si1) + 
        s2k*(32*s01*s02 + (32*s01 - 32*s12)*si0 + 
           32*s01*si2) + 
        s0k*(-64*s12*si0 + 64*s02*si1 + 64*s01*si2) - 
        64*s01*s02*sik) + 
     Q2inv*(si2*(s01*(-16*s01 + 16*s02 - 32*s12) + 
           16*s01*si2) + 
        si1*(16*s01*s02 + s02*(-16*s02 - 32*s12) + 
           16*s02*si1 + (-16*s01 - 16*s02)*si2) + 
        si0*(16*s01*s12 + 16*s02*s12 - 32*s12*si0 - 
           16*s12*si1 - 16*s12*si2 + 32*POWER2(s12))) + 
     si0*(-32*s01 + 32*si1)*POWER2(s2k))/
   ((-s01 + si0 + si1)*(-s02 + si0 + si2)*POWER2(si0))
);
term[4][2] = (
(s0k*s1k*(64*s02*si0 + 64*s02*si2) + 
     s1k*s2k*(64*s02*si0 + 64*s02*si2) - 
     64*s1k*sik*POWER2(s02) + 
     Q2inv*((32*s01*s02 + 32*s02*s12)*si0 + 
        (32*s01*s02 + 32*s02*s12)*si2 - 32*si1*POWER2(s02)))
    /(si0*si2*POWER2(-s02 + si0 + si2))
);
term[5][2] = (
(s1k*(-32*s02*s1k*si2 + s2k*(-32*s12*si0 + 32*s01*si2)) + 
     s0k*(32*s0k*si1*si2 + 
        s1k*(-32*s12*si0 - 32*s02*si1 + 32*s01*si2) + 
        s2k*(-32*si0*si1 + 32*si1*si2)) + 
     sik*(s0k*(32*s12*si0 - 32*s01*si2) + 
        s2k*(32*s12*si0 + 32*s02*si1 - 32*s01*si2) + 
        s1k*(32*s02*s12 + 32*s02*si2) - 32*s02*s12*sik) + 
     Q2inv*(si2*(s01*(16*s01 + 16*s12) - 16*s01*si2) + 
        si1*(-16*s01*s02 + 16*s02*s12 - 16*s02*si2) + 
        si0*(-16*s01*s12 + 16*s12*si0 - 16*s02*si1 + 
           (-16*s01 + 16*s12)*si2 - 16*POWER2(s12))) - 
     32*si0*si1*POWER2(s2k))/
   (si0*si2*(-s02 + si0 + si2)*(-s12 + si1 + si2))
);
term[6][2] = (
(s1k*s2k*(32*s02*si0 + 32*s02*si2) + 
     Q2inv*((32*s01*s02 + 32*s02*s12)*si0 + 
        (32*s01*s02 + 32*s02*s12)*si2) + 
     s0k*(-32*s02*s0k*si1 - 64*s02*s2k*si1 + 
        s1k*(32*s02*si0 + 32*s02*si2)) + 
     sik*(s0k*(32*s01*s02 + 32*s02*s12) + 
        (32*s01*s02 + 32*s02*s12)*s2k - 64*s1k*POWER2(s02))\
      - 32*s02*si1*POWER2(s2k))/
   (s01*(s01 + s02 + s12)*si0*(s02 - si0 - si2))
);
term[7][2] = (
(s1k*(-32*s02*s1k*si0 + s2k*
         (32*s12*si0 - 32*s02*si1 - 32*s01*si2)) + 
     s0k*(-32*s0k*si1*si2 + s1k*(32*s12*si0 - 32*s01*si2) + 
        s2k*(32*si0*si1 - 32*si1*si2)) + 
     sik*(s1k*(32*s01*s02 + 32*s02*si0) + 
        s2k*(-32*s12*si0 + 32*s01*si2) + 
        s0k*(-32*s12*si0 + 32*s02*si1 + 32*s01*si2) - 
        32*s01*s02*sik) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s12) + 16*s01*si2) + 
        si1*(16*s01*s02 - 16*s02*s12 - 16*s02*si2) + 
        si0*(16*s01*s12 - 16*s12*si0 - 16*s02*si1 + 
           (16*s01 - 16*s12)*si2 + 16*POWER2(s12))) + 
     32*si0*si1*POWER2(s2k))/
   (s01*si0*(s01 - si0 - si1)*(-s02 + si0 + si2))
);
term[3][3] = (
(s1k*s2k*(si0*(64*s01 - 64*si1) + 128*s01*si1) + 
     s0k*s2k*(64*s01*si0 + si1*(64*s01 + 64*si1)) + 
     s2k*sik*(-64*s01*si1 + 64*si0*si1 - 64*POWER2(s01)) + 
     Q2inv*(si1*(s01*(32*s02 + 64*s12) + 32*s02*si1 - 
           32*s01*si2) + 
        si0*(s01*(32*s02 + 32*s12) + 
           si1*(-32*s12 + 32*si2)) - 32*si2*POWER2(s01)))/
   (POWER2(si0)*POWER2(-s01 + si0 + si1))
);
term[4][3] = (
(s1k*(-32*s02*s1k*si0 + s2k*
         (32*s12*si0 - 32*s02*si1 - 32*s01*si2)) + 
     s0k*(-32*s0k*si1*si2 + s1k*(32*s12*si0 - 32*s01*si2) + 
        s2k*(32*si0*si1 - 32*si1*si2)) + 
     sik*(s1k*(32*s01*s02 + 32*s02*si0) + 
        s2k*(-32*s12*si0 + 32*s01*si2) + 
        s0k*(-32*s12*si0 + 32*s02*si1 + 32*s01*si2) - 
        32*s01*s02*sik) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s12) + 16*s01*si2) + 
        si1*(16*s01*s02 - 16*s02*s12 - 16*s02*si2) + 
        si0*(16*s01*s12 - 16*s12*si0 - 16*s02*si1 + 
           (16*s01 - 16*s12)*si2 + 16*POWER2(s12))) + 
     32*si0*si1*POWER2(s2k))/
   (si0*(-s01 + si0 + si1)*si2*(-s02 + si0 + si2))
);
term[5][3] = (
(s1k*(32*s02*s1k*si1 + s2k*(-32*s01*si1 + 32*si0*si1)) + 
     Q2inv*(si0*si1*(32*s12 - 32*si2) + 
        si1*(-32*s01*s12 + 32*s01*si2)) + 
     sik*(-64*s02*s1k*si1 + s2k*(32*s01*si1 - 32*si0*si1) + 
        s0k*si1*(32*s12 - 32*si2) + 32*s02*si1*sik) + 
     s0k*(s1k*si1*(-32*s12 + 32*si2) - 64*s2k*POWER2(si1)))/
   (si0*(-s01 + si0 + si1)*si2*(-s12 + si1 + si2))
);
term[6][3] = (
(s1k*(-32*s02*s1k*si0 + s2k*(32*s12*si0 - 32*s01*si2)) + 
     s0k*(-32*s0k*s12*si1 + (32*s01 - 32*s12)*s2k*si1 + 
        s1k*(32*s12*si0 + 32*s02*si1 - 32*s01*si2)) + 
     sik*(s1k*(32*s01*s02 + 32*s02*si0) + 
        s0k*(-32*s12*si0 + 32*s01*si2) + 
        s2k*(-32*s12*si0 - 32*s02*si1 + 32*s01*si2) - 
        32*s01*s02*sik) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s12) + 16*s01*si2) + 
        si1*(-16*s01*s02 - 16*s02*s12 - 16*s02*si2) + 
        si0*(16*s01*s12 - 16*s12*si0 + 16*s02*si1 + 
           (16*s01 - 16*s12)*si2 + 16*POWER2(s12))) + 
     32*s01*si1*POWER2(s2k))/
   (s01*(s01 + s02 + s12)*si0*(s01 - si0 - si1))
);
term[7][3] = (
(s1k*s2k*(-64*s01*si1 + 64*si0*si1) + 
     Q2inv*(si0*si1*(32*s12 - 32*si2) + 
        si1*(-32*s01*s12 - 32*s02*si1 + 32*s01*si2)) + 
     s2k*(64*s01*si1 - 64*si0*si1)*sik - 
     64*s0k*s2k*POWER2(si1))/
   (s01*si0*POWER2(s01 - si0 - si1))
);
term[4][4] = (
(s1k*s2k*(si0*(64*s02 + 64*si0) + 64*s02*si2) + 
     s0k*s1k*(si0*(128*s02 - 64*si2) + 64*s02*si2) + 
     s1k*sik*(si0*(-64*s02 + 64*si2) - 64*POWER2(s02)) + 
     Q2inv*((32*s01*s02 + 32*s02*s12)*si2 + 
        si0*(64*s01*s02 + 32*s02*s12 + 32*s12*si0 - 
           32*s01*si2 + si1*(-32*s02 + 32*si2)) - 
        32*si1*POWER2(s02)))/
   (POWER2(s02 - si0 - si2)*POWER2(si2))
);
term[5][4] = (
(s1k*(s2k*(si0*(32*s12 + 32*si1 - 32*si2) - 32*s01*si2) + 
        s1k*(32*s02*si2 - 32*si0*si2)) + 
     s0k*(s1k*(64*s12*si0 + 64*s02*si1 - 64*s01*si2) + 
        s2k*(32*si0*si1 + si1*(32*s02 - 32*si2) - 
           32*s01*si2) + s0k*(32*s12*si2 - 32*si1*si2)) + 
     sik*(s2k*(-64*s12*si0 - 64*s02*si1 + 64*s01*si2) + 
        s1k*(-32*s02*s12 - 32*s02*si1 + 
           (32*s01 - 32*s02)*si2) + 
        s0k*(-32*s02*s12 - 32*s12*si0 + 
           (32*s01 - 32*s12)*si2) + 64*s02*s12*sik) + 
     Q2inv*(si1*(32*s01*s02 + s02*(16*s02 - 16*s12) - 
           16*s02*si1 + 16*s01*si2) + 
        si2*(s01*(-32*s01 - 16*s02 - 16*s12) + 
           32*s01*si2) + 
        si0*(32*s01*s12 - 16*s02*s12 - 16*s12*si0 + 
           (16*s02 + 16*s12)*si1 + 16*s01*si2 + 
           16*POWER2(s12))) + 64*si0*si1*POWER2(s2k))/
   ((s02 - si0 - si2)*(-s12 + si1 + si2)*POWER2(si2))
);
term[6][4] = (
(s1k*(-32*s02*s1k*si0 + s2k*(32*s01*si0 + 32*s02*si2)) + 
     s0k*(s2k*(-32*s12*si0 - 32*s02*si1 - 32*s01*si2) + 
        s1k*(32*s02*si0 + 32*s02*si1 + 
           (-32*s01 + 32*s02)*si2) + 
        s0k*(-32*s02*si1 + 32*s12*si2)) + 
     sik*(s2k*(32*s01*s02 + 32*s02*s12 - 32*s02*si1 + 
           32*s01*si2) + s0k*(32*s01*s02 - 32*s12*si2) + 
        32*s02*s12*sik + 
        s1k*(32*s12*si0 - 32*s02*si1 + 32*s01*si2 - 
           64*POWER2(s02))) + 
     Q2inv*(si2*(s01*(-16*s01 + 16*s02) + 32*s02*s12 + 
           16*s01*si2) + 
        si0*(48*s02*s12 + s01*(32*s02 + 16*s12) + 
           16*s12*si1 - 16*s12*si2) + 
        si1*(16*s01*s02 - 16*s02*si1 + 
           (16*s01 - 16*s02)*si2 + 16*POWER2(s02))) + 
     (32*s01*si0 - 32*s02*si1)*POWER2(s2k))/
   (s01*(s01 + s02 + s12)*(s02 - si0 - si2)*si2)
);
term[7][4] = (
(s1k*(-32*s02*s1k*si0 + s2k*
         (si0*(32*s12 + 64*si0) - 32*s02*si1 - 32*s01*si2))\
      + s0k*(s1k*(si0*(32*s02 + 32*s12 - 32*si2) - 
           32*s01*si2) + s2k*(32*s01*si0 - 32*si1*si2) + 
        s0k*(-32*s12*si0 - 32*si1*si2)) + 
     sik*(s0k*(32*s12*si0 + 32*s02*si1 + 32*s01*si2) + 
        s2k*(si0*(-32*s01 - 32*s12 + 32*si1) + 
           32*s01*si2) + s1k*(32*s01*s02 + 32*si0*si2) + 
        (-32*s01*s02 - 32*s12*si0)*sik) + 
     Q2inv*(si2*(s01*(-16*s01 - 16*s12) + 16*s01*si2) + 
        si1*(16*s01*s02 - 16*s02*s12 - 16*s02*si2) + 
        si0*(s01*(32*s02 + 16*s12) - 16*s12*si0 + 
           (-16*s01 - 16*s12)*si2 + 
           si1*(-48*s02 + 32*si2) + 16*POWER2(s12))) + 
     32*si0*si1*POWER2(s2k))/
   (s01*(s01 - si0 - si1)*si2*(-s02 + si0 + si2))
);
term[5][5] = (
(s0k*(s2k*(si1*(64*s12 + 64*si1) + 64*s12*si2) + 
        s1k*(si1*(128*s12 - 64*si2) + 64*s12*si2)) + 
     Q2inv*((32*s01*s12 + 32*s02*s12)*si2 + 
        si1*(64*s01*s12 + 32*s02*s12 + 32*s02*si1 - 
           32*s01*si2) + 
        si0*(si1*(-32*s12 + 32*si2) - 32*POWER2(s12))) + 
     s0k*sik*(si1*(-64*s12 + 64*si2) - 64*POWER2(s12)))/
   (POWER2(s12 - si1 - si2)*POWER2(si2))
);
term[6][5] = (
(s1k*(s2k*(32*s12*si0 + 32*s02*si1 + 32*s01*si2) + 
        s1k*(32*s12*si0 - 32*s02*si2)) + 
     s0k*(32*s0k*s12*si1 + 
        s1k*(-32*s12*si0 - 32*s12*si1 + 
           (32*s01 - 32*s12)*si2) + 
        s2k*(-32*s01*si1 - 32*s12*si2)) + 
     Q2inv*(si2*(s01*(16*s01 - 16*s12) - 32*s02*s12 - 
           16*s01*si2) + 
        si1*(s01*(-16*s02 - 32*s12) - 48*s02*s12 + 
           16*s02*si2) + 
        si0*(-16*s01*s12 + 16*s12*si0 - 16*s02*si1 + 
           (-16*s01 + 16*s12)*si2 - 16*POWER2(s12))) + 
     sik*(s2k*(-32*s01*s12 - 32*s02*s12 + 32*s12*si0 - 
           32*s01*si2) + s1k*(-32*s01*s12 + 32*s02*si2) - 
        32*s02*s12*sik + 
        s0k*(32*s12*si0 - 32*s02*si1 - 32*s01*si2 + 
           64*POWER2(s12))) + 
     (32*s12*si0 - 32*s01*si1)*POWER2(s2k))/
   (s01*(s01 + s02 + s12)*(s12 - si1 - si2)*si2)
);
term[7][5] = (
(s1k*(s2k*(-32*s01*si1 + 32*si0*si2) + 
        s1k*(32*s02*si1 + 32*si0*si2)) + 
     s0k*(32*s0k*s12*si1 + 
        s2k*(32*s12*si0 + (-32*s02 - 64*si1)*si1 + 
           32*s01*si2) + 
        s1k*(32*s01*si2 + si1*(-32*s02 - 32*s12 + 32*si2)))\
      + sik*(s1k*(-32*s12*si0 - 32*s02*si1 - 32*s01*si2) + 
        s2k*((32*s01 + 32*s02)*si1 - 32*si0*si1 - 
           32*s01*si2) + s0k*(-32*s01*s12 - 32*si1*si2) + 
        (32*s01*s12 + 32*s02*si1)*sik) + 
     Q2inv*(si2*(s01*(16*s01 + 16*s02) - 16*s01*si2) + 
        si0*(-16*s01*s12 + 16*s02*s12 + 
           si1*(48*s12 - 32*si2) + 16*s12*si2) + 
        si1*(s01*(-16*s02 - 32*s12) + 16*s02*si1 + 
           (16*s01 + 16*s02)*si2 - 16*POWER2(s02))) - 
     32*si0*si1*POWER2(s2k))/
   (s01*(s01 - si0 - si1)*si2*(-s12 + si1 + si2))
);
term[6][6] = (
(sik*((64*s01*s12 + s02*(-64*s02 + 64*s12))*s1k + 
        (128*s02*s12 + s01*(64*s02 + 64*s12))*s2k + 
        s0k*(64*s01*s02 + 64*s02*s12 - 64*POWER2(s12))) + 
     Q2inv*((32*s01*s12 + s02*(-32*s02 + 32*s12))*si1 + 
        (64*s02*s12 + s01*(32*s02 + 32*s12))*si2 + 
        si0*(32*s01*s02 + 32*s02*s12 - 32*POWER2(s12))))/
   (POWER2(s01)*POWER2(s01 + s02 + s12))
);
term[7][6] = (
(s1k*(64*s02*s1k*si0 + s2k*
         ((-32*s01 + 32*s02 - 32*s12)*si0 + 32*s01*si2)) + 
     s0k*(64*s0k*s12*si1 + 
        s2k*((-32*s01 - 32*s02 + 32*s12)*si1 + 
           32*s01*si2) + 
        s1k*(-64*s12*si0 - 64*s02*si1 + 64*s01*si2)) + 
     sik*(s2k*(64*s12*si0 + 64*s02*si1 - 64*s01*si2) + 
        s1k*(-32*s01*s02 - 32*s02*si0 + 32*s02*si1 - 
           32*s01*si2) + 
        s0k*(-32*s01*s12 + 32*s12*si0 - 32*s12*si1 - 
           32*s01*si2) + s01*(32*s02 + 32*s12)*sik) + 
     Q2inv*(si2*(s01*(32*s01 + 16*s02 + 16*s12) - 
           32*s01*si2) + 
        si1*(s02*(-16*s02 + 16*s12) + 16*s02*si1 + 
           (-16*s01 + 32*s02)*si2) + 
        si0*(16*s02*s12 + 16*s12*si0 + 
           (-16*s02 - 16*s12)*si1 + 
           (-16*s01 + 32*s12)*si2 - 16*POWER2(s12))) + 
     (-32*s01*si0 - 32*s01*si1)*POWER2(s2k))/
   ((s01 + s02 + s12)*(s01 - si0 - si1)*POWER2(s01))
);
term[7][7] = (
(s1k*s2k*(si0*(64*si0 - 64*si1) + 64*s01*si1) + 
     Q2inv*(si1*(32*s01*s12 + 32*s02*si1 - 32*s01*si2) + 
        si0*(32*s01*s02 + 32*s12*si0 - 32*s01*si2 + 
           si1*(-32*s02 - 32*s12 + 64*si2))) + 
     s2k*(-64*s01*si1 + si0*(-64*s01 + 128*si1))*sik + 
     s0k*s2k*(si0*(64*s01 - 64*si1) + 64*POWER2(si1)))/
   (POWER2(s01)*POWER2(s01 - si0 - si1))
);
}
// DIS subtraction terms
// Created on {1997, 12, 16, 13, 32, 0}
#include <math.h>
#include <stdio.h>
#include "mathdefs.h"
#include "global.h"
#include "proc.h"
#include "me.h"
double DISMatrixElement::GetLorentzSubtraction(
   int list, LimitType sub_type, int i, int j
                       ) {
   double term_L=0;
   switch (list) {
      case 1:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
(16*vi0*(2*s02*s0k*s2k - 2*Q2inv*s02*si0 - 2*s0k*s2k*si0 + 
       2*Q2inv*si0*si2 - 2*s02*s0k*sik + 4*s0k*si0*sik + 
       2*s2k*si0*sik + 2*s0k*si2*sik + 2*s2k*si2*sik + 
       Q2inv*POWER2(s02) + 2*si2*POWER2(s0k) + 
       2*Q2inv*POWER2(si0) + Q2inv*POWER2(si2) - 
       2*s02*POWER2(sik)))/(s02*si2*v01*vi1)
                     );
                     break;
                  case 2:
                     term_L=(
(16*vi0*(2*s01*s0k*s1k - 2*Q2inv*s01*si0 - 2*s0k*s1k*si0 + 
       2*Q2inv*si0*si1 - 2*s01*s0k*sik + 4*s0k*si0*sik + 
       2*s1k*si0*sik + 2*s0k*si1*sik + 2*s1k*si1*sik + 
       Q2inv*POWER2(s01) + 2*si1*POWER2(s0k) + 
       2*Q2inv*POWER2(si0) + Q2inv*POWER2(si1) - 
       2*s01*POWER2(sik)))/(s01*si1*v02*vi2)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
(-32*(1 + POWER2(lambda))*
     (-2*Q2inv*s2h*sih + 2*Q2inv*si2*sih + 2*s2k*si2*sik + 
       2*s2k*sih*sik + 2*s2h*s2k*skh - 2*s2k*sih*skh - 
       2*s2h*sik*skh + 2*si2*sik*skh + 4*sih*sik*skh + 
       Q2inv*POWER2(s2h) + Q2inv*POWER2(si2) + 
       2*Q2inv*POWER2(sih) - 2*s2h*POWER2(sik) + 
       2*si2*POWER2(skh)))/((-1 + lambda)*s2h*si2)
                           );
                           break;
                        case 2:
                           term_L=(
(-32*(1 + POWER2(lambda))*
     (-2*Q2inv*s1h*sih + 2*Q2inv*si1*sih + 2*s1k*si1*sik + 
       2*s1k*sih*sik + 2*s1h*s1k*skh - 2*s1k*sih*skh - 
       2*s1h*sik*skh + 2*si1*sik*skh + 4*sih*sik*skh + 
       Q2inv*POWER2(s1h) + Q2inv*POWER2(si1) + 
       2*Q2inv*POWER2(sih) - 2*s1h*POWER2(sik) + 
       2*si1*POWER2(skh)))/((-1 + lambda)*s1h*si1)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
(-32*(1 + POWER2(eta))*(-2*Q2inv*s02*s0h + 
       2*Q2inv*s0h*s2h + 2*s02*s0k*s2k - 2*s0h*s0k*s2k - 
       2*s02*s0k*skh + 4*s0h*s0k*skh + 2*s0k*s2h*skh + 
       2*s0h*s2k*skh + 2*s2h*s2k*skh + Q2inv*POWER2(s02) + 
       2*Q2inv*POWER2(s0h) + 2*s2h*POWER2(s0k) + 
       Q2inv*POWER2(s2h) - 2*s02*POWER2(skh)))/
   ((-1 + eta)*eta*s02*s2h)
                           );
                           break;
                        case 0:
                           term_L=(
(32*(2 - 2*lambda + POWER2(lambda))*
     (-2*Q2inv*s2h*sih + 2*Q2inv*si2*sih + 2*s2k*si2*sik + 
       2*s2k*sih*sik + 2*s2h*s2k*skh - 2*s2k*sih*skh - 
       2*s2h*sik*skh + 2*si2*sik*skh + 4*sih*sik*skh + 
       Q2inv*POWER2(s2h) + Q2inv*POWER2(si2) + 
       2*Q2inv*POWER2(sih) - 2*s2h*POWER2(sik) + 
       2*si2*POWER2(skh)))/(lambda*s2h*si2)
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
(-32*(1 + POWER2(eta))*(-2*Q2inv*s01*s0h + 
       2*Q2inv*s0h*s1h + 2*s01*s0k*s1k - 2*s0h*s0k*s1k - 
       2*s01*s0k*skh + 4*s0h*s0k*skh + 2*s0k*s1h*skh + 
       2*s0h*s1k*skh + 2*s1h*s1k*skh + Q2inv*POWER2(s01) + 
       2*Q2inv*POWER2(s0h) + 2*s1h*POWER2(s0k) + 
       Q2inv*POWER2(s1h) - 2*s01*POWER2(skh)))/
   ((-1 + eta)*eta*s01*s1h)
                           );
                           break;
                        case 0:
                           term_L=(
(32*(2 - 2*lambda + POWER2(lambda))*
     (-2*Q2inv*s1h*sih + 2*Q2inv*si1*sih + 2*s1k*si1*sik + 
       2*s1k*sih*sik + 2*s1h*s1k*skh - 2*s1k*sih*skh - 
       2*s1h*sik*skh + 2*si1*sik*skh + 4*sih*sik*skh + 
       Q2inv*POWER2(s1h) + Q2inv*POWER2(si1) + 
       2*Q2inv*POWER2(sih) - 2*s1h*POWER2(sik) + 
       2*si1*POWER2(skh)))/(lambda*s1h*si1)
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
(64*(-2*Q2inv*s02*s0h + 2*Q2inv*s0h*s2h + 2*s02*s0k*s2k - 
       2*s0h*s0k*s2k - 2*s02*s0k*skh + 4*s0h*s0k*skh + 
       2*s0k*s2h*skh + 2*s0h*s2k*skh + 2*s2h*s2k*skh + 
       Q2inv*POWER2(s02) + 2*Q2inv*POWER2(s0h) + 
       2*s2h*POWER2(s0k) + Q2inv*POWER2(s2h) - 
       2*s02*POWER2(skh)))/(s02*s2h)
                           );
                           break;
                        case 0:
                           term_L=(
(64*(-2*Q2inv*s2h*sih + 2*Q2inv*si2*sih + 2*s2k*si2*sik + 
       2*s2k*sih*sik + 2*s2h*s2k*skh - 2*s2k*sih*skh - 
       2*s2h*sik*skh + 2*si2*sik*skh + 4*sih*sik*skh + 
       Q2inv*POWER2(s2h) + Q2inv*POWER2(si2) + 
       2*Q2inv*POWER2(sih) - 2*s2h*POWER2(sik) + 
       2*si2*POWER2(skh)))/(s2h*si2)
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
(64*(-2*Q2inv*s01*s0h + 2*Q2inv*s0h*s1h + 2*s01*s0k*s1k - 
       2*s0h*s0k*s1k - 2*s01*s0k*skh + 4*s0h*s0k*skh + 
       2*s0k*s1h*skh + 2*s0h*s1k*skh + 2*s1h*s1k*skh + 
       Q2inv*POWER2(s01) + 2*Q2inv*POWER2(s0h) + 
       2*s1h*POWER2(s0k) + Q2inv*POWER2(s1h) - 
       2*s01*POWER2(skh)))/(s01*s1h)
                           );
                           break;
                        case 0:
                           term_L=(
(64*(-2*Q2inv*s1h*sih + 2*Q2inv*si1*sih + 2*s1k*si1*sik + 
       2*s1k*sih*sik + 2*s1h*s1k*skh - 2*s1k*sih*skh - 
       2*s1h*sik*skh + 2*si1*sik*skh + 4*sih*sik*skh + 
       Q2inv*POWER2(s1h) + Q2inv*POWER2(si1) + 
       2*Q2inv*POWER2(sih) - 2*s1h*POWER2(sik) + 
       2*si1*POWER2(skh)))/(s1h*si1)
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 2:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
(16*(v12*vi0 - v02*vi1 - v01*vi2)*
     (2*s02*s0k*s2k - 2*Q2inv*s02*si0 - 2*s0k*s2k*si0 + 
       2*Q2inv*si0*si2 - 2*s02*s0k*sik + 4*s0k*si0*sik + 
       2*s2k*si0*sik + 2*s0k*si2*sik + 2*s2k*si2*sik + 
       Q2inv*POWER2(s02) + 2*si2*POWER2(s0k) + 
       2*Q2inv*POWER2(si0) + Q2inv*POWER2(si2) - 
       2*s02*POWER2(sik)))/(s02*si2*v01*v12*vi1)
                     );
                     break;
                  case 2:
                     term_L=(
(16*(v12*vi0 - v02*vi1 - v01*vi2)*
     (2*s01*s0k*s1k - 2*Q2inv*s01*si0 - 2*s0k*s1k*si0 + 
       2*Q2inv*si0*si1 - 2*s01*s0k*sik + 4*s0k*si0*sik + 
       2*s1k*si0*sik + 2*s0k*si1*sik + 2*s1k*si1*sik + 
       Q2inv*POWER2(s01) + 2*si1*POWER2(s0k) + 
       2*Q2inv*POWER2(si0) + Q2inv*POWER2(si1) - 
       2*s01*POWER2(sik)))/(s01*si1*v02*v12*vi2)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
(128*(Q2inv*(2*si0*(-s0h + sih) + POWER2(s0h) + 
            2*POWER2(si0) + POWER2(sih)) + 
         2*(2*s0k*si0*sik + s0k*sih*sik - s0k*si0*skh + 
            si0*sik*skh + sih*sik*skh + sih*POWER2(s0k) - 
            s0h*(s0k*sik - s0k*skh + POWER2(sik)))) - 
      256*lambda*(Q2inv*(2*si0*(-s0h + sih) + POWER2(s0h) + 
            2*POWER2(si0) + POWER2(sih)) + 
         2*(2*s0k*si0*sik + s0k*sih*sik - s0k*si0*skh + 
            si0*sik*skh + sih*sik*skh + sih*POWER2(s0k) - 
            s0h*(s0k*sik - s0k*skh + POWER2(sik)))) + 
      128*POWER2(lambda)*
       (Q2inv*(4*si0*(-s0h + sih) + 4*POWER2(si0) + 
            3*(POWER2(s0h) + POWER2(sih))) + 
         2*(4*s0k*si0*sik + 2*s0k*sih*sik - 2*s0k*si0*skh - 
            s0k*sih*skh + 2*si0*sik*skh + 3*sih*sik*skh + 
            2*sih*POWER2(s0k) - 
            s0h*(2*s0k*sik - 3*s0k*skh + sik*skh + 
               2*POWER2(sik)) + si0*POWER2(skh))) - 
      256*(2*skh*(s0h*s0k - s0h*sik + sih*(-s0k + sik) + 
            si0*skh) + Q2inv*(POWER2(s0h) + POWER2(sih)))*
       POWER3(lambda) + 128*
       (2*skh*(s0h*s0k - s0h*sik + sih*(-s0k + sik) + 
            si0*skh) + Q2inv*(POWER2(s0h) + POWER2(sih)))*
       POWER4(lambda))/((-1 + lambda)*lambda*s0h*sih) + 
   512*(-1 + lambda)*lambda*
    (-2*dt0*dti*Q2inv*s0h*si0*sih - 
      4*dt0*dti*s0h*s0k*sih*sik + 
      2*dt0*dti*s0h*s0k*sih*skh - 
      2*dt0*dti*s0h*sih*sik*skh + 
      2*dt0*dti*Q2inv*sih*POWER2(s0h) + 
      2*dti*dtk*s0k*sih*POWER2(s0h) + 
      2*dti*dtk*sih*sik*POWER2(s0h) - 
      2*dti*dtk*sih*skh*POWER2(s0h) + 
      Q2inv*si0*POWER2(dti)*POWER2(s0h) + 
      Q2inv*sih*POWER2(dti)*POWER2(s0h) + 
      2*s0k*sik*POWER2(dti)*POWER2(s0h) - 
      2*s0k*skh*POWER2(dti)*POWER2(s0h) - 
      2*dt0*dti*Q2inv*s0h*POWER2(sih) - 
      2*dt0*dtk*s0h*s0k*POWER2(sih) - 
      2*dt0*dtk*s0h*sik*POWER2(sih) - 
      2*dt0*dtk*s0h*skh*POWER2(sih) - 
      Q2inv*s0h*POWER2(dt0)*POWER2(sih) + 
      Q2inv*si0*POWER2(dt0)*POWER2(sih) + 
      2*s0k*sik*POWER2(dt0)*POWER2(sih) + 
      2*sik*skh*POWER2(dt0)*POWER2(sih) + 
      2*POWER2(dtk)*POWER2(s0h)*POWER2(sih) + 
      2*dt0*dti*s0h*sih*POWER2(skh) - 
      Q2inv*POWER2(dti)*POWER3(s0h) + 
      Q2inv*POWER2(dt0)*POWER3(sih))*POWERM2(s0h)*
    POWERM2(sih)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
(128*(Q2inv*(2*si0*(-s0h + sih) + POWER2(s0h) + 
            2*POWER2(si0) + POWER2(sih)) + 
         2*(2*s0k*si0*sik + s0k*sih*sik - s0k*si0*skh + 
            si0*sik*skh + sih*sik*skh + sih*POWER2(s0k) - 
            s0h*(s0k*sik - s0k*skh + POWER2(sik)))) - 
      256*lambda*(Q2inv*(2*si0*(-s0h + sih) + POWER2(s0h) + 
            2*POWER2(si0) + POWER2(sih)) + 
         2*(2*s0k*si0*sik + s0k*sih*sik - s0k*si0*skh + 
            si0*sik*skh + sih*sik*skh + sih*POWER2(s0k) - 
            s0h*(s0k*sik - s0k*skh + POWER2(sik)))) + 
      128*POWER2(lambda)*
       (Q2inv*(4*si0*(-s0h + sih) + 4*POWER2(si0) + 
            3*(POWER2(s0h) + POWER2(sih))) + 
         2*(4*s0k*si0*sik + 2*s0k*sih*sik - 2*s0k*si0*skh - 
            s0k*sih*skh + 2*si0*sik*skh + 3*sih*sik*skh + 
            2*sih*POWER2(s0k) - 
            s0h*(2*s0k*sik - 3*s0k*skh + sik*skh + 
               2*POWER2(sik)) + si0*POWER2(skh))) - 
      256*(2*skh*(s0h*s0k - s0h*sik + sih*(-s0k + sik) + 
            si0*skh) + Q2inv*(POWER2(s0h) + POWER2(sih)))*
       POWER3(lambda) + 128*
       (2*skh*(s0h*s0k - s0h*sik + sih*(-s0k + sik) + 
            si0*skh) + Q2inv*(POWER2(s0h) + POWER2(sih)))*
       POWER4(lambda))/((-1 + lambda)*lambda*s0h*sih) + 
   512*(-1 + lambda)*lambda*
    (-2*dt0*dti*Q2inv*s0h*si0*sih - 
      4*dt0*dti*s0h*s0k*sih*sik + 
      2*dt0*dti*s0h*s0k*sih*skh - 
      2*dt0*dti*s0h*sih*sik*skh + 
      2*dt0*dti*Q2inv*sih*POWER2(s0h) + 
      2*dti*dtk*s0k*sih*POWER2(s0h) + 
      2*dti*dtk*sih*sik*POWER2(s0h) - 
      2*dti*dtk*sih*skh*POWER2(s0h) + 
      Q2inv*si0*POWER2(dti)*POWER2(s0h) + 
      Q2inv*sih*POWER2(dti)*POWER2(s0h) + 
      2*s0k*sik*POWER2(dti)*POWER2(s0h) - 
      2*s0k*skh*POWER2(dti)*POWER2(s0h) - 
      2*dt0*dti*Q2inv*s0h*POWER2(sih) - 
      2*dt0*dtk*s0h*s0k*POWER2(sih) - 
      2*dt0*dtk*s0h*sik*POWER2(sih) - 
      2*dt0*dtk*s0h*skh*POWER2(sih) - 
      Q2inv*s0h*POWER2(dt0)*POWER2(sih) + 
      Q2inv*si0*POWER2(dt0)*POWER2(sih) + 
      2*s0k*sik*POWER2(dt0)*POWER2(sih) + 
      2*sik*skh*POWER2(dt0)*POWER2(sih) + 
      2*POWER2(dtk)*POWER2(s0h)*POWER2(sih) + 
      2*dt0*dti*s0h*sih*POWER2(skh) - 
      Q2inv*POWER2(dti)*POWER3(s0h) + 
      Q2inv*POWER2(dt0)*POWER3(sih))*POWERM2(s0h)*
    POWERM2(sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
(-128*(-2*Q2inv*s0h*si0 + 2*Q2inv*si0*sih - 2*s0h*s0k*sik + 
       4*s0k*si0*sik + 2*s0k*sih*sik + 2*s0h*s0k*skh - 
       2*s0k*si0*skh + 2*si0*sik*skh + 2*sih*sik*skh + 
       Q2inv*POWER2(s0h) + 2*sih*POWER2(s0k) + 
       2*Q2inv*POWER2(si0) + Q2inv*POWER2(sih) - 
       2*s0h*POWER2(sik)))/(s0h*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
(-128*(-2*Q2inv*s0h*si0 + 2*Q2inv*si0*sih - 2*s0h*s0k*sik + 
       4*s0k*si0*sik + 2*s0k*sih*sik + 2*s0h*s0k*skh - 
       2*s0k*si0*skh + 2*si0*sik*skh + 2*sih*sik*skh + 
       Q2inv*POWER2(s0h) + 2*sih*POWER2(s0k) + 
       2*Q2inv*POWER2(si0) + Q2inv*POWER2(sih) - 
       2*s0h*POWER2(sik)))/(s0h*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 3:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
0
                     );
                     break;
                  case 2:
                     term_L=(
(16*v01*(4*s01*s0k*s1k - 2*Q2inv*s01*si0 - 2*s0k*s1k*si0 - 
       2*Q2inv*s01*si1 - 2*s0k*s1k*si1 - 2*s01*s0k*sik - 
       2*s01*s1k*sik + 2*s0k*si0*sik + 2*s1k*si1*sik + 
       2*Q2inv*POWER2(s01) + 2*si1*POWER2(s0k) + 
       2*si0*POWER2(s1k) + Q2inv*POWER2(si0) + 
       Q2inv*POWER2(si1)))/(si0*si1*v02*v12)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
(32*(1 - 2*eta + 2*POWER2(eta))*
     (-2*Q2inv*s12*s1h + 2*Q2inv*s1h*s2h + 2*s12*s1k*s2k - 
       2*s1h*s1k*s2k - 2*s12*s1k*skh + 4*s1h*s1k*skh + 
       2*s1k*s2h*skh + 2*s1h*s2k*skh + 2*s2h*s2k*skh + 
       Q2inv*POWER2(s12) + 2*Q2inv*POWER2(s1h) + 
       2*s2h*POWER2(s1k) + Q2inv*POWER2(s2h) - 
       2*s12*POWER2(skh)))/(eta*s12*s2h)
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
(-32*(1 + POWER2(lambda))*
     (-2*Q2inv*s1h*si1 - 2*Q2inv*s1h*sih - 2*s1h*s1k*sik + 
       2*s1k*si1*sik + 4*s1h*s1k*skh - 2*s1k*si1*skh - 
       2*s1k*sih*skh - 2*s1h*sik*skh + 2*sih*sik*skh + 
       2*Q2inv*POWER2(s1h) + 2*sih*POWER2(s1k) + 
       Q2inv*POWER2(si1) + Q2inv*POWER2(sih) + 
       2*si1*POWER2(skh)))/((-1 + lambda)*si1*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
(32*(1 - 2*eta + 2*POWER2(eta))*
     (-2*Q2inv*s02*s0h + 2*Q2inv*s0h*s2h + 2*s02*s0k*s2k - 
       2*s0h*s0k*s2k - 2*s02*s0k*skh + 4*s0h*s0k*skh + 
       2*s0k*s2h*skh + 2*s0h*s2k*skh + 2*s2h*s2k*skh + 
       Q2inv*POWER2(s02) + 2*Q2inv*POWER2(s0h) + 
       2*s2h*POWER2(s0k) + Q2inv*POWER2(s2h) - 
       2*s02*POWER2(skh)))/(eta*s02*s2h)
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
(-32*(1 + POWER2(lambda))*
     (-2*Q2inv*s0h*si0 - 2*Q2inv*s0h*sih - 2*s0h*s0k*sik + 
       2*s0k*si0*sik + 4*s0h*s0k*skh - 2*s0k*si0*skh - 
       2*s0k*sih*skh - 2*s0h*sik*skh + 2*sih*sik*skh + 
       2*Q2inv*POWER2(s0h) + 2*sih*POWER2(s0k) + 
       Q2inv*POWER2(si0) + Q2inv*POWER2(sih) + 
       2*si0*POWER2(skh)))/((-1 + lambda)*si0*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
(32*(2 - 2*lambda + POWER2(lambda))*
     (-2*Q2inv*s1h*si1 - 2*Q2inv*s1h*sih - 2*s1h*s1k*sik + 
       2*s1k*si1*sik + 4*s1h*s1k*skh - 2*s1k*si1*skh - 
       2*s1k*sih*skh - 2*s1h*sik*skh + 2*sih*sik*skh + 
       2*Q2inv*POWER2(s1h) + 2*sih*POWER2(s1k) + 
       Q2inv*POWER2(si1) + Q2inv*POWER2(sih) + 
       2*si1*POWER2(skh)))/(lambda*si1*sih)
                           );
                           break;
                        case 1:
                           term_L=(
(32*(2 - 2*lambda + POWER2(lambda))*
     (-2*Q2inv*s0h*si0 - 2*Q2inv*s0h*sih - 2*s0h*s0k*sik + 
       2*s0k*si0*sik + 4*s0h*s0k*skh - 2*s0k*si0*skh - 
       2*s0k*sih*skh - 2*s0h*sik*skh + 2*sih*sik*skh + 
       2*Q2inv*POWER2(s0h) + 2*sih*POWER2(s0k) + 
       Q2inv*POWER2(si0) + Q2inv*POWER2(sih) + 
       2*si0*POWER2(skh)))/(lambda*si0*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
(64*(-2*Q2inv*s1h*si1 - 2*Q2inv*s1h*sih - 2*s1h*s1k*sik + 
       2*s1k*si1*sik + 4*s1h*s1k*skh - 2*s1k*si1*skh - 
       2*s1k*sih*skh - 2*s1h*sik*skh + 2*sih*sik*skh + 
       2*Q2inv*POWER2(s1h) + 2*sih*POWER2(s1k) + 
       Q2inv*POWER2(si1) + Q2inv*POWER2(sih) + 
       2*si1*POWER2(skh)))/(si1*sih)
                           );
                           break;
                        case 1:
                           term_L=(
(64*(-2*Q2inv*s0h*si0 - 2*Q2inv*s0h*sih - 2*s0h*s0k*sik + 
       2*s0k*si0*sik + 4*s0h*s0k*skh - 2*s0k*si0*skh - 
       2*s0k*sih*skh - 2*s0h*sik*skh + 2*sih*sik*skh + 
       2*Q2inv*POWER2(s0h) + 2*sih*POWER2(s0k) + 
       Q2inv*POWER2(si0) + Q2inv*POWER2(sih) + 
       2*si0*POWER2(skh)))/(si0*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 4:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
0
                     );
                     break;
                  case 2:
                     term_L=(
(16*(-(v12*vi0) - v02*vi1 + v01*vi2)*
     (4*s01*s0k*s1k - 2*Q2inv*s01*si0 - 2*s0k*s1k*si0 - 
       2*Q2inv*s01*si1 - 2*s0k*s1k*si1 - 2*s01*s0k*sik - 
       2*s01*s1k*sik + 2*s0k*si0*sik + 2*s1k*si1*sik + 
       2*Q2inv*POWER2(s01) + 2*si1*POWER2(s0k) + 
       2*si0*POWER2(s1k) + Q2inv*POWER2(si0) + 
       Q2inv*POWER2(si1)))/(si0*si1*v02*v12*vi2)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
((128*(2*skh*(s0k*s1h + s1h*s1k + s0h*(s0k + s1k) - 
              s01*skh) + Q2inv*(POWER2(s0h) + POWER2(s1h)))\
         - 256*eta*(2*skh*
            (s0k*s1h + s1h*s1k + s0h*(s0k + s1k) - s01*skh)\
            + Q2inv*(POWER2(s0h) + POWER2(s1h))) + 
        128*POWER2(eta)*(Q2inv*
            (-4*s01*(s0h + s1h) + 4*POWER2(s01) + 
              3*(POWER2(s0h) + POWER2(s1h))) + 
           2*(4*s01*s0k*s1k - 2*s0h*s0k*s1k - 
              2*s01*s0k*skh + 3*s0h*s0k*skh - 
              2*s01*s1k*skh + s0h*s1k*skh + 
              s1h*(-2*s0k*s1k + s0k*skh + 3*s1k*skh + 
                 2*POWER2(s0k)) + 2*s0h*POWER2(s1k) - 
              s01*POWER2(skh))) - 
        256*(Q2inv*(-2*s01*(s0h + s1h) + 2*POWER2(s01) + 
              POWER2(s0h) + POWER2(s1h)) + 
           2*(2*s01*s0k*s1k - s0h*s0k*s1k - s01*s0k*skh + 
              s0h*s0k*skh - s01*s1k*skh + 
              s1h*(-(s0k*s1k) + s1k*skh + POWER2(s0k)) + 
              s0h*POWER2(s1k)))*POWER3(eta) + 
        128*(Q2inv*(-2*s01*(s0h + s1h) + 2*POWER2(s01) + 
              POWER2(s0h) + POWER2(s1h)) + 
           2*(2*s01*s0k*s1k - s0h*s0k*s1k - s01*s0k*skh + 
              s0h*s0k*skh - s01*s1k*skh + 
              s1h*(-(s0k*s1k) + s1k*skh + POWER2(s0k)) + 
              s0h*POWER2(s1k)))*POWER4(eta))*POWERM2(eta))/
    ((-1 + eta)*s0h*s1h) - 
   512*(-1 + eta)*(2*dt0*dt1*Q2inv*s01*s0h*s1h + 
      4*dt0*dt1*s0h*s0k*s1h*s1k - 
      2*dt0*dt1*s0h*s0k*s1h*skh - 
      2*dt0*dt1*s0h*s1h*s1k*skh - 
      2*dt0*dt1*Q2inv*s1h*POWER2(s0h) - 
      2*dt1*dtk*s0k*s1h*POWER2(s0h) + 
      2*dt1*dtk*s1h*s1k*POWER2(s0h) - 
      2*dt1*dtk*s1h*skh*POWER2(s0h) - 
      Q2inv*s01*POWER2(dt1)*POWER2(s0h) + 
      Q2inv*s1h*POWER2(dt1)*POWER2(s0h) - 
      2*s0k*s1k*POWER2(dt1)*POWER2(s0h) + 
      2*s0k*skh*POWER2(dt1)*POWER2(s0h) - 
      2*dt0*dt1*Q2inv*s0h*POWER2(s1h) + 
      2*dt0*dtk*s0h*s0k*POWER2(s1h) - 
      2*dt0*dtk*s0h*s1k*POWER2(s1h) - 
      2*dt0*dtk*s0h*skh*POWER2(s1h) - 
      Q2inv*s01*POWER2(dt0)*POWER2(s1h) + 
      Q2inv*s0h*POWER2(dt0)*POWER2(s1h) - 
      2*s0k*s1k*POWER2(dt0)*POWER2(s1h) + 
      2*s1k*skh*POWER2(dt0)*POWER2(s1h) + 
      2*POWER2(dtk)*POWER2(s0h)*POWER2(s1h) + 
      2*dt0*dt1*s0h*s1h*POWER2(skh) + 
      Q2inv*POWER2(dt1)*POWER3(s0h) + 
      Q2inv*POWER2(dt0)*POWER3(s1h))*POWERM2(eta)*
    POWERM2(s0h)*POWERM2(s1h)
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
(128*(2*Q2inv*s01*s0h + 2*Q2inv*s01*s1h - 4*s01*s0k*s1k + 
       2*s0h*s0k*s1k + 2*s0k*s1h*s1k + 2*s01*s0k*skh - 
       2*s0h*s0k*skh + 2*s01*s1k*skh - 2*s1h*s1k*skh - 
       2*Q2inv*POWER2(s01) - Q2inv*POWER2(s0h) - 
       2*s1h*POWER2(s0k) - Q2inv*POWER2(s1h) - 
       2*s0h*POWER2(s1k)))/(s0h*s1h)
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 5:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
0
                     );
                     break;
                  case 2:
                     term_L=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
(-64*lambda*(2*skh*(s0h*s0k - s0h*sik + sih*(-s0k + sik) + 
            si0*skh) + Q2inv*(POWER2(s0h) + POWER2(sih))) + 
      64*POWER2(lambda)*(2*skh*
          (s0h*s0k - s0h*sik + sih*(-s0k + sik) + si0*skh)\
          + Q2inv*(POWER2(s0h) + POWER2(sih))) + 
      32*(Q2inv*(2*si0*(-s0h + sih) + POWER2(s0h) + 
            2*POWER2(si0) + POWER2(sih)) + 
         2*(2*s0k*si0*sik + s0k*sih*sik - s0k*si0*skh + 
            si0*sik*skh + sih*sik*skh + sih*POWER2(s0k) - 
            s0h*(s0k*sik - s0k*skh + POWER2(sik)))))/
    (s0h*sih) + 256*(-1 + lambda)*lambda*
    (-2*dt0*dti*Q2inv*s0h*si0*sih - 
      4*dt0*dti*s0h*s0k*sih*sik + 
      2*dt0*dti*s0h*s0k*sih*skh - 
      2*dt0*dti*s0h*sih*sik*skh + 
      2*dt0*dti*Q2inv*sih*POWER2(s0h) + 
      2*dti*dtk*s0k*sih*POWER2(s0h) + 
      2*dti*dtk*sih*sik*POWER2(s0h) - 
      2*dti*dtk*sih*skh*POWER2(s0h) + 
      Q2inv*si0*POWER2(dti)*POWER2(s0h) + 
      Q2inv*sih*POWER2(dti)*POWER2(s0h) + 
      2*s0k*sik*POWER2(dti)*POWER2(s0h) - 
      2*s0k*skh*POWER2(dti)*POWER2(s0h) - 
      2*dt0*dti*Q2inv*s0h*POWER2(sih) - 
      2*dt0*dtk*s0h*s0k*POWER2(sih) - 
      2*dt0*dtk*s0h*sik*POWER2(sih) - 
      2*dt0*dtk*s0h*skh*POWER2(sih) - 
      Q2inv*s0h*POWER2(dt0)*POWER2(sih) + 
      Q2inv*si0*POWER2(dt0)*POWER2(sih) + 
      2*s0k*sik*POWER2(dt0)*POWER2(sih) + 
      2*sik*skh*POWER2(dt0)*POWER2(sih) + 
      2*POWER2(dtk)*POWER2(s0h)*POWER2(sih) + 
      2*dt0*dti*s0h*sih*POWER2(skh) - 
      Q2inv*POWER2(dti)*POWER3(s0h) + 
      Q2inv*POWER2(dt0)*POWER3(sih))*POWERM2(s0h)*
    POWERM2(sih)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
(-64*lambda*(2*skh*(s0h*s0k - s0h*sik + sih*(-s0k + sik) + 
            si0*skh) + Q2inv*(POWER2(s0h) + POWER2(sih))) + 
      64*POWER2(lambda)*(2*skh*
          (s0h*s0k - s0h*sik + sih*(-s0k + sik) + si0*skh)\
          + Q2inv*(POWER2(s0h) + POWER2(sih))) + 
      32*(Q2inv*(2*si0*(-s0h + sih) + POWER2(s0h) + 
            2*POWER2(si0) + POWER2(sih)) + 
         2*(2*s0k*si0*sik + s0k*sih*sik - s0k*si0*skh + 
            si0*sik*skh + sih*sik*skh + sih*POWER2(s0k) - 
            s0h*(s0k*sik - s0k*skh + POWER2(sik)))))/
    (s0h*sih) + 256*(-1 + lambda)*lambda*
    (-2*dt0*dti*Q2inv*s0h*si0*sih - 
      4*dt0*dti*s0h*s0k*sih*sik + 
      2*dt0*dti*s0h*s0k*sih*skh - 
      2*dt0*dti*s0h*sih*sik*skh + 
      2*dt0*dti*Q2inv*sih*POWER2(s0h) + 
      2*dti*dtk*s0k*sih*POWER2(s0h) + 
      2*dti*dtk*sih*sik*POWER2(s0h) - 
      2*dti*dtk*sih*skh*POWER2(s0h) + 
      Q2inv*si0*POWER2(dti)*POWER2(s0h) + 
      Q2inv*sih*POWER2(dti)*POWER2(s0h) + 
      2*s0k*sik*POWER2(dti)*POWER2(s0h) - 
      2*s0k*skh*POWER2(dti)*POWER2(s0h) - 
      2*dt0*dti*Q2inv*s0h*POWER2(sih) - 
      2*dt0*dtk*s0h*s0k*POWER2(sih) - 
      2*dt0*dtk*s0h*sik*POWER2(sih) - 
      2*dt0*dtk*s0h*skh*POWER2(sih) - 
      Q2inv*s0h*POWER2(dt0)*POWER2(sih) + 
      Q2inv*si0*POWER2(dt0)*POWER2(sih) + 
      2*s0k*sik*POWER2(dt0)*POWER2(sih) + 
      2*sik*skh*POWER2(dt0)*POWER2(sih) + 
      2*POWER2(dtk)*POWER2(s0h)*POWER2(sih) + 
      2*dt0*dti*s0h*sih*POWER2(skh) - 
      Q2inv*POWER2(dti)*POWER3(s0h) + 
      Q2inv*POWER2(dt0)*POWER3(sih))*POWERM2(s0h)*
    POWERM2(sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 6:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
0
                     );
                     break;
                  case 2:
                     term_L=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
((64*(2*skh*(s1k*s2h + s2h*s2k + s1h*(s1k + s2k) - 
              s12*skh) + Q2inv*(POWER2(s1h) + POWER2(s2h)))\
         - 64*eta*(2*skh*
            (s1k*s2h + s2h*s2k + s1h*(s1k + s2k) - s12*skh)\
            + Q2inv*(POWER2(s1h) + POWER2(s2h))) + 
        32*POWER2(eta)*(Q2inv*
            (-2*s12*(s1h + s2h) + 2*POWER2(s12) + 
              POWER2(s1h) + POWER2(s2h)) + 
           2*(2*s12*s1k*s2k - s1h*s1k*s2k - s12*s1k*skh + 
              s1h*s1k*skh - s12*s2k*skh + 
              s2h*(-(s1k*s2k) + s2k*skh + POWER2(s1k)) + 
              s1h*POWER2(s2k))))*POWERM2(eta))/(s1h*s2h) + 
   256*(-1 + eta)*(2*dt1*dt2*Q2inv*s12*s1h*s2h + 
      4*dt1*dt2*s1h*s1k*s2h*s2k - 
      2*dt1*dt2*s1h*s1k*s2h*skh - 
      2*dt1*dt2*s1h*s2h*s2k*skh - 
      2*dt1*dt2*Q2inv*s2h*POWER2(s1h) - 
      2*dt2*dtk*s1k*s2h*POWER2(s1h) + 
      2*dt2*dtk*s2h*s2k*POWER2(s1h) - 
      2*dt2*dtk*s2h*skh*POWER2(s1h) - 
      Q2inv*s12*POWER2(dt2)*POWER2(s1h) + 
      Q2inv*s2h*POWER2(dt2)*POWER2(s1h) - 
      2*s1k*s2k*POWER2(dt2)*POWER2(s1h) + 
      2*s1k*skh*POWER2(dt2)*POWER2(s1h) - 
      2*dt1*dt2*Q2inv*s1h*POWER2(s2h) + 
      2*dt1*dtk*s1h*s1k*POWER2(s2h) - 
      2*dt1*dtk*s1h*s2k*POWER2(s2h) - 
      2*dt1*dtk*s1h*skh*POWER2(s2h) - 
      Q2inv*s12*POWER2(dt1)*POWER2(s2h) + 
      Q2inv*s1h*POWER2(dt1)*POWER2(s2h) - 
      2*s1k*s2k*POWER2(dt1)*POWER2(s2h) + 
      2*s2k*skh*POWER2(dt1)*POWER2(s2h) + 
      2*POWER2(dtk)*POWER2(s1h)*POWER2(s2h) + 
      2*dt1*dt2*s1h*s2h*POWER2(skh) + 
      Q2inv*POWER2(dt2)*POWER3(s1h) + 
      Q2inv*POWER2(dt1)*POWER3(s2h))*POWERM2(eta)*
    POWERM2(s1h)*POWERM2(s2h)
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 7:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
0
                     );
                     break;
                  case 2:
                     term_L=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 8:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
0
                     );
                     break;
                  case 2:
                     term_L=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
((64*(2*skh*(s1k*s2h + s2h*s2k + s1h*(s1k + s2k) - 
              s12*skh) + Q2inv*(POWER2(s1h) + POWER2(s2h)))\
         - 64*eta*(2*skh*
            (s1k*s2h + s2h*s2k + s1h*(s1k + s2k) - s12*skh)\
            + Q2inv*(POWER2(s1h) + POWER2(s2h))) + 
        32*POWER2(eta)*(Q2inv*
            (-2*s12*(s1h + s2h) + 2*POWER2(s12) + 
              POWER2(s1h) + POWER2(s2h)) + 
           2*(2*s12*s1k*s2k - s1h*s1k*s2k - s12*s1k*skh + 
              s1h*s1k*skh - s12*s2k*skh + 
              s2h*(-(s1k*s2k) + s2k*skh + POWER2(s1k)) + 
              s1h*POWER2(s2k))))*POWERM2(eta))/(s1h*s2h) + 
   256*(-1 + eta)*(2*dt1*dt2*Q2inv*s12*s1h*s2h + 
      4*dt1*dt2*s1h*s1k*s2h*s2k - 
      2*dt1*dt2*s1h*s1k*s2h*skh - 
      2*dt1*dt2*s1h*s2h*s2k*skh - 
      2*dt1*dt2*Q2inv*s2h*POWER2(s1h) - 
      2*dt2*dtk*s1k*s2h*POWER2(s1h) + 
      2*dt2*dtk*s2h*s2k*POWER2(s1h) - 
      2*dt2*dtk*s2h*skh*POWER2(s1h) - 
      Q2inv*s12*POWER2(dt2)*POWER2(s1h) + 
      Q2inv*s2h*POWER2(dt2)*POWER2(s1h) - 
      2*s1k*s2k*POWER2(dt2)*POWER2(s1h) + 
      2*s1k*skh*POWER2(dt2)*POWER2(s1h) - 
      2*dt1*dt2*Q2inv*s1h*POWER2(s2h) + 
      2*dt1*dtk*s1h*s1k*POWER2(s2h) - 
      2*dt1*dtk*s1h*s2k*POWER2(s2h) - 
      2*dt1*dtk*s1h*skh*POWER2(s2h) - 
      Q2inv*s12*POWER2(dt1)*POWER2(s2h) + 
      Q2inv*s1h*POWER2(dt1)*POWER2(s2h) - 
      2*s1k*s2k*POWER2(dt1)*POWER2(s2h) + 
      2*s2k*skh*POWER2(dt1)*POWER2(s2h) + 
      2*POWER2(dtk)*POWER2(s1h)*POWER2(s2h) + 
      2*dt1*dt2*s1h*s2h*POWER2(skh) + 
      Q2inv*POWER2(dt2)*POWER3(s1h) + 
      Q2inv*POWER2(dt1)*POWER3(s2h))*POWERM2(eta)*
    POWERM2(s1h)*POWERM2(s2h)
                           );
                           break;
                        case 1:
                           term_L=(
(-64*lambda*(2*skh*(s2h*s2k - s2h*sik + sih*(-s2k + sik) + 
            si2*skh) + Q2inv*(POWER2(s2h) + POWER2(sih))) + 
      64*POWER2(lambda)*(2*skh*
          (s2h*s2k - s2h*sik + sih*(-s2k + sik) + si2*skh)\
          + Q2inv*(POWER2(s2h) + POWER2(sih))) + 
      32*(Q2inv*(2*si2*(-s2h + sih) + POWER2(s2h) + 
            2*POWER2(si2) + POWER2(sih)) + 
         2*(2*s2k*si2*sik + s2k*sih*sik - s2k*si2*skh + 
            si2*sik*skh + sih*sik*skh + sih*POWER2(s2k) - 
            s2h*(s2k*sik - s2k*skh + POWER2(sik)))))/
    (s2h*sih) + 256*(-1 + lambda)*lambda*
    (-2*dt2*dti*Q2inv*s2h*si2*sih - 
      4*dt2*dti*s2h*s2k*sih*sik + 
      2*dt2*dti*s2h*s2k*sih*skh - 
      2*dt2*dti*s2h*sih*sik*skh + 
      2*dt2*dti*Q2inv*sih*POWER2(s2h) + 
      2*dti*dtk*s2k*sih*POWER2(s2h) + 
      2*dti*dtk*sih*sik*POWER2(s2h) - 
      2*dti*dtk*sih*skh*POWER2(s2h) + 
      Q2inv*si2*POWER2(dti)*POWER2(s2h) + 
      Q2inv*sih*POWER2(dti)*POWER2(s2h) + 
      2*s2k*sik*POWER2(dti)*POWER2(s2h) - 
      2*s2k*skh*POWER2(dti)*POWER2(s2h) - 
      2*dt2*dti*Q2inv*s2h*POWER2(sih) - 
      2*dt2*dtk*s2h*s2k*POWER2(sih) - 
      2*dt2*dtk*s2h*sik*POWER2(sih) - 
      2*dt2*dtk*s2h*skh*POWER2(sih) - 
      Q2inv*s2h*POWER2(dt2)*POWER2(sih) + 
      Q2inv*si2*POWER2(dt2)*POWER2(sih) + 
      2*s2k*sik*POWER2(dt2)*POWER2(sih) + 
      2*sik*skh*POWER2(dt2)*POWER2(sih) + 
      2*POWER2(dtk)*POWER2(s2h)*POWER2(sih) + 
      2*dt2*dti*s2h*sih*POWER2(skh) - 
      Q2inv*POWER2(dti)*POWER3(s2h) + 
      Q2inv*POWER2(dt2)*POWER3(sih))*POWERM2(s2h)*
    POWERM2(sih)
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
(-64*lambda*(2*skh*(s2h*s2k - s2h*sik + sih*(-s2k + sik) + 
            si2*skh) + Q2inv*(POWER2(s2h) + POWER2(sih))) + 
      64*POWER2(lambda)*(2*skh*
          (s2h*s2k - s2h*sik + sih*(-s2k + sik) + si2*skh)\
          + Q2inv*(POWER2(s2h) + POWER2(sih))) + 
      32*(Q2inv*(2*si2*(-s2h + sih) + POWER2(s2h) + 
            2*POWER2(si2) + POWER2(sih)) + 
         2*(2*s2k*si2*sik + s2k*sih*sik - s2k*si2*skh + 
            si2*sik*skh + sih*sik*skh + sih*POWER2(s2k) - 
            s2h*(s2k*sik - s2k*skh + POWER2(sik)))))/
    (s2h*sih) + 256*(-1 + lambda)*lambda*
    (-2*dt2*dti*Q2inv*s2h*si2*sih - 
      4*dt2*dti*s2h*s2k*sih*sik + 
      2*dt2*dti*s2h*s2k*sih*skh - 
      2*dt2*dti*s2h*sih*sik*skh + 
      2*dt2*dti*Q2inv*sih*POWER2(s2h) + 
      2*dti*dtk*s2k*sih*POWER2(s2h) + 
      2*dti*dtk*sih*sik*POWER2(s2h) - 
      2*dti*dtk*sih*skh*POWER2(s2h) + 
      Q2inv*si2*POWER2(dti)*POWER2(s2h) + 
      Q2inv*sih*POWER2(dti)*POWER2(s2h) + 
      2*s2k*sik*POWER2(dti)*POWER2(s2h) - 
      2*s2k*skh*POWER2(dti)*POWER2(s2h) - 
      2*dt2*dti*Q2inv*s2h*POWER2(sih) - 
      2*dt2*dtk*s2h*s2k*POWER2(sih) - 
      2*dt2*dtk*s2h*sik*POWER2(sih) - 
      2*dt2*dtk*s2h*skh*POWER2(sih) - 
      Q2inv*s2h*POWER2(dt2)*POWER2(sih) + 
      Q2inv*si2*POWER2(dt2)*POWER2(sih) + 
      2*s2k*sik*POWER2(dt2)*POWER2(sih) + 
      2*sik*skh*POWER2(dt2)*POWER2(sih) + 
      2*POWER2(dtk)*POWER2(s2h)*POWER2(sih) + 
      2*dt2*dti*s2h*sih*POWER2(skh) - 
      Q2inv*POWER2(dti)*POWER3(s2h) + 
      Q2inv*POWER2(dt2)*POWER3(sih))*POWERM2(s2h)*
    POWERM2(sih)
                           );
                           break;
                        case 2:
                           term_L=(
(-64*lambda*(2*skh*(s0h*s0k - s0h*sik + sih*(-s0k + sik) + 
            si0*skh) + Q2inv*(POWER2(s0h) + POWER2(sih))) + 
      64*POWER2(lambda)*(2*skh*
          (s0h*s0k - s0h*sik + sih*(-s0k + sik) + si0*skh)\
          + Q2inv*(POWER2(s0h) + POWER2(sih))) + 
      32*(Q2inv*(2*si0*(-s0h + sih) + POWER2(s0h) + 
            2*POWER2(si0) + POWER2(sih)) + 
         2*(2*s0k*si0*sik + s0k*sih*sik - s0k*si0*skh + 
            si0*sik*skh + sih*sik*skh + sih*POWER2(s0k) - 
            s0h*(s0k*sik - s0k*skh + POWER2(sik)))))/
    (s0h*sih) + 256*(-1 + lambda)*lambda*
    (-2*dt0*dti*Q2inv*s0h*si0*sih - 
      4*dt0*dti*s0h*s0k*sih*sik + 
      2*dt0*dti*s0h*s0k*sih*skh - 
      2*dt0*dti*s0h*sih*sik*skh + 
      2*dt0*dti*Q2inv*sih*POWER2(s0h) + 
      2*dti*dtk*s0k*sih*POWER2(s0h) + 
      2*dti*dtk*sih*sik*POWER2(s0h) - 
      2*dti*dtk*sih*skh*POWER2(s0h) + 
      Q2inv*si0*POWER2(dti)*POWER2(s0h) + 
      Q2inv*sih*POWER2(dti)*POWER2(s0h) + 
      2*s0k*sik*POWER2(dti)*POWER2(s0h) - 
      2*s0k*skh*POWER2(dti)*POWER2(s0h) - 
      2*dt0*dti*Q2inv*s0h*POWER2(sih) - 
      2*dt0*dtk*s0h*s0k*POWER2(sih) - 
      2*dt0*dtk*s0h*sik*POWER2(sih) - 
      2*dt0*dtk*s0h*skh*POWER2(sih) - 
      Q2inv*s0h*POWER2(dt0)*POWER2(sih) + 
      Q2inv*si0*POWER2(dt0)*POWER2(sih) + 
      2*s0k*sik*POWER2(dt0)*POWER2(sih) + 
      2*sik*skh*POWER2(dt0)*POWER2(sih) + 
      2*POWER2(dtk)*POWER2(s0h)*POWER2(sih) + 
      2*dt0*dti*s0h*sih*POWER2(skh) - 
      Q2inv*POWER2(dti)*POWER3(s0h) + 
      Q2inv*POWER2(dt0)*POWER3(sih))*POWERM2(s0h)*
    POWERM2(sih)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
((64*(2*skh*(s0k*s1h + s1h*s1k + s0h*(s0k + s1k) - 
              s01*skh) + Q2inv*(POWER2(s0h) + POWER2(s1h)))\
         - 64*eta*(2*skh*
            (s0k*s1h + s1h*s1k + s0h*(s0k + s1k) - s01*skh)\
            + Q2inv*(POWER2(s0h) + POWER2(s1h))) + 
        32*POWER2(eta)*(Q2inv*
            (-2*s01*(s0h + s1h) + 2*POWER2(s01) + 
              POWER2(s0h) + POWER2(s1h)) + 
           2*(2*s01*s0k*s1k - s0h*s0k*s1k - s01*s0k*skh + 
              s0h*s0k*skh - s01*s1k*skh + 
              s1h*(-(s0k*s1k) + s1k*skh + POWER2(s0k)) + 
              s0h*POWER2(s1k))))*POWERM2(eta))/(s0h*s1h) + 
   256*(-1 + eta)*(2*dt0*dt1*Q2inv*s01*s0h*s1h + 
      4*dt0*dt1*s0h*s0k*s1h*s1k - 
      2*dt0*dt1*s0h*s0k*s1h*skh - 
      2*dt0*dt1*s0h*s1h*s1k*skh - 
      2*dt0*dt1*Q2inv*s1h*POWER2(s0h) - 
      2*dt1*dtk*s0k*s1h*POWER2(s0h) + 
      2*dt1*dtk*s1h*s1k*POWER2(s0h) - 
      2*dt1*dtk*s1h*skh*POWER2(s0h) - 
      Q2inv*s01*POWER2(dt1)*POWER2(s0h) + 
      Q2inv*s1h*POWER2(dt1)*POWER2(s0h) - 
      2*s0k*s1k*POWER2(dt1)*POWER2(s0h) + 
      2*s0k*skh*POWER2(dt1)*POWER2(s0h) - 
      2*dt0*dt1*Q2inv*s0h*POWER2(s1h) + 
      2*dt0*dtk*s0h*s0k*POWER2(s1h) - 
      2*dt0*dtk*s0h*s1k*POWER2(s1h) - 
      2*dt0*dtk*s0h*skh*POWER2(s1h) - 
      Q2inv*s01*POWER2(dt0)*POWER2(s1h) + 
      Q2inv*s0h*POWER2(dt0)*POWER2(s1h) - 
      2*s0k*s1k*POWER2(dt0)*POWER2(s1h) + 
      2*s1k*skh*POWER2(dt0)*POWER2(s1h) + 
      2*POWER2(dtk)*POWER2(s0h)*POWER2(s1h) + 
      2*dt0*dt1*s0h*s1h*POWER2(skh) + 
      Q2inv*POWER2(dt1)*POWER3(s0h) + 
      Q2inv*POWER2(dt0)*POWER3(s1h))*POWERM2(eta)*
    POWERM2(s0h)*POWERM2(s1h)
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
(-64*lambda*(2*skh*(s0h*s0k - s0h*sik + sih*(-s0k + sik) + 
            si0*skh) + Q2inv*(POWER2(s0h) + POWER2(sih))) + 
      64*POWER2(lambda)*(2*skh*
          (s0h*s0k - s0h*sik + sih*(-s0k + sik) + si0*skh)\
          + Q2inv*(POWER2(s0h) + POWER2(sih))) + 
      32*(Q2inv*(2*si0*(-s0h + sih) + POWER2(s0h) + 
            2*POWER2(si0) + POWER2(sih)) + 
         2*(2*s0k*si0*sik + s0k*sih*sik - s0k*si0*skh + 
            si0*sik*skh + sih*sik*skh + sih*POWER2(s0k) - 
            s0h*(s0k*sik - s0k*skh + POWER2(sik)))))/
    (s0h*sih) + 256*(-1 + lambda)*lambda*
    (-2*dt0*dti*Q2inv*s0h*si0*sih - 
      4*dt0*dti*s0h*s0k*sih*sik + 
      2*dt0*dti*s0h*s0k*sih*skh - 
      2*dt0*dti*s0h*sih*sik*skh + 
      2*dt0*dti*Q2inv*sih*POWER2(s0h) + 
      2*dti*dtk*s0k*sih*POWER2(s0h) + 
      2*dti*dtk*sih*sik*POWER2(s0h) - 
      2*dti*dtk*sih*skh*POWER2(s0h) + 
      Q2inv*si0*POWER2(dti)*POWER2(s0h) + 
      Q2inv*sih*POWER2(dti)*POWER2(s0h) + 
      2*s0k*sik*POWER2(dti)*POWER2(s0h) - 
      2*s0k*skh*POWER2(dti)*POWER2(s0h) - 
      2*dt0*dti*Q2inv*s0h*POWER2(sih) - 
      2*dt0*dtk*s0h*s0k*POWER2(sih) - 
      2*dt0*dtk*s0h*sik*POWER2(sih) - 
      2*dt0*dtk*s0h*skh*POWER2(sih) - 
      Q2inv*s0h*POWER2(dt0)*POWER2(sih) + 
      Q2inv*si0*POWER2(dt0)*POWER2(sih) + 
      2*s0k*sik*POWER2(dt0)*POWER2(sih) + 
      2*sik*skh*POWER2(dt0)*POWER2(sih) + 
      2*POWER2(dtk)*POWER2(s0h)*POWER2(sih) + 
      2*dt0*dti*s0h*sih*POWER2(skh) - 
      Q2inv*POWER2(dti)*POWER3(s0h) + 
      Q2inv*POWER2(dt0)*POWER3(sih))*POWERM2(s0h)*
    POWERM2(sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 9:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
0
                     );
                     break;
                  case 2:
                     term_L=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 2:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 10:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 11:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
(8*(Q2inv*si0 + 2*s0k*sik)*vi0)/(v01*vi1)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
(-16*(Q2inv*sih + 2*sik*skh)*(1 + POWER2(lambda)))/
   (-1 + lambda)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
(-16*(Q2inv*s0h + 2*s0k*skh)*(1 + POWER2(eta)))/
   (-eta + POWER2(eta))
                           );
                           break;
                        case 0:
                           term_L=(
(16*(Q2inv*sih + 2*sik*skh)*
     (2 - 2*lambda + POWER2(lambda)))/lambda
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
32*(Q2inv*s0h + 2*s0k*skh)
                           );
                           break;
                        case 0:
                           term_L=(
32*(Q2inv*sih + 2*sik*skh)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 12:
         switch (sub_type) {
            case SOFT:
               switch (i) {
                  case 0:
                     term_L=(
0
                     );
                     break;
                  case 1:
                     term_L=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
(16*(Q2inv*s1h + 2*s1k*skh)*(1 - 2*eta + 2*POWER2(eta)))/eta
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
(16*(Q2inv*s0h + 2*s0k*skh)*(1 - 2*eta + 2*POWER2(eta)))/eta
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            case SOFT_AND_COLLINEAR:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 1:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           term_L=(
0
                           );
                           break;
                        case 0:
                           term_L=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzSubtraction: 1");
         } // of switch (sub_type)
         break;
   } // of switch (list)
return term_L;
} // of function
// DIS added subtraction terms
// Created on {1997, 12, 16, 13, 39, 7}
#include <math.h>
#include <stdio.h>
#include "mathdefs.h"
#include "global.h"
#include "proc.h"
#include "me.h"
void DISMatrixElement::GetLorentzAddedSubtraction(
   int list, LimitType sub_type, int i, int j
                ) {
   switch (list) {
      case 1:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     // list=1 soft=1 alpha=-1 beta=0
                     softContributes[0]=(TRUE);
                     softCoefficient[0]=(1);
                     softAngle[0]=(vi0);
                     index1[0]=(-1);
                     index2[0]=(0);
                     nCoefficient=(1);
                     bornTerm=(
(16*(Q2inv*(-2*s02*si0 + 2*si0*si2 + POWER2(s02) + 
          2*POWER2(si0) + POWER2(si2)) + 
       2*(-(s0k*s2k*si0) + 2*s0k*si0*sik + s2k*si0*sik + 
          s0k*si2*sik + s2k*si2*sik + si2*POWER2(s0k) + 
          s02*(s0k*s2k - s0k*sik - POWER2(sik)))))/(s02*si2)
                     );
                     break;
                  case 2:
                     // list=1 soft=2 alpha=-1 beta=0
                     softContributes[0]=(TRUE);
                     softCoefficient[0]=(1);
                     softAngle[0]=(vi0);
                     index1[0]=(-1);
                     index2[0]=(0);
                     nCoefficient=(1);
                     bornTerm=(
(16*(Q2inv*(-2*s01*si0 + 2*si0*si1 + POWER2(s01) + 
          2*POWER2(si0) + POWER2(si1)) + 
       2*(-(s0k*s1k*si0) + 2*s0k*si0*sik + s1k*si0*sik + 
          s0k*si1*sik + s1k*si1*sik + si1*POWER2(s0k) + 
          s01*(s0k*s1k - s0k*sik - POWER2(sik)))))/(s01*si1)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=1 i=0 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=1 i=1 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_q);
                           bornTerm=(
(16*(2*(skh*(s0h*s2k + s2h*s2k - s02*skh) + 
          s0k*(-(s0h*s2k) + s02*(s2k - skh) + 2*s0h*skh + 
             s2h*skh) + s2h*POWER2(s0k)) + 
       Q2inv*(-2*s02*s0h + 2*s0h*s2h + POWER2(s02) + 
          2*POWER2(s0h) + POWER2(s2h))))/(s02*s2h)
                           );
                           break;
                        case 0:
                           // list=1 i=1 j=0
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_G_from_q);
                           bornTerm=(
(16*(Q2inv*(POWER2(s2h) + POWER2(si2) + 
          2*(-(s2h*sih) + si2*sih + POWER2(sih))) + 
       2*(-(s2h*sik*skh) + si2*sik*skh + 2*sih*sik*skh + 
          s2k*(si2*sik + sih*(sik - skh) + s2h*skh) - 
          s2h*POWER2(sik) + si2*POWER2(skh))))/(s2h*si2)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=1 i=2 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_q);
                           bornTerm=(
(16*(2*(skh*(s0h*s1k + s1h*s1k - s01*skh) + 
          s0k*(-(s0h*s1k) + s01*(s1k - skh) + 2*s0h*skh + 
             s1h*skh) + s1h*POWER2(s0k)) + 
       Q2inv*(-2*s01*s0h + 2*s0h*s1h + POWER2(s01) + 
          2*POWER2(s0h) + POWER2(s1h))))/(s01*s1h)
                           );
                           break;
                        case 0:
                           // list=1 i=2 j=0
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_G_from_q);
                           bornTerm=(
(16*(Q2inv*(POWER2(s1h) + POWER2(si1) + 
          2*(-(s1h*sih) + si1*sih + POWER2(sih))) + 
       2*(-(s1h*sik*skh) + si1*sik*skh + 2*sih*sik*skh + 
          s1k*(si1*sik + sih*(sik - skh) + s1h*skh) - 
          s1h*POWER2(sik) + si1*POWER2(skh))))/(s1h*si1)
                           );
                           break;
                        case 1:
                           // list=1 i=2 j=1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 2:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     // list=2 soft=1 alpha=-1 beta=0
                     softContributes[0]=(TRUE);
                     softCoefficient[0]=(1);
                     softAngle[0]=(vi0);
                     index1[0]=(-1);
                     index2[0]=(0);
                     // list=2 soft=1 alpha=-1 beta=2
                     softContributes[1]=(TRUE);
                     softCoefficient[1]=(-1);
                     softAngle[1]=(vi2);
                     index1[1]=(-1);
                     index2[1]=(2);
                     // list=2 soft=1 alpha=0 beta=2
                     softContributes[2]=(TRUE);
                     softCoefficient[2]=(-1);
                     softAngle[2]=(v02);
                     index1[2]=(0);
                     index2[2]=(2);
                     nCoefficient=(3);
                     bornTerm=(
(16*(Q2inv*(-2*s02*si0 + 2*si0*si2 + POWER2(s02) + 
          2*POWER2(si0) + POWER2(si2)) + 
       2*(-(s0k*s2k*si0) + 2*s0k*si0*sik + s2k*si0*sik + 
          s0k*si2*sik + s2k*si2*sik + si2*POWER2(s0k) + 
          s02*(s0k*s2k - s0k*sik - POWER2(sik)))))/(s02*si2)
                     );
                     break;
                  case 2:
                     // list=2 soft=2 alpha=-1 beta=0
                     softContributes[0]=(TRUE);
                     softCoefficient[0]=(1);
                     softAngle[0]=(vi0);
                     index1[0]=(-1);
                     index2[0]=(0);
                     // list=2 soft=2 alpha=-1 beta=1
                     softContributes[1]=(TRUE);
                     softCoefficient[1]=(-1);
                     softAngle[1]=(vi1);
                     index1[1]=(-1);
                     index2[1]=(1);
                     // list=2 soft=2 alpha=0 beta=1
                     softContributes[2]=(TRUE);
                     softCoefficient[2]=(-1);
                     softAngle[2]=(v01);
                     index1[2]=(0);
                     index2[2]=(1);
                     nCoefficient=(3);
                     bornTerm=(
(16*(Q2inv*(-2*s01*si0 + 2*si0*si1 + POWER2(s01) + 
          2*POWER2(si0) + POWER2(si1)) + 
       2*(-(s0k*s1k*si0) + 2*s0k*si0*sik + s1k*si0*sik + 
          s0k*si1*sik + s1k*si1*sik + si1*POWER2(s0k) + 
          s01*(s0k*s1k - s0k*sik - POWER2(sik)))))/(s01*si1)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=2 i=0 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=2 i=1 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=2 i=1 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=2 i=2 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=2 i=2 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 1:
                           // list=2 i=2 j=1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(-8);
                           splittingFunction=(Q_G_from_G);
                           bornTerm=(
(16*(2*(sik*(-(s0h*sik) + (si0 + sih)*skh) + 
          s0k*(2*si0*sik + sih*sik - si0*skh + 
             s0h*(-sik + skh)) + sih*POWER2(s0k)) + 
       Q2inv*(-2*s0h*si0 + 2*si0*sih + POWER2(s0h) + 
          2*POWER2(si0) + POWER2(sih))))/(s0h*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 3:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 2:
                     // list=3 soft=2 alpha=0 beta=1
                     softContributes[0]=(TRUE);
                     softCoefficient[0]=(1);
                     softAngle[0]=(v01);
                     index1[0]=(0);
                     index2[0]=(1);
                     nCoefficient=(1);
                     bornTerm=(
(16*(2*(-(s0k*s1k*si0) - s0k*s1k*si1 + s0k*si0*sik + 
          s1k*si1*sik + s01*
           (2*s0k*s1k - s0k*sik - s1k*sik) + 
          si1*POWER2(s0k) + si0*POWER2(s1k)) + 
       Q2inv*(-2*s01*(si0 + si1) + 2*POWER2(s01) + 
          POWER2(si0) + POWER2(si1))))/(si0*si1)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=3 i=0 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_G);
                           bornTerm=(
(16*(2*(skh*(s1h*s2k + s2h*s2k - s12*skh) + 
          s1k*(-(s1h*s2k) + s12*(s2k - skh) + 2*s1h*skh + 
             s2h*skh) + s2h*POWER2(s1k)) + 
       Q2inv*(-2*s12*s1h + 2*s1h*s2h + POWER2(s12) + 
          2*POWER2(s1h) + POWER2(s2h))))/(s12*s2h)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=3 i=1 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_G);
                           bornTerm=(
(16*(2*(skh*(s0h*s2k + s2h*s2k - s02*skh) + 
          s0k*(-(s0h*s2k) + s02*(s2k - skh) + 2*s0h*skh + 
             s2h*skh) + s2h*POWER2(s0k)) + 
       Q2inv*(-2*s02*s0h + 2*s0h*s2h + POWER2(s02) + 
          2*POWER2(s0h) + POWER2(s2h))))/(s02*s2h)
                           );
                           break;
                        case 0:
                           // list=3 i=1 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=3 i=2 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=3 i=2 j=0
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_G_from_q);
                           bornTerm=(
(16*(-2*(-(skh*(-(s1h*sik) + sih*sik + si1*skh)) + 
          s1k*(s1h*(sik - 2*skh) + sih*skh + 
             si1*(-sik + skh)) - sih*POWER2(s1k)) + 
       Q2inv*(-2*s1h*(si1 + sih) + 2*POWER2(s1h) + 
          POWER2(si1) + POWER2(sih))))/(si1*sih)
                           );
                           break;
                        case 1:
                           // list=3 i=2 j=1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_G_from_q);
                           bornTerm=(
(16*(-2*(-(skh*(-(s0h*sik) + sih*sik + si0*skh)) + 
          s0k*(s0h*(sik - 2*skh) + sih*skh + 
             si0*(-sik + skh)) - sih*POWER2(s0k)) + 
       Q2inv*(-2*s0h*(si0 + sih) + 2*POWER2(s0h) + 
          POWER2(si0) + POWER2(sih))))/(si0*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 4:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 2:
                     // list=4 soft=2 alpha=-1 beta=0
                     softContributes[0]=(TRUE);
                     softCoefficient[0]=(-1);
                     softAngle[0]=(vi0);
                     index1[0]=(-1);
                     index2[0]=(0);
                     // list=4 soft=2 alpha=-1 beta=1
                     softContributes[1]=(TRUE);
                     softCoefficient[1]=(-1);
                     softAngle[1]=(vi1);
                     index1[1]=(-1);
                     index2[1]=(1);
                     // list=4 soft=2 alpha=0 beta=1
                     softContributes[2]=(TRUE);
                     softCoefficient[2]=(1);
                     softAngle[2]=(v01);
                     index1[2]=(0);
                     index2[2]=(1);
                     nCoefficient=(3);
                     bornTerm=(
(16*(2*(-(s0k*s1k*si0) - s0k*s1k*si1 + s0k*si0*sik + 
          s1k*si1*sik + s01*
           (2*s0k*s1k - s0k*sik - s1k*sik) + 
          si1*POWER2(s0k) + si0*POWER2(s1k)) + 
       Q2inv*(-2*s01*(si0 + si1) + 2*POWER2(s01) + 
          POWER2(si0) + POWER2(si1))))/(si0*si1)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=4 i=0 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=4 i=1 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=4 i=1 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=4 i=2 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(-8);
                           splittingFunction=(Q_G_from_G);
                           bornTerm=(
(16*(2*(s0k*(2*s01*s1k - s0h*s1k - s1h*s1k - s01*skh + 
             s0h*skh) + s1k*(s0h*s1k + (-s01 + s1h)*skh) + 
          s1h*POWER2(s0k)) + 
       Q2inv*(-2*s01*(s0h + s1h) + 2*POWER2(s01) + 
          POWER2(s0h) + POWER2(s1h))))/(s0h*s1h)
                           );
                           break;
                        case 0:
                           // list=4 i=2 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 1:
                           // list=4 i=2 j=1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 5:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 2:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=5 i=0 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=5 i=1 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=5 i=1 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=5 i=2 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=5 i=2 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 1:
                           // list=5 i=2 j=1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_G);
                           bornTerm=(
(16*(2*(sik*(-(s0h*sik) + (si0 + sih)*skh) + 
          s0k*(2*si0*sik + sih*sik - si0*skh + 
             s0h*(-sik + skh)) + sih*POWER2(s0k)) + 
       Q2inv*(-2*s0h*si0 + 2*si0*sih + POWER2(s0h) + 
          2*POWER2(si0) + POWER2(sih))))/(s0h*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 6:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 2:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=6 i=0 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_G_from_q);
                           bornTerm=(
(16*(2*(s1k*(2*s12*s2k - s1h*s2k - s2h*s2k - s12*skh + 
             s1h*skh) + s2k*(s1h*s2k + (-s12 + s2h)*skh) + 
          s2h*POWER2(s1k)) + 
       Q2inv*(-2*s12*(s1h + s2h) + 2*POWER2(s12) + 
          POWER2(s1h) + POWER2(s2h))))/(s1h*s2h)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=6 i=1 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=6 i=1 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=6 i=2 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=6 i=2 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 1:
                           // list=6 i=2 j=1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 7:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 2:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=7 i=0 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=7 i=1 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=7 i=1 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=7 i=2 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=7 i=2 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 1:
                           // list=7 i=2 j=1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 8:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 2:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=8 i=0 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_G_from_q);
                           bornTerm=(
(16*(2*(s1k*(2*s12*s2k - s1h*s2k - s2h*s2k - s12*skh + 
             s1h*skh) + s2k*(s1h*s2k + (-s12 + s2h)*skh) + 
          s2h*POWER2(s1k)) + 
       Q2inv*(-2*s12*(s1h + s2h) + 2*POWER2(s12) + 
          POWER2(s1h) + POWER2(s2h))))/(s1h*s2h)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=8 i=1 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=8 i=1 j=0
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_G);
                           bornTerm=(
(16*(2*(sik*(-(s2h*sik) + (si2 + sih)*skh) + 
          s2k*(2*si2*sik + sih*sik - si2*skh + 
             s2h*(-sik + skh)) + sih*POWER2(s2k)) + 
       Q2inv*(-2*s2h*si2 + 2*si2*sih + POWER2(s2h) + 
          2*POWER2(si2) + POWER2(sih))))/(s2h*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=8 i=2 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_G_from_q);
                           bornTerm=(
(16*(2*(s0k*(2*s01*s1k - s0h*s1k - s1h*s1k - s01*skh + 
             s0h*skh) + s1k*(s0h*s1k + (-s01 + s1h)*skh) + 
          s1h*POWER2(s0k)) + 
       Q2inv*(-2*s01*(s0h + s1h) + 2*POWER2(s01) + 
          POWER2(s0h) + POWER2(s1h))))/(s0h*s1h)
                           );
                           break;
                        case 0:
                           // list=8 i=2 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 1:
                           // list=8 i=2 j=1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_G);
                           bornTerm=(
(16*(2*(sik*(-(s0h*sik) + (si0 + sih)*skh) + 
          s0k*(2*si0*sik + sih*sik - si0*skh + 
             s0h*(-sik + skh)) + sih*POWER2(s0k)) + 
       Q2inv*(-2*s0h*si0 + 2*si0*sih + POWER2(s0h) + 
          2*POWER2(si0) + POWER2(sih))))/(s0h*sih)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 9:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 2:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=9 i=0 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=9 i=1 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=9 i=1 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 2:
                     switch (j) {
                        case -1:
                           // list=9 i=2 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 0:
                           // list=9 i=2 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                        case 1:
                           // list=9 i=2 j=1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 11:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     // list=11 soft=1 alpha=-1 beta=0
                     softContributes[0]=(TRUE);
                     softCoefficient[0]=(1);
                     softAngle[0]=(vi0);
                     index1[0]=(-1);
                     index2[0]=(0);
                     nCoefficient=(1);
                     bornTerm=(
8*(Q2inv*si0 + 2*s0k*sik)
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=11 i=0 j=-1
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=11 i=1 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_q);
                           bornTerm=(
8*(Q2inv*s0h + 2*s0k*skh)
                           );
                           break;
                        case 0:
                           // list=11 i=1 j=0
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_G_from_q);
                           bornTerm=(
8*(Q2inv*sih + 2*sik*skh)
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
      case 12:
         switch (sub_type) {
            case SOFT_ADDED:
               switch (i) {
                  case 0:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
                  case 1:
                     nCoefficient=(0);
                     bornTerm=(
0
                     );
                     break;
               } // of switch (i)
               break;
            case COLLINEAR_ADDED:
               switch (i) {
                  case 0:
                     switch (j) {
                        case -1:
                           // list=12 i=0 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_G);
                           bornTerm=(
8*(Q2inv*s1h + 2*s1k*skh)
                           );
                           break;
                     } // of switch (j)
                     break;
                  case 1:
                     switch (j) {
                        case -1:
                           // list=12 i=1 j=-1
                           collinearContributes=(TRUE);
                           collinearCoefficient=(2);
                           splittingFunction=(Q_q_from_G);
                           bornTerm=(
8*(Q2inv*s0h + 2*s0k*skh)
                           );
                           break;
                        case 0:
                           // list=12 i=1 j=0
                           collinearContributes=(FALSE);
                           collinearCoefficient=(0);
                           splittingFunction=(NoSplittingFunction);
                           bornTerm=(
0
                           );
                           break;
                     } // of switch (j)
                     break;
               } // of switch (i)
               break;
            default:
               errf(-1, "DISMatrixElement::GetLorentzAddedSubtraction: 1");
         } // of switch (sub_type)
         break;
   } // of switch (list)
} // of function
// DIS Born terms
// Created on {1997, 12, 16, 13, 41, 52}
#include <math.h>
#include <stdio.h>
#include "mathdefs.h"
#include "global.h"
#include "proc.h"
#include "me.h"
double DISMatrixElement::GetBornTerm(
   DISSubProcess iborn
   ) {
   double born = 0;
   switch (iborn) {
      case qXq:
         born=(
8*(s0l*skh + s0k*slh)
);
         break;
      case qXqg:
         born=(
(16*(-(s0h*s0l*s1k) + 2*s0h*s0l*skh + s0l*s1h*skh + 
       s0h*s1l*skh + s1h*s1l*skh + s0h*s1k*slh + 
       s1h*s1k*slh + s0k*
        (2*s0l*s1h + s01*s1l - s0h*s1l - s01*slh + 
          2*s0h*slh + s1h*slh) + 
       s01*(s0l*s1k - s0l*skh - 2*skh*slh)))/(s01*s1h)
);
         break;
      case gXqq:
         born=(
(-16*(s0h*s0l*s1k + s0l*s1h*s1k - 2*s0h*s1k*s1l - 
       s0h*s0l*skh - s1h*s1l*skh - s1h*s1k*slh + 
       s0k*(-2*s0l*s1h - 2*s01*s1l + s0h*s1l + s1h*s1l + 
          s01*slh - s0h*slh) + 
       s01*(-2*s0l*s1k + s0l*skh + s1l*skh + s1k*slh)))/
   (s0h*s1h)
);
         break;
   default:
      errf(-1, "DISMatrixElement::GetBornTerm: not known");
   } // of switch
   return born;
} // of function;
// ============================================================================
//
// --> Interface of PDFLIB to C++
//
// file:              pdflib_cc.cc
// created:           16.04.1997
// last modification: 11.12.1997
//
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fortran.h"

#include "pdflib.h"

// ----------------------------------------------------------------------------

// wrapper for FORTRAN function to set pdf parameters

   PROTOCCALLSFSUB9(PDFSETPAR, pdfsetpar,
                    INT, INT,
                    PINT, PINT, PDOUBLE, PDOUBLE, PDOUBLE, PDOUBLE, PDOUBLE)

#define pdfsetparMacro(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) \
        CCALLSFSUB9(PDFSETPAR, pdfsetpar, \
                    INT, INT, \
                    PINT, PINT, PDOUBLE, PDOUBLE, PDOUBLE,PDOUBLE, PDOUBLE, \
                    arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)

// the following is the function (with C linkage!) to be called 

void pdfsetparC(
        int ipd, int ipdset,
        int& nfl, int& norder,
        double& dlambda4,
        double& xmin, double& xmax, double& scale_min, double& scale_max
     ) {

   double q2min, q2max;

   pdfsetparMacro(
      ipd, ipdset,
      nfl, norder,
      dlambda4,
      xmin, xmax, q2min, q2max
   );

   scale_min = sqrt(q2min);
   scale_max = sqrt(q2max);
}

// ----------------------------------------------------------------------------

// wrapper for FORTRAN function to get PDF distribution

   PROTOCCALLSFSUB3(PDFGET, pdfget,
                    DOUBLE, DOUBLE, DOUBLEV)

#define pdfgetMacro(arg1, arg2, arg3) \
        CCALLSFSUB3(PDFGET, pdfget, \
                    DOUBLE, DOUBLE, DOUBLEV, \
                    arg1, arg2, arg3)

// the following is the function (with C linkage!) to be called.
// note that -6 = tbar, -5 = bbar, ..., 0 = glue, 1 = d, ... 6 = t
// it is assumed that the array range is suitably defined!

void pdfgetC(
        double x, double scale, double *pdf
     ) {

   double *pdfloc;
   pdfloc = pdf - 6;

   pdfgetMacro(x, scale, pdfloc);

   double one_over_x = 1. / x;
   register int i;
   // this loop destroys the value of pdfloc
   for (i = -6; i <= 6; ++i)
       *pdfloc++ *= one_over_x;   
}

// ----------------------------------------------------------------------------

// wrapper for FORTRAN function to set alpha_s parameters (PDFLIB)

   PROTOCCALLSFSUB2(PDFASSETPAR, pdfassetpar,
                    INT,
                    PDOUBLE
                   )

#define pdfassetparMacro( \
           arg1, \
           arg2  \
        ) \
        CCALLSFSUB2(PDFASSETPAR, pdfassetpar, \
                    INT,     \
                    PDOUBLE, \
                    arg1,    \
                    arg2     \
                   )

// the following is the function (with C linkage!) to be called 

void pdfassetparC(
        int ilo,
        double* dlambda4
     ) {

   double dlambda4_R;

   pdfassetparMacro(
      ilo, 
      dlambda4_R
   );
  
   *dlambda4 = dlambda4_R;
}

// ----------------------------------------------------------------------------

// wrapper for FORTRAN function to calculate alpha_s (PDFLIB)

   PROTOCCALLSFFUN1(
      DOUBLE,
      PDFAS, pdfas,
      DOUBLE
   )

#define pdfasMacro(                  \
           arg1                      \
        )                            \
        CCALLSFFUN1(                 \
           PDFAS, pdfas,             \
           DOUBLE,                   \
           arg1                      \
        )

// the following is the function (with C linkage!) to be called 

double pdfasC(
        double scale
     ) {

   return pdfasMacro(scale);
}

// ----------------------------------------------------------------------------

// dummy to avoid warning messages (static functions...)

void dummyfunctionPdflib() {
   static char c[]="a";
   c2fstrv(c,c,0,0);
   f2cstrv(c,c,0,0);
   kill_trailing(c,'a');
   vkill_trailing(c,0,0,'a');
   num_elem(c,0,0,0);
}

// ============================================================================
//        
// --> End of file.
//      
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

// ============================================================================
//
// --> Interface of old FORTRAN routines to C++
//
// file:              fold_cc.cc
// created:           04.05.1997
// last modification: 11.12.1997
//
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fortran.h"
#include "fold.h"

#include "user.h"

// ----------------------------------------------------------------------------

// wrapper for FORTRAN function to get virtual contributions

   PROTOCCALLSFSUB6(EVALVIRT, evalvirt,
                    DOUBLEV,
                    INT, INT, INT,
                    PDOUBLE, PDOUBLE
   )

#define evalvirtMacro(arg1, arg2, arg3, arg4, arg5, arg6) \
        CCALLSFSUB6(EVALVIRT, evalvirt, \
                    DOUBLEV, INT, INT, INT, PDOUBLE, PDOUBLE, \
                    arg1, arg2, arg3, arg4, arg5, arg6 \
        )

// the following is the function (with C linkage!) to be called.

void evalvirtC(
        double *par,
        int i0loop, int i1loop, int ihelicity,
        double& xq, double& xg
     ) {
   evalvirtMacro(par, i0loop, i1loop, ihelicity, xq, xg);
}

// ----------------------------------------------------------------------------

// wrapper for FORTRAN function to get the Spence function

   PROTOCCALLSFFUN1(DOUBLE, SPENCE, spence,
                    DOUBLE
   )

#define RealSpenceMacro(arg1) \
        CCALLSFFUN1(SPENCE, spence, \
                    DOUBLE, \
                    arg1 \
        )

// the following is the function (with C linkage!) to be called.

double RealSpence(double x) {
   return RealSpenceMacro(x);
}

// ----------------------------------------------------------------------------

// wrapper for FORTRAN function for jet clustering

   PROTOCCALLSFSUB7(CIARRAY, ciarray,
                    DOUBLEV, INT, DOUBLE, PINT, DOUBLE, PINT, INT
   )

#define ciarrayMacro(arg1, arg2, arg3, arg4, arg5, arg6, arg7) \
        CCALLSFSUB7(CIARRAY, ciarray, \
                    DOUBLEV, INT, DOUBLE, PINT, DOUBLE, PINT, INT, \
                    arg1, arg2, arg3, arg4, arg5, arg6, arg7 \
        )

// the following is the function (with C linkage!) to be called.

void ciarrayC(
        double *parray, int npa, double dcscale2, 
        int& njetsout, double etain, int& icut, int iflags
     ) {
   ciarrayMacro(parray, npa, dcscale2, njetsout, etain, icut, iflags);
}

// ----------------------------------------------------------------------------
//
// wrapper for FORTRAN functions of the KTCLUS package
//
// ----------------------------------------------------------------------------

   PROTOCCALLSFFUN5(
      INT, KTCLUS, ktclus,
      INT,
      DOUBLEV,
      INT,
      FLOAT, 
      FLOATV
   )

#define ktclus_Macro(      \
           arg1,           \
           arg2,           \
           arg3,           \
           arg4,           \
           arg5            \
        )                  \
        CCALLSFFUN5(       \
           KTCLUS, ktclus, \
           INT,            \
           DOUBLEV,         \
           INT,            \
           FLOAT,          \
           FLOATV,         \
           arg1,           \
           arg2,           \
           arg3,           \
           arg4,           \
           arg5            \
        )

// the following is the function (with C linkage!) to be called.
int ktclus_C(
       int imode,
       double* pp,
       int n,
       float ecut,   
       float* y
    ) {

   return 
      ktclus_Macro(imode, pp, n, ecut, y);
}

// ----------------------------------------------------------------------------

   PROTOCCALLSFFUN10(
      INT, KTRECO, ktreco,
      INT,
      DOUBLEV,
      INT,
      FLOAT, 
      FLOAT, 
      FLOAT, 
      DOUBLEV,
      INTV,
      PINT,
      PINT
   )

#define ktreco_Macro(      \
           arg1,           \
           arg2,           \
           arg3,           \
           arg4,           \
           arg5,           \
           arg6,           \
           arg7,           \
           arg8,           \
           arg9,           \
           arg10           \
        )                  \
        CCALLSFFUN10(      \
           KTRECO, ktreco, \
           INT,            \
           DOUBLEV,        \
           INT,            \
           FLOAT,          \
           FLOAT,          \
           FLOAT,          \
           DOUBLEV,        \
           INTV,           \
           PINT,           \
           PINT,           \
           arg1,           \
           arg2,           \
           arg3,           \
           arg4,           \
           arg5,           \
           arg6,           \
           arg7,           \
           arg8,           \
           arg9,           \
           arg10           \
        )

// the following is the function (with C linkage!) to be called.
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
    ) {

   return 
      ktreco_Macro(
         reco,
         pp,
         nn,
         ecut,
         ycut,
         ymac,
         pjet,
         jet,
         njet,
         nsub
      );
}

// ----------------------------------------------------------------------------
//
// wrapper for FORTRAN functions of the user interface
//
// ----------------------------------------------------------------------------

   PROTOCCALLSFSUB1(USER1, user1,
                    INT
   )

#define user1_Macro(arg1) \
        CCALLSFSUB1(USER1, user1, \
                    INT, \
                    arg1 \
        )

// the following is the function (with C linkage!) to be called.

void user1_C(
        int iaction
     ) {
   user1_Macro(iaction);
}

// ----------------------------------------------------------------------------

   PROTOCCALLSFFUN12(
      DOUBLE, USER2, user2,
      INT,
      INT,
      DOUBLEV, 
      INTV,
      DOUBLEV, 
      DOUBLEV, 
      DOUBLEV, 
      DOUBLEV,
      INTV,
      INTV,
      INTV,
      INTV             
   )

#define user2_Macro(     \
           arg1,         \
           arg2,         \
           arg3,         \
           arg4,         \
           arg5,         \
           arg6,         \
           arg7,         \
           arg8,         \
           arg9,         \
           arg10,        \
           arg11,        \
           arg12         \
        )                \
        CCALLSFFUN12(    \
           USER2, user2, \
           INT,          \
           INT,          \
           DOUBLEV,      \
           INTV,         \
           DOUBLEV,      \
           DOUBLEV,      \
           DOUBLEV,      \
           DOUBLEV,      \
           INTV,         \
           INTV,         \
           INTV,         \
           INTV,         \
           arg1,         \
           arg2,         \
           arg3,         \
           arg4,         \
           arg5,         \
           arg6,         \
           arg7,         \
           arg8,         \
           arg9,         \
           arg10,        \
           arg11,        \
           arg12         \
        )

// the following is the function (with C linkage!) to be called.
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
       ) {
   return
      user2_Macro(
         nr,
         nl,
         fvect,
         npartons,
         xb,  
         q2,  
         xi,  
         weight,
         irps,
         ialphas,
         ialphaem,
         lognf
      );
}

// ----------------------------------------------------------------------------

   PROTOCCALLSFSUB2(USER3, user3,
                    DOUBLE, DOUBLE
   )

#define user3_Macro(arg1, arg2)   \
        CCALLSFSUB2(USER3, user3, \
                    DOUBLE, \
                    DOUBLE, \
                    arg1,   \
                    arg2    \
        )

// the following is the function (with C linkage!) to be called.

void user3_C(
        double value, 
        double error
     ) {
   user3_Macro(value, error);
}

// ----------------------------------------------------------------------------
//
// wrapper for the FORTRAN function to print a string of characters
//
// ----------------------------------------------------------------------------

   PROTOCCALLSFSUB3(
      Underscore_FORTRAN_Name(PRINT_F), 
      Underscore_FORTRAN_Name(print_f),
      STRING,
      INT, 
      INT
   )

#define print_f_Macro(arg1, arg2, arg3)       \
        CCALLSFSUB3(                          \
           Underscore_FORTRAN_Name(PRINT_F),  \
           Underscore_FORTRAN_Name(print_f),  \
           STRING,                            \
           INT,                               \
           INT,                               \
           arg1,                              \
           arg2,                              \
           arg3                               \
        )

// the following is the function (with C linkage!) to be called.

void print_f_C(
        const char* string, 
        int length,
        int unit
     ) {

//   const char prefix[] = "FORTRAN >> ";
   const char prefix[] = "";
   int prefixLength = strlen(prefix);

   // have to copy in order to be able to have the const qualifier

   char* duplicate = new char[prefixLength + length + 1];
   register int i;
   
   for (i = 0; i < prefixLength; ++i)
       duplicate[i] = prefix[i];

   for (i = 0; i < length; ++i)
       duplicate[prefixLength + i] = string[i];

   duplicate[prefixLength + i] = '\0';

#if 0
   // for the time being, terminate after max. 80 characters
   if (prefixLength + length > 80) {
      duplicate[80] = '\0';
      print_f_Macro(duplicate, 80, unit);
   }
   else 
#endif
      print_f_Macro(duplicate, prefixLength + length, unit);

   delete duplicate;
}

void strFORT(const char* str, int length, int unit) {
   printf("FORTRAN [%3d,%3d] >> %s\n", length, unit, str);
   fflush(stdout);
}

void strFORT_x(
        const char* string, 
        int length,
        int unit
     ) {

   unit = unit;  // :_TAWM_:

//   const char prefix[] = "FORTRAN >> ";
   const char prefix[] = "";
   int prefixLength = strlen(prefix);

   // have to copy in order to be able to have the const qualifier

   char* duplicate = new char[prefixLength + length + 1];
   register int i;
   
   for (i = 0; i < prefixLength; ++i)
       duplicate[i] = prefix[i];

   for (i = 0; i < length; ++i)
       duplicate[prefixLength + i] = string[i];

   duplicate[prefixLength + i] = '\0';

#if 0
   // for the time being, terminate after max. 80 characters
   if (prefixLength + length > 80) {
      duplicate[80] = '\0';
      print_f_Macro(duplicate, 80, unit);
   }
   else 
#endif
//   printf("FORTRAN [%3d,%3d] >> %s\n", length, unit, duplicate);

   printf("F %s\n", duplicate);
   fflush(stdout);

   delete duplicate;
}

// ----------------------------------------------------------------------------
//
// wrapper for the FORTRAN function to provoke a floating point exception
//
// ----------------------------------------------------------------------------

   PROTOCCALLSFSUB0(
      Underscore_FORTRAN_Name(PFPX_F), \
      Underscore_FORTRAN_Name(pfpx_f)  \
   )

#define pfpx_f_Macro()                      \
        CCALLSFSUB0(                        \
           Underscore_FORTRAN_Name(PFPX_F), \
           Underscore_FORTRAN_Name(pfpx_f)  \
        )

// the following is the function (with C linkage!) to be called.

void pfpx_f_C() {

   pfpx_f_Macro();
}

// ----------------------------------------------------------------------------
//
// wrapper for C functions of the user interface
//
// ----------------------------------------------------------------------------

// to start DISASTER++

   FCALLSCSUB0(
      disaster_c_from_f_start, 
      Underscore_FORTRAN_Name(DISASTER_CA), 
      Underscore_FORTRAN_Name(disaster_ca)
   )

// ----------------------------------------------------------------------------

// to end DISASTER++

   FCALLSCSUB0(
      disaster_c_from_f_stop, 
      Underscore_FORTRAN_Name(DISASTER_CB), 
      Underscore_FORTRAN_Name(disaster_cb)
   )

// ----------------------------------------------------------------------------

// to set an integer parameter

   FCALLSCSUB2(
      disaster_c_from_f_setInt, 
      Underscore_FORTRAN_Name(DISASTER_CI), 
      Underscore_FORTRAN_Name(disaster_ci),
      STRING, \
      INT \
   )

// ----------------------------------------------------------------------------

// to set a double parameter

   FCALLSCSUB2(
      disaster_c_from_f_setDouble,  
      Underscore_FORTRAN_Name(DISASTER_CD), 
      Underscore_FORTRAN_Name(disaster_cd),
      STRING, \
      DOUBLE \
   )

// ----------------------------------------------------------------------------

// to execute a command

   FCALLSCSUB1(
      disaster_c_from_f_executeCommand, 
      Underscore_FORTRAN_Name(DISASTER_CC), 
      Underscore_FORTRAN_Name(disaster_cc),
      STRING
   )

// ----------------------------------------------------------------------------

// to calculate the adaptation observable

   FCALLSCFUN7(
      DOUBLE,
      disaster_c_from_f_adaptationObservable, 
      Underscore_FORTRAN_Name(DISASTER_CAO), 
      Underscore_FORTRAN_Name(disaster_cao),
      INT, 
      INT, 
      INT, 
      INT, 
      INT, 
      DOUBLE,
      INT
   )

// ----------------------------------------------------------------------------

// to calculate the adaptation observable

   FCALLSCSUB10(
      disaster_c_from_f_FlavourFactors, 
      Underscore_FORTRAN_Name(DISASTER_CFF), 
      Underscore_FORTRAN_Name(disaster_cff),
      INT, 
      INT, 
      INT, 
      INT, 
      INT, 
      DOUBLE,
      INT,
      INT,
      DOUBLEV,
      DOUBLEV
   )

// ----------------------------------------------------------------------------

// dummy to avoid warning messages (static functions...)

void dummyfunctionFold() {
   static char c[]="a";
   c2fstrv(c,c,0,0);
   f2cstrv(c,c,0,0);
   kill_trailing(c,'a');
   vkill_trailing(c,0,0,'a');
   num_elem(c,0,0,0);
}

// ============================================================================
//        
// --> End of file.
//      
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

// ============================================================================
//
// --> container templates
//
// file:              container_t.cc
// created:           10.06.1997
// last modification: 27.11.1997
//
// ============================================================================
  
// load only once...   
#ifndef CONTAINER_T
#define CONTAINER_T

// ============================================================================

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "strng.h"

#include "container.h"

// ============================================================================

void errf(int level, const char* text);

// ============================================================================

// ----------------------------------------------------------------------------
//
// --> Arrays with arbitrary index range
//
// ----------------------------------------------------------------------------

template <class T>
Array <T> :: Array() {

   data=NULL;
}

// ----------------------------------------------------------------------------

template <class T>
Array <T> :: Array(int rminIn, int rmaxIn) {

   data=NULL;
   Define(rminIn, rmaxIn);
}

// ----------------------------------------------------------------------------

template <class T>
Array <T> :: Array(int length) {

   data=NULL;
   Define(0, length-1);
}

// ----------------------------------------------------------------------------

template <class T>
Array <T> :: ~Array() {

   Delete();
}

// ----------------------------------------------------------------------------

template <class T>
void Array <T> :: Define(int rminIn, int rmaxIn) {

   if (data !=NULL)
      errf(-1, "Array::Define: not empty!");
   else {
      rmin=rminIn;
      rmax=rmaxIn;
      data = (new T [rmax-rmin+1]) - rmin;
   }
}

// ----------------------------------------------------------------------------

template <class T>
void Array <T> :: Delete() {

   if (data != NULL) {
      delete [] (data+rmin);
      data=NULL;
   }
}

// ----------------------------------------------------------------------------

template <class T>
void Array <T> :: CreateSame(Array <T> *new_array) {

   new_array = new Array <T>;
   new_array->Define(rmin,rmax);
}

// ----------------------------------------------------------------------------

template <class T>
void Array <T> :: CopyInto(Array <T> * into) {
   
   if (rmin != into -> GetRmin() || rmax != into -> GetRmax())
      errf(-1, "Array<T>::CopyInto: boundaries do not match.");

   for (register int i = rmin; i <= rmax; ++i)
       into -> SetData(i, data[i]);
}

// ----------------------------------------------------------------------------
//
// --> Array of pointers with arbitrary index range
//
// ----------------------------------------------------------------------------

template <class T>
ArrayOfPointers <T> :: ArrayOfPointers(int rminIn, int rmaxIn)
   : Array <T> (rminIn, rmaxIn) {
}

// ----------------------------------------------------------------------------

template <class T>
void ArrayOfPointers <T> :: DeleteAllEntries() {

   for (register int i = rmin; i <= rmax; ++i) {
       delete data[i];
       data[i] = NULL;
   }
}

// ----------------------------------------------------------------------------

template <class T>
void ArrayOfPointers <T> :: SetToNull() {

   for (register int i = rmin; i <= rmax; ++i) {
       data[i] = NULL;
   }
}

// ----------------------------------------------------------------------------
//
// --> Two-dimensional arrays with arbitrary index range
//
// ----------------------------------------------------------------------------

template <class T>
TwoD_Array <T> :: TwoD_Array() {

   data=NULL;
}

// ----------------------------------------------------------------------------

template <class T>
TwoD_Array<T> :: TwoD_Array(int r0minIn, int r0maxIn,
                            int r1minIn, int r1maxIn) {
   data=NULL;
   Define(r0minIn, r0maxIn,
          r1minIn, r1maxIn);
}

// ----------------------------------------------------------------------------

template <class T>
TwoD_Array <T> :: ~TwoD_Array() {

   Delete();
}

// ----------------------------------------------------------------------------

template <class T>
void TwoD_Array<T> :: Define(int r0minIn, int r0maxIn,
                             int r1minIn, int r1maxIn) {

   if (data !=NULL)
      errf(-1, "TwoD_Array::Define: not empty!");
   else {
      r0min=r0minIn;
      r0max=r0maxIn;
      r1min=r1minIn;
      r1max=r1maxIn;
 
      // create only if not empty; else data is the NULL pointer
      if ( r0max >= r0min && r1max >= r1min ) {

         data = (new T* [r0max-r0min+1]) - r0min;

         register int i;
         for (i=r0min; i<=r0max; ++i)
             data[i] = (new T [r1max-r1min+1]) - r1min;
      }
   }
}

// ----------------------------------------------------------------------------

template <class T>
void TwoD_Array <T> :: Delete() {

   if (data != NULL) {
      register int i;
      for (i=r0min; i<=r0max; ++i)
          delete (data[i] + r1min);
      delete [] (data + r0min);
      data=NULL;
   }
}

// ----------------------------------------------------------------------------

template <class T>
T **TwoD_Array <T> :: GetDataPtrWithNullWarning() {

  if (data==NULL)
     errf(-1, "TwoD_Array::GetDataPtrWithNullWarning: data ptr is NULL"); 

  return data;
}
   
// ----------------------------------------------------------------------------
//
// --> two-dimensional array of pointers with arbitrary index range
//
// ----------------------------------------------------------------------------

template <class T>
void TwoD_ArrayOfPointers <T> :: DeleteAllEntries() {

   if (data != NULL) {
      for (register int i = r0min; i <= r0max; ++i)
          for (register int j = r1min; j <= r1max; ++j) {
              delete data[i][j];
              data[i][j] = NULL;
          }
   }
}

// ----------------------------------------------------------------------------
//
// --> Free list container
//
// ----------------------------------------------------------------------------

extern long coCont;

// ----------------------------------------------------------------------------

template <class T>
FreeListContainer <T> :: FreeListContainer() {

   ident = coCont++;
   data = NULL;
}

// ----------------------------------------------------------------------------

template <class T>
FreeListContainer <T> :: ~FreeListContainer() {
}

// ----------------------------------------------------------------------------
//
// --> Simple free list
//
// ----------------------------------------------------------------------------

template <class T>
SimpleFreeList <T> :: SimpleFreeList(Global* globalIn) {

   defined = FALSE;
   Define(globalIn);
}

// ----------------------------------------------------------------------------

template <class T>
SimpleFreeList <T> :: ~SimpleFreeList() {

   Delete();
}

// ----------------------------------------------------------------------------
               
// :_TAWM_: this is defined to avoid a warning message from gcc

template <class T>
SimpleFreeList <T> & SimpleFreeList <T>
                        :: operator=(const SimpleFreeList <T> & object) {
   
   const SimpleFreeList <T> * dummy = &object;
   dummy = dummy;
         
   errf(-1, "SimpleFreeList <T>: assignment operator not defined!");
         
   return *this;  // :_TAWM_:
}
      
// ----------------------------------------------------------------------------

template <class T>
void SimpleFreeList <T> :: Define(Global* globalIn) {

   if (defined) 
      errf(-1, "SimpleFreeList <T> :: Define: error");

   Set_global(globalIn);

   defined = FALSE;
   Reset();
   defined = TRUE;
}

// ----------------------------------------------------------------------------

template <class T>
void SimpleFreeList <T> :: Delete() {

   global -> To_status() << FS("SimpleFreeList Delete: %s\n", GetNameSafe());
   // check numbers
   if (   listFullLength != createCount 
       || listEmptyLength+listFullLength != containerCount 
       || TRUE //&& FALSE
      ) {
      global -> To_status() 
         << FS("# of requested objects:     %8ld\n", getCount);
      global -> To_status() 
         << FS("# of returned objects:      %8ld\n", returnCount);
      global -> To_status() 
         << FS("# of created objects:       %8ld\n", createCount);
      global -> To_status() 
         << FS("# of objects in Full list:  %8ld\n", listFullLength);
      global -> To_status() 
         << FS("  mismatch:                 %8ld\n", 
               listFullLength - createCount);
      global -> To_status() 
         << FS("# of created containers:    %8ld\n", containerCount);
      global -> To_status() 
         << FS("# of objects in Empty list: %8ld\n", listEmptyLength);
      global -> To_status() 
         << FS("  mismatch:                 %8ld\n",
               listEmptyLength + listFullLength - containerCount);
   } 

   // delete lists
   FreeListContainer <T> * ptemp;
   long co;
   co = 0;
   while (listFull != NULL) {
      ++co;
      ptemp = listFull;
      listFull = listFull -> GetNext();
      delete ptemp -> GetDataPtr();
      delete ptemp;
   }
   if (listFullLength != co)
      global -> To_status() << FS("DELETE Full: %ld ", listFullLength)
                            << FS("%ld\n", co);
   co = 0;
   while (listEmpty != NULL) {
      ++co;
      ptemp = listEmpty;
      listEmpty = listEmpty -> GetNext();
      delete ptemp;
   }
   if (listEmptyLength != co)
      global -> To_status() << FS("DELETE Empty: %ld ", listEmptyLength)
                            << FS("%ld\n", co);
}

// ----------------------------------------------------------------------------

template <class T>
void SimpleFreeList <T> :: Reset() {

   if (defined) 
      Delete();

   listFull = NULL;
   listEmpty = NULL;

   createCount = 0;
   containerCount = 0;

   listFullLength = 0;
   listEmptyLength = 0;

   getCount = 0;
   returnCount = 0;
}

// ----------------------------------------------------------------------------

template <class T>
int SimpleFreeList <T> :: GetObject(T* &object) {

   int ifNewObject;
   FreeListContainer <T> *ptemp;

   if (listFull==NULL) {

      createCount++;

      object = new T;     

      ifNewObject = TRUE;

      if (trapCreate >= 0 && trapCreate == object -> GetIdent()) {
         global -> To_error() 
            << FS(":%s:\n", GetNameSafe());
         global -> To_error() 
            << FS("FreeList (Create): trap found %ld\n", trapCreate);
         provoke_segmentation_fault();
      }

   } else {
     
      ++listEmptyLength;
      --listFullLength;

      object = listFull -> GetDataPtr();
      listFull -> SetDataPtr(NULL);

      ptemp = listEmpty;
      listEmpty = listFull;
      listFull = listFull -> GetNext();
      listEmpty -> SetNext(ptemp);

      ifNewObject = FALSE;

      if (trapGetNoCreate >= 0 && trapGetNoCreate == object -> GetIdent()) {
         global -> To_error() 
            << FS(":%s:\n", GetNameSafe());
         global -> To_error()
            << FS("FreeList (GetNoCreate): trap found %ld\n", trapGetNoCreate);
         provoke_segmentation_fault();
      }
   }

   ++ getCount;

   return ifNewObject;
}

// ----------------------------------------------------------------------------

template <class T>
void SimpleFreeList <T> :: ReturnObject(T *object) {

   object -> ReturnToFreeList();

   if (trapReturn >= 0 && trapReturn == object -> GetIdent()) {
      global -> To_error() 
         << FS(":%s:\n", GetNameSafe());
      global -> To_error() 
         << FS("FreeList (Return): trap found %ld\n", trapReturn);
      provoke_segmentation_fault();
   }

   ++listFullLength;
   
   FreeListContainer <T> *ptemp;

   // get a container
   if (listEmpty == NULL) { 
      ptemp = new FreeListContainer <T>;
      ++ containerCount;
   } else {
      --listEmptyLength;
      ptemp = listEmpty;
      listEmpty = listEmpty -> GetNext();
   }

   // assign and include in free list
   ptemp -> SetDataPtr(object);
   ptemp -> SetNext(listFull);
   listFull=ptemp;
 
   ++ returnCount;
}

// ----------------------------------------------------------------------------

template <class T>
void SimpleFreeList <T> :: Print(FILE *out) {

   Mstream s_out;
   s_out.Associate(out);

   FreeListContainer <T> * p;

   s_out << "\n-------------------------------------\n";
   s_out << FS("objects    created: %8ld\n", createCount);
   s_out << FS("containers created: %8ld\n", containerCount);
   s_out << "-------------------------------------\n";
   s_out << FS("Full (%ld):\n", listFullLength);

   p = listFull;
   while (p != NULL) {
      s_out << FS("%ld --> ", p -> GetIdent())
            << FS("%ld\n", p -> GetDataPtr() -> GetIdent());
      p = p -> GetNext();
   }   

   s_out << "-------------------------------------\n";
   s_out << FS("Empty (%ld):\n", listEmptyLength);

   p = listEmpty;
   while (p != NULL) {
      s_out << FS("%ld --> NULL\n", p -> GetIdent());
      p = p -> GetNext();
   }   

   s_out << "-------------------------------------\n\n";
}

// ----------------------------------------------------------------------------
//
// --> Free list
//
// ----------------------------------------------------------------------------

template <class T>
FreeList <T> :: FreeList(Global* globalIn)
   : SimpleFreeList <T> (globalIn)
   {
}

// ----------------------------------------------------------------------------

template <class T>
FreeList <T> :: ~FreeList() {
}

// ----------------------------------------------------------------------------

template <class T>
void FreeList <T> :: ReturnObjectsFromUseList(
                        HorizontalList <T> * (T :: *hf)(),
                        UseList <T> * ul
                     ) {

   ul -> TraverseReset(hf);
   T* t = ul -> TraverseNext();
   while (t != NULL) {  
      ReturnObject(t); 
      t = ul -> TraverseNext(); 
   }
   ul -> DestructContainerList();
}

// ----------------------------------------------------------------------------

template <class T>
void FreeList <T> :: ReturnObjectsFromHorizontalList(
                        HorizontalList <T> * (T :: *hf)(),
                        T* hl
                     ) {

   T* hlnext;
   if (hl != NULL) {
      hlnext = (hl ->* hf)() -> GetHorizontalObject();
      while (hl != NULL) {  
         ReturnObject(hl); 
         hl = hlnext;
         if (hl != NULL) 
            hlnext = (hl ->* hf)() -> GetHorizontalObject();
      }
   }
}

// ----------------------------------------------------------------------------
//
// --> Free list: no parameters
//
// ----------------------------------------------------------------------------

template <class T>
FreeList_0 <T> :: FreeList_0(Global* globalIn) 
   : FreeList <T> (globalIn)
   {
}

// ----------------------------------------------------------------------------

template <class T>
FreeList_0 <T> :: ~FreeList_0() {
}

// ----------------------------------------------------------------------------

template <class T>
T *FreeList_0 <T> :: GetDefinedObject() {

   T *object;

   if (GetObject(object))
      object -> Define();
   else
      object -> Reset();

   return object;
}

// ----------------------------------------------------------------------------
//
// --> Free list: one parameter
//
// ----------------------------------------------------------------------------

template <class T, class U>
FreeList_1 <T, U> :: FreeList_1(Global* globalIn) 
   : FreeList <T> (globalIn)
   {
}

// ----------------------------------------------------------------------------

template <class T, class U>
FreeList_1 <T, U> :: ~FreeList_1() {
}

// ----------------------------------------------------------------------------

template <class T, class U>
T *FreeList_1 <T, U> :: GetDefinedObject(U parIn) {

   T *object;

   if (GetObject(object))
      object -> Define(parIn);
   else
      object -> Reset(parIn);

   return object;
}

// ----------------------------------------------------------------------------
//
// --> Array of free lists
//
// ----------------------------------------------------------------------------

template <class T>
FreeListArray <T> :: FreeListArray(Global* globalIn) {

   Set_global(globalIn);
   defined = FALSE;
}

// ----------------------------------------------------------------------------

template <class T>
FreeListArray <T> :: FreeListArray(Global* globalIn, int rlength) {

   Set_global(globalIn);
   defined = FALSE;
   Define(0, rlength-1);
}

// ----------------------------------------------------------------------------

template <class T>
FreeListArray <T> :: FreeListArray(Global* globalIn, int rminIn, int rmaxIn) {

   Set_global(globalIn);
   defined = FALSE;
   Define(rminIn, rmaxIn);
}

// ----------------------------------------------------------------------------

template <class T>
FreeListArray <T> :: ~FreeListArray() {

   Delete();
}

// ----------------------------------------------------------------------------

// :_TAWM_: this is defined to avoid a warning message from gcc

template <class T>
FreeListArray <T> & FreeListArray <T> 
                        :: operator=(const FreeListArray <T> & object) {

   const FreeListArray <T> * dummy = &object;
   dummy = dummy;

   errf(-1, "FreeListArray <T>: assignment operator not defined!");

   return *this;  // :_TAWM_:
}

// ----------------------------------------------------------------------------

template <class T>
void FreeListArray <T> :: Define(int rminIn, int rmaxIn) {

   if (defined)
      errf(-1, "FreeListArray <T> :: Define: already defined!");

   defined = TRUE;

   rmin = rminIn;
   rmax = rmaxIn;

   // flist is an array of pointers to potential free lists.
   // Initialized to NULL.
   flist = new ArrayOfPointers < FreeList <T> *> (rmin, rmax);
   flptr = flist -> GetDataPtr();
   register int i;
   for (i = rmin; i <= rmax; ++i)
       flptr[i] = NULL;
}

// ----------------------------------------------------------------------------

template <class T>
void FreeListArray <T> :: Delete() {

   if (!defined)
      errf(-1, "FreeListArray <T> :: Delete: not defined!");

   global -> To_status() << FS("FreeListArray Delete: %s\n", GetNameSafe());
   flist -> DeleteAllEntries();
   delete flist;
}

// ----------------------------------------------------------------------------
   
template <class T>
FreeList <T> * FreeListArray <T> :: CreateNewFreeList(int parIn) {

   FreeList <T> * flp;
   flp = new FreeList <T> (global);
   flp -> DefineName(String(GetNameSafe(), String(" (array element)")));
   flp -> GetTrapsFrom(this);
   global -> To_status() 
      << FS("Created new entry (%d)", parIn)
      << FS("\nof FreeListArray :%s:\n", flp -> GetNameSafe());
   return flp;
}

// ----------------------------------------------------------------------------
   
template <class T>
int FreeListArray <T> :: GetObject(T* &object, int parIn) {

   if (parIn < rmin || parIn > rmax) {
      global -> To_error() << FS("%d ", parIn)
                           << FS("%d ", rmin)
                           << FS("%d\n", rmax);
      errf(-1, "FreeListArray <T> :: GetDefinedObject: out of range");
   }

   // create a new free list if required
   if (flptr[parIn] == NULL)
      flptr[parIn] = CreateNewFreeList(parIn);

   // now get an object from the free list
   int ifNewObject;
   ifNewObject = flptr[parIn] -> GetObject(object);

   return ifNewObject;
}

// ----------------------------------------------------------------------------
   
template <class T>
void FreeListArray <T> :: ReturnObject(T* object) {

   int par = object -> GetDefineParameter();

   if (par < rmin || par > rmax) {
      global -> To_error() << FS("%d ", par)
                           << FS("%d ", rmin)
                           << FS("%d\n", rmax);
      errf(-1, "FreeListArray <T> :: ReturnObject: out of range");
   }

   // create a new free list if required
   if (flptr[par] == NULL)
      flptr[par] = CreateNewFreeList(par);

   flptr[par] -> ReturnObject(object);
}

// ----------------------------------------------------------------------------
   
template <class T>
void FreeListArray <T> :: Print(FILE* out) {

   Mstream s_out;
   s_out.Associate(out);

   int par; 

   for (par = rmin; par <= rmax; ++par) {
       s_out << FS("\npar=%d\n", par);
       if (flptr[par] == NULL)
          s_out << "no free list allocated\n";
       else
          flptr[par] -> Print(out); 
   }
}

// ----------------------------------------------------------------------------

template <class T>
void FreeListArray <T> :: ReturnObjectsFromUseList(
                             HorizontalList <T> * (T :: *hf)(),
                             UseList <T> * ul
                          ) {

   ul -> TraverseReset(hf);
   T* t = ul -> TraverseNext();
   while (t != NULL) {  
      ReturnObject(t); 
      t = ul -> TraverseNext(); 
   }
   ul -> DestructContainerList();
}

// ----------------------------------------------------------------------------

template <class T>
void FreeListArray <T> :: ReturnObjectsFromHorizontalList(
                             HorizontalList <T> * (T :: *hf)(),
                             T* hl
                          ) {

   T* hlnext;
   if (hl != NULL) {
      hlnext = (hl ->* hf)() -> GetHorizontalObject();
      while (hl != NULL) {  
         ReturnObject(hl); 
         hl = hlnext;
         if (hl != NULL) 
            hlnext = (hl ->* hf)() -> GetHorizontalObject();
      }
   }
}

// ----------------------------------------------------------------------------
//
// --> Two-dimensional array of free lists
//
// ----------------------------------------------------------------------------

template <class T>
FreeListTwoD_Array <T> :: FreeListTwoD_Array(Global* globalIn) {

   Set_global(globalIn);
   defined = FALSE;
}

// ----------------------------------------------------------------------------

template <class T>
FreeListTwoD_Array <T> :: FreeListTwoD_Array(
                             Global* globalIn,
                             int r0minIn, int r0maxIn,
                             int r1minIn, int r1maxIn
                          ) {

   Set_global(globalIn);
   defined = FALSE;
   Define(r0minIn, r0maxIn,
          r1minIn, r1maxIn
   );
}

// ----------------------------------------------------------------------------

template <class T>
FreeListTwoD_Array <T> :: FreeListTwoD_Array(
                             Global* globalIn,
                             int r0length, int r1length
                          ) {

   Set_global(globalIn);
   defined = FALSE;
   Define(0, r0length-1, 0, r1length-1);
}

// ----------------------------------------------------------------------------

template <class T>
FreeListTwoD_Array <T> :: ~FreeListTwoD_Array() {

   Delete();
}

// ----------------------------------------------------------------------------

// :_TAWM_: this is defined to avoid a warning message from gcc

template <class T>
FreeListTwoD_Array <T> & FreeListTwoD_Array <T> 
                        :: operator=(const FreeListTwoD_Array <T> & object) {

   const FreeListTwoD_Array <T> * dummy = &object;
   dummy = dummy;

   errf(-1, "FreeListTwoD_Array <T>: assignment operator not defined!");

   return *this;  // :_TAWM_:
}

// ----------------------------------------------------------------------------

template <class T>
void FreeListTwoD_Array <T> :: Define(int r0minIn, int r0maxIn,
                                      int r1minIn, int r1maxIn
                               ) {

   if (defined)
      errf(-1, "FreeListTwoD_Array <T> :: Define: already defined!");

   defined = TRUE;

   r0min = r0minIn;
   r0max = r0maxIn;
   r1min = r1minIn;
   r1max = r1maxIn;

   // flist is an array of pointers to potential free lists.
   // Initialized to 0.
   flist = new TwoD_ArrayOfPointers < FreeList <T> *> (r0min, r0max,
                                                       r1min, r1max);
   flptr = flist -> GetDataPtr();
   register int i, j;
   for (i=r0min; i<=r0max; ++i)
       for (j=r1min; j<=r1max; ++j)
           flptr[i][j]=NULL;
}

// ----------------------------------------------------------------------------

template <class T>
void FreeListTwoD_Array <T> :: Delete() {

   if (!defined)
      errf(-1, "FreeListTwoD_Array <T> :: Delete: not defined!");

   global -> To_status() 
      << FS("FreeListTwoD_Array Delete: %s\n", GetNameSafe());
   flist -> DeleteAllEntries();
   delete flist;
}

// ----------------------------------------------------------------------------
   
template <class T>
FreeList <T> * FreeListTwoD_Array <T> :: CreateNewFreeList(
                                            int par0In, int par1In
                                         ) {
   FreeList <T> * flp;
   flp = new FreeList <T> (global);
   flp -> DefineName(String(GetNameSafe(), String(" (array element)")));
   flp -> GetTrapsFrom(this);
   global -> To_status() 
      << FS("Created new entry (%d, ", par0In)
      << FS("%d)\n", par1In)
      << FS("of FreeListArray :%s:\n", flp -> GetNameSafe());
   return flp;
}

// ----------------------------------------------------------------------------
   
template <class T>
int FreeListTwoD_Array <T> :: GetObject(T* &object, int par0In, int par1In) {

   if (   par0In < r0min || par0In > r0max
       || par1In < r1min || par1In > r1max
      ) {
      global -> To_error() 
         << FS("%d ", par0In)
         << FS("%d ", r0min)
         << FS("%d ", r0max)
         << FS("%d ", par1In)
         << FS("%d ", r1min)
         << FS("%d \n", r1max);
      errf(-1, "FreeListTwoD_Array <T> :: GetDefinedObject: out of range");
   }

   // create a new free list if required
   if (flptr[par0In][par1In] == NULL)
      flptr[par0In][par1In] = CreateNewFreeList(par0In, par1In);

   // now get an object from the free list
   int ifNewObject;
   ifNewObject = flptr[par0In][par1In] -> GetObject(object);

   return ifNewObject;
}

// ----------------------------------------------------------------------------
   
template <class T>
void FreeListTwoD_Array <T> :: ReturnObject(T* object) {

   int par0, par1;
   object -> GetDefineParameter(par0, par1);

   if (   par0 < r0min || par0 > r0max
       || par1 < r1min || par1 > r1max
      ) {
      global -> To_error() 
         << FS("%d ", par0)
         << FS("%d ", r0min)
         << FS("%d ", r0max)
         << FS("%d ", par1)
         << FS("%d ", r1min)
         << FS("%d \n", r1max);
      errf(-1, "FreeListTwoD_Array <T> :: ReturnObject: out of range");
   }

   // create a new free list if required
   if (flptr[par0][par1] == NULL)
      flptr[par0][par1] = CreateNewFreeList(par0, par1);

   flptr[par0][par1] -> ReturnObject(object);
}

// ----------------------------------------------------------------------------
   
template <class T>
void FreeListTwoD_Array <T> :: Print(FILE* out) {

   Mstream s_out;
   s_out.Associate(out);

   int par0, par1; 

   for (par0 = r0min; par0 <= r0max; ++par0) {
       for (par1 = r1min; par1 <= r1max; ++par1) {
           s_out << FS("\npar0=%d ", par0)
                 << FS("par1=%d\n", par1);
           if (flptr[par0][par1] == NULL)
              s_out << "no free list allocated\n";
           else
              flptr[par0][par1] -> Print(out); 
      }
   }
}

// ----------------------------------------------------------------------------

template <class T>
void FreeListTwoD_Array <T> :: ReturnObjectsFromUseList(
                                  HorizontalList <T> * (T :: * hf)(),
                                  UseList <T> * ul
                               ) {

   ul -> TraverseReset(hf);
   T* t = ul -> TraverseNext();
   while (t != NULL) {  
      ReturnObject(t); 
      t = ul -> TraverseNext(); 
   }
   ul -> DestructContainerList();
}

// ----------------------------------------------------------------------------

template <class T>
void FreeListTwoD_Array <T> :: ReturnObjectsFromHorizontalList(
                                  HorizontalList <T> * (T :: *hf)(),
                                  T* hl
                               ) {

   T* hlnext;
   if (hl != NULL) {
      hlnext = (hl ->* hf)() -> GetHorizontalObject();
      while (hl != NULL) {  
         ReturnObject(hl); 
         hl = hlnext;
         if (hl != NULL) 
            hlnext = (hl ->* hf)() -> GetHorizontalObject();
      }
   }
}

// ----------------------------------------------------------------------------
//
// --> Array of free lists: no explicit parameter
//
// ----------------------------------------------------------------------------

template <class T>
FreeListArray_0 <T> :: FreeListArray_0(
                          Global* globalIn, 
                          int rminIn, int rmaxIn
                       ) 
   : FreeListArray <T> (globalIn, rminIn, rmaxIn) {
}

// ----------------------------------------------------------------------------

template <class T>
FreeListArray_0 <T> :: ~FreeListArray_0() {
}

// ----------------------------------------------------------------------------

template <class T>
T *FreeListArray_0 <T> :: GetDefinedObject(int parIn) {

   T *object;

   if ( GetObject(object, parIn) )
      object -> Define(parIn);
   else
      object -> Reset();

   return object;
}

// ----------------------------------------------------------------------------
//
// --> Array of free lists: one explicit parameter
//
// ----------------------------------------------------------------------------

template <class T, class U>
FreeListArray_1 <T, U> :: FreeListArray_1(Global* globalIn, int rlength) 
   : FreeListArray <T> (globalIn, rlength) {
}

// ----------------------------------------------------------------------------

template <class T, class U>
FreeListArray_1 <T, U> :: FreeListArray_1(
                             Global* globalIn,
                             int rminIn, int rmaxIn
                          ) 
   : FreeListArray <T> (globalIn, rminIn, rmaxIn) {
}

// ----------------------------------------------------------------------------

template <class T, class U>
FreeListArray_1 <T, U> :: ~FreeListArray_1() {
}

// ----------------------------------------------------------------------------

template <class T, class U>
T *FreeListArray_1 <T, U> :: GetDefinedObject(U explicitParIn, int parIn) {

   T *object;

   if ( GetObject(object, parIn) )
      object -> Define(explicitParIn, parIn);
   else
      object -> Reset(explicitParIn);

   return object;
}

// ----------------------------------------------------------------------------
//
// --> Two-dimensional array of free lists: one explicit parameter
//
// ----------------------------------------------------------------------------

template <class T, class U>
FreeListTwoD_Array_1 <T, U> :: FreeListTwoD_Array_1(Global* globalIn)
   : FreeListTwoD_Array <T> (globalIn)
   {
}

// ----------------------------------------------------------------------------

template <class T, class U>
FreeListTwoD_Array_1 <T, U> :: FreeListTwoD_Array_1(
                                  Global* globalIn,
                                  int r0minIn, int r0maxIn,
                                  int r1minIn, int r1maxIn
                               ) 
   : FreeListTwoD_Array <T> (globalIn,
                             r0minIn, r0maxIn,  
                             r1minIn, r1maxIn) {
}

// ----------------------------------------------------------------------------

template <class T, class U>
FreeListTwoD_Array_1 <T, U> 
   :: FreeListTwoD_Array_1(Global* globalIn, int r0length, int r1length)
    : FreeListTwoD_Array <T> (globalIn, r0length, r1length) {
}

// ----------------------------------------------------------------------------

template <class T, class U>
FreeListTwoD_Array_1 <T, U> :: ~FreeListTwoD_Array_1() {
}

// ----------------------------------------------------------------------------

template <class T, class U>
T *FreeListTwoD_Array_1 <T, U> 
     :: GetDefinedObject(U explicitParIn, int par0In, int par1In) {

   T *object;

   if ( GetObject(object, par0In, par1In) )
      object -> Define(explicitParIn, par0In, par1In);
   else
      object -> Reset(explicitParIn);

   return object;
}

// ----------------------------------------------------------------------------
//
// --> UseList container
//
// ----------------------------------------------------------------------------

template <class T>
UseListContainer <T> :: UseListContainer() {

   ident = coCont++;
   data = NULL;
}

// ----------------------------------------------------------------------------

template <class T>
UseListContainer <T> :: ~UseListContainer() {
}

// ----------------------------------------------------------------------------

template <class T>
void UseListContainer <T> :: ReturnToFreeList() {
}

// ----------------------------------------------------------------------------
//
// --> list of items (UseList)
//
// ----------------------------------------------------------------------------

template <class T>
UseList <T> :: UseList(Global* globalIn) {

   Set_global(globalIn);

   global -> To_status() << "UseList: Constructor called\n";

   list = NULL;

   // free list for containers
   freeList = new SimpleFreeList < UseListContainer <T> > (global);  
}

// ----------------------------------------------------------------------------

template <class T>
UseList <T> :: ~UseList() {

   global -> To_status() << FS("UseList destructor: :%s:\n", GetNameSafe());

#if 0
   // :_CHECK_: this led to a problem (is this a compiler bug?)
   // (cf. dmalloc: Error: tried to free pointer 
   // which is already freed (err 61))
   freeList -> DefineName(String((String) String(GetNameSafe()),
                                 (String) String(": UseList container list")
                          )
               );
#else
   String s1(GetNameSafe());
   String s2(": UseList container list");
   freeList -> DefineName(String(s1, s2));
#endif

   delete freeList;
}

// ----------------------------------------------------------------------------

template <class T>
void UseList <T> :: InsertObject(T* object) {

   UseListContainer <T> *container;
   freeList -> GetObject(container);

   container -> SetDataPtr(object);   

   container -> SetNext(list);   
   list = container;
}

// ----------------------------------------------------------------------------

template <class T>
void UseList <T> :: AppendObject(T* object) {

   // find a pointer to the last container in the list
   UseListContainer <T> *last;
   last = list;
   while (last != NULL && last -> GetNext() != NULL)
         last = last -> GetNext();

   UseListContainer <T> *container;
   freeList -> GetObject(container);

   container -> SetDataPtr(object);   
   container -> SetNext(NULL);   

   if (list == NULL)
      list = container;
   else
      last -> SetNext(container);
}

// ----------------------------------------------------------------------------

template <class T>
T* UseList <T> :: GetFirstObject() {

   if (list == NULL)
      errf(-1, "UseList <T> :: GetFirstObject(): error");

   UseListContainer <T> *container;
   container = list;
   list = container -> GetNext();

   T *data;
   data = container -> GetDataPtr();

   freeList -> ReturnObject(container);

   return data;
}

// ----------------------------------------------------------------------------

template <class T>
T* UseList <T> :: GetLastObject() {

   if (list == NULL)
      errf(-1, "UseList <T> :: GetLastObject(): error");

   // find a pointer to the ``next-to-last'' container in the list
   UseListContainer <T> *last;
   last = list;
   UseListContainer <T> *before_last;
   before_last = NULL;
   while (last != NULL && last -> GetNext() != NULL) {
         before_last = last;
         last = last -> GetNext();
   }

   T *data;
   data = last -> GetDataPtr();

   if (before_last == NULL)
      list = NULL;
   else
      before_last -> SetNext(NULL);

   freeList -> ReturnObject(last);

   return data;
}

// ----------------------------------------------------------------------------

// this destructs the linked list of containers; the appended objects
// of type T are not touched

template <class T>
void UseList <T> :: DestructContainerList() {

   UseListContainer <T> *p;
   UseListContainer <T> *q;
   p = list;
   while (p!=NULL) {
      q = p;
      p = p -> GetNext();
      freeList -> ReturnObject(q);
   } 
   list = NULL;  
}

// ----------------------------------------------------------------------------

template <class T>
void UseList <T> :: Print(FILE* out) {

   Mstream s_out;
   s_out.Associate(out);

   UseListContainer <T> *p;
   s_out << "-------------------------------------\n";
   s_out << "UseList print:\n";
   p = list;
   while (p!=NULL) {
      p -> GetDataPtr() -> Print(out);
      p = p -> GetNext();
   }   
   s_out << "-------------------------------------\n";
}

// ----------------------------------------------------------------------------

template <class T>
void UseList <T> :: PrintIdent(
                       FILE* out,
                       HorizontalList <T> * (T :: *hf)()
                    ) {

   Mstream s_out;
   s_out.Associate(out);

   const char str0[] = "";
   const char str1[] = "   ";

   int count;
   UseListContainer <T> * actual = list;
   HorizontalList <T> * loc;

   s_out << "UseList :: Print\n";

   while (actual != NULL) {
      count = 0;
      loc = (actual -> GetDataPtr() ->* hf)();
      while (loc != NULL) {
         s_out << FS("--- %s ", (count == 0) ? str0 : str1)
               << FS("%8ld\n", loc -> GetObject() -> GetIdent());
         ++ count;
         loc = loc -> GetHorizontalPtr();
      }
      actual = actual -> GetNext();
   }   
}

// ----------------------------------------------------------------------------

template <class T>
T* UseList <T> :: FindSame(T* prot) {

   UseListContainer <T> * p = list;

   while ( p != NULL && ! p -> GetDataPtr() -> UseListIsSame(prot) ) 
      p = p -> GetNext();

   if (p == NULL) 
      return NULL;
   else
      return p -> GetDataPtr();
}

// ----------------------------------------------------------------------------

template <class T>
T* UseList <T> :: FindIdent(
                     HorizontalList <T> * (T :: *hf)(),
                     long ident
                  ) {

   TraverseReset(hf);
   T* t = TraverseNext();
   while (t != NULL && t -> GetIdent() != ident)
      t = TraverseNext();

   return t;
}

// ----------------------------------------------------------------------------

template <class T>
void UseList <T> :: TraverseReset(
                       HorizontalList <T> * (T :: *hFunctionIn)()
                    ) {

   hFunction = hFunctionIn;

   tVertical = list;
   if (tVertical != NULL) {
      tHorizontal = tVertical -> GetDataPtr();
      tVerticalNext = tVertical -> GetNext();
   }
   else
      tHorizontal = NULL;
}

// ----------------------------------------------------------------------------

// this method is safe in the sense that the object referred to by the
// pointer may be destructed without interfering with the method
// (implemented via look-ahead)

template <class T>
T* UseList <T> :: TraverseNext() {

   T* ret;

   if (tVertical == NULL)
      ret = NULL;
   else {
      if (tHorizontal != NULL) {
         ret = tHorizontal;
         tHorizontal = (tHorizontal ->* hFunction)() 
                          -> GetHorizontalObject();
      } else {
         tVertical = tVertical -> GetNext();
         tVertical = tVerticalNext;
         if (tVertical != NULL) {
            tVerticalNext = tVertical -> GetNext();
            ret = tVertical -> GetDataPtr();
            tHorizontal = (ret ->* hFunction)() 
                             -> GetHorizontalObject();
         } else
            ret = NULL;
            // nothing else, because for tVertical == NULL
            // we never return here!
      }
   }

   return ret;
}

// ----------------------------------------------------------------------------
//
// HorizontalList
//
// ----------------------------------------------------------------------------

template <class T>
HorizontalList <T> :: HorizontalList(
                         T* objectIn, 
                         HorizontalList <T> * (T :: * getHIn)()
                      ) {
   object = objectIn;
   getH = getHIn;
   hor = NULL;
}

// ----------------------------------------------------------------------------

template <class T>
HorizontalList <T> :: ~HorizontalList() {
}

// ----------------------------------------------------------------------------

template <class T>
void HorizontalList <T> :: Reset() {

   hor = NULL;
}

// ----------------------------------------------------------------------------

template <class T>
void HorizontalList <T> :: InsertAfterFirstElement(HorizontalList <T> * link) {

   link -> SetHorizontalPtr(GetHorizontalPtr());
   SetHorizontalPtr(link);
}

// ----------------------------------------------------------------------------
//
// Cache container
//
// ----------------------------------------------------------------------------

template <class T>
CacheContainer <T> :: CacheContainer() {
}

// ----------------------------------------------------------------------------

template <class T>
CacheContainer <T> :: ~CacheContainer() {
   Delete();
}

// ----------------------------------------------------------------------------

template <class T>
void CacheContainer <T> :: Define(int nIntIn, int nDoubleIn) {

    if (IsDefined())
       errf(-1, "CacheContainer <T> :: Define: already defined");

    nInt       = nIntIn;
    nDouble    = nDoubleIn;

    intKey    = new Array_int(nInt);
    doubleKey = new Array_double(nDouble);

    intKeyData = intKey -> GetDataPtr();
    doubleKeyData = doubleKey -> GetDataPtr();

    Reset();
}

// ----------------------------------------------------------------------------

template <class T>
void CacheContainer <T> :: Delete() {

   delete intKey;
   delete doubleKey;
}

// ----------------------------------------------------------------------------

template <class T>
void CacheContainer <T> :: Reset() {

   SetOccupied(FALSE);
   SetReused(0);
}

// ----------------------------------------------------------------------------

template <class T>
void CacheContainer <T> :: CopyKeyInto(CacheContainer <T> * into) {

   register int i;

   for (i = 0; i < nInt; ++i)
       into -> intKeyData[i] = intKeyData[i];

   for (i = 0; i < nDouble; ++i)
       into -> doubleKeyData[i] = doubleKeyData[i];
}

// ----------------------------------------------------------------------------

template <class T>
int CacheContainer <T> :: CompareHashValues(
                             double hashIn, 
                             int doubleCrit, 
                             double doubleRes
                          ) {

   return compareDouble(hashIn, hashValue, doubleRes, doubleCrit);

#if 0
   int equal = TRUE; 

   switch (doubleCrit) {
      case 0:
         if (  fabs(hashIn - hashValue)
             > fabs(hashValue * doubleRes))
            equal = FALSE;
         break;
      default:
        errf(-1, "CacheContainer <T> :: CompareHashValues: crit. not known");
   }

   return equal;
#endif
}

// ----------------------------------------------------------------------------

template <class T>
int CacheContainer <T> :: CompareKeys(int*    tIntKey, 
                                      double* tDoubleKey, 
                                      int doubleCrit, 
                                      double doubleRes) {

   register int i;
   register int equal = TRUE;

   int* intKeyDataPtr = intKeyData;
   int* tIntKeyPtr    = tIntKey;
   for (i = 0; i < nInt && equal; ++i) {
       if (*(intKeyDataPtr++) != *(tIntKeyPtr++))
          equal = FALSE;
   }

   double* doubleKeyDataPtr = doubleKeyData;
   double* tDoubleKeyPtr    = tDoubleKey;
   for (i = 0; i < nDouble && equal; ++i) {
       equal = equal && compareDouble(
                           *tDoubleKeyPtr, 
                           *doubleKeyDataPtr,
                           doubleRes, 
                           doubleCrit
                        ); 
       ++tDoubleKeyPtr;
#if 0
       switch (doubleCrit) {
          case 0:
             if (  fabs(*(tDoubleKeyPtr++) - *doubleKeyDataPtr) 
                 > fabs(*doubleKeyDataPtr * doubleRes))
                equal = FALSE;
             break;
          default:
            errf(-1, "CacheContainer <T> :: CompareKeys: crit. not known");
       }
#endif
       ++doubleKeyDataPtr;
   }

   return equal;
}

// ----------------------------------------------------------------------------
//
// Cache
//
// ----------------------------------------------------------------------------

template <class T>
Cache <T> :: Cache(
                Global* globalIn, 
                int nCacheSizeIn, 
                int nIntIn, 
                int nDoubleIn
             ) {

   Set_global(globalIn);
   Define(nCacheSizeIn, nIntIn, nDoubleIn);
}

// ----------------------------------------------------------------------------

template <class T>
Cache <T> :: ~Cache() {

//   global -> To_status() << "Cache <T> :: ~Cache\n";
   Delete();
}

// ----------------------------------------------------------------------------

template <class T>
void Cache <T> :: Define(int nCacheSizeIn, int nIntIn, int nDoubleIn) {

   if (IsDefined())
      errf(-1, "Cache <T> :: Define: already defined");

   nCacheSize = nCacheSizeIn;
   nInt       = nIntIn;
   nDouble    = nDoubleIn;

   logBase = 10.;
   logBaseLog = log(logBase);
   oneOverLogBaseLog = 1. / logBaseLog;
   integerFactor = nCacheSize * 0.1;

   critHash = 0;
   accuracyHash = 1.e-8;

   critDouble = 0;
   accuracyDouble = 1.e-8;

   global -> To_status() << FS("defining Cache: %d ", nCacheSize)
                         << FS("%d ", nInt)
                         << FS("%d\n", nDouble);

   register int i;
   store = new CacheContainer <T> [nCacheSize];
   for (i = 0; i < nCacheSize; ++i)
       (store+i) -> Define(nInt, nDouble);

   freeList = NULL;

   Reset();
}

// ----------------------------------------------------------------------------

// :_CAUTION_: have to take care of data still attached to container
//             elements --> return mechanism!

template <class T>
void Cache <T> :: Delete() {

   global -> To_status() << FS("Cache Delete: %s\n", GetNameSafe());

   CacheContainer <T> * object;
   int hash;

   // statistics and destruct
   for (hash = 0; hash < nCacheSize; ++hash) {
        gl -> To_status().Flush();
        object = store + hash;
        if (object -> GetOccupied()) {
           StatisticsForAnObject(object);
           if (freeList != NULL) {
              freeList -> ReturnObject(object -> GetDataPtr());
           }
           else {
              delete object -> GetDataPtr();
           }
        }
   }

   global -> To_status() 
      << FS("# of different stored objects:                 %16ld\n", 
            numberOfObjects);
   if (numberOfObjects > 0) {
      global -> To_status() 
         << FS("average survival time:                         %16.6e\n",
            survivalTime / (double) numberOfObjects);
      global -> To_status() 
         << FS("average number of reuses:                      %16.6e\n",
               reused / (double) numberOfObjects);
   }
   global -> To_status() 
      << FS("# of accesses:                                 %16ld\n", 
            numberOfAccesses);
   if (numberOfAccesses > 0) {
      global -> To_status() 
         << FS("fraction of unsuccessful accesses:            %16.6e\n",
               unsuccessful / (double) numberOfAccesses);
   }

   delete [] store;  // :_CAUTION_: should take care of the T* part, too!
}

// ----------------------------------------------------------------------------

template <class T>
void Cache <T> :: Reset() {

   clock = 0;
   survivalTime = 0;
   reused = 0;
   numberOfObjects = 0;
   unsuccessful = 0;
   numberOfAccesses = 0;
}

// ----------------------------------------------------------------------------

template <class T>
void Cache <T> :: StatisticsForAnObject(CacheContainer <T> * object) {

   // statistics
   numberOfObjects++;
   survivalTime += (clock - object -> GetTime());
   reused += object -> GetReused();
//   object -> GetDataPtr() -> Print(gl -> To_status());
}

// ----------------------------------------------------------------------------

template <class T>
void Cache <T> :: StoreInCache(Array_int* tIntArray, 
                               Array_double* tDoubleArray, 
                               T* data) {

   ++clock;

   // get the hash value
   double fullHash;
   int hash = GetHashValue(tIntArray -> GetDataPtr(),
                           tDoubleArray -> GetDataPtr(), 
                           &fullHash);

   CacheContainer <T> * object = store + hash;

   if (object -> GetOccupied()) {

      StatisticsForAnObject(object);
      if (freeList != NULL)
         freeList -> ReturnObject(object -> GetDataPtr());
      else
         delete object -> GetDataPtr();

   }

   object -> SetDataPtr(data);
   tIntArray    -> CopyInto(object -> GetIntKeyPtr());
   tDoubleArray -> CopyInto(object -> GetDoubleKeyPtr());
   object -> SetHashValue(fullHash);
   object -> SetTime(clock);
   object -> SetOccupied(TRUE);
   object -> SetReused(0);
   object -> IncrementFreq();
}

// ----------------------------------------------------------------------------

template <class T>
T* Cache <T> :: FindInCache(Array_int* tIntArray, 
                            Array_double* tDoubleArray) {

   ++ numberOfAccesses;

   T* ret;

   double fullHash;
   int hash = GetHashValue(tIntArray -> GetDataPtr(),
                           tDoubleArray -> GetDataPtr(), 
                           &fullHash);

   
   CacheContainer <T> * object = store + hash;
   
   if (   object -> GetOccupied()
       && object -> CompareHashValues(fullHash, critHash, accuracyHash)
       && object -> CompareKeys(tIntArray -> GetDataPtr(), 
                                tDoubleArray -> GetDataPtr(),
                                critDouble, 
                                accuracyDouble)
      ) {
      object -> IncrementReused();
      ret = object -> GetDataPtr();
   } else {
      ++unsuccessful;
      ret = NULL;
   }

   return ret;
}

// ----------------------------------------------------------------------------

template <class T>
int Cache <T> :: GetHashValue(int*    tIntKey, 
                              double* tDoubleKey, 
                              double* fullHash) {

   int hash;

   register int i;

   double prod = 1.;
   for (i = 0; i < nDouble; ++i) 
       prod *= tDoubleKey[i];
   prod = fabs(prod);

   int iSum = 0;
   for (i = 0; i < nInt; ++i) 
       iSum += tIntKey[i];
   iSum = abs(iSum);  // to guarantee that hash will be positive

   double dHash; 
   dHash = log(prod) * oneOverLogBaseLog;
   dHash = dHash - floor(dHash);
   dHash *= nCacheSize;

   double iHash;
   iHash = iSum * integerFactor;

   double fh = dHash + iHash;
   *fullHash = fh; 

   hash = labs(rclv(floor(fh)));
   hash %= nCacheSize;
   
   return hash;
}

// ============================================================================
  
#endif // of CONTAINER_T
  
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
// ============================================================================
//
// --> REST
//
// file:              rest.cc
// created:           09.12.1996
// last modification: 15.12.1997
//
// ============================================================================

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <errno.h>

#include "switch.h"
#include "global.h"
#include "mth.h"
#include "strng.h"
#include "user.h"
#include "mcintg.h"
#include "cmb.h"
#include "qcd.h"
#include "disaster.h"

#include "rest.h"

// ============================================================================
//
// definitions
//
// ============================================================================

// PSI workstation specific?
#define PSW_MACHINE 0

// ============================================================================
//
// enumerations
//
// ============================================================================

// ============================================================================
//
// global variables
//
// ============================================================================

#if PSW_MACHINE
    extern unsigned long __noshrink;
    extern size_t __minshrink;
    extern double __minshrinkfactor;
    extern size_t __mingrow;
    extern double __mingrowfactor;
    extern unsigned long __madvisor;
    extern unsigned long __small_buff;
    extern int __fast_free_max;
    extern unsigned long __sbrk_override;
    extern unsigned long __taso_mode;
#endif
    
// ============================================================================
//
// --> classes
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> Test MatrixElements
//
// ----------------------------------------------------------------------------

TestMatrixElement::TestMatrixElement(int mel_indexIn) {

   const char tname[] = "TestMatrixElement ";

   mel_index=mel_indexIn;

   gl -> To_status() << FS("TestMatrixElement: constructor. %d\n", mel_index);

   char number[20];
   sprintf(number,"%8d",mel_index);

   name->Define(strlen(tname)+strlen(number));
   name->CopyFromChar(tname);
   name->AppendFromChar(number);

   Set_nIntVar(3*GetOut()-4);
}

// ----------------------------------------------------------------------------

TestMatrixElement::~TestMatrixElement() {
   gl -> To_status() << "TestMatrixElement: destructor called.\n";
}

// ----------------------------------------------------------------------------

int TestMatrixElement::GetIn() {
   int ret;
   switch(mel_index) {
      case 9999:
         ret = 0;
         break;
      default:
         ret = (abs(mel_index)/100)%10;
   }
   return ret;
}

// ----------------------------------------------------------------------------

int TestMatrixElement::GetOut() {
   int ret;
   switch(mel_index) {
      case 9999:
         ret = 0;
         break;
      default:
         ret = (abs(mel_index)/1000);
   }
   return ret;
}

// ----------------------------------------------------------------------------

// assumes that invariants have already been calculated in ``e''

void TestMatrixElement::Evaluate(
                            Event *e,
                            Contribution *con,
                            double factor
                          ) {

   int ifix = 0;
   int jfix = 0; 
   LimitType limittype = e -> GetLimitType();

   if (limittype != NO_LIMIT) {

      ifix = e->GetLimit1();
      jfix = e->GetLimit2();

      if (ifix == jfix || ifix < 0) {
         errf(-1,"TestMatrixElement::Evaluate: ifix==jfix");
      }
   }

   Event *limit;

   if (ifix<0 && (limittype==SOFT || limittype==SOFT_AND_COLLINEAR))
      errf(-1,"TestMatrixElement::Evaluate: soft<0?");

   // determine index of matrix element
   int me_index;
   switch(mel_index) {
      case 9999:
         me_index=0;
         break;
      default:
         me_index=lsign(mel_index)*(abs(mel_index)%100);
   }
   
   double ecm;
   ecm=e->GetEcm();

   // kinematical variables
   CopyPartonInvariantsFrom(e);

//   Invariant::Print(stdout);

   double yy,vv;
   yy=e0/(0.5*ecm);
   vv=v01;
      
   // order the indices (for the collinear limit only)
   int io,jo;
   if (limittype==COLLINEAR) {
      if (ifix<jfix) {
         io=ifix;
         jo=jfix;
      } else {
         io=jfix;
         jo=ifix;
      }
   } else {
      io=ifix;
      jo=jfix;
   }

   double term;
   term=0;

   int i0, i1, i2, i3, i2p, i3p;

   switch (me_index) {

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case 1000:

         switch(limittype) {
            case NO_LIMIT:
               term=0.;
               break;
            case SOFT:
               switch (ifix) {
                  case 0:
                     term=0.;
                     break;
                  case 1:
                     term=0.;
                     break;
                  case 2:
                     term=0.;
                     break;
               }
               break;
            case COLLINEAR:
               switch (io) {
                  case -1:
                     switch (jo) {
                        case 0:
                           term=0.;
                           break;
                        case 1:
                           term=0.;
                           break;
                        case 2:
                           term=0.;
                           break;
                     }
                     break;
                  case 0:
                     switch (jo) {
                        case 1:
                           term=0.;
                           break;
                        case 2:
                           term=0.;
                           break;
                     }
                     break;
                  case 1:
                     switch (jo) {
                        case 2:
                           term=0.;
                           break;
                     }
                     break;
               }
               break;
            case SOFT_AND_COLLINEAR:
               switch (io) {
                  case 0:
                     switch (jo) {
                        case -1:
                           term=0.;
                           break;
                        case 1:
                           term=0.;
                           break;
                        case 2:
                           term=0.;
                           break;
                     }
                     break;
                  case 1:
                     switch (jo) {
                        case -1:
                           term=0.;
                           break;
                        case 0:
                           term=0.;
                           break;
                        case 2:
                           term=0.;
                           break;
                     }
                     break;
                  case 2:
                     switch (jo) {
                        case -1:
                           term=0.;
                           break;
                        case 0:
                           term=0.;
                           break;
                        case 1:
                           term=0.;
                           break;
                     }
                     break;
               }
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1.;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case 5:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(ei*e0*f1*f2*vi1*vi2*v01*v02);
               break;
            case SOFT:
               switch (ifix) {
                  case 0:
                     term=0.;
                     break;
                  case 1:
                     term=1./(ei*e0*f2*vi1*vi2*v01*v02);
                     break;
                  case 2:
                     term=1./(ei*e0*f1*vi1*vi2*v01*v02);
                     break;
               }
               break;
            case COLLINEAR:
               switch (io) {
                  case -1:
                     switch (jo) {
                        case 0:
                           term=0.;
                           break;
                        case 1:
                           term=1./(ei*e0*f1*f2*vi2*v01*v02);
                           break;
                        case 2:
                           term=1./(ei*e0*f1*f2*vi1*v01*v02);
                           break;
                     }
                     break;
                  case 0:
                     switch (jo) {
                        case 1:
                           term=1./(ei*e0*f1*f2*vi1*vi2*v02);
                           break;
                        case 2:
                           term=1./(ei*e0*f1*f2*vi1*vi2*v01);
                           break;
                     }
                     break;
                  case 1:
                     switch (jo) {
                        case 2:
                           term=0.;
                           break;
                     }
                     break;
               }
               break;
            case SOFT_AND_COLLINEAR:
               switch (io) {
                  case 0:
                     switch (jo) {
                        case -1:
                           term=0.;
                           break;
                        case 1:
                           term=0.;
                           break;
                        case 2:
                           term=0.;
                           break;
                     }
                     break;
                  case 1:
                     switch (jo) {
                        case -1:
                           term=1./(ei*e0*f2*vi2*v01*v02);
                           break;
                        case 0:
                           term=1./(ei*e0*f2*vi1*vi2*v02);
                           break;
                        case 2:
                           term=0.;
                           break;
                     }
                     break;
                  case 2:
                     switch (jo) {
                        case -1:
                           term=1./(ei*e0*f1*vi1*v01*v02);
                           break;
                        case 0:
                           term=1./(ei*e0*f1*vi1*vi2*v01);;
                           break;
                        case 1:
                           term=0.;
                           break;
                     }
                     break;
               }
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,6);
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case 4:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(e0*e1*f2*v02*v12);
               break;
            case SOFT:
               switch (ifix) {
                  case 2:
                     term=1./(e0*e1*v02*v12);
                     break;
               }
               break;
            case COLLINEAR:
               switch (io) {
                  case 0:
                     switch (jo) {
                        case 2:
                           term=1./(e0*e1*f2*v12);
                           break;
                     }
                     break;
                  case 1:
                     switch (jo) {
                        case 2:
                           term=1./(e0*e1*f2*v02);
                           break;
                     }
                     break;
               }
               break;
            case SOFT_AND_COLLINEAR:
               switch (io) {
                  case 2:
                     switch (jo) {
                        case 0:
                           term=1./(e0*e1*v12);
                           break;
                        case 1:
                           term=1./(e0*e1*v02);
                           break;
                     }
                     break;
               }
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,4);
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case 0:   // matrix element is 0
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -1:

         switch(limittype) {
            case NO_LIMIT:
               term=1.;
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (io==0 && jo==1)
                  term=0.;
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -2:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(f0*sqrt(v01));
               break;
            case SOFT:
               if (ifix==0)
                  term=1./sqrt(v01);
               break;
            case COLLINEAR:
               if (io==0 && jo==1)
                  term=0;
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,2);
         switch (e->npartons) {
            case 2:
               term/=1.3445;
               break;
            case 3:
               term/=-6.1807;
               break;
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -3:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(sqrt(f0)*v01);
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (io==0 && jo==1)
                  term=1./(sqrt(f0));
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,1);
         switch (e->npartons) {
            case 2:
               term/=2.9999;
               break;
            case 3:
               term/=3.0;
               break;
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -4:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(f0*v01);
               break;
            case SOFT:
               if (ifix==0)
                  term=1./v01;
               break;
            case COLLINEAR:
               if (io==0 && jo==1)
                  term=1./f0;
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=1.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,2);
         switch (e->npartons) {
            case 2:
               term/=11.739;
               break;
            case 3:
               term/=13.1592;
               break;
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -5:

         switch(limittype) {
            case NO_LIMIT:
               term=1./sqrt(v01);
               break;
            case SOFT:
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         term/=1.4442;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -6:

         switch(limittype) {
            case NO_LIMIT:
               term=4.*(-1+2*vv-yy*vv*vv)/yy/(1-yy)/sqrt(vv);
               break;
            case SOFT:
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         term/=-6.1807;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -7:

         switch(limittype) {
            case NO_LIMIT:
               term=(-1+2*vv-yy*vv*vv)/yy/sqrt(vv);
               break;
            case SOFT:
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         term/=-0.72686;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -8:

         switch(limittype) {
            case NO_LIMIT:
               term=(-1+2*vv-yy*vv*vv)/sqrt(vv);
               break;
            case SOFT:
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         term/=-0.2546;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -9:

         switch(limittype) {
            case NO_LIMIT:
               term=(2*vv-yy*vv*vv)/sqrt(vv);
               break;
            case SOFT:
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         term/=1.1816;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -10:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(f0*sqrt(vi0));
               break;
            case SOFT:
               if (ifix==0)
                  term=1./sqrt(vi0);
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,2);
         switch (e->npartons) {
            case 2:
               term/=4.90965;
               break;
            case 3:
               term/=-6.1807;
               break;
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -11:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(sqrt(f0)*vi0);
               break;
            case SOFT:
               break;
            case COLLINEAR:
               if (io==-1 && jo==0)
                  term=1./(sqrt(f0));
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,1);
         switch (e->npartons) {
            case 2:
               term/=2.0;
               break;
            case 3:
               term/=3.0;
               break;
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -12:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(f0*vi0);
               break;
            case SOFT:
               if (ifix==0)
                  term=1./vi0;
               break;
            case COLLINEAR:
               if (io==-1 && jo==0)
                  term=1./f0;
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==-1)
                  term=1.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,2);
         switch (e->npartons) {
            case 2:
               term/=6.580;
               break;
            case 3:
               term/=13.1592;
               break;
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -13:

         switch(limittype) {
            case NO_LIMIT:
               term=1./sqrt(vi0);
               break;
            case SOFT:
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         term/=2.0;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -14:

         switch(limittype) {
            case NO_LIMIT:
               term=pow((vi0*v12-v02*vi1)/sqrt(v01),2)/v01;
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (io==0 && jo==1)
                  term=pow(e->wlist[-1][2],2);
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -15:

         switch(limittype) {
            case NO_LIMIT:
               limit = new Event;
               double jrel,en,va;
               e->Limit1(COLLINEAR,FS_REFERENCE,limit,&jrel,&en,&va);
//               limit->CalculateWList();
               limit->CalculatePartonInvariantsLocal();
               double a, b;
               a = pow((vi0*v12-v02*vi1)/sqrt(v01),2)/v01;
               b = jrel*pow(limit->wlist[-1][2],2)/va;
               term=a-b;
               delete limit;
               break;
            case SOFT:
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -16:

         switch(limittype) {
            case NO_LIMIT:
               term=pow((vi0*v12-v02*vi1)/sqrt(v01),2)/v01;
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=pow(e->wlist[-1][2],2);
               if (ifix==1 && jfix==0)
                  term=pow(e->wlist[-1][2],2);
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -17:

         switch(limittype) {
            case NO_LIMIT:
               switch (global->itest1) {
                  case 0:
                     term=pow(e->mapping[global->itest3]->GetCosPhi(),
                              global->itest2);
                     break;
                  case 1:
                     term=pow(e->mapping[global->itest3]->GetSinPhi(),
                              global->itest2);
                     break;
               }
               if (global->itest2 % 2 == 0)
                  term/=(8.*Pi/TwoPi*pow(2,-global->itest2)/global->itest2
                        *dfactorial(global->itest2-1)
                        /pow(dfactorial(global->itest2/2-1),2));
//               term=1.;
               break;
            case SOFT:
            case COLLINEAR:
            case SOFT_AND_COLLINEAR:
               break;      
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -18:

         switch(limittype) {
            case NO_LIMIT:
               double vl,afac,eps1;
               int iref;
               switch (global->itest3) {
                  case 0:
                     vl=V(e->mapping[0], e->mapping[1]);
                     iref=0;
                     eps1=2.*e->mapping[iref]->GetEnergy()/e->W;
                     afac=log(e->W/(e->W-2.*e->mapping[iref]->GetEnergy()))
                         /2./e->mapping[iref]->GetEnergy()
                         *(e->W-2.*e->mapping[iref]->GetEnergy()*vl)
                         /(1.-eps1)*(1.-eps1*(2.-eps1)*vl);
                     break;
                  case 1:
                     vl=V(e->mapping[-1], e->mapping[0]);
                     iref=0;
                     afac=log(e->W/(e->W-2.*e->mapping[iref]->GetEnergy()))
                         /2./e->mapping[iref]->GetEnergy()
                         *(e->W-2.*e->mapping[iref]->GetEnergy()*vl);
                     break;
                  case 2:
                     vl=V(e->mapping[-1], e->mapping[1]);
                     iref=1;
                     afac=log(e->W/(e->W-2.*e->mapping[iref]->GetEnergy()))
                         /2./e->mapping[iref]->GetEnergy()
                         *(e->W-2.*e->mapping[iref]->GetEnergy()*vl);
                     break;
                  default:
                     // TAWM
                     afac=0.;
                     vl=0.;
               }
               term=pow(vl,global->itest2)
                   *(global->itest2+1)
                   *afac;
               break;
            case SOFT:
            case COLLINEAR:
            case SOFT_AND_COLLINEAR:
               break;      
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -19:

         switch(limittype) {
            case NO_LIMIT:
               term=1./vi0;
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (io==-1 && jo==0)
                  term=1.;
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term/=3./4.;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -20:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(2.*e0/e->W)/sqrt(vi0);
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (io==-1 && jo==0)
                  term=0.;
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term/=2.7726;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -21:

         switch(limittype) {
            case NO_LIMIT:
               limit = new Event;
               double jrel,en,va;
               e->Limit1(COLLINEAR,IS_REFERENCE,limit,&jrel,&en,&va);
               limit->CalculatePartonInvariantsLocal();
//               limit->CalculateWList();
               double a, b;
               a = pow((vi1*v02-vi2*v01)/sqrt(vi0),2)/vi0;
               b = jrel*pow(limit->wlist[1][2],2)/va;
               term=a-b;
               gl -> To_status()
                  << FS("ydiff=%16.6e\n",pow((vi1*v02-vi2*v01)/sqrt(vi0),2));
               gl -> To_status()
                  << FS("jrel= %16.6e\n",jrel);
               gl -> To_status()
                  << FS("a, b, a-b= %16.6e", a)
                  << FS(" %16.6e", b)
                  << FS(" %16.6e\n", term);
               gl -> To_status()
                  << FS("vi0, va = %16.6e", vi0)
                  << FS(" %16.6e\n", va);
               gl -> To_status()
                  << FS("NN01=%16.6e", v01)
                  << FS(" L %16.6e\n", limit->v01);
               gl -> To_status()
                  << FS("NN02=%16.6e", v02)
                  << FS(" L %16.6e\n", limit->v02);
               gl -> To_status()
                  << FS("NNi1=%16.6e", vi1)
                  << FS(" L %16.6e\n", limit->vi1);
               gl -> To_status()
                  << FS("NNi2=%16.6e", vi2)
                  << FS(" L %16.6e\n", limit->vi2);
               delete limit;
               break;
            case SOFT:
               break;
            case COLLINEAR:
               break;
            case SOFT_AND_COLLINEAR:
               break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -22:

         switch(limittype) {
            case NO_LIMIT:
               term=pow((vi0*v12-vi2*v01)/sqrt(vi1),2)/vi1;
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (ifix==1 && jfix==-1)
                  term=pow(e->wlist[0][1],2);
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -23:

         switch(limittype) {
            case NO_LIMIT:
               term=pow((si0*s12-s02*si1),2)/s01/s01;
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (io==0 && jo==1) 
                  term=1./eh/eh*pow(sih * dt2 - s2h * dti, 2);
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -24:

         i0  =  global->i0;
         i1  =  global->i1;
         i2  =  global->i2;
         i3  =  global->i3;
         i2p =  global->i2p;
         i3p =  global->i3p;

         switch(limittype) {
            case NO_LIMIT:
               term=(e->DotMap(i0,i2)*e->DotMap(i1,i3)
                    -e->DotMap(i0,i3)*e->DotMap(i1,i2))
                   *(e->DotMap(i0,i2p)*e->DotMap(i1,i3p)
                    -e->DotMap(i0,i3p)*e->DotMap(i1,i2p))
                   /pow(e->DotMap(i0,i1),2);
               break;
            case SOFT:
               if (ifix==0)
                  term=0.;
               break;
            case COLLINEAR:
               if (io==i0 && jo==i1 || io==i1 && jo==i0) 
                  term=1./eh/eh
                      *(e->DotMap(-3,i2) * e->dtlist[i3]
                       -e->DotMap(-3,i3) * e->dtlist[i2])
                      *(e->DotMap(-3,i2p) * e->dtlist[i3p]
                       -e->DotMap(-3,i3p) * e->dtlist[i2p]);
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=0.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=1;
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -25:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(f0*v01)+1./(e0*sqrt(v01));
               break;
            case SOFT:
               if (ifix==0)
                  term=1./v01;
               break;
            case COLLINEAR:
               if (io==0 && jo==1)
                  term=1./f0;
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==0 && jfix==1)
                  term=1.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,2);
         switch (e->npartons) {
            case 2:
               term/=11.739;
               break;
            case 3:
               term/=13.1592;
               break;
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      case -40:

         switch(limittype) {
            case NO_LIMIT:
               term=1./(f1*v01);
               break;
            case SOFT:
               if (ifix==1)
                  term=1./v01;
               break;
            case COLLINEAR:
               if (io==0 && jo==1)
                  term=1./f1;
               break;
            case SOFT_AND_COLLINEAR:
               if (ifix==1 && jfix==0)
                  term=1.;
                break;
            default:
               errf(-1, "limittype not available");
         }
         term*=pow(ecm,2);
         switch (e->npartons) {
            case 2:
               term/=11.739;
               break;
            case 3:
               term/=13.1592;
               break;
         }
         break;

// _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

      default:
         global -> To_error()
            << FS("%d\n", me_index);
         errf(-1,"TestMatrixElement::Evaluate: matrix element not known");
   }

   // normalize by the phase space volume

   int nfinal;
   nfinal=GetOut();

   double ecmfull;
   ecmfull=e->GetEcmFull();

   double volume;
   if (GetIn()==1 && nfinal==1) {
      volume=1./8./Pi; // both for fixed and integrated lepton
   } else {
      volume=pow(TwoPi,(double)4-3*nfinal)*pow(Pi/2.,(double)nfinal-1)
            /dfactorial(nfinal-1)/dfactorial(nfinal-2)
            *pow(ecmfull,(double)2*nfinal-4);

      // modify if we are dealing with DIS
      if (GetIn()==1) {
         volume*=1./16./Pi/Pi*e->SH;
         if (!global->disaster->Get_xi_integration()) {
            if (!global->disaster->Get_lepton_integration()) {
               volume*=pow(e->y,nfinal-1)*pow(e->xi-e->xB,nfinal-2);
            } else {
               volume*=1./((double)(nfinal-1)*nfinal);
            }
         } else {
            if (!global->disaster->Get_lepton_integration()) {
               volume*=pow(e->y,nfinal-1)*pow(1.-e->xB,nfinal-1)
                          /((double)(nfinal-1));
            } else {
               volume*=1./((double)nfinal*nfinal*(nfinal-1));
            }
         }
      }
   }

   term/=volume;
   term *= factor;

   con -> GetData() [0] += term * global -> VEGAS_weight_current;
}

// ----------------------------------------------------------------------------
//
// --> Test process _1
//
//     integrates a given matrix element (mel_index)
//
// ----------------------------------------------------------------------------

TestProcess_1::TestProcess_1() {
   gl -> To_status() << "TestProcess_1: constructor.\n";
}

// ----------------------------------------------------------------------------

TestProcess_1::~TestProcess_1() {
   gl -> To_status() << "TestProcess_1: destructor called.\n";
}

// ----------------------------------------------------------------------------

void TestProcess_1::Define(Global *globalIn) {

   global=globalIn;
}

// ----------------------------------------------------------------------------

void TestProcess_1::AssignMatrixElements() {

   // delete old definition, if any
   Delete();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 1;

   // *** fill in number of components here
   CreateNME(1);

   // *** fill in number of matrix elements per component here
   n_me             [0] = 1;
   to_be_subtracted [0] = Subtraction;

   CreateOtherArrays();

   // *** fill in matrix elements here
   me_list [0][0] = new TestMatrixElement(global->mel_index); 

   AssignInformationToMatrixElements();

   // frame to be used for phase space
   if (GetIn() >= 1 && GetMaxOut() <= 2)
      frame=hCMS;
   else
      frame=pCMS;
}

// ----------------------------------------------------------------------------
//
// --> Test process _2
//
//     integrates a given matrix element, but adds a constant of 1
//     (required for cases when the original integral is 0, for 
//      MC stability reasons)
//
// ----------------------------------------------------------------------------

TestProcess_2::TestProcess_2() {
   gl -> To_status() << "TestProcess_2: constructor.\n";
}

// ----------------------------------------------------------------------------

TestProcess_2::~TestProcess_2() {
   gl -> To_status() << "TestProcess_2: destructor called.\n";
}

// ----------------------------------------------------------------------------

void TestProcess_2::Define(Global *globalIn) {

   global=globalIn;
}

// ----------------------------------------------------------------------------

void TestProcess_2::AssignMatrixElements() {

   // delete old definition, if any
   Delete();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 1;

   // *** fill in number of components here
   switch (global->mel_index) {
      case -3110:
      case -3111:
      case -3112:
         CreateNME(2);
         break;
      case -13110:
      case -13111:
      case -13112:
         CreateNME(1);
         break;
      default:
         errf(-1,"TestProcess_2::AssignMatrixElements: ME not known");
   }

   // *** fill in number of matrix elements per component here
   switch (global->mel_index) {
      case -3110:
      case -3111:
      case -3112:
         n_me             [0] = 1;
         to_be_subtracted [0] = NoSubtraction;
         n_me             [1] = 1;
         to_be_subtracted [1] = Subtraction;
         break;
      case -13110:
      case -13111:
      case -13112:
         n_me             [0] = 2;
         to_be_subtracted [0] = Subtraction;
         break;
   }

   CreateOtherArrays();

   // *** fill in matrix elements here
   switch (global->mel_index) {
      case -3110:
      case -3111:
      case -3112:
         me_list          [0][0] = new TestMatrixElement(-3101); 
         me_list          [1][0] = new TestMatrixElement(global->mel_index); 
         break;
      case -13110:
      case -13111:
      case -13112:
         me_list          [0][0] = new TestMatrixElement(-3101); 
         me_list          [0][1] = new TestMatrixElement(
                                          sign(global->mel_index)
                                         *(abs(global->mel_index)-10000)
                                       ); 
         break;
   }

   AssignInformationToMatrixElements();

   // frame to be used for phase space
   if (GetIn() >= 1 && GetMaxOut() <= 2)
      frame=hCMS;
   else
      frame=pCMS;
}

// ----------------------------------------------------------------------------
//
// --> Test process _3
//
//     integrates a given specific matrix element unsubtracted
//
// ----------------------------------------------------------------------------

TestProcess_3::TestProcess_3() {
   gl -> To_status() << "TestProcess_3: constructor.\n";
}

// ----------------------------------------------------------------------------

TestProcess_3::~TestProcess_3() {
   gl -> To_status() << "TestProcess_3: destructor called.\n";
}

// ----------------------------------------------------------------------------

void TestProcess_3::Define(Global *globalIn) {

   global=globalIn;
}

// ----------------------------------------------------------------------------

void TestProcess_3::AssignMatrixElements() {

   // delete old definition, if any
   Delete();

   // set up matrix elements, etc.

   // (default) order of the strong coupling constant
   alpha_s_Order = 1;

   // *** fill in number of components here
   CreateNME(1);

   // *** fill in number of matrix elements per component here
   n_me             [0] = 1;
   to_be_subtracted [0] = NoSubtraction;

   CreateOtherArrays();

   // *** fill in matrix elements here
   me_list          [0][0] = new TestMatrixElement(global->mel_index); 

   AssignInformationToMatrixElements();

   // frame to be used for phase space
   if (GetIn() >= 1 && GetMaxOut() <= 2)
      frame=hCMS;
   else
      frame=pCMS;
}

// ============================================================================

// test routine

double testps(double *x,double *wgt) {

      x=x;
      wgt=wgt;
      return 0;
}

// ----------------------------------------------------------------------------
//
// --> test function
//
// ----------------------------------------------------------------------------

double testmel(Global *global, Event *event) {
   global=global;
//   event->Print(stdout);
//   event->Invariant::Print(stdout);
   double v01=event->v01;
   double v02=event->v02;
   double v12=event->v12;
   double vi0=event->vi0;
   double vi1=event->vi1;
//   double vi2=event->vi2;
   double a,b;
   a=(vi0*v12-v02*vi1)/sqrt(v01);
   b=event->wlist[-1][2];
   global -> To_status() 
      << FS("a=%16.6e", a)
      << FS("   wlist=%16.6e", b)
      << FS(", a/wlist=%16.6e\n", a/b);
   return (a*a-b*b)/v01;
}

// ----------------------------------------------------------------------------

double testmeli(Global *global, Event *event) {
   global=global;
//   event->Print(stdout);
//   event->Invariant::Print(stdout);
   double v01=event->v01;
   double v02=event->v02;
//   double v12=event->v12;
   double vi0=event->vi0;
   double vi1=event->vi1;
   double vi2=event->vi2;
   double a,b;
   a=(v01*vi2-v02*vi1)/sqrt(vi0);
   b=event->wlist[1][2];
   global -> To_status() 
      << FS("i a=%16.6e", a)
      << FS("   wlist=%16.6e", b)
      << FS(", a/wlist=%16.6e\n", a/b);
   return (a*a-b*b)/vi0;
}

// ----------------------------------------------------------------------------
//
// --> rescaled subtracted term
//
// ----------------------------------------------------------------------------

double rescaled(Global *global,Event *ein,
                LimitType ltype,double frac,
                double (*pmel)(Global*,Event*)
               ) {

   pmel=pmel;

//   gl->v_delta=1.0;
//   gl->v_delta1=1.0;

   Event event;
   event.Define(20,global);
   double u[100];

   ein->CopyInto(&event);

   event.MapPSToUnit(u);

   double x0_lim;
   x0_lim=pow(frac,1.0);   

   if (ltype==SOFT) {
      // approximative soft limit
      u[0]=u[0]*frac;
   } else if (ltype==COLLINEAR) {
      // approximative collinear limit
      u[1]=u[1]*frac;
   } else if (ltype==SOFT_AND_COLLINEAR) {
      // approximative soft and limit
      u[0]=u[0]*frac;
      u[1]=u[1]*frac;
   } else
      errf(-1,"ltype not known.");

   event.MapUnitToPS(u);

   if (event.GetFrame()==hCMS)
      event.EcmAddRemnantAndIncident();

//   event.Print(stdout);

   event.DetermineHLambdaEta();
   event.CalculatePartonInvariantsLocal();
//   event.CalculateWList();

   double term=0.;

   term=pmel(global,&event);

   // check for NaN
   int nan;
   nan=fpnan->DoubleStatusReset(term);
   if (nan) {
      global -> To_status() << FS("NaN term=%e\n", term);
   }

   return term;
}

// ----------------------------------------------------------------------------
//
// --> rescaled subtracted term, use sumoverdiff
//
// ----------------------------------------------------------------------------

double rescalednew(Global *global,Event *ein,
                LimitType ltype,double frac,
                int variant, Process *process, int ic,
                double jacobian
               ) {

   Event event;
   event.Define(ein->npartons,global);
   double u[100];

   ein->CopyInto(&event);

   event.MapPSToUnit(u);

   double x0_lim;
   x0_lim=pow(frac,1.0);   

   if (ltype==SOFT) {
      // approximative soft limit
      u[0]=u[0]*frac;
   } else if (ltype==COLLINEAR) {
      // approximative collinear limit
      u[1]=u[1]*frac;
   } else if (ltype==SOFT_AND_COLLINEAR) {
      // approximative soft and limit
      u[0]=u[0]*frac;
      u[1]=u[1]*frac;
   } else
      errf(-1,"ltype not known.");

   event.MapUnitToPS(u);

   if (event.GetFrame()==hCMS)
      event.EcmAddRemnantAndIncident();

//   event.Print(stdout);

   event.DetermineHLambdaEta();
   event.CalculatePartonInvariantsLocal();
//   event.CalculateWList();

   double term=0.;

   switch(variant) {

      case 0:
               
         event.DoSumOverSubtractionTerms(
                                   process,ic,
                                   jacobian,
                                   global->process->GetCA()
                                );

         break;

      case 1:

         event.DoSumOverDiffSubtractionTerms(
                                   process,ic,
                                   jacobian,
                                   global->process->GetCA()
                                );

         break;
   }

   // check for NaN
   int nan;
   nan=fpnan->DoubleStatusReset(term);
   if (nan) {
      global -> To_status() << FS("NaN term=%e\n", term);
   }

   term=-99.;
   return term;
}

// ----------------------------------------------------------------------------
//
// --> determine the exponent
//
// ----------------------------------------------------------------------------

#if 0
double get_exponent(Global *global,Event *ein,
                    LimitType ltype,
                    double (*pmel)(Global*,LimitType,int,int,Event*) 
                   ) {
   double rho0,rho0_old;
   double rhoa,rhob,rhomed;
   double rho_min=1.e-10;
   double rho_diffmin=1.e-5;
   double rho_reduc=0.5;
   double frac1=1.e-9;  // these two values should already be 
   double frac2=1.e-10; // in the scaling region
   double frac3=1.e-12; // final check
   double val1,val2,val3;
   int deca,decb,decmed;   
   
   rho0=1./rho_reduc;
   do {
      rho0_old=rho0;
      rho0=rho_reduc*rho0;
      val1=pow(frac1,1.0-rho0)*rescaled(global,ein,ltype,frac1,pmel);
      val2=pow(frac2,1.0-rho0)*rescaled(global,ein,ltype,frac2,pmel);
   } while (rho0>rho_min && (   fabs(val1)<fabs(val2)
                             || (val1>0 && !val2>0)
                             || (!val1>0 && val2>0)
                            )
           );

   if (rho0_old<2.0 && rho0>rho_min) {
      // divide interval; rhoa<=rhob
      rhoa=rho0;
      rhob=rho0_old;
      do {
         rhomed=0.5*(rhoa+rhob);
         val1=pow(frac1,1.0-rhoa)*rescaled(global,ein,ltype,frac1,pmel);
         val2=pow(frac2,1.0-rhoa)*rescaled(global,ein,ltype,frac2,pmel);
         deca=(fabs(val1)>fabs(val2));
         val1=pow(frac1,1.0-rhob)*rescaled(global,ein,ltype,frac1,pmel);
         val2=pow(frac2,1.0-rhob)*rescaled(global,ein,ltype,frac2,pmel);
         decb=(fabs(val1)>fabs(val2));
         val1=pow(frac1,1.0-rhomed)*rescaled(global,ein,ltype,frac1,pmel);
         val2=pow(frac2,1.0-rhomed)*rescaled(global,ein,ltype,frac2,pmel);
         decmed=(fabs(val1)>fabs(val2));
         if (decmed)
            rhoa=rhomed;
         else
            rhob=rhomed;
      } while (fabs(rhob-rhoa)>rho_diffmin);
      rho0=rhomed;
   }

   val1=pow(frac1,1.0-rho0)*rescaled(global,ein,ltype,frac1,pmel);
   val2=pow(frac2,1.0-rho0)*rescaled(global,ein,ltype,frac2,pmel);
   val3=pow(frac3,1.0-rho0)*rescaled(global,ein,ltype,frac3,pmel);

   gl -> To_status() 
      << FS("rho0=%13.6e", rho0)
      << FS(" val1=%13.6e", val1)
      << FS(" val2=%13.6e", val2)
      << FS(" val3=%13.6e\n", val3);

   if (sign(val1)!=sign(val2))
      gl -> To_status() << "DIFFERENT SIGNS! ILLEGAL!\n";

   return rho0;   
}
#endif

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890

#if 0
// ++++++++++++++++++++++++++ GARBAGE +++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++ END OF GARBAGE +++++++++++++++++++++++++++++++++++++++++
#endif
