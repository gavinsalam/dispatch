// ============================================================================
//
// --> include file for user routines
//
// file:              user.h
// created:           29.03.1997
// last modification: 02.12.1997
//
// ============================================================================

#ifndef USERH
#define USERH

// ============================================================================

#include "switch.h"
#include "container.h"
#include "qcd.h"

// ============================================================================
//
// --> classes 
//
// ============================================================================

class Disaster;
class Library;
class Book;
class Process;
class ContributionArray;
class Contribution;
class Event;
class PartonDensityParametrization;
class AlphaS;
class MCUserForDisaster;
class MCUserInterface;

// ----------------------------------------------------------------------------
//
// --> products of parton densities and charges (== flavour factors)
//
// ----------------------------------------------------------------------------
 
class PDFC
   : public DefinedObject
   {

private:

   Disaster* disaster;

   double ff[11];
   double alphaS;
   double alphaEM;

protected:

public:

   inline double* GetFlavourFactorPtr() { return ff; }
   inline double* GetAlphaSPtr()  { return &alphaS; }
   inline double* GetAlphaEMPtr() { return &alphaEM; }

   PDFC();
  ~PDFC();

   void Define(Disaster*);
   void Reset(Disaster*);
   void ReturnToFreeList();
   void Print(Mstream&);
   void Print(FILE*);
   int UseListIsSame(PDFC*);
   long GetIdent();
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
   typedef FreeList_1 <PDFC, Disaster*> PDFCFreeList;
  
// ----------------------------------------------------------------------------
//
// --> flavour factors
//
// ----------------------------------------------------------------------------
 
class FlavourFactors
   : public DefinedObject
   {

private:

   Disaster* disaster;

   PDFCFreeList* pdfcfl;

   Cache <PDFC> * cache;
   Array_int* intKeyArray;
   int* intKey;
   Array_double* doubleKeyArray;   
   double* doubleKey;

   PartonDensityParametrization* pdf_Server;
   AlphaS* alphaS_Server;

   AlphaEM* alphaEM_Server;
   
protected:

public:

   void SetPDFCFL(PDFCFreeList* in) { pdfcfl = in; }
   
   inline void Set_pdf_Server(PartonDensityParametrization* in)
                  { pdf_Server = in; }

   inline void Set_alphaS_Server(AlphaS* in)
                  { alphaS_Server = in; }

   inline void Set_alphaEM_Server(AlphaEM* in)
                  { alphaEM_Server = in; }

   FlavourFactors();
  ~FlavourFactors();

   void Define(Disaster*);
   void Delete();

   void CalculateFlavourFactors(
           PDFC* pdfc,
           double xi, 
           double muF,
           double muR,
           int nflavour_pden, 
           int nflavour_out,
           int pdf_collection,
           int pdf_parametrization,
           int pdf_set, 
           int alphaSVariant,
           int alphaSOrder, 
           double alphasLambdaQCD4,
           int alphaEMVariant,
           double alphaEMScale
        );

   PDFC* FindOrCreateAndFillAndCacheFlavourFactors(
            double xi, 
            double muF,
            double muR,
            int nflavour_pden, 
            int nflavour_out,
            int pdf_collection,
            int pdf_parametrization,
            int pdf_set,
            int alphaSVariant,
            int alphaSOrder, 
            double alphasLambdaQCD4,
            int alphaEMVariant,
            double alphaEMScale
         );

   void CreateCache(int);
   void DeleteCache();
};

// ----------------------------------------------------------------------------
//
// --> MC User routine
//
// ----------------------------------------------------------------------------
 
class MCUser 
   {

private:

   DefinedObject definedObject;

protected:

   Disaster* global;

public:

   MCUser();
   virtual ~MCUser();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void Start() = 0;
   virtual void End() = 0;
   virtual int  SetParameter(char* pname, char* parameter) = 0;
   virtual void StartIntegration1() = 0;
   virtual void StartIntegration2() = 0;
   virtual void EndIntegration(double average,double error) = 0;
   virtual void BeginEvent(Process* process) = 0;
   virtual void EndOfEvent(Process* process) = 0;
   virtual double EndEvent(
                     Process* process,
                     ContributionArray* ca,
                     int ifFinalRun
                  ) = 0;
   virtual void DropEvent(Process* process) = 0;
   virtual void AcceptEvent(Process* process) = 0;
   virtual void PhaseSpace(
                   Process* process,
                   Event* event,
                   int& flag,
                   Contribution* con, 
                   int a = 0)
                   { process = process;
                     event = event;
                     flag = flag;
                     con = con; 
                     a = a; 
                     errf(-1, "PhaseSpace from user.h ...");
                   }
   virtual void IncidentPhaseSpace(
                     Process* process,
                     Event* event,
                     int& flag) = 0;
};

// ----------------------------------------------------------------------------
//
// --> User Interface: called from Disaster
//
// ----------------------------------------------------------------------------
 
class MCUserForDisaster 
   : public MCUser
   {

private:

   DefinedObject definedObject;

protected:

   MCUserInterface* mcUserInterface;

   FlavourFactors* flavourFactors;

   PDFCFreeList* pdfcFree;

   PartonDensityFreeList* pdFree;
   AlphaS* alphaS_Server;
   PartonDensityParametrization* pdf_Server;

   AlphaEM* alphaEM_Server;

   int variant_alpha_s;
   int user_defined_alpha_s;
   int order_alpha_s;
   double lambda_4_alpha_s;
   int pdf_collection;
   int pdf_parametrization;
   int pdf_set;
   int variant_alpha_em;

   int alphaSVariant;
   int alphaSOrder;
   double alphasLambdaQCD4;

   int alphaEMVariant;

   int nOutFlavours;      // number of flavours in outgoing q qbar - pairs
   int nPdenFlavours;     // number of flavours taken from the parton densities

   int testCase;
   int n_events;

   // flags for observables
   int universalObservableChoice;

   // parameters for universal observable
   double uoEnergyPower;
   double uoAnglePower;

   // flags for components
   int cFlagIndex[100];
   int cFlagValue[100];
   int nFlag;

   double sum[100], sum_2[100], err[100], relative_err[100], contr[100];
   int nused;

   int observablesD;

   Library* globalLibrary;

   Book* adaptationObservable;
   Book** unitCubeHistogram;

   // scales and numbers
   double event_sum;
   double event_sum_q;
   double event_sum_g;
   double tsum[200];
   double tsumq[200];
   double tsumg[200];
   double event_ff[11];

   double QArray[200];
   int noutArray[200];
   int npdenArray[200];
   int nfSwitchArray[200];
   int nRSArray[10][200];
   int nFSArray[10][200];
   double rsArray[10][200];
   double fsArray[10][200];

   double univObsStore[200];

   Event* breit_event;

public:

   FlavourFactors* Get_flavourFactors() { return flavourFactors; }

   MCUserForDisaster();
   virtual ~MCUserForDisaster();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void Start();
   virtual void End();
   virtual int  SetParameter(char* pname, char* parameter);
   virtual void StartIntegration1();
   virtual void StartIntegration2();
   virtual void EndIntegration(double average,double error);
   virtual void BeginEvent(Process* process);
   virtual void EndOfEvent(Process* process);
   virtual double EndEvent(
                     Process* process,
                     ContributionArray* ca,
                     int ifFinalRun
                  );
   virtual void ConvoluteWithDistributions(        
                   Process* process,
                   ContributionArray* ca,
                   double*, 
                   double*, 
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
                );
   virtual void DropEvent(Process* process);
   virtual void AcceptEvent(Process* process);
   virtual double UniversalObservable(Event* e);
   virtual double CalculateAdaptationObservable(
                     Process* process,
                     ContributionArray* ca,
                     int pdf_collectionIn,   
                     int pdf_parametrizationIn,
                     int pdf_setIn,
                     int alphaSVariant,
                     int alphaSOrder,
                     double alphasLambdaQCD4,
                     int alphaEMVariant
                  );
   virtual void PhaseSpace(
                   Process* process,
                   Event* event,
                   int& flag,
                   Contribution* con, 
                   int a = 0);
   virtual void BoostToBreitFrame(
                   Process* process, 
                   Event* event,
                   ContributionArray* ca, 
                   int ifFinalRun);
   virtual void IncidentPhaseSpace(
                     Process* process,
                     Event* event,
                     int& flag);

   inline virtual void             Set_mcUserInterface(MCUserInterface* in)
                                      { mcUserInterface = in; }
   inline virtual MCUserInterface* Get_mcUserInterface()
                                      { return mcUserInterface; }
};

// ----------------------------------------------------------------------------
//
// --> methods that should be implemented by the user routine
//
// ----------------------------------------------------------------------------

class MCUserCalls {

private:

   Disaster* glob;
   DefinedObject definedObject;

protected:

public:

   MCUserCalls();
   virtual ~MCUserCalls();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void Start_User() = 0;
   virtual void End_User() = 0;
   virtual int SetParameter_User(char* pname, char* parameter) = 0;
   virtual void StartIntegration1_User() = 0;
   virtual void StartIntegration2_User() = 0;
   virtual void EndIntegration_User(double average,double error) = 0;
   virtual void BeginEvent_User(Process* process) = 0;
   virtual void EndOfEvent_User(Process* process) = 0;
   virtual double EndEvent_User(
                     Process* process,
                     ContributionArray* ca,
                     int ifFinalRun
                  ) = 0;
   virtual void DropEvent_User(Process* process) = 0;
   virtual void AcceptEvent_User(Process* process) = 0;
   virtual void IncidentPhaseSpace_User(
                     Process* process,
                     Event* event,
                     int& flag) = 0;
};

// ----------------------------------------------------------------------------
//
// --> 
//
// ----------------------------------------------------------------------------

class MCUserCallsOther {

private:

   Disaster* glob;
   DefinedObject definedObject;

protected:

public:

   MCUserCallsOther();
   virtual ~MCUserCallsOther();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void ConvoluteWithDistributions_User(        
                   Process* process,
                   ContributionArray* ca,
                   double*, 
                   double*, 
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
                ) = 0;
   virtual double CalculateAdaptationObservable_User(
                     Process* process,
                     ContributionArray* ca,
                     int pdf_collectionIn,   
                     int pdf_parametrizationIn,
                     int pdf_setIn,
                     int alphaSVariant,
                     int alphaSOrder,
                     double alphasLambdaQCD4,
                     int alphaEMVariant
                  ) = 0;
};

// ----------------------------------------------------------------------------
//
// --> MC User routine: called from DISASTER++
//
// ----------------------------------------------------------------------------
 
class MCUserCalled
   : public MCUserForDisaster, 
     public MCUserCalls
   {

private:

   DefinedObject definedObject;

protected:

public:

   MCUserCalled();
   virtual ~MCUserCalled();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void ConvoluteWithDistributions_User(        
                   Process* process,
                   ContributionArray* ca,
                   double*, 
                   double*, 
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
                );
   virtual double CalculateAdaptationObservable_User(
                     Process* process,
                     ContributionArray* ca,
                     int pdf_collectionIn,   
                     int pdf_parametrizationIn,
                     int pdf_setIn,
                     int alphaSVariant,
                     int alphaSOrder,
                     double alphasLambdaQCD4,
                     int alphaEMVariant
                  );
};

// ----------------------------------------------------------------------------
//
// --> an interface class to switch handlers
//
// ----------------------------------------------------------------------------
 
class MCUserInterface
   : public MCUser, 
     public MCUserCallsOther
   {

private:

   DefinedObject definedObject;

protected:

   MCUserCalled*      mcUserCalled;
   MCUserForDisaster* mcUserForDisaster;

public:

   MCUserInterface();
   virtual ~MCUserInterface();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void Start();
   virtual void End();
   virtual int SetParameter(char* pname, char* parameter);
   virtual void StartIntegration1();
   virtual void StartIntegration2();
   virtual void EndIntegration(double average, double error);
   virtual void BeginEvent(Process* process);
   virtual void EndOfEvent(Process* process);
   virtual double EndEvent(
                     Process* process,
                     ContributionArray* ca,
                     int ifFinalRun
                  );
   virtual void DropEvent(Process* process);
   virtual void AcceptEvent(Process* process);
   virtual void IncidentPhaseSpace(
                     Process* process,
                     Event* event,
                     int& flag);
   virtual void ConvoluteWithDistributions_User(        
                   Process* process,
                   ContributionArray* ca,
                   double*, 
                   double*, 
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
                );
   virtual double CalculateAdaptationObservable_User(
                     Process* process,
                     ContributionArray* ca,
                     int pdf_collectionIn,   
                     int pdf_parametrizationIn,
                     int pdf_setIn,
                     int alphaSVariant,
                     int alphaSOrder,
                     double alphasLambdaQCD4,
                     int alphaEMVariant
                  );

   inline virtual void Set_mcUserCalled(MCUserCalled* in)
                          { mcUserCalled = in; }
   inline virtual void Set_mcUserForDisaster(MCUserForDisaster* in)
                          { mcUserForDisaster = in; }
};

// ----------------------------------------------------------------------------
//
// --> an interface class to switch handlers for the FORTRAN interface
//
// ----------------------------------------------------------------------------
 
class MCUserInterface_FORTRAN
   : public MCUserInterface 
   {

private:

   DefinedObject definedObject;

   Counter changed_print_counter;

protected:

public:

   MCUserInterface_FORTRAN();
   virtual ~MCUserInterface_FORTRAN();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void Start();
   virtual void End();
   virtual int SetParameter(char* pname, char* parameter);
   virtual void StartIntegration1();
   virtual void StartIntegration2();
   virtual void EndIntegration(double average, double error);
   virtual void BeginEvent(Process* process);
   virtual void EndOfEvent(Process* process);
   virtual double EndEvent(
                     Process* process,
                     ContributionArray* ca,
                     int ifFinalRun
                  );
   virtual void DropEvent(Process* process);
   virtual void AcceptEvent(Process* process);
   virtual void IncidentPhaseSpace(
                     Process* process,
                     Event* event,
                     int& flag);
   virtual void ConvoluteWithDistributions_User(        
                   Process* process,
                   ContributionArray* ca,
                   double*, 
                   double*, 
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
                );
   virtual double CalculateAdaptationObservable_User(
                     Process* process,
                     ContributionArray* ca,
                     int pdf_collectionIn,   
                     int pdf_parametrizationIn,
                     int pdf_setIn,
                     int alphaSVariant,
                     int alphaSOrder,
                     double alphasLambdaQCD4,
                     int alphaEMVariant
                  );
   
   void Control_To_FORTRAN();
   void Control_To_C();
};

// ----------------------------------------------------------------------------
//
// --> my own user class
//
// ----------------------------------------------------------------------------
 
class MyUserClass
      : public MCUserCalled 
      {

private:

   DefinedObject definedObject;

   int cluster_flag, min_combine, max_combine;
   ClusterAlgorithm cluster_algorithm;
   int cluster_scale_choice;
   double cluster_cut;
   int ren_scale_choice;
   int fact_scale_choice;
   int decisionMode;
   int minClusterNumber;
   int maxClusterNumber;

   ClusterAlgorithmType cat;
   RecombinationType rt;

   double ratioProtonElectronInLaboratory;

   int testObservables;
   int observablesA;
   int observablesB;
   int observablesC;
   int observablesE;

   int obsSet;
   int minScaleJets, maxScaleJets;
   int minJets, maxJets;
   int minPtJets, maxPtJets;
   int rsfsscale;

   Library* library;

   // booking arrays

   Book* total;

   Book* t_total;   
   Book* t_q;   
   Book* t_g;   

   Book* t_cut_total;   
   Book* t_cut_q;   
   Book* t_cut_g;   

   Book* t_cut_old_total;   
   Book* t_cut_old_q;   
   Book* t_cut_old_g;   

   Book* three_q;
   Book* three_g;
   Book* direct_sum;

   Book* JADE_E_hCMS_cut;
   Book* JADE_E0_hCMS_cut;
   Book* JADE_P_hCMS_cut;
   Book* JADE_JADE_hCMS_cut;
   Book* kT_E_Breit_cut;
   Book* Thrust_Breit;

   Book* JADE_EBreit_scale;
   Book* JADE_Q_scale;
   Book* JADE_Q_rs_scale;
   Book* JADE_Q_fs_scale;
   Book* kT_Q_scale;
   Book* kT_Q_rs_scale;
   Book* kT_Q_fs_scale;
 
   Book* pT_hCMS;
   Book* kT_Breit;
   Book* energy_Breit;

   Book* test2History;
   Book* test2Other;
   Book* test2Old;
   Book* test3History;
   Book* test3Other;
   Book* test3Old;
   Book* test4History;
   Book* test4Other;
   Book* test4Old;

   Book* kt_cmp;
   Book* kt_cut_cmp;

   double EBreitArray[200];

   Event* hCMS_event;
   Event* lab_event;

protected:

public:

   MyUserClass();
   virtual ~MyUserClass();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void Start_User();
   virtual void End_User();
   virtual int SetParameter_User(char* pname, char* parameter);
   virtual void StartIntegration1_User();
   virtual void StartIntegration2_User();
   virtual void EndIntegration_User(double average,double error);
   virtual void BeginEvent_User(Process* process);
   virtual void EndOfEvent_User(Process* process);
   virtual double EndEvent_User(
                     Process* process,
                     ContributionArray* ca,
                     int ifFinalRun
                  );
   virtual void DropEvent_User(Process* process);
   virtual void AcceptEvent_User(Process* process);
   virtual void CalculateObservablesForEvent_User(
                   Process* process, 
                   Event* event,
                   ContributionArray* ca, 
                   int ifFinalRun);
   virtual void IncidentPhaseSpace_User(
                   Process* process,
                   Event* event,
                   int& flag);
};

// ----------------------------------------------------------------------------
//
// --> the FORTRAN interface user class
//
// ----------------------------------------------------------------------------
 
class FORTRAN_User
      : public MCUserCalled
      {

private:

   DefinedObject definedObject;

   ContributionArray* caSave;
   Process* processSave;

protected:

public:

   ContributionArray* Get_caSave() { return caSave; }
   Process* Get_processSave() { return processSave; }

   FORTRAN_User();
   virtual ~FORTRAN_User();

   virtual void Define(Disaster*);
   virtual void Delete();

   virtual void Start_User();
   virtual void End_User();
   virtual int SetParameter_User(char* pname, char* parameter);
   virtual void StartIntegration1_User();
   virtual void StartIntegration2_User();
   virtual void EndIntegration_User(double average,double error);
   virtual void BeginEvent_User(Process* process);
   virtual void EndOfEvent_User(Process* process);
   virtual double EndEvent_User(
                     Process* process,
                     ContributionArray* ca,
                     int ifFinalRun
                  );
   virtual void DropEvent_User(Process* process);
   virtual void AcceptEvent_User(Process* process);
   virtual void CalculateObservablesForEvent_User(
                   Process* process, 
                   Event* event,
                   ContributionArray* ca, 
                   int ifFinalRun);
   virtual void IncidentPhaseSpace_User(
                     Process* process,
                     Event* event,
                     int& flag);
};

// ----------------------------------------------------------------------------
//
// --> individual routines called from FORTRAN
//
// ----------------------------------------------------------------------------
 
extern "C" {
void disaster_c_from_f_start();
}

// ----------------------------------------------------------------------------
 
extern "C" {
void disaster_c_from_f_stop();
}

// ----------------------------------------------------------------------------
 
extern "C" {
void disaster_c_from_f_setInt(void*, int);
}

// ----------------------------------------------------------------------------
 
extern "C" {
void disaster_c_from_f_setDouble(void*, double);
}

// ----------------------------------------------------------------------------
 
extern "C" {
void disaster_c_from_f_executeCommand(void*);
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
       );
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
     );
}

// ============================================================================ 

#endif // of USERH

// ============================================================================ 
// 
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
