// ============================================================================
// 
// --> include file for processes
//
// file:              proc.h
// created:           29.03.1997
// last modification: 24.11.1997
//
// ============================================================================

#include "switch.h"

// ============================================================================

// load only once...
#ifndef PROC_H
#define PROC_H

// ============================================================================

#include "container.h"
#include "mth.h"

// ============================================================================
//      
// symbols
//
// ============================================================================

// check consistency of the event record

#ifndef CHECK_EVENT_RECORD_FLAG
#define CHECK_EVENT_RECORD_FLAG 0
#endif

#if CHECK_EVENT_RECORD_FLAG
#define CHECK_EVENT_RECORD(A) (A) -> CER();
#else
#define CHECK_EVENT_RECORD(A) /**/
#endif

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

class Global;
class Process;
class Event;
class String;
class PartonDensity;
class Contribution;
class ContributionArray;
class Array_double;
class Array_int;
class TheoryEvaluation;
class Book;
class VariableTransformationParameters;

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

class Particle {

private:

   long ident;

   TheoryEvaluation* tEval;

   Event* event;

   int label;
   char type;
   int location;

   double energy;
   double momentum;
   double v;
   double sv;
   double sinPhi;
   double cosPhi;

   double sinPhiRef;
   double cosPhiRef;

   double fvect[4];
   double nvect[4];

public:

   // everything public for efficiency reasons
   double vz, cdp, sq;

   inline long GetIdent() { return ident; }

   inline double* GetFvect() { return fvect; }
   inline void    SetFvect(int i, double d) { fvect[i] = d; }
   inline double  GetFvect(int i) { return fvect[i]; }

   inline double* GetNvect() { return nvect; }
   inline double  GetNvect(int i) { return nvect[i]; }

   inline void   SetEnergy(double in) { energy = in; }
   inline double GetEnergy() { return energy; }
   inline void   MultiplyEnergyBy(double in) { energy *= in; }
   inline void   DivideEnergyBy(double in)   { energy /= in; }

   inline void   SetMomentum (double in) { momentum = in; }
   inline double GetMomentum() { return momentum; }
   inline void   MultiplyMomentumBy(double in) { momentum *= in; }
   inline void   DivideMomentumBy(double in)   { momentum /= in; }

   void CalculateSV();

   inline void   SetV(double in) 
                    { v = in; 
                      CalculateSV(); }

   inline double GetV() { return v; }

   inline void   SetSV(double in) { sv = in; };
   inline double GetSV() { return sv; };

   inline void   SetVAndSV(double vIn, double svIn)
                    { v  = vIn;
                      sv = svIn;
                    };

   inline void   SetVAndSV(Particle* in)
                    { v  = in -> v; 
                      sv = in -> sv;
                    };

   inline void   SetSinPhi(double in) { sinPhi = in; };
   inline double GetSinPhi() { return sinPhi; };

   inline void   SetCosPhi(double in) { cosPhi = in; };
   inline double GetCosPhi() { return cosPhi; };

   inline void   SetSinCosPhi(double sinIn, double cosIn) 
                    { sinPhi = sinIn;
                      cosPhi = cosIn; };

   inline void   GetSinCosPhi(double* sinOut, double* cosOut)
                    { *sinOut = sinPhi; *cosOut = cosPhi; }
 
   inline void   SetSinCosPhi(Particle* in)
                    { sinPhi = in -> sinPhi; 
                      cosPhi = in -> cosPhi;
                    };

   inline void   SetSinPhiRef(double in) { sinPhiRef = in; };
   inline double GetSinPhiRef() { return sinPhiRef; };

   inline void   SetCosPhiRef(double in) { cosPhiRef = in; };
   inline double GetCosPhiRef() { return cosPhiRef; };

   inline void   SetSinCosPhiRef(double sinIn, double cosIn) 
                    { sinPhiRef = sinIn;
                      cosPhiRef = cosIn; };
   inline void   GetSinCosPhiRef(double* sinOut, double* cosOut)
                    { *sinOut = sinPhiRef; *cosOut = cosPhiRef; }

   inline double CalculatePhi() { return atan2ZeroTwoPi(sinPhi, cosPhi); }

   inline void   SetPhi(double in) 
                    { sinPhi = sin(in); 
                      cosPhi = cos(in); }

   inline void SetLabel(int in) { label = in; } 
   inline int  GetLabel() { return label; }

   inline void SetType(char in) { type = in; }
   inline char GetType() { return type; }

   inline void SetEvent(Event *in) { event = in; }
   inline Event *GetEvent() { return event; }

   inline void SetLocation(int in) { location = in; }
   inline int GetLocation() { return location; }

   INLINE Particle();
   INLINE ~Particle();

   Particle& operator=(const Particle&);

   void Define(TheoryEvaluation*);
   void Delete();
   void Reset(TheoryEvaluation*);
   void ReturnToFreeList();

   void CrossFromCartesian(Particle*, Particle*);
   double DotIntoCrossFromCartesian(Particle*, Particle*);
   void SetToZ();
   void SetToX();
   void SetToY();
   void SetToE();
   void CalculateRelativeAzimuth(
           Particle* zaxis, 
           Particle* xaxis,
           double* sinPhi,
           double* cosPhi
        );
   void CopyInto(Particle*);
   void Add(Particle*, Particle*);
   void AddCartesian(Particle*, Particle*);
   void Add_E0Cartesian(Particle*, Particle*);
   void Add_PCartesian(Particle*, Particle*);
   void Subtract(Particle*, Particle*);
   void ScaleBy(double);
   void Divide(Particle*, Particle*);
   void TransferLabelsInto(Particle*);
   void Print(FILE*);
   inline void Print() { Print(stdout); }
   void SetToZeroMomentum();

   void ScaleByCartesian(double);
   void LinearCombinationCartesian(
           double, 
           Particle*, 
           double, 
           Particle*
        );
   void LinearCombination(
           double, 
           Particle*, 
           double, 
           Particle*
        );
   INLINE double Get_pT();

   void EpsilonCartesian(Particle* a, Particle* b, Particle* c);
   void Epsilon(Particle* a, Particle* b, Particle* c);

   double CalculateMass();
   double CalculateImaginaryMass();

   double CalculateMassCartesian();
   double CalculateImaginaryMassCartesian();

   // --> change of frame
   void NewComponentsCartesian(
           Particle* old, 
           Particle* nE, 
           Particle* nX, 
           Particle* nY, 
           Particle* nZ);

   double Calculate_pT();
   double Calculate_PseudoRapidity();

   int UseListIsSame(Particle*);

   void CalculateNormalFromPolar();
   void CalculateCartesianFromPolar();
   void CalculateCartesianAndNormalFromPolar();
   void CalculatePolarFromCartesian();

   friend double Dot(Particle*, Particle*);
   friend double DotOld(Particle*, Particle*);
   friend double DotCartesian(Particle*, Particle*);
   friend double InvariantMassCartesian(Particle*, Particle*);
   friend double DotSpaceCartesian(Particle*, Particle*);
   friend double V(Particle*, Particle*);
   friend double VCartesian(Particle*, Particle*);
   friend double VNormal(Particle*, Particle*);
   friend void   DotAndVMasslessNormal(
                    Particle*, 
                    Particle*, 
                    double* d, 
                    double* v
                 );
   friend void   InvariantAndVMasslessNormal(
                    Particle*, 
                    Particle*, 
                    double* d, 
                    double* v
                 );
};

// ----------------------------------------------------------------------------

// the Lorentz product (p.q)
double Dot(Particle*, Particle*); 
double DotOld(Particle*, Particle*); 

// the Cartesian Lorentz product (p.q)
double DotCartesian(Particle*, Particle*); 

// the Cartesian invariant mass (p+q)^2
double InvariantMassCartesian(Particle*, Particle*); 

// the Cartesian product of the space components (vect(p).vect(q))
double DotSpaceCartesian(Particle*, Particle*); 

// the variable 0.5*(1-cos(theta))
double V(Particle*, Particle*); 

// the variable 0.5*(1-cos(theta))
double VCartesian(Particle*, Particle*); 

// the variable 0.5*(1-cos(theta))
double VNormal(Particle*, Particle*); 

// calculate both the inner product and the v-variable
void DotAndVMasslessNormal(
        Particle*, 
        Particle*, 
        double* d, 
        double* v
     );

// calculate both the invariant mass and the v-variable
void InvariantAndVMasslessNormal(
        Particle*, 
        Particle*, 
        double* d, 
        double* v
     );

// ----------------------------------------------------------------------------

typedef FreeList_1 <Particle, TheoryEvaluation*> ParticleFreeList;

// ----------------------------------------------------------------------------

class Invariant {

private:

   int n_in_partons_local, npartons_local;
   int clusterDefined;
   int dtTermsDefined;

protected:

public:

   // public for efficiency...

   // kinematical variables

   double Q2inv;

   double e0,e1,e2;
   double f0,f1,f2;
   double ei,fi;

   double e3,f3,v03,v13,v23;

   double vi0,vi1,vi2,v01,v02,v12;
   double vkl,vik,vil,v0k,v0l,v1k,v1l,v2k,v2l;

   double vi3;

   double si0,si1,si2,s01,s02,s12;
   double skl,sik,sil,s0k,s0l,s1k,s1l,s2k,s2l;

   double lambda, eta;
   double eh, fh;
   double dti, dt0, dt1, dt2, dt3, dtk, dtl;
   double sih, s0h, s1h, s2h, s3h, skh, slh;

   Invariant();
   ~Invariant();

   inline void SetLambda(double in) { lambda = in; }
   inline void SetEta(double in)    { eta    = in; }

   INLINE void CalculatePartonInvariants(Event* e);
   INLINE void CopyPartonInvariantsFrom(Invariant*);
   void Print(FILE*);
};

// ----------------------------------------------------------------------------

class Event 
   : public Invariant, 
     public DefinedObject
   {

private:

   HorizontalList <Event>* hl;

   TheoryEvaluation* tEval;

   int maxparticles,nparticles,addentries;

   Global* global;

   Disaster* disaster;

   Frame frame;
   P0_Reference p0ref;

   int limit1, limit2;
   LimitType limitType;

   long ident;
   long id;

   Frame userFrame;

   int isConstructed;
   SubtractionType toBeSubtracted;

   EventType eventType;

   double xJacobian;
   int varpos; 
   double ustore;

   double fsJacobian;

   Contribution *contribution;

   int alreadyInUseList;

   VariableTransformationParameters* vtp;
   
public:

   // public for efficiency reasons

   double ecm,ecmfull,SH,proton_energy;
   double xB,y,xi,Q,W,Q2,W2;

   double u, ujacobian;
   double xi_hard;

   double xi_int, xi_pden;

   int n_in_partons;
   Particle **particle; 
   Particle **mapping; 
   Particle **mappingSave; 
   int npartons;
   int nvlist;
   double *vlist,*flist;
   int **mapij;
   int *imap, *jmap;
   double **wlist;
   double *dtlist;
   int *iwlist,iwentries;

   inline int** GetMapij() { return mapij; }

   inline void SetNpartons(int in) { npartons = in; }
   inline int GetNpartons() { return npartons; }

   inline void SetP0ref(P0_Reference in) { p0ref = in; }
   inline P0_Reference GetP0ref() { return p0ref; }

   inline int GetN_in_partons() { return n_in_partons; }
   inline Frame GetFrame() { return frame; }

   inline void SetXJacobian(double in) { xJacobian = in; }
   inline double GetXJacobian() { return xJacobian; }

   inline void SetVarpos(int in) {varpos = in; }
   inline int GetVarpos() { return varpos; }

   inline void SetUstore(double in) { ustore = in; }
   inline double GetUstore() { return ustore; }

   inline void SetFSJacobian(double in) { fsJacobian = in; }
   inline double GetFSJacobian() { return fsJacobian; }

   inline long GetIdent() { return ident; }
   inline long GetId() { return id; }

   inline void SetLimit1(int in) { limit1  = in; }
   inline void SetLimit2(int in) { limit2  = in; }
   inline void SetLimitType(LimitType in)  { limitType  = in; }
   inline int GetLimit1() { return limit1; }
   inline int GetLimit2() { return limit2; }
   inline LimitType GetLimitType() { return limitType; }

   inline Event* GetHorizontal() { return NULL; }

   inline void       SetMapping(int i, Particle *p) { mapping[i] = p; }
   inline Particle*  GetMapping(int i) { return mapping[i]; }
   inline Particle** GetMapping() { return mapping; }

   inline int GetDefineParameter() { return maxparticles; }

   inline HorizontalList <Event> * GetH() { return hl; }  

   inline void SetUserFrame(Frame in) { userFrame = in; }
   inline Frame GetUserFrame() { return userFrame; }

   inline void SetIsConstructed (int in) { isConstructed = in; }
   inline int GetIsConstructed () { return isConstructed; }

   inline void SetToBeSubtracted (SubtractionType in) { toBeSubtracted = in; }
   inline SubtractionType GetToBeSubtracted () { return toBeSubtracted; }

   inline void SetEventType(EventType in) {eventType = in; }
   inline EventType GetEventType() { return eventType; }

   inline void SetContribution(Contribution* in) { contribution = in; }
   inline Contribution* GetContribution() { return contribution; }

   inline void SetAlreadyInUseList(int in) { alreadyInUseList = in; }
   inline int  GetAlreadyInUseList() { return alreadyInUseList; }

   inline void Set_vtp(VariableTransformationParameters* in)
                  { vtp = in; }
   inline VariableTransformationParameters* Get_vtp() { return vtp; }

   inline double Get_Q2() { return Q2; }
   inline double Get_xB() { return xB; }
   
   // --> class administration
   Event();
   ~Event();
   void Define(int, Global*);
   void DefineX(int, Global*);
   void Define(TheoryEvaluation*, int);
   void Reset();
   void Reset(TheoryEvaluation*);
   void ReturnToFreeList();

   // --> utilities
   void Print(Mstream&);
   void Print(FILE*);
   inline void Print() { Print(stdout); }
   void PrintMapping(FILE*);

   inline double Get_xi_hard() { return xi_hard; }

   // GPS modification for compilation with gcc 2.95.2
   int GetMaxNumber();
   void SetPartonNumber(int);
   int GetPartonNumber(); // :_MOD_: redundant
   void SetIn(int);
   int GetIn();           // :_MOD_: redundant
   Particle **GetParticle();
   inline Global* GetGlobal() { return global; };

   // --> copy operations
   void CopyPSInformationInto(Event*);
   void CopyInto(Event*);
   void CopyAdditionalParticlesInto(Event*);

   // --> permutations and reorganization
   void Permute(int, ...);
   void PermuteTwoTo01(int,int);
   void PermuteOneTo0(int);
   void PermuteBackTwoTo01(int,int);
   void RenamePartons();
   void SwapLabelAndType(int, int);
   void SwapLabelAndTypeTwoTo01(int, int);
   void SwapLabelAndTypeOneTo0(int);
   void Rename(int oldlabel, int newlabel);
   void TransferPartonLabelsInto(Event*);
   void PermuteLabels(int co_max, int* list);

   // --> event record administration
   Particle* AddParticle();
   Particle* AddParticle(int label, char type);

   void SchemeCombineCartesianMap(RecombinationType rt, int, int);
   void Divide(Event *a, Event *b);
   void CreateMapping();
   void SetMappingToNULL();
   void CER();  // Check Event Record (structural integrity)
   void MakeVList();   
   void MakeEList();   
   void PrintVELists(FILE*);   
   double SumPartialV(int ifixed, int jfixed);
   double SumPartialE(int ifixed);
   void DoSumOverSubtractionTerms(
           Process*,
           int ic,
           double weight_factor,
           ContributionArray* ca
        );
   void DoSumOverDiffSubtractionTerms(
           Process*,
           int ic,
           double weight_factor,
           ContributionArray* ca
        );
   void DoSumOverAddedSubtractionTerms(
           Process*, 
           int ic,
           double weight_factor,
           ContributionArray* ca
        );
   void DoSumOverCollinearInitialTerms(
           Process*, 
           int ic,
           double weight_factor,
           ContributionArray* ca
        );

   void DefineEcm();
   double GetEcm();
   void SetEcmFull(double);
   double GetEcmFull();
   void SetLeptonVariables(double xB, double y);
   double GetxB();
   double Gety();
   void Setxi(double xi);
   double Getxi();

   // -> parton phase space parametrizations
   void MapUnitToPSGeneral(
           PSMode, 
           P0_Reference,
           double* unit,
           double* jac1, 
           double* jacrest
        );
   void MapUnitToPSGeneralPCMS(
           PSMode, 
           P0_Reference,
           double* unit,
           double* jac1, 
           double* jacrest 
        );
   void MapUnitToPSGeneralHCMS(
           PSMode, 
           P0_Reference,
           double* unit,
           double* jac1, 
           double* jacrest
        );
   double MapUnitToPSChoice(
             P0_Reference,
             double* unit
          );
   double MapUnitToPS(double* unit);
   void MapPSToUnitGeneral(
           P0_Reference, 
           double* unit
        );
   void MapPSToUnitGeneralPCMS(
           P0_Reference, 
           double* unit
        );
   void MapPSToUnitGeneralHCMS(
           P0_Reference,
           double* unit
        );
   void MapPSToUnit(
           double* unit
        );

   // --> limit operations on phase space parametrizations
   void Limit1(
           LimitType,
           P0_Reference,
           Event* limit,
           double* jrelative,
           double* energy,
           double* v
        );
   void Limit(
           LimitType,
           int i,
           int j,
           Event* limit,
           double* jrelative,
           double* energy,
           double* v
        );

   // --> DIS lepton phase space
   double MapUnitToxBy(double zx, double zq);
   double LeptonJacobian();

   // --> initial state parton momentum parametrizations
   double MapUnitToXi(double zxi);
   void   MapUnitTo_u(double zu);

   // --> methods to add additional particles
   void EEAddParticles();
   void DISAddParticles();
   void DISAddRemnantAndIncident();
   void EcmAddRemnantAndIncident();

   void SetFrame(Frame);
   void CalculatePartonInvariantsLocal();
   void CalculateWList();
   void SetKinematicalDefaults();
   void DetermineHLambdaEta();
   int CheckTechnicalCut(Process*, char* location);
   Event* CreateAdded(Process*, LimitType lt, int iangle, int jangle);
   void AddedRenumbering(int iangle, int jangle);

   // --> Lorentz boosts
   void ZBoost(double alpha);
   void BoostTo_hCMS();
   void BoostTo_Breit();
   void BoostTo_DIS_Lab(double ratio_Ep_Ee);
   void BoostWithCartesianVector(
           Particle* nE, 
           Particle* nX, 
           Particle* nY,
           Particle* nZ
        );
   void BoostEZX(
           Particle* vE, 
           Particle* vZ, 
           Particle* vX
        );

   // --> observables: jets and event shapes
   void SmallestInvariant(
           ClusterAlgorithmType cat,
           RecombinationType rt,
           double& minscale2,
           int& index1, 
           int& index2
        );
   void SmallestInvariantCartesian(
           ClusterAlgorithmType cat,
           RecombinationType rt,
           double& minscale2,
           int& index1, 
           int& index2
        );
   void ClusterGeneral(
           ClusterAlgorithmType cat,
           RecombinationType rt,
           double scale2,
           int& startCluster, 
           int& finalCluster
        );
   void ClusterGeneralHistory(
           ClusterAlgorithmType cat,
           RecombinationType rt,
           Array_double* scale2,
           Array_int* nCluster,
           int& nEntries
        );
   void ClusterGeneralOneStep(
           ClusterAlgorithmType cat,
           RecombinationType rt
        );
   void ClusterGeneralWithEvent(
           Event* copy,
           ClusterAlgorithmType cat,
           RecombinationType rt,
           double scale2,
           int& startCluster, 
           int& finalCluster
        );
   void ClusterGeneralHistoryWithEvent(
           Event* copy,
           ClusterAlgorithmType cat,
           RecombinationType rt,
           Array_double* scale2,
           Array_int* nCluster,
           int& nEntries
       );
   int ClusterInvariantOld(
          ClusterAlgorithm cluster_algorithm,
          double scale,
          int& startCluster, 
          int& etaCut, 
          int iflags
       );
   double CurrentThrust();

   // --> phase space cuts
   int Cut_Jets_pT_PseudoRapidity(
          double pT_Min, 
          double pT_Max,
          double eta_Min,
          double eta_Max
       );
   int Cut_Lepton_pT_PseudoRapidity(
          double pT_Min, 
          double pT_Max,
          double eta_Min,
          double eta_Max  
       );

   int UseListIsSame(Event*);
   Event* UserFrameFindEvent(Frame);

   // --> support to prepare histograms  
   // :_MOD_: should be the other way round? 
   //         ... derive a class from Book...
   void PrestoreGraph(Book* book,
                      double x, 
                      int entry
        );
   void PrestoreHistogram(
        Book* book,
        double x
        );
   void PrestoreHistogram(
           Book* book,
           double x, 
           double weight          
        );
   
   // --> inner products of entries in the event record
   double DotMap(int, int);
   double DotCartesianMap(int, int);
   double InvariantMassCartesianMap(int, int);
   double VNormalMap(int, int);
   double VMap(int, int);
   void   DotAndVMasslessNormalMap(
             int, 
             int,
             double* d, 
             double* v
          );
   void   InvariantAndVMasslessNormalMap(
             int, 
             int,
             double* d, 
             double* v
          );

   // --> calculation of alternative sets of phase space variables
   void CalculateNormalFromPolarRangeMap(int imin, int imax);
   void CalculateCartesianFromPolarRangeMap(int imin, int imax);
   void CalculateCartesianAndNormalFromPolarRangeMap(int imin, int imax);
   void CalculatePolarFromCartesianRangeMap(int imin, int imax);
};

// ----------------------------------------------------------------------------

typedef FreeListArray_1 <Event, TheoryEvaluation*> EventFreeListArray;

typedef UseListContainer <Event> EventUseListContainer;
typedef UseList <Event> EventUseList;

// ----------------------------------------------------------------------------

class VariableTransformationParameters {
   
private:

   Array_int* offset;
   int* offsetData;
   int n_offset;
   
   Array_int* mapFunction;
   int* mapFunctionData;
   int n_mapFunction;

   Array_double* mapParameter;
   double* mapParameterData;
   int n_mapParameter;

protected:
  
public:

   inline int    Get_offset(int i) { return offsetData[i]; }
   inline int    Get_mapFunction(int i) { return mapFunctionData[i]; }
   inline double Get_mapParameter(int i) { return mapParameterData[i]; }

   VariableTransformationParameters();
   ~VariableTransformationParameters();

   void CopyInto(VariableTransformationParameters*);

   void Set_offset(      int pos, int    in);
   void Set_mapFunction( int pos, int    in);
   void Set_mapParameter(int pos, double in);

   void PrintParameters(Mstream&);
};

// ----------------------------------------------------------------------------

class MatrixElement : public Invariant {

private:

protected:

   Global* global;
   Process* process;
   String* name;
   String* particleAssignment;

   int nIntVar;

public:

   MatrixElement();
   virtual ~MatrixElement();

   inline int Get_nIntVar() { return nIntVar; }
   inline void Set_nIntVar(int in) { nIntVar = in; }

   void AssignGlobal(Global*);
   void AssignProcess(Process*);
   virtual void Define();
   Process *GetProcess();
   String *GetName();
   char *GetParticleAssignment();
   virtual int GetIn() = 0;
   virtual int GetOut() = 0;
   virtual void Evaluate(
                   Event*, 
                   Contribution* con,
                   double factor
                ) = 0;
};

// ----------------------------------------------------------------------------

class Process {

private:

   int minout;
   int maxout;
   int max_nIntVar;
   
   ContributionArray* ca;

   EventUseList* eventList;

protected:

   Global* global;
   Disaster* disaster;

   MatrixElement*** me_list;  // this is also the flag indicating whether
                              // the matrix elements have been defined
   SubtractionType* to_be_subtracted;
   int n_components;
   int* n_me;

   int* componentFlags;

   VariableTransformationParameters* vtp;

   int alpha_s_Order;
   Frame frame;

public:

   inline Global* GetGlobal() { return global; }
   inline int GetAlphaS_Order() { return alpha_s_Order; }
   inline ContributionArray* GetCA() { return ca; }
   inline Frame GetFrame() { return frame; }   
   inline void  SetFrame(Frame in) { frame = in; }   
   inline EventUseList* GetEventList() { return eventList; }

   // --> variable transformations
   inline VariableTransformationParameters* Get_vtp(int i)
                                               { return vtp + i; }

   inline int Get_componentFlag(int i) { return componentFlags[i]; }

   Process();
   virtual ~Process();

   void Delete();
   virtual void Define(Global*) = 0;
   void DefineProcess(Global*);

   void DeleteMatrixElements();
   virtual void AssignMatrixElements() = 0;
   void CreateNME(int n_componentsI);
   void CreateOtherArrays();
   void AssignInformationToMatrixElements();
   int GetIn();
   int GetOut(int);
   int GetMaxOut();
   inline int GetMinOut() { return minout; };
   void SetMinMaxOut();
   int GetNComponents();
   int GetNMe(int);
   int Get_nIntVar(int);
   void Set_Max_nIntVar();
   int Get_Max_nIntVar();
   MatrixElement* GetPMatrixElement(int,int);
   SubtractionType GetToBeSubtracted(int);
   void Evaluate(int flag, int ic, LimitType, 
                 int ifix, int jfix, Event*, double factor,
                        int lstore);
   void CalculateWeight(
           double* unit,
           ContributionArray* ca
        );

   // copy variable transformation parameters
   void Copy_vtp(int to, int from);

   // --> management of components
   void Set_componentFlag(int i, int val);
   void Print_componentFlags(Mstream&);
   inline void Print_componentFlags() { 
                  Print_componentFlags(global -> To_status()); 
               }
};

// ============================================================================

#endif // of PROC_H

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
