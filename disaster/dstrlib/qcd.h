// ============================================================================
// 
// --> QCD include file: 
//     class definitions etc.
//
// file:              qcd.h
// created:           24.06.1997
// last modification: 15.12.1997
//
// ============================================================================

// load only once...
#ifndef QCDH
#define QCDH

// ============================================================================

#include "switch.h"
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

// ============================================================================
//
// --> classes
//
// ============================================================================

// function that determines the jet multiplicity
// :_MOD_: should have a class Cluster...

int clusterNumber(Array_double *scale2,
                  Array_int *nCluster,
                  int nEntries,
                  double testscale2
);

// ----------------------------------------------------------------------------

// booking device

class Book {

private:

   Global* global;

   int defined;

   BookMode bookMode; 
   ScaleType scaleType;
   int size; 
   char *name;
   char *id;

   int as_counter;
   int as_counter_max;

   Array_double *xarray; 
   double *xdata;

   Array_double *yarray; 
   double *ydata;

   Array_double *y2array; 
   double *y2data;

   Array_double *errors; 
   double *errorsData;

   Array_double *relerrors; 
   double *relerrorsData;

   Array_double *zarray; 
   double *zdata;

   Array_double* prestore;
   double* psData;

   Array_int* prestoreFlag;
   int* psFData;

   Array_double* prestoreFactor;
   double* psFactData;

   TwoD_Array_double* prestore2;
   double** ps2Data;
   
   int prestoreSize;

   int counter; 

protected:

public:

   Book(Global*);
   ~Book();

   inline double GetX(int i) { return xdata[i]; }
   inline int    GetSize() { return size; }
   inline double Get_prestore(int i) { return psData[i]; }
   inline double Get_prestore2(int event, int entry) 
                    { return ps2Data[event][entry]; }
   inline double GetBeforeSum(int i) { return zdata[i]; }

   void Define(
           BookMode, 
           ScaleType, 
           int size, 
           double smin, 
           double smax,
           char* name, 
           char* id, 
           int prestoreSize
        );
   void Delete();
   void ResetStorage();
   void CalculateErrors();

   void New_X_Entry(double);

   void PrepareForEvent();
   void AddUpEvent();

   void StoreHistogram(double, double);
   void StoreGraph    (int,    double);

   void PrestoreHistogram(double x, int);
   void PrestoreHistogram(double x, int, double weight);
   void PrestoreGraph(double x, int event, int entry);

   void SPS_Histogram(int event, double);
   void SPS_Graph(int event, int entry, double);

   void Print(Mstream&);
   void PrintMathForm(FILE*, int numid);

   void GetEntry(
           int i, 
           double* average, 
           double* error, 
           double* relError
        );
};

// ----------------------------------------------------------------------------

// manager for booking devices

class Library {

private:

   Global* global;

   int defined;

   ArrayOfPointers <Book*> *books;
   Book **bookData;
   int nBooks;
   int maxBooks;

protected:

public:

   Library(Global*, int);
   ~Library();

   void Define(int);
   void Delete();

   Book *CreateNewBook();   
   void Print(Mstream&);   
   void PrintMathForm(FILE *);   
   void ResetStorage();
   void CalculateErrors();
   void PrepareForEvent();
   void AddUpEvent();
};

// ----------------------------------------------------------------------------

class Contribution {

private:

   int defined;
   Global *global;
   Array <double> *clist;
   double *data;
   int length;
   int orderAlphaS;
   ScaleLogarithmFactor slf;
   int orderAlphaEM;
   int matrixElement;
   int location;
   Event* event;
   Contribution* sameEvent;

   // defined in and used by the user routine
   int nflavour_out_U;
   int nflavour_pden_U;
   double ren_scale_U;
   double fact_scale_U;
   double xi_U;
   double Q_U;
   double EBreit_U;

protected:

public:

   Contribution();
   Contribution(Global*, int);
   ~Contribution();

   inline int  GetDefineParameter() { return length; }
   inline long GetIdent()           { return -98; }

   inline int  GetLocation() { return location; }
   inline void SetLocation(int in) { location = in; }

   inline Array <double> *GetClist() { return clist; }
   inline double *GetData() { return data; }

   inline int GetLength() { return length; }

   inline int  GetOrderAlphaS() { return orderAlphaS; }
   inline void SetOrderAlphaS(int in) { orderAlphaS = in; }

   inline ScaleLogarithmFactor GetSLF() { return slf; }
   inline void                 SetSLF(ScaleLogarithmFactor in) { slf = in; }

   inline int  Get_orderAlphaEM() { return orderAlphaEM; }
   inline void Set_orderAlphaEM(int in) { orderAlphaEM = in; }

   inline int  GetMatrixElement() { return matrixElement; }
   inline void SetMatrixElement(int in) { matrixElement = in; }
  
   inline void   SetEvent(Event* in) { event = in; }
   inline Event* GetEvent() { return event; }

   inline void   SetSameEvent(Contribution* in) { sameEvent = in; }
   inline Contribution* GetSameEvent() { return sameEvent; }

   inline int  Get_nflavour_out_U() { return nflavour_out_U; }
   inline void Set_nflavour_out_U(int in) { nflavour_out_U = in; }

   inline int  Get_nflavour_pden_U() { return nflavour_pden_U; }
   inline void Set_nflavour_pden_U(int in) { nflavour_pden_U = in; }

   inline double Get_ren_scale_U() { return ren_scale_U; }
   inline void   Set_ren_scale_U(double in) { ren_scale_U = in; }

   inline double Get_fact_scale_U() { return fact_scale_U; }
   inline void   Set_fact_scale_U(double in) { fact_scale_U = in; }

   inline double Get_xi_U() { return xi_U; }
   inline void   Set_xi_U(double in) { xi_U = in; }

   inline double Get_Q_U() { return Q_U; }
   inline void   Set_Q_U(double in) { Q_U = in; }

   inline double Get_EBreit_U() { return EBreit_U; }
   inline void   Set_EBreit_U(double in) { EBreit_U = in; }

   void Define(Global*, int length);
   void Delete();

   void Reset(Global*);
   inline void ReturnToFreeList() {}
   INLINE void ClearEntry();

   void Print(FILE*, int mode);
   void Print(FILE*);

   int UseListIsSame(Contribution*);
};

// ----------------------------------------------------------------------------

class ContributionArray {

private:

   int defined;
   Global* global;
   Array <Contribution> * entry;
   Contribution* cptr;
   int length;
   int maxContributions;
   int actual;
    
protected:

public:

   ContributionArray();
   ContributionArray(Global*, int maxContributions, int length);
   ~ContributionArray();

   inline void GetDefineParameter(int &maxContributionsOut, int &lengthOut) { 
                  maxContributionsOut=maxContributions;
                  lengthOut=length;
               };

   inline long GetIdent() { return -999999; };

   inline int  GetActual() { return actual; }
   inline void IncrementActual() { actual++; }
   inline void ResetActual() { actual = -1; }
   inline void SetActual(int in) { actual = in; }

   inline Array <Contribution> *GetContributionList() { return entry; }
   inline Contribution *GetContribution(int i) { return cptr+i; }
   inline int GetMaxContributions() { return maxContributions; }
   inline int GetLength() { return length; }
   inline Array <double> *GetCList(int i) { return cptr[i].GetClist(); }
   inline double *GetData(int i) { return cptr[i].GetData(); }

   void Define(Global*, int maxContributions, int length);
   void Delete();
   void Reset(Global *);
   void ReturnToFreeList() {}
   void Print(FILE *, int mode);
   void Print(FILE *);

   int GetActual_SetToZero_Increment();
   int Increment_SetToZero_GetActual();
   Contribution *NextContribution();

   int UseListIsSame(ContributionArray*);
};

// ----------------------------------------------------------------------------

class PartonList 
   {

private:

   int defined;
   Array_double *data;
   int nflavours;
    
protected:

   Global *global;

public:

   PartonList();
   PartonList(Global*, int);
   ~PartonList();

   inline int GetDefineParameter() { return nflavours; };
   inline long GetIdent() { return -999999; };

   inline double *GetDataPtr()
                     { return data -> GetDataPtr(); }

   inline int     GetNflavours() 
                     { return nflavours; }

   void Define(Global *, int);
   void Delete();
   void Reset(Global *);

   void Print(Mstream&);
   void Print(FILE*);
   void Print() { Print(stdout); }

   int UseListIsSame(PartonList *);
   
   void ReturnToFreeList();
};

// ----------------------------------------------------------------------------

class PartonListArray {

private:

   int defined;
   Array <PartonList> *data;
   int length;
   PartonList *ptr;
    
protected:

   Global *global;

public:

   PartonListArray();
   PartonListArray(Global*, int);
   ~PartonListArray();

   inline int GetDefineParameter() { return length; }
   inline long GetIdent() { return -999999; };

   inline PartonList *GetDataPtr()
                         { return data -> GetDataPtr(); }

   inline int GetLength() 
                 { return length; }

   void Define(Global *, int);
   void Delete();
   void Reset(Global *);
   inline void ReturnToFreeList() {}
};

// ----------------------------------------------------------------------------

   typedef FreeListArray_1 <PartonList, Global*> PartonListFreeList;
//   typedef FreeListArray_1 <PartonListArray, Global*> PartonListArrayFreeList;
   
   typedef FreeListArray_1 <Contribution, Global*>
              ContributionFreeList;
   typedef FreeListTwoD_Array_1 <ContributionArray, Global*>
              ContributionArrayFreeList;

// ----------------------------------------------------------------------------

class PartonDensity 
  : public PartonList
   {

private:

   int defined;
   int nflavours;

protected:

public:

   PartonDensity();
   ~PartonDensity();

   inline int GetDefineParameter() { return nflavours; }

   void Define(Global*, int);
   void Delete();
   void Reset(Global *);

   void Print(Mstream&);
   void Print(FILE*);
   void Print() { Print(stdout); }

   int UseListIsSame(PartonDensity *);
   void ReturnToFreeList();
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   typedef FreeListArray_1 <PartonDensity, Global*> PartonDensityFreeList;

// ----------------------------------------------------------------------------

class PartonDensityParametrization 
   : public GlobalObject
   {

private:

   int defined;
   int collection, parametrization, set;
   int max_flavours;

   int nflavours, order;
   double lambda_4;
   double x_min, x_max, scale_min, scale_max;

   PartonDensityFreeList* pdfl;

   Cache <PartonDensity> * cache;
   Array_int* intKeyArray;
   int* intKey;
   Array_double* doubleKeyArray;
   double* doubleKey;

protected:

public:

   void SetPDFL(PartonDensityFreeList* in) { pdfl = in; }

   PartonDensityParametrization(Global*);
   ~PartonDensityParametrization();

   inline void SetPDFCollection(int cIn)
                  { collection = cIn; }

   inline int GetPDFCollection()
                  { return collection; }

   inline void SetPDFParametrization(int pIn)
                  { parametrization = pIn; }

   inline int GetPDFParametrization()
                  { return parametrization; }

   inline void SetPDFSet(int sIn)
                  { set = sIn; }

   inline int GetPDFSet()
                  { return set; }

   inline double GetPDFLambda_4()
                  { return lambda_4; }

   void Define();
   void Delete();
   void FillPartonDensity(
           PartonDensity*, 
           double x, 
           double scale, 
           int nquarks
        );
   void FillPartonDensity_0(
           PartonDensity*,
           double x, 
           double scale, 
           int nquarks
        );
   void FillPartonDensityPDFLIB(
           PartonDensity*,
           double x, 
           double scale, 
           int nquarks
        );
   void FillPartonDensity(
           PartonDensity*, 
           double x, 
           double scale, 
           int nquarks,      
           int collection,
           int parametrization,
           int set
        );
   PartonDensity* CreateAndFillPartonDensity(
                     double x, 
                     double scale, 
                     int nquarks,      
                     int collection,
                     int parametrization,
                     int set
                  );
   PartonDensity* FindOrCreateAndFillAndCachePartonDensity(
                     double x, 
                     double scale,  
                     int nquarks,      
                     int collection,
                     int parametrization,
                     int set
                  );
   void InitializeCollection();
   void InitializeParametrization(
           int collection, 
           int parametrization, 
           int set
        );
   void CreateCache(int);
   void DeleteCache();
};

// ----------------------------------------------------------------------------

class CouplingConstantEvaluation 
   : public GlobalObject
   {

private:

   int defined;  // :_MOD_: make this protected, and save the GetDefined()!

protected:

   int variant;
   int order;

   int outputFlag;

public:

   CouplingConstantEvaluation(Global*);
   virtual ~CouplingConstantEvaluation();

   inline int GetDefined()
                  { return defined; }

   inline void SetVariant(int vIn)
                  { variant = vIn; }

   inline int GetVariant()
                  { return variant; }

   inline void SetOrder(int oIn)
                  { order = oIn; }

   inline int GetOrder()
                  { return order; }

   inline void Set_outputFlag(int in) { outputFlag = in; }
   inline int  Get_outputFlag()       { return outputFlag; }

   void Define();
   void Delete();
   double CouplingConstantEvaluation::EvaluateCoupling(double,int);
};

// ----------------------------------------------------------------------------

// :_MOD_: should split this into a class HeavyQuarks and AlphaS.
// :_MOD_: also should put lambda into the QCD class.

class AlphaS
   : public CouplingConstantEvaluation
   {

private:

   double mCharm;
   double mBottom;
   double mTop;

   double lambdaQCD[6+1];

protected:

public:

   AlphaS(Global*);
   ~AlphaS();

   void Define();
   void Delete();

   inline void SetMcharm (double m) { mCharm  = m; }
   inline void SetMbottom(double m) { mBottom = m; }
   inline void SetMtop   (double m) { mTop    = m; }

   inline double GetMCharm()  { return mCharm;  };
   inline double GetMBottom() { return mBottom; };
   inline double GetMTop()    { return mTop;    };

   inline void   SetLambdaQCD_4 (double lambdaIn) { lambdaQCD[4] = lambdaIn; }
   inline double GetLambdaQCD_4 ()                { return lambdaQCD[4]; }

   double EvaluateAlphaS(double scale, int nquarks);
   double EvaluateAlphaSAutomaticNquarks(double scale);
   double EvaluateAlphaSAutomaticNquarks(
             double scale,
             int variantIn,
             int orderIn, 
             double lambdaQCD4In
          );

   void   CalculateLambdaValues();
   double EvaluateRunningCoupling(double scale, int nquarks);
   int    FlavourSwitch(double scale);
   double EvaluateRunningCouplingAutomaticNquarks(double scale);
   void   InitializeAlphaS(
             int variantIn,
             int orderIn, 
             double lambdaQCD4In
          );
};

// ----------------------------------------------------------------------------

class AlphaEM
   : public CouplingConstantEvaluation
   {

private:

   double fine_structure_constant;

protected:

public:

   AlphaEM(Global*);
   ~AlphaEM();

   void Define();
   void Delete();

   double EvaluateAlphaEM(double scale, int variantIn);
   double EvaluateAlphaEM(double scale);
};

// ----------------------------------------------------------------------------

class QCD 
   : public virtual Global
   {

private:

   int defined;

   Global* global;

   PartonList* quark_charge; 
   int nflavours;

   Array_double* quarkChargeSquaredArray;
   double* quarkChargeSquared;

   TwoD_Array_double* quarkChargeProductArray;
   double** quarkChargeProduct;

   Array_double* quarkChargeSquaredPartialSumSingleArray;
   double* quarkChargeSquaredPartialSumSingle;

   TwoD_Array_double* quarkChargePartialSumDoubleArray;
   double** quarkChargePartialSumDouble;

   TwoD_Array_double* quarkChargeSquaredPartialSumDoubleArray;
   double** quarkChargeSquaredPartialSumDouble;

   PartonDensityFreeList* pdFree;

   AlphaS* alphaS_Server;  

   PartonDensityParametrization* pdf_Server;

protected:

public:

   double ncolour, cf, quarkfactor, gluonfactor;

   inline PartonDensityFreeList* 
             Get_pdFree() { return pdFree; }

   inline AlphaS* 
             Get_alphaS_Server() { return alphaS_Server; }

   inline PartonDensityParametrization* 
             Get_pdf_Server() { return pdf_Server; }
             
   QCD();
   ~QCD();

   inline PartonList* GetQuarkCharge()
                         { return quark_charge; }

   inline Array_double* GetQuarkChargeSquared()
                           { return quarkChargeSquaredArray; }

   inline TwoD_Array_double* GetQuarkChargeProduct()
                                { return quarkChargeProductArray; }

   inline Array_double* GetQuarkChargeSquaredPartialSumSingle()
                           { return quarkChargeSquaredPartialSumSingleArray; }

   inline TwoD_Array_double* GetQuarkChargePartialSumDouble()
                                { return 
                                     quarkChargePartialSumDoubleArray; }

   inline TwoD_Array_double* GetQuarkChargeSquaredPartialSumDouble()
                                { return 
                                     quarkChargeSquaredPartialSumDoubleArray; }

   void Define(Global*);
   void Delete();

   double QSplitRegularEps0(SplittingFunction, double);
   double QSplitRegularEps1(SplittingFunction, double);
   double QSplitHalfRegularEps0(SplittingFunction, double);
   double QSplitHalfRegularEps1(SplittingFunction, double);
   double QSplitAt0(SplittingFunction);
   double QSplitAt1(SplittingFunction);
   double QSplitRegularIntEps0(SplittingFunction);
   double QSplitRegularIntEps1(SplittingFunction);

   double PSplitCounterTermDelta(
             SplittingFunction, 
             double lower
          );
   double PSplitCounterTermSubtracted( 
             SplittingFunction, 
             double u
          );
   double PSplitCounterTermRegular(
             SplittingFunction, 
             double u
          );
   double PSplitCounterTermDeltaNf(
             SplittingFunction
          );
};

// ----------------------------------------------------------------------------

class ElectroWeak
   : public virtual Global
   {

private:

   int defined;
   Global* global;

   AlphaEM* alphaEM_Server;  

protected:

public:

   inline AlphaEM* 
             Get_alphaEM_Server() { return alphaEM_Server; }

   ElectroWeak();
   ~ElectroWeak();

   void Define(Global *);
   void Delete();
};

// ----------------------------------------------------------------------------

class TheoryEvaluation
   : public virtual Global,
     public QCD,
     public ElectroWeak
   {

private:

   int defined;
   double GeV_2IntoPb;

   ParticleFreeList   *pFree;
   EventFreeListArray *eFree;

protected:

   ContributionFreeList      *cFree;
   ContributionArrayFreeList *caFree;

public:

   TheoryEvaluation();
   ~TheoryEvaluation();

   inline double GetGeV_2IntoPb() 
                    { return GeV_2IntoPb; }

   inline ContributionFreeList*      GetCFree()  { return cFree; }
   inline ContributionArrayFreeList* GetCAFree() { return caFree; }

   inline ParticleFreeList*          GetPFree()  { return pFree; }
   inline EventFreeListArray*        GetEFree()  { return eFree; }

   void Define();
   void Delete();
};

// ============================================================================

#endif // of QCDH

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
