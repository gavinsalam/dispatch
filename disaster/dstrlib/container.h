// ============================================================================
// 
// --> container classes
//
// file:              container.h
// created:           10.06.1997
// last modification: 16.11.1997
//
// ============================================================================

// load only once...
#ifndef CONTAINER_H
#define CONTAINER_H

// ============================================================================

#include "global.h"
#include "switch.h"

// ============================================================================
//
// --> global variables
//
// ============================================================================

extern long uCount;

extern long coTest;

// ============================================================================
//
// --> classes, templates
//
// ============================================================================

// ----------------------------------------------------------------------------
//
// --> array template
//
// ----------------------------------------------------------------------------

template <class T>
class Array {

private:

protected:

   T *data;
   int rmin, rmax;

public:

   Array();
   Array(int rmin, int rmax);
   Array(int length);
   ~Array();

   inline T   *GetDataPtr() 
                  { return data; }

   inline T    GetData(int i) 
                  { return data[i]; }

   inline void SetData(int i, T dataIn)  // :_MOD_: should be made safe
                  { data[i]=dataIn; }

   inline int  GetRmin() 
                  { return rmin; }

   inline int  GetRmax() 
                  { return rmax; }

   INLINE void Define(int, int);
   inline void Define(int size) { Define(0, size-1); }

   INLINE void Delete();
   INLINE void CreateSame(Array <T> *);
   INLINE void CopyInto(Array <T> *);
   INLINE void DeleteAllEntries();
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_Array(T) \
           instantiate_Template_1(Array, T)

// ----------------------------------------------------------------------------
//
// --> array of pointers
//
// ----------------------------------------------------------------------------

template <class T>
class ArrayOfPointers : public Array <T> {

public:

   ArrayOfPointers() {};
   ArrayOfPointers(int rmin, int rmax);
   ~ArrayOfPointers() {};

   INLINE void DeleteAllEntries();
   INLINE void SetToNull();
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_ArrayOfPointers(T) \
           I_Array(T) \
           instantiate_Template_1(ArrayOfPointers, T)

// ----------------------------------------------------------------------------
//
// --> array of double
//
// ----------------------------------------------------------------------------

class Array_double : public Array <double> {

public:

   Array_double();
   Array_double(int rsize);
   Array_double(int rmin, int rmax);
   ~Array_double();

   void SetToZero();
   void LinearScale(double smin, double smax);
   void LogScale(double smin, double smax);
   int FindBinLocation(double, int&);
   void Print(Mstream&, int, int);
   void Print(FILE*, int, int);
   inline void Print(Mstream& s) { Print(s, GetRmin(), GetRmax()); }
   inline void Print(FILE* f) { Print(f, GetRmin(), GetRmax()); }
   inline void Print() { Print(stdout); }
};

// ----------------------------------------------------------------------------
//
// --> array of ints
//
// ----------------------------------------------------------------------------

class Array_int : public Array <int> {

public:

   Array_int();
   Array_int(int rsize);
   Array_int(int rmin, int rmax);
   ~Array_int();

   void SetToZero();
   void Print(Mstream&, int, int);
   void Print(FILE*, int, int);
   inline void Print(Mstream& s) { Print(s, GetRmin(), GetRmax()); }
   inline void Print(FILE *f) { Print(f, GetRmin(), GetRmax()); }
   inline void Print() { Print(stdout); }
};

// ----------------------------------------------------------------------------
//
// --> two-dimensional array
//
// ----------------------------------------------------------------------------

template <class T>
class TwoD_Array {

private:

protected:

   T **data;

   int r0min, r0max;
   int r1min, r1max;

public:

   TwoD_Array();
   TwoD_Array(int r0minIn, int r0maxIn,
              int r1minIn, int r1maxIn);
   ~TwoD_Array();

   inline T   **GetDataPtr() 
                  { return data; }

   inline T    GetData(int i, int j) 
                  { return data[i][j]; }

   inline void SetData(int i, int j, T dataIn) 
                  { data[i][j]=dataIn; }

   inline void GetBounds(int &r0minOut, int &r0maxOut,
                         int &r1minOut, int &r1maxOut) 
                  { r0minOut = r0min; r0maxOut = r0max; 
                    r1minOut = r1min; r1maxOut = r1max; }

   INLINE void Define(int r0minIn, int r0maxIn,
                      int r1minIn, int r1maxIn);

   INLINE void Define(int r0size,
                      int r1size) { Define(0, r0size-1, 0, r1size-1); }

   INLINE void Delete();

   INLINE T **GetDataPtrWithNullWarning(); 
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_TwoD_Array(T) \
           instantiate_Template_1(TwoD_Array, T)

// ----------------------------------------------------------------------------
//
// two-dimensional array of pointers
//
// ----------------------------------------------------------------------------

template <class T>
class TwoD_ArrayOfPointers : public TwoD_Array <T> {

public:

   TwoD_ArrayOfPointers() {};
   TwoD_ArrayOfPointers(int r0min_In, int r0max_In,
                        int r1min_In, int r1max_In)
      : TwoD_Array <T> (r0min_In, r0max_In, 
                        r1min_In, r1max_In) {};
   ~TwoD_ArrayOfPointers() {};

   INLINE void DeleteAllEntries();
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_TwoD_ArrayOfPointers(T) \
           I_TwoD_Array(T) \
           instantiate_Template_1(TwoD_ArrayOfPointers, T)

// ----------------------------------------------------------------------------
//
// --> two-dimensional array of double
//
// ----------------------------------------------------------------------------

class TwoD_Array_double : public TwoD_Array <double> {

public:

   TwoD_Array_double();
   TwoD_Array_double(int r0min, int r0max,
                     int r1min, int r1max
                    );
   ~TwoD_Array_double();

   INLINE void SetToZero();
};

// ----------------------------------------------------------------------------
//
// --> Free list container
//
// ----------------------------------------------------------------------------

template <class T>
class FreeListContainer {

private:

   T *data;
   FreeListContainer <T>  *next;

   long ident;

protected:

public:

   FreeListContainer();
  ~FreeListContainer();

   inline T *GetDataPtr() { return data; }
   inline void SetDataPtr( T *dataIn ) { data = dataIn; }

   inline void SetNext( FreeListContainer <T> *nextIn ) { next = nextIn; }
   inline FreeListContainer <T> *GetNext() { return next; }

   inline long GetIdent() { return ident; }
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_FreeListContainer(T) \
           instantiate_Template_1(FreeListContainer, T)

// ----------------------------------------------------------------------------
//
// --> Base class for free lists: common ``return object''
//
// ----------------------------------------------------------------------------

template <class T>
class CommonFreeList {

private:

protected:

public:

   CommonFreeList() {}
   virtual ~CommonFreeList() {}

   virtual void ReturnObject(T*) = 0;  // this is to be replaced appropriately
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_CommonFreeList(T) \
           instantiate_Template_1(CommonFreeList, T)

// ----------------------------------------------------------------------------
//
// --> Simple free list (no UseList can be destructed)
//
// ----------------------------------------------------------------------------

template <class T>
class SimpleFreeList
   : public CommonFreeList <T>,
     public virtual GlobalObject,
     public NamedObject, 
     public TrapObject
   {

private:

   int defined; 

   FreeListContainer <T> *listFull;
   FreeListContainer <T> *listEmpty;

   long createCount;
   long containerCount;

   long listFullLength;
   long listEmptyLength;

   long getCount;
   long returnCount;

   void Delete();

protected:

   void Define(Global*);

public:

   SimpleFreeList(Global*);
  ~SimpleFreeList();

   SimpleFreeList <T> & operator=(const SimpleFreeList <T> &);

   void Reset();

   int GetObject(T* &);
   void ReturnObject(T*);

   void Print(FILE*);
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_SimpleFreeList(T) \
           I_FreeListContainer(T) \
           I_CommonFreeList(T) \
           instantiate_Template_1(SimpleFreeList, T)

// ----------------------------------------------------------------------------
//
// --> UseList container
//
// ----------------------------------------------------------------------------

template <class T>
class UseListContainer {
   
private:

   T *data;
   UseListContainer <T> *next;

   long ident;

protected:

public:
   
   UseListContainer();
  ~UseListContainer();

   inline T *GetDataPtr() { return data; }
   inline void SetDataPtr( T *dataIn ) { data = dataIn; }
   
   inline void SetNext( UseListContainer <T> *nextIn ) { next = nextIn; }
   inline UseListContainer <T> *GetNext() { return next; }

   inline long GetIdent() { return ident; }
   void ReturnToFreeList();
}; 
   
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_UseListContainer(T) \
           instantiate_Template_1(UseListContainer, T)

// ----------------------------------------------------------------------------
//
// --> HorizontalList
//
// ----------------------------------------------------------------------------

template <class T>
class HorizontalList {

private:

   HorizontalList <T> * hor;
   T* object;
   HorizontalList <T> * (T :: *getH)(); 

protected:

public:

   HorizontalList(T*, HorizontalList <T> * (T :: *)());
  ~HorizontalList();

   void Reset();

   inline void SetHorizontalPtr(HorizontalList <T> * in) { hor = in; }
   inline HorizontalList <T> * GetHorizontalPtr() { return hor; }
   inline HorizontalList <T> * (T :: *GetGetH())()
             {  return getH; } 

   inline void SetHorizontalObject(T* tIn) { 
                     hor = (tIn ->* getH)();
               }

   inline T* GetHorizontalObject() { return (hor == NULL) 
                                               ? (T*) NULL 
                                               : hor -> GetObject(); }
   inline T* GetObject() { return object; }

   void InsertAfterFirstElement(HorizontalList <T> *);   
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_HorizontalList(T) \
           instantiate_Template_1(HorizontalList, T)

// ----------------------------------------------------------------------------
//
// --> UseList
//
// ----------------------------------------------------------------------------

template <class T>
class UseList 
   : public virtual GlobalObject,
     public NamedObject 
   {

private:

   UseListContainer <T> * list;

   UseListContainer <T> * tVertical;
   UseListContainer <T> * tVerticalNext;
   T *tHorizontal;

   HorizontalList <T> * (T :: *hFunction)();

   SimpleFreeList < UseListContainer <T> > *freeList;

protected:

public:

   UseListContainer <T> * GetPtrToFirstContainer() { return list; }

   UseList(Global*);
   ~UseList();

   int IsNotEmpty() { return (list == NULL) ? 0 : 1; }

   void InsertObject(T*);  // insert in the beginning
   void AppendObject(T*);  // append at the end
   T* GetFirstObject();    // removes the first object and returns it
   T* GetLastObject();     // removes the last object and returns it

   // return containers to the free list, do not touch data
   void DestructContainerList();         

   void Print(FILE*);
   void PrintIdent(
           FILE*, 
           HorizontalList <T> * (T :: *hFunctionIn)()
        );

   // find an object with the same parameters (defined by T::UseListIsEqual())
   // and return its pointer; if none exists, return NULL; 
   // scans the vertical list only
   T* FindSame(T*);

   // find an object with the given ``ident''
   // and return its pointer; if none exists, return NULL;
   // also scans the horizontal list; 
   // uses the ``Traverse...'' methods
   T* FindIdent(
         HorizontalList <T> * (T :: *hFunctionIn)(),
         long
      );

   // methods to transverse globally all elements of a UseList, 
   // including the horizontal part
   void TraverseReset(HorizontalList <T> * (T :: *hFunctionIn)());
   T* TraverseNext();
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_UseList(T) \
           I_UseListContainer(T) \
           I_SimpleFreeList(UseListContainer <T>) \
           I_HorizontalList(T) \
           instantiate_Template_1(UseList, T)

// ----------------------------------------------------------------------------
//
// --> Free list
//
// ----------------------------------------------------------------------------

template <class T>
class FreeList
   : public SimpleFreeList <T>
   {

private:

protected:

public:

   FreeList(Global*);
  ~FreeList();

   void ReturnObjectsFromUseList(
           HorizontalList <T> * (T :: *hFunctionIn)(),
           UseList <T> *
        );

   void ReturnObjectsFromHorizontalList(
           HorizontalList <T> * (T :: *hFunctionIn)(),
           T*
        );
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_FreeList(T) \
           I_SimpleFreeList(T) \
           I_UseList(T) \
           instantiate_Template_1(FreeList, T)

// ----------------------------------------------------------------------------
//
// Free list: no explicit parameter
//
// ----------------------------------------------------------------------------

template <class T>
class FreeList_0 : public FreeList <T> {

private:

protected:

public:

   FreeList_0(Global*);
  ~FreeList_0();

   INLINE T *GetDefinedObject();
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_FreeList_0(T) \
           I_FreeList(T) \
           instantiate_Template_1(FreeList_0, T)

// ----------------------------------------------------------------------------
//
// Free list: one explicit parameter
//
// ----------------------------------------------------------------------------

template <class T, class U>
class FreeList_1 : public FreeList <T> {

private:

protected:

public:

   FreeList_1(Global*);
  ~FreeList_1();

   INLINE T *GetDefinedObject(U);
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_FreeList_1(T, U) \
           I_FreeList(T) \
           instantiate_Template_2(FreeList_1, T, U)

// ----------------------------------------------------------------------------
//
// Array of free lists
//
// ----------------------------------------------------------------------------

template <class T>
class FreeListArray
   : public CommonFreeList <T>,
     public virtual GlobalObject,
     public NamedObject,
     public TrapObject
   {

private:

   int defined;

   ArrayOfPointers <FreeList <T> *> *flist;
   FreeList <T> **flptr;

   int rmin, rmax;

protected:

public:

   FreeListArray(Global*);
   FreeListArray(Global*, int rlength);
   FreeListArray(Global*, int rmin, int rmax);
  ~FreeListArray();

   FreeListArray <T> & operator=(const FreeListArray <T> &);

   void Define(int rmin, int rmax);
   void Delete();

   FreeList <T> *CreateNewFreeList(int);
   int  GetObject(T* &, int par);
   void ReturnObject(T*);

   void Print(FILE *);

   void ReturnObjectsFromUseList(
           HorizontalList <T> * (T :: *hFunctionIn)(),
           UseList <T> *
        );

   void ReturnObjectsFromHorizontalList(
           HorizontalList <T> * (T :: *hFunctionIn)(),
           T* 
        );
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_FreeListArray(T) \
           I_FreeList(T) \
           I_ArrayOfPointers(FreeList <T> *) \
           instantiate_Template_1(FreeListArray, T)

// ----------------------------------------------------------------------------
//
// Two-dimensional array of free lists
//
// ----------------------------------------------------------------------------

template <class T>
class FreeListTwoD_Array
   : public CommonFreeList <T>,
     public virtual GlobalObject,
     public NamedObject,
     public TrapObject
   {

private:

   int defined;

   TwoD_ArrayOfPointers <FreeList <T> *> *flist;
   FreeList <T> ***flptr;

   int r0min, r0max;
   int r1min, r1max;

protected:

public:

   FreeListTwoD_Array(Global*);
   FreeListTwoD_Array(
      Global*,
      int r0min, int r0max,
      int r1min, int r1max
   );
   FreeListTwoD_Array(
      Global*,
      int r0length, int r1length
   );
  ~FreeListTwoD_Array();

   FreeListTwoD_Array <T> & operator=(const FreeListTwoD_Array <T> &);

   void Define(int r0min, int r0max,
               int r1min, int r1max
        );
   void Delete();

   FreeList <T> *CreateNewFreeList(int, int);
   int  GetObject(T* &, int par0, int par1);
   void ReturnObject(T*);

   void Print(FILE *);

   void ReturnObjectsFromUseList(
           HorizontalList <T> * (T :: *hFunctionIn)(),
           UseList <T> *
        );

   void ReturnObjectsFromHorizontalList(
           HorizontalList <T> * (T :: *hFunctionIn)(),
           T* 
        );
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define I_FreeListTwoD_Array(T) \
           I_FreeList(T) \
           I_TwoD_ArrayOfPointers(FreeList < T > *) \
           instantiate_Template_1(FreeListTwoD_Array, T)

// ----------------------------------------------------------------------------
//
// Array of free lists: no explicit parameter
//
// ----------------------------------------------------------------------------

template <class T>
class FreeListArray_0 : public FreeListArray <T> {

private:

protected:

public:

   FreeListArray_0(Global*, int rminIn, int rmaxIn);
  ~FreeListArray_0();

   INLINE T *GetDefinedObject(int);
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
#define I_FreeListArray_0(T) \
           I_FreeListArray(T) \
           instantiate_Template_1(FreeListArray_0, T)
           
// ----------------------------------------------------------------------------
//
// Array of free lists: one explicit parameter of type U
//
// ----------------------------------------------------------------------------

template <class T, class U>
class FreeListArray_1 : public FreeListArray <T> {

private:

protected:

public:

   FreeListArray_1(Global*, int rlength);
   FreeListArray_1(Global*, int rminIn, int rmaxIn);
  ~FreeListArray_1();

   INLINE T *GetDefinedObject(U, int);
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
#define I_FreeListArray_1(T, U) \
           I_FreeListArray(T) \
           instantiate_Template_2(FreeListArray_1, T, U)
           
// ----------------------------------------------------------------------------
//
// Two-dimensional array of free lists: one explicit parameter of type U
//
// ----------------------------------------------------------------------------

template <class T, class U>
class FreeListTwoD_Array_1 
   : public FreeListTwoD_Array <T> 
   {

private:

protected:

public:

   FreeListTwoD_Array_1(Global*);
   FreeListTwoD_Array_1(
      Global*, 
      int r0minIn, int r0maxIn,
      int r1minIn, int r1maxIn
   );
   FreeListTwoD_Array_1(
      Global*,
      int r0length, int r1length
   );
  ~FreeListTwoD_Array_1();

   INLINE T *GetDefinedObject(U, int par1, int par2);
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
#define I_FreeListTwoD_Array_1(T, U) \
           I_FreeListTwoD_Array(T) \
           instantiate_Template_2(FreeListTwoD_Array_1, T, U)
           
// ----------------------------------------------------------------------------
//
// Cache container
//
// ----------------------------------------------------------------------------

template <class T>
class CacheContainer
   : public DefinedObject 
   {

private:

   T* data;

   int nInt; 
   int nDouble;

   Array_int* intKey;
   Array_double* doubleKey;
   
   int *intKeyData;
   double *doubleKeyData;

   double hashValue;

   int occupied;

   long time;
   long reused;

   long freq;

protected:

public:

   void SetDataPtr(T* in) { data = in; }
   T* GetDataPtr() { return data; }

   void SetHashValue(double in) { hashValue = in; }
   double GetHashValue() { return hashValue; }

   void SetOccupied(int in) { occupied = in; }
   int GetOccupied() { return occupied; }

   void SetTime(long in) { time = in; }
   long GetTime() { return time; }

   void SetReused(long in) { reused = in; }
   long GetReused() { return reused; }
   void IncrementReused() { reused++; }

   void SetFreq(long in) { freq = in; }
   long GetFreq() { return freq; }
   void IncrementFreq() { freq++; }

   Array_int* GetIntKeyPtr() { return intKey; }
   Array_double* GetDoubleKeyPtr() { return doubleKey; }

   CacheContainer();
   ~CacheContainer();

   void Define(int nInt, int nDouble);
   void Delete();
   void Reset();

   void CopyKeyInto(CacheContainer <T> *);
   
   int CompareHashValues(double, int, double);
   int CompareKeys(int*, double*, int, double);
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
#define I_CacheContainer(T) \
           instantiate_Template_1(CacheContainer, T)
           
// ----------------------------------------------------------------------------
//
// Cache
//
// ----------------------------------------------------------------------------

template <class T>
class Cache 
   : public GlobalObject,
     public DefinedObject, 
     public NamedObject 
   {

private:
 
   CacheContainer <T> * store;

   long clock;
   long survivalTime;
   long reused;
   long numberOfObjects;
   long unsuccessful; 
   long numberOfAccesses;

   int nCacheSize;
   int nInt; 
   int nDouble;

   Array_int* intKey;
   Array_double* doubleKey;

   double logBase;
   double logBaseLog;
   double oneOverLogBaseLog;
   double integerFactor;

   int critHash;
   double accuracyHash;   
 
   int critDouble;
   double accuracyDouble;   

   HorizontalList <T> * (T :: *getH)();
   CommonFreeList <T> * freeList;

protected:

public:

   void SetFreeList(CommonFreeList <T> * in) { freeList = in; } 

   Cache(Global*, int nCacheSize, int nInt, int nDouble);
   ~Cache();

   void Define(int nCacheSize, int nInt, int nDouble);
   void Delete();
   void Reset();

   void StatisticsForAnObject(CacheContainer <T> *);

   void StoreInCache(Array_int*, Array_double*, T*);
   T* FindInCache(Array_int*, Array_double*);

   int GetHashValue(int*, double*, double*);
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
#define I_Cache(T) \
           I_CacheContainer(T) \
           instantiate_Template_1(Cache, T)
           
// ----------------------------------------------------------------------------

#if 0
// multi-dimensional array
template <class T>
class MDArray {

private:

   MDArray <T> *nextArray;
   T *data;
   int rmin, rmax;

protected:

public:

   MDArray();
   ~MDArray();

   inline T   *GetDataPtr() 
                  { return data; }

   inline T    GetData(int i) 
                  { return data[i]; }

   inline void SetData(int i, T dataIn) 
                  { data[i]=dataIn; }

   inline int  GetRmin() 
                  { return rmin; }

   inline int  GetRmax() 
                  { return rmax; }

   INLINE void Define(int, int);
   INLINE void Delete();
   INLINE void CreateSame(Array <T> *);
   INLINE void CopyInto(Array <T> *);
};
#endif

// ============================================================================

// also include the templates
#include "container_t.cc"

// ============================================================================

#endif // of CONTAINER_H

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
