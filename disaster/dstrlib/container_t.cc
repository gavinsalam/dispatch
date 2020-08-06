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
