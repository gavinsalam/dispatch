// ============================================================================
// 
// --> global include file: 
//     class definitions etc.
//
// file:              global.h
// created:           29.03.1997
// last modification: 02.12.1997
//
// ============================================================================
/*
   Note: ``Global'' is supposed to be a general class for supplying 
   true global objects, such as FILE pointers to output and error
   files, handling of exceptions (floating point exceptions)
   and other run-time errors, and for help in debugging (traps).

   For historical reasons the ``Global'' class at present still contains things
   that do not belong to it, such as Monte Carlo parameters, etc.
   They will later be put into the appropriate derived class
   (in this case, ``Disaster'').
*/
// ============================================================================

// load only once...
#ifndef GLOBALH
#define GLOBALH

// ============================================================================
//      
// --> symbols
//
// ============================================================================

#include <signal.h>
#include <iostream.h>
#include <fstream.h>

#include "standard.h"

#include "switch.h"
#include "enum.h"
#include "cmb.h"

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

class Permute;
class MapSquareLinear;
class Process;
class TrivialList;
class FPNaN;
class Global;
class TestProcess;
class String;
class MCFunction;
class MCUser;

// ============================================================================
//
// --> global variables
//
// ============================================================================

// global objects
extern Global* gl;

// floating point error handling
extern FPNaN *fpnan;

// ============================================================================
//
// --> functions
//
// ============================================================================

void errf(int level, const char* text);
void provoke_segmentation_fault();
void provoke_floating_point_exception();

// to catch floating point exceptions
int signal_fp_status();
void signal_fp_receive(int);
void signal_fp_install(void);   

// ============================================================================
//
// --> classes
//
// ============================================================================

class Disaster;

// ============================================================================

class GlobalError {

private:

public:

   GlobalError();
   ~GlobalError();
};

// ============================================================================

class FatalGlobalError
      : public GlobalError {

private:

   String* message;

public:

   FatalGlobalError(const String* message);
   FatalGlobalError(const char*   message);
   ~FatalGlobalError();

   void CreateAndSetMessage(const String* message);

   String* GetMessageString();
   char*   GetMessageChar();
};

// ----------------------------------------------------------------------------
//
// --> base class to check for already defined objects
//
// ----------------------------------------------------------------------------

class DefinedObject {

private:

   int isDefined;

   void SetDefined();

protected:

public:

   DefinedObject();
   ~DefinedObject();

   void ResetDefined();
   int IsDefined();
   int Status() { return isDefined; }
};

// ----------------------------------------------------------------------------
//
// --> named object
//
//     to be used as a base class for objects that carry a name
//     (for bookkeeping reasons)
//
// ----------------------------------------------------------------------------

class NamedObject {

private:

   String *name;
   char *noName;

protected:

public:

   NamedObject();
   ~NamedObject();

   void DefineName(char *nameIn); 
   void DefineName(const String *nameIn);
   void DefineName(const String &nameIn);
   char *GetName();
   char *GetNameSafe();
};

// ----------------------------------------------------------------------------

// class for Formatted Strings (used for printing formatted data via stream io)

class FS {

private:

   char *string;

public:

   // default constructor
   FS();
   ~FS();
   FS(const FS& fs);
   FS& operator=(const FS& fs);
   FS(char *format, const double d);
   FS(char *format, const int d);
   FS(char *format, const long d);
   FS(char *format, const unsigned long d);
   FS(char *format, const char *d);
   FS(char *format, const char d);
   void Print(FILE *);
   void CreateLargeString();
   char *GetString();
};

// ----------------------------------------------------------------------------

// operator overloading of the stream io output operator
// specify how to print a formatted quantity

ostream& operator << (ostream& s, FS fs);

// ----------------------------------------------------------------------------
//
// --> output stream for multiple destinations: can be set to
//
//     FILE*, ostream*, char*, void chout(char*)
//
//     required to have a synchronous output in conjunction
//     with FORTRAN programs
// 
// ----------------------------------------------------------------------------

enum Mstream_target {
        Mstream_None       = 0,
        Mstream_file       = 1,
        Mstream_ostream    = 2,
        Mstream_charString = 3,
        Mstream_function   = 4,
        Mstream_FORTRAN    = 5,
        Mstream_Mstream    = 6
};

// ----------------------------------------------------------------------------

class Mstream
   : public DefinedObject,
     public NamedObject
   {

private:

   Mstream_target target;

   FILE* fileOut;

   ostream* ostreamOut;

   char* charStringOut;
   int maxLength;  // this is the size of the string excluding the '\0'
   int position;

   void (*functionOut)(const char*);

   void (*FORTRAN_Out)(const char*, int, int);
   int unitNumber;
   int bufferSize;
   int bufferPos;
   char* buffer;

   Mstream* mstreamOut;

protected:

public:

   Mstream();
   ~Mstream();

   Mstream* Get_mstreamOut() { return mstreamOut; }

   void Associate(FILE* out);
   void Associate(ostream* out);
   void Associate(char* out, int maxLength_I); 
   void Associate(void (*out)(const char*));
   void Associate(void (*out)(const char*, int, int), int unitNumber_I);
   void Associate(Mstream* out);

   void CharStringOut(const char*);
   void Flush();
};

// ----------------------------------------------------------------------------

// operator overloading of the output operator
Mstream& operator << (Mstream& m, const char* s);

// specify how to print a formatted quantity
Mstream& operator << (Mstream& s, FS fs);


// ----------------------------------------------------------------------------
//
// --> base class for objects that have a pointer to Global
//
// ----------------------------------------------------------------------------

class GlobalObject {

private:

protected:

   Global* global;

public:

   inline void    Set_global(Global* in) { global = in; }
   inline Global* Get_global()           { return global; }

   GlobalObject();
   GlobalObject(Global* globalIn) { global = globalIn; }
   virtual ~GlobalObject();
};

// ----------------------------------------------------------------------------

class TrivialList
   : public GlobalObject,
     public NamedObject
   {

private:

   int identifier,maximum,actual;
   long *index,unique,trap;
   long *count;
   int is_soft_trap;
   long trap_check_count;

public:

   TrivialList(Global*);
   ~TrivialList();

   void Define(int,int);
   void Clear();
   long Unique();
   long Mark();
   int Check(long id, int situation);
   void Remove(long);
   void Print(FILE*);
   void SetTrap(long trap, int is_soft_trap, long trap_check_count);
};

// ----------------------------------------------------------------------------
//
// --> trap object
//
//     to be used as a base class for objects that can trap other objects
//
// ----------------------------------------------------------------------------

class TrapObject {

private:

protected:

   long trapCreate;
   long trapGetNoCreate;
   long trapReturn;

public:

   inline void SetTrapCreate(long in) { trapCreate = in; 
                                        Notify(in); }
   inline long GetTrapCreate() { return trapCreate; }

   inline void SetTrapGetNoCreate(long in) { trapGetNoCreate = in;
                                             Notify(in); }
   inline long GetTrapGetNoCreate() { return trapGetNoCreate; }

   inline void SetTrapReturn(long in) { trapReturn = in;
                                        Notify(in); }
   inline long GetTrapReturn() { return trapReturn; }

   TrapObject();
   ~TrapObject();

   void GetTrapsFrom(TrapObject*);
   void Notify(long);
};

// ----------------------------------------------------------------------------
//
// --> counter
//
// ----------------------------------------------------------------------------

class Counter {

private:

   long counter;
   long maximum;

protected:

public:

   Counter();
  ~Counter();

   inline void Set_maximum(long in) { maximum = in; }

   void Reset();
   Boolean Check();
   void    Increment();
   Boolean IncrementAndCheck();
   Boolean BorderlineCase();
};

// ----------------------------------------------------------------------------

class FPSignalHandler
   : public GlobalObject
   {

private:

   int flag;

   struct sigaction originalHandler_action;
   struct sigaction newHandler_action;

   void (*newHandler_function)(int);

   int counter;
   
public:

   int maxprint, nprint;
 
   FPSignalHandler(Global*);
   ~FPSignalHandler();

   int Status();
   void Set(int);
   void Reset();

   void Install(void (*)(int));
   void Reinstall_old();
   void Reinstall_new();

   void Increment_counter() { ++ counter; }
   int  Get_counter()       { return counter; }
   int  Get_flag()          { return flag; }
   struct sigaction * Get_p_originalHandler_action()
                         { return &originalHandler_action; }
};

// ----------------------------------------------------------------------------

class FPNaN_info {

private:

   char* no_message_errno; 
   char* no_message_mth; 

protected:

   int flag_SIGFPE;
   int counter_SIGFPE;

   int flag_errno;
   int value_errno;

   int flag_NaN;
   int counter_NaN;

   int flag_mth;
   int error_mth;
   int counter_mth;
   char* text_mth;

public:

   FPNaN_info();
   ~FPNaN_info();

   void Reset();

   void Copy_text_mth_From(const char*);

   const char* Get_msg_errno(int flag, int errno_in); 
   const char* Get_msg_mth(); 
   void PrintStatus(Mstream&);
   void CopyInto(FPNaN_info*);

   friend Boolean operator==(const FPNaN_info&, const FPNaN_info&);
};

Boolean operator==(const FPNaN_info&, const FPNaN_info&);

// ----------------------------------------------------------------------------

class FPNaN 
   : public GlobalObject,
     public FPNaN_info
   {

private:

   Boolean onHold;
   FPNaN_info saveStatus;

public:
 
   FPSignalHandler* fpsh;

   FPNaN(Global*, Boolean);
   ~FPNaN();

   void Modify_fpHandler(Boolean);

   inline int Get_flag_NaN()    {return flag_NaN; }

   inline int  Get_flag_mth()   { return flag_mth; }

   int   Get_flag_SIGFPE();
   int   Get_counter_SIGFPE();

   int   Get_flag_errno(); 
   int   Get_value_errno(); 

   void Set_flag_mth(int in, const char* text);   

   void Reset();
   void LogDoubleStatus(double);
   Boolean Get_Error();
   void PrintStatus(Mstream&);
   void PrintSaveStatus(Mstream& s) { saveStatus.PrintStatus(s); }

   void    Hold();
   Boolean Continue();

   int DoubleStatus(double);
   int DoubleStatusReset(double);
};

// ----------------------------------------------------------------------------
//
// --> global class: streams for status and error output
//
// ----------------------------------------------------------------------------

class Global {

private:

   Mstream statusStream;
   Mstream errorStream;
   Mstream logStream;

protected:

public:

   // --> I/O

   Mstream& To_status() { return statusStream; }
   Mstream& To_error()  { return errorStream; }
   Mstream& To_log()    { return logStream; }

   // --> floating point handler

   Boolean if_fpHandler;

   // public for efficiency reasons:

   // pointer to the function to be integrated by the Monte Carlo method
   MCFunction* mc_function;

   double VEGAS_weight_current; // VEGAS weight of an event
   int    VEGAS_pts_in_last_iteration; 

   // global objects 
   MapSquareLinear *pMapSquareLinear;
   Permute *pPermute;

   // variables

   int finalrun;
   Process *process;
   MCUser *user;

   Disaster *disaster;

   // for reading trap data
   int trap_set;
   int is_soft_trap;
   long trap_check_count;

   // input parameters for test purposes
   double rho_lim;        // regularizing factor for check of limits
   double alpha_lim;      // exponent for test
   int ia,ib;             // for limit checks
   int itest1,itest2,itest3; // integers for test purposes
   int mel_index;         // matrix element to be integrated
   int process_index;     // process under consideration, partly fixes MEs
   int numberOfFinalStatePartonsInBornTerm; // the number of final-state
                                            // partons of the Born term
   int subtraction_type;  // type of subtraction procedure
   int icatch, jcatch, tcatch;
   int i0,i1,i2,i3,i2p,i3p;

   // variables for tests
   int istore, jstore, tstore;
   LimitType lstore;

   TrivialList** tlist;
   TrivialList*  tEvent;
   TrivialList*  tParticle;

   TrivialList* idEvent; // to identify individual event records

   Global();
   virtual ~Global();
};

// ============================================================================

#endif // of GLOBALH

// ============================================================================
//
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
