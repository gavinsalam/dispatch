// ============================================================================
// 
// --> string class include file
//
// file:              strng.h
// created:           29.03.1997
// last modification: 11.12.1997
//
// ============================================================================

// load only once...
#ifndef STRNG_H
#define STRNG_H

// ============================================================================

#include "switch.h"
//#include "global.h"

//#include "container.h"

// ============================================================================
//
// --> definitions
//
// ============================================================================

class Global;
class String;

// ============================================================================
//
// --> functions
//
// ============================================================================

long   ReadLong(Global *global,char *pname,char *ins);
int    ReadInt(Global *global,char *pname,char *ins);
int    ReadBool(Global *global,char *pname,char *ins);
double ReadDouble(Global *global,char *pname,char *ins);
void   ReadString(Global *global,char *pname,char *ins,String *str);
int    SkipWhiteSpace(FILE *in);

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

class String {

private:

protected:

public:

   char *data;
   int maxlength;

   String(int);
   String (const char*); // copy constructor char -> String
   String (const String*, const String*); // ``concatenation'' constructor
   String (const String&, const String&); // ``concatenation'' constructor
   ~String();

   char *GetData() const { return data; };

   void Define(int);
   void Define(const String*, const String*);
   void Define(const String&, const String&);
   void Delete();

   void Print(FILE*);
   int  DetermineLength() const;   
   void CopyFromChar(const char *);
   void AppendFromChar(const char *);
   void DefineAndCopyFromChar(const char *);
   void CopyFromString(const String *);
   void AppendFromString(const String *);
   void DefineAndCopyFromString(const String *);
   void DefineAndCopyFromString(const String &);
   void BackslashToNewline();
   void ReadBoundedFromFile(FILE *);
   void ReadNonWhiteFromFile(FILE *);
};

// ============================================================================

#endif // of STRNG_H

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
