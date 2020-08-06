// ============================================================================
// 
// --> enumerations
//
// file:              enum.h
// created:           29.03.1997
// last modification: 17.09.1997
//
// ============================================================================

#ifndef ENUMH
#define ENUMH

// ============================================================================

enum particles {
        Undefined='u',
        Quark='q',
        Antiquark='a',
        Gluon='g'  
     };

// ----------------------------------------------------------------------------

enum PSMode {
        ALLUNIT=0,
        FIRSTEXPLICIT=1
     };

// ----------------------------------------------------------------------------

enum P0_Reference {
        IS_REFERENCE=0,
        FS_REFERENCE=1
     };

extern const char* P0_ReferenceName[];

// ----------------------------------------------------------------------------

enum LimitType {
        NO_LIMIT = 0,
        SOFT = 1,
        COLLINEAR = 2,
        SOFT_AND_COLLINEAR = 3,
        SOFT_ADDED = 4,
        COLLINEAR_ADDED = 5,
        SOFT_AND_COLLINEAR_ADDED = 6
     };

extern const char* LimitTypeName[];

// ----------------------------------------------------------------------------

enum Frame {
        none = 0,
        original = 1,
        pCMS=2,
        hCMS=3,
        breit = 4,
        lab = 5
     };

extern const char* FrameName[];

// ----------------------------------------------------------------------------

enum ClusterAlgorithm {
        NONE=0,
        LORENTZ_INVARIANT_1=1,
        LORENTZ_INVARIANT_OLD_1=2
     };

// ----------------------------------------------------------------------------

enum ClusterAlgorithmType {
        JADE_algorithm=1,
        kT_algorithm=2
     };

// ----------------------------------------------------------------------------

enum RecombinationType {
        JADE_Scheme=1,
        E_Scheme=2,
        E0_Scheme=3,
        P_Scheme=4, 
        P0_Scheme=5
     };

// ----------------------------------------------------------------------------

enum BookMode {
        Histogram = 1,
        Graph = 2
     };

// ----------------------------------------------------------------------------

enum ScaleType {
        Arbitrary = 0,
        Linear = 1,
        Log = 2
     };

// ----------------------------------------------------------------------------

enum DISSubProcess {
        NoDISSubProcess=0,
        qXq    =  1,
        qXqg   =  2,
        gXqq   =  3,
        qXqgg  =  4,
        gXqqg  =  5,
        qXqqq  =  6,
        qXqgV  =  7,
        gXqqV  =  8,
        qXqV   =  9,
        qXqgA  = 10,
        gXqqA  = 11,
        qXqggA = 12,
        gXqqgA = 13,
        qXqqqA = 14,
        qXqgRS = 15,
        gXqqRS = 16,
        qXqFS  = 17,
        qXqgFS = 18,
        gXqqFS = 19,
        qXqgRSNF = 20,
        gXqqRSNF = 21,
        gXqqFSNF = 22,
        testDIS = 100
     };

// ----------------------------------------------------------------------------

enum ScaleLogarithmFactor {
        NoScaleLogarithm=0,
        RenormalizationScaleLogarithm=1,
        RenormalizationScaleLogarithmWithNf=2,
        FactorizationScaleLogarithm=3,
        FactorizationScaleLogarithmWithNf=4
     };

// ----------------------------------------------------------------------------

enum SubtractionType {
        NoSubtraction = 0,
        Subtraction = 1,
        AddedSubtraction = 2, 
        CollinearInitial = 3, 
        Collinear = 4
     };

extern const char* SubtractionTypeName[];   

// ----------------------------------------------------------------------------

enum EventType {
        EventTypeNone = 0,
        Additional = 1, 
        Final = 2
     };

extern const char* EventTypeName[];   

// ----------------------------------------------------------------------------

enum SplittingFunction {
        NoSplittingFunction = 0,
        Q_q_from_q = 1,
        Q_G_from_q = 2,
        Q_q_from_G = 3,
        Q_G_from_G = 4
     };

extern const char *QName[];

// ============================================================================

#endif // of ENUMH

// ============================================================================
//
// --> End of file.
//
// ============================================================================

// 01234567890123456789012345678901234567890123456789012345678901234567890
