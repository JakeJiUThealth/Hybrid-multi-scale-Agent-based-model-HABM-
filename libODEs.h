//
// MATLAB Compiler: 4.17 (R2012a)
// Date: Wed Dec 17 11:27:12 2014
// Arguments: "-B" "macro_default" "-W" "cpplib:libODEs" "-T" "link:lib"
// "Suv_Adh_MIC.m" "cal_stiffness.m" 
//

#ifndef __libODEs_h
#define __libODEs_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#include "mclcppclass.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libODEs
#define PUBLIC_libODEs_C_API __global
#else
#define PUBLIC_libODEs_C_API /* No import statement needed. */
#endif

#define LIB_libODEs_C_API PUBLIC_libODEs_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libODEs
#define PUBLIC_libODEs_C_API __declspec(dllexport)
#else
#define PUBLIC_libODEs_C_API __declspec(dllimport)
#endif

#define LIB_libODEs_C_API PUBLIC_libODEs_C_API


#else

#define LIB_libODEs_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libODEs_C_API 
#define LIB_libODEs_C_API /* No special import/export declaration */
#endif

extern LIB_libODEs_C_API 
bool MW_CALL_CONV libODEsInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_libODEs_C_API 
bool MW_CALL_CONV libODEsInitialize(void);

extern LIB_libODEs_C_API 
void MW_CALL_CONV libODEsTerminate(void);



extern LIB_libODEs_C_API 
void MW_CALL_CONV libODEsPrintStackTrace(void);

extern LIB_libODEs_C_API 
bool MW_CALL_CONV mlxSuv_Adh_MIC(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

extern LIB_libODEs_C_API 
bool MW_CALL_CONV mlxCal_stiffness(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__BORLANDC__)

#ifdef EXPORTING_libODEs
#define PUBLIC_libODEs_CPP_API __declspec(dllexport)
#else
#define PUBLIC_libODEs_CPP_API __declspec(dllimport)
#endif

#define LIB_libODEs_CPP_API PUBLIC_libODEs_CPP_API

#else

#if !defined(LIB_libODEs_CPP_API)
#if defined(LIB_libODEs_C_API)
#define LIB_libODEs_CPP_API LIB_libODEs_C_API
#else
#define LIB_libODEs_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_libODEs_CPP_API void MW_CALL_CONV Suv_Adh_MIC(int nargout, mwArray& MIC_Suv_rate, mwArray& MIC_Adh_rate, const mwArray& st, const mwArray& BTZ);

extern LIB_libODEs_CPP_API void MW_CALL_CONV cal_stiffness(int nargout, mwArray& st, const mwArray& types, const mwArray& SDF1);

#endif
#endif
