#ifndef PAR_H_
#define PAR_H_

#include <string>
using namespace std;

namespace par {
  extern int iMaxTime;
  extern double iA;
  extern double iB;
  extern double iC;
  extern string sParFile;
  extern string sSeedFile;
  extern string sOutDir;
  extern bool bPrintCyt;
  extern string sM;
  extern int iT;
  extern bool bP; 
  extern int iP_Start; 
  extern int iP_End; 
  extern double basicstiffness;
  extern double inistiffness;
  extern double maxstiffness;
  extern double iniSDF1;
  extern double basicSDF1;

  extern int sw_mic;
  extern int sw_mm;
  extern int sw_cd8;
  extern int sw_treg;
}

#endif /*PAR_H_*/
