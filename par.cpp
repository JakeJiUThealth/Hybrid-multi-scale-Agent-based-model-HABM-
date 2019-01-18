#include "par.hpp"

using namespace std;

namespace par {
  int iMaxTime = 10;
  double iA = 0;  //BTZ
  double iB = 0;    //LEN
  double iC = 0;    //Third drug
  string sParFile = "";
  string sSeedFile = "";
  string sOutDir = "Data";
  bool bPrintCyt = false;
  string sM = "M";
  int iT = -1;
  bool bP = false;
  int iP_Start = 100; 
  int iP_End = 200; 

  double iniSDF1=0.00535; 
  double inistiffness=400;
  
  double basicSDF1=0.0012;//0.00535;
  double basicstiffness=250;  //default for mBMSC
  double maxstiffness=530;

  
  int sw_mic=1;
  int sw_mm=1;
  int sw_cd8=0;
  int sw_treg=0;
  
}

