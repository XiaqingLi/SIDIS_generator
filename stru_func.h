#ifndef _STRU_FUNC_H
#define _STRU_FUNC_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"
#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/SpecFuncMathCore.h"

#include "LHAPDF/LHAPDF.h"

#include "constant.h"

namespace StructureFunction{
  LHAPDF::PDF * pdfs = LHAPDF::mkPDF("CJ15lo", 0);
  LHAPDF::PDF * ffs;
  
  double getFUUT_p_TL(const double x, const double z, const double Q2, const double Pt, const int hid){
    if(hid >= 0)
      ffs = LHAPDF::mkPDF("DSSFFlo", hid);
    else
      ffs = LHAPDF::mkPDF("DSSFFlo", -hid+1000);
    double f1[6], D1[6];
    f1[0] = pdfs->xfxQ2(2, x, Q2) / x;
    f1[1] = pdfs->xfxQ2(1, x, Q2) / x;
    f1[2] = pdfs->xfxQ2(3, x, Q2) / x;
    f1[3] = pdfs->xfxQ2(-2, x, Q2) / x;
    f1[4] = pdfs->xfxQ2(-1, x, Q2) / x;
    f1[5] = pdfs->xfxQ2(-3, x, Q2) / x; 
    D1[0] = ffs->xfxQ2(2, z, Q2) / z;//u
    D1[1] = ffs->xfxQ2(1, z, Q2) / z;//d
    D1[2] = ffs->xfxQ2(3, z, Q2) / z;//s
    D1[3] = ffs->xfxQ2(-2, z, Q2) / z;//ubar
    D1[4] = ffs->xfxQ2(-1, z, Q2) / z;//dbar
    D1[5] = ffs->xfxQ2(-3, z, Q2) / z;//sbar
    double PhT2 = z * z * 0.57 + 0.12;
    return x * (eu2 * (f1[0] * D1[0] + f1[3] * D1[3]) + ed2 * (f1[1] * D1[1] + f1[2] * D1[2] + f1[4] * D1[4] + f1[5] * D1[5])) * exp(-Pt * Pt / PhT2) / (PI * PhT2);    
  }

  double getFUUT_n_TL(const double x, const double z, const double Q2, const double Pt, const int hid){
    if(hid >= 0)
      LHAPDF::PDF * ffs = LHAPDF::mkPDF("DSSFFlo", hid);
    else
      LHAPDF::PDF * ffs = LHAPDF::mkPDF("DSSFFlo", -hid+1000);
    double f1[6], D1[6];
    f1[0] = pdfs->xfxQ2(2, x, Q2) / x;
    f1[1] = pdfs->xfxQ2(1, x, Q2) / x;
    f1[2] = pdfs->xfxQ2(3, x, Q2) / x;
    f1[3] = pdfs->xfxQ2(-2, x, Q2) / x;
    f1[4] = pdfs->xfxQ2(-1, x, Q2) / x;
    f1[5] = pdfs->xfxQ2(-3, x, Q2) / x; 
    D1[0] = ffs->xfxQ2(2, z, Q2) / z;//u
    D1[1] = ffs->xfxQ2(1, z, Q2) / z;//d
    D1[2] = ffs->xfxQ2(3, z, Q2) / z;//s
    D1[3] = ffs->xfxQ2(-2, z, Q2) / z;//ubar
    D1[4] = ffs->xfxQ2(-1, z, Q2) / z;//dbar
    D1[5] = ffs->xfxQ2(-3, z, Q2) / z;//sbar
    double PhT2 = z * z * 0.57 + 0.12;
    return x * (eu2 * (f1[1] * D1[0] + f1[4] * D1[3]) + ed2 * (f1[0] * D1[1] + f1[2] * D1[2] + f1[3] * D1[4] + f1[5] * D1[5])) * exp(-Pt * Pt / PhT2) / (PI * PhT2);
  }

  double getFUUT_p_FD(const double x, const double z, const double Q2, const double Pt){
    return x*eu2*getIntegral_u()+x*ed2*getIntegral_d();
  }
  
  double getFUUT_n_FD(const double x, const double z, const double Q2, const double Pt){
     return x*eu2*getIntegral_u()+x*ed2*getIntegral_d();
  }
  
  TFile *file = new TFile("grids_hist.root");
  TH3D *h_f1u = (TH3D*)file->Get("f1u");
  TH3D *h_D1u = (TH3D*)file->Get("D1u");  
  
  double Tint(const double x, const double z, const double Q2, const double Pt, const double * par){
    double pt = sqrt(Pt * Pt + z * z * par[0] * par[0] - 2 * z * Pt * par[0] * cos(par[1]));
    return par[0] * h_f1u->Interpolate(x, par[0], Q2) * h_D1u->Interpolate(z, pt, Q2);
  }
  
  double getIntegral_u(const double x, const double z, const double Q2, const double Pt){
    double xl[2] = {0, -PI};//integration lower bound
    double xu[2] = {1.0, PI};//integration upper bound
    ROOT::Math::WrappedParamFunction<> wf(&Tint, 2, 2);
    wf.SetParameters(par);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001, 100000);
    ig.SetFunction(wf);
    double val = ig.Integral(xl, xu);
    return val; 
  }
  
  
}

#endif
