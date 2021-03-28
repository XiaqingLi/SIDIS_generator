#ifndef _STRU_FUNC_H
#define _STRU_FUNC_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"
#include "TMath.h"
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

/////////////////////////////////  Mode 0  ///////////////////////////////
////////////////// use Gaussian ansatz for PDF and FF ////////////////////

  LHAPDF::PDF * pdfs;
  LHAPDF::PDF * ffs;
  bool st_pdf = false;
  bool st_ff = false;
  int hid;

  int setPDF(){
    pdfs = LHAPDF::mkPDF("CJ15lo", 0);
    //pdfs = LHAPDF::mkPDF("MMHT2014lo68cl", 0);
    st_pdf = true;
    return 0;
  }  
  
  int setFF(const int hid){
    if(hid >= 0)
      ffs = LHAPDF::mkPDF("DSSFFlo", hid);
    else
      ffs = LHAPDF::mkPDF("DSSFFlo", -hid+1000); 
    st_ff = true;
    return 0;
  }

  double getFUUT_p_Gauss(const double x, const double z, const double Q2, const double Pt){
    if(!st_pdf) setPDF();
    if(!st_ff) setFF(hid);
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
    //double PhT2 = z * z * 0.4 + 0.25; //for rho production
    return x * (eu2 * (f1[0] * D1[0] + f1[3] * D1[3]) + ed2 * (f1[1] * D1[1] + f1[2] * D1[2] + f1[4] * D1[4] + f1[5] * D1[5])) * exp(-Pt * Pt / PhT2) / (PI * PhT2);    
  }

  double getFUUT_n_Gauss(const double x, const double z, const double Q2, const double Pt){
    if(!st_pdf) setPDF();
    if(!st_ff) setFF(hid);
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



/////////////////////////////////  Mode 1  ///////////////////////////////
////////////// use TMD grids and interpolation for PDF and FF ////////////
//////// notation: x=x_bj; Q2=Q^2 (GeV); kt=k_perp; pt=p_perp; Pt=PhT; 
  TFile *file = new TFile("/Users/xiaqingli/Projects/EIC/SIDIS_generator/grids/grids_hist_MMHT_old.root");
  TH3D *h_f1u = (TH3D*)file->Get("f1u");
  TH3D *h_D1u = (TH3D*)file->Get("D1u");  
  TH3D *h_f1d = (TH3D*)file->Get("f1d");
  TH3D *h_D1d = (TH3D*)file->Get("D1d");  
  TH3D *h_f1s = (TH3D*)file->Get("f1s");
  TH3D *h_D1s = (TH3D*)file->Get("D1s");  
  TH3D *h_f1ubar = (TH3D*)file->Get("f1ubar");
  TH3D *h_D1ubar = (TH3D*)file->Get("D1ubar");  
  TH3D *h_f1dbar = (TH3D*)file->Get("f1dbar");
  TH3D *h_D1dbar = (TH3D*)file->Get("D1dbar");  
  TH3D *h_f1sbar = (TH3D*)file->Get("f1sbar");
  TH3D *h_D1sbar = (TH3D*)file->Get("D1sbar");  
  const double xmin = 0;//h_f1u->GetXaxis()->GetBinLowEdge(2);
  const double xmax = 0.99;
  const double zmin = 0.2;//h_D1u->GetXaxis()->GetBinLowEdge(2);
  const double zmax = 0.8;
  const double Q2min = 1;//h_f1u->GetYaxis()->GetBinLowEdge(2);
  const double Q2max = 500.0;
  const double ktmin = 0;//h_f1u->GetZaxis()->GetBinLowEdge(2);
  const double ktmax = 1;
  const double ptmin = 0;//h_D1u->GetZaxis()->GetBinLowEdge(2);
  const double ptmax = 1;
  
  // int Nmax_x = h_f1u->GetXaxis()->GetNbins();
  // double xlow = pow(10, h_f1u->GetXaxis()->GetBinLowEdge(2));
  // double xhigh = pow(10, h_f1u->GetXaxis()->GetBinLowEdge(Nmax_x-1));
  // int Nmax_Q2 = h_f1u->GetYaxis()->GetNbins();
  // double Q2low = pow(10, h_f1u->GetYaxis()->GetBinLowEdge(2));
  // double Q2high = pow(10, h_f1u->GetYaxis()->GetBinLowEdge(Nmax_Q2-1));
  // int Nmax_kt = h_f1u->GetZaxis()->GetNbins();
  // double ktlow = h_f1u->GetZaxis()->GetBinLowEdge(2);
  // double kthigh = h_f1u->GetZaxis()->GetBinLowEdge(Nmax_kt-1);
  // int Nmax_z = h_D1u->GetXaxis()->GetNbins();
  // double zlow = h_D1u->GetXaxis()->GetBinLowEdge(2);
  // double zhigh = h_D1u->GetXaxis()->GetBinLowEdge(Nmax_z-1);
  // int Nmax_pt = h_D1u->GetZaxis()->GetNbins();
  // double ptlow = h_D1u->GetZaxis()->GetBinLowEdge(2);
  // double pthigh = h_D1u->GetZaxis()->GetBinLowEdge(Nmax_pt-1);

  double Integrand_p(const double *var, const double *par){
    //////////var[2]={kt, angle}; par[4]={x, z, Q2, PhT}
    double _Q = sqrt(par[2]);
    double pt = sqrt(par[3] * par[3] + par[1] * par[1] * var[0] * var[0] - 2 * par[1] * par[3] * var[0] * cos(var[1]));    //p_perp^2 = (z*k_perp - PhT)^2
    /* if(var[0]<ktmin || var[0]>ktmax || pt<ptmin || pt>ptmax){ */
    if(var[0]>_Q*1 || pt>_Q/par[1]){ //set k_perp and p_perp cut off values in convolution in momentum space for FUUT
      return 0;
    // }else if(par[0]>xhigh || par[0]<xlow || par[2]>Q2high || par[2]<Q2low || var[0]/_Q<ktlow || var[0]/_Q>kthigh || par[1]>zhigh || par[1]<zlow || pt/(_Q/par[1])<ptlow || pt/(_Q/par[1])>pthigh){
    //   return 0;
    }else{
      double val_u = eu2 * h_f1u->Interpolate(log10(par[0]), log10(par[2]), var[0]/_Q) * h_D1u->Interpolate(par[1], log10(par[2]), pt/(_Q/par[1]));
      double val_ubar = eu2 * h_f1ubar->Interpolate(log10(par[0]), log10(par[2]), var[0]/_Q) * h_D1ubar->Interpolate(par[1], log10(par[2]), pt/(_Q/par[1]));
      double val_d = ed2 * h_f1d->Interpolate(log10(par[0]), log10(par[2]), var[0]/_Q) * h_D1d->Interpolate(par[1], log10(par[2]), pt/(_Q/par[1]));
      double val_dbar = ed2 * h_f1dbar->Interpolate(log10(par[0]), log10(par[2]), var[0]/_Q) * h_D1dbar->Interpolate(par[1], log10(par[2]), pt/(_Q/par[1]));
      double val_s = ed2 * h_f1s->Interpolate(log10(par[0]), log10(par[2]), var[0]/_Q) * h_D1s->Interpolate(par[1], log10(par[2]), pt/(_Q/par[1]));
      double val_sbar = ed2 * h_f1sbar->Interpolate(log10(par[0]), log10(par[2]), var[0]/_Q) * h_D1sbar->Interpolate(par[1], log10(par[2]), pt/(_Q/par[1]));
      if(val_u < 0)  val_u = 0;
      if(val_ubar < 0)  val_ubar = 0;
      if(val_d < 0)  val_d = 0;
      if(val_dbar < 0)  val_dbar = 0;
      if(val_s < 0)  val_s = 0;
      if(val_sbar < 0)  val_sbar = 0;
      return var[0] * (val_u + val_ubar + val_d + val_dbar + val_s + val_sbar);
      //return val_d / ed2;
    }
  }
  
  double getIntegral_p(const double x, const double z, const double Q2, const double Pt){
    double par[4] = {x, z, Q2, Pt};
    double xl[2] = {0.0, -PI};//integration lower bound {k_t, phi}
    double xu[2] = {sqrt(Q2), PI};//integration upper bound {k_t, phi}
    ROOT::Math::WrappedParamFunction<> wf(&Integrand_p, 2, 4);
    wf.SetParameters(par);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.01, 1000000);
    ig.SetFunction(wf);
    return ig.Integral(xl, xu);
  }

  double getFUUT_p_grid(const double x, const double z, const double Q2, const double Pt){
    if(x<xmin || x>xmax || z<zmin || z>zmax || Q2<Q2min || Q2>Q2max){
      perror("Kinematics out of interpolation range!");
      return 0;     
    }else
      return x * getIntegral_p(x, z, Q2, Pt);
  }


  
  double Integrand_n(const double *var, const double *par){
    //////////var[2]={kt, angle}; par[4]={x, z, Q2, Pt}
    double pt = sqrt(par[3] * par[3] + par[1] * par[1] * var[0] * var[0] - 2 * par[1] * par[3] * var[0] * cos(var[1]));
    if(var[0]<=ktmin || var[0]>=ktmax || pt<=ptmin || pt>=ptmax){
      return 0;
    }else{
      double val_u = ed2 * h_f1d->Interpolate(log10(par[0]), log10(par[2]), var[0]) * h_D1u->Interpolate(par[1], log10(par[2]), pt);
      double val_ubar = ed2 * h_f1dbar->Interpolate(log10(par[0]), log10(par[2]), var[0]) * h_D1ubar->Interpolate(par[1], log10(par[2]), pt);
      double val_d = eu2 * h_f1u->Interpolate(log10(par[0]), log10(par[2]), var[0]) * h_D1d->Interpolate(par[1], log10(par[2]), pt);
      double val_dbar = eu2 * h_f1ubar->Interpolate(log10(par[0]), log10(par[2]), var[0]) * h_D1dbar->Interpolate(par[1], log10(par[2]), pt);
      double val_s = ed2 * h_f1s->Interpolate(log10(par[0]), log10(par[2]), var[0]) * h_D1s->Interpolate(par[1], log10(par[2]), pt);
      double val_sbar = ed2 * h_f1sbar->Interpolate(log10(par[0]), log10(par[2]), var[0]) * h_D1sbar->Interpolate(par[1], log10(par[2]), pt);
      if(val_u < 0)  val_u = 0;
      if(val_ubar < 0)  val_ubar = 0;
      if(val_d < 0)  val_d = 0;
      if(val_dbar < 0)  val_dbar = 0;
      if(val_s < 0)  val_s = 0;
      if(val_sbar < 0)  val_sbar = 0;
      return var[0] * (val_u + val_ubar + val_d + val_dbar + val_s + val_sbar);
    }
  }

  double getIntegral_n(const double x, const double z, const double Q2, const double Pt){
    double par[4] = {x, z, Q2, Pt};
    double xl[2] = {0.0, -PI};//integration lower bound
    double xu[2] = {2.0, PI};//integration upper bound
    ROOT::Math::WrappedParamFunction<> wf(&Integrand_n, 2, 4);
    wf.SetParameters(par);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.00001, 10000000);
    ig.SetFunction(wf);
    double val = ig.Integral(xl, xu);
    return val; 
  }

  double getFUUT_n_grid(const double x, const double z, const double Q2, const double Pt){
    if(x<xmin || x>xmax || z<zmin || z>zmax || Q2<Q2min || Q2>Q2max){
      perror("Kinematics out of interpolation range!");
      return 0;     
    }else
      return x * getIntegral_n(x, z, Q2, Pt);
  }
  
  
}

#endif
