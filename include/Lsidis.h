//This is a class for sidis generator or cross section calculation
//Tested with gcc v4.4.7
//LHAPDF is required. Developed with v6.1
//ROOT is required. Developed with v5.34.21
//Potential error with other versions not tested
//Last update date 9 April 2019 version 3.0

#ifndef _LSIDIS_H_
#define _LSIDIS_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "TLorentzVector.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TH1D.h"
#include "TF1.h"

#include "constant.h"
#include "stru_func.h"

class Lsidis{
 protected:

 private:
  bool st_nucleus;
  bool st_pdf;
  bool st_ff;
  bool st_hadron;
  bool st_l;
  bool st_P;
  bool st_lp;
  bool st_Ph;
  bool st_gibbs;
  double Np;//number of protons in the nucleus
  double Nn;//number of neutrons in the nucleus
  double Polp;//effective polarization of the proton in the nucleus
  double Poln;//effective polarization of the neutron in the nucleus
  double Mh;//final hadron mass
  double MXminp;//minimum of the invariant mass of the final X from proton
  double MXminn;//minimum of the invariant mass of the final X from neutron
  char hname[10];//final hadron name
  int lid;//incoming lepton id
  int lic;//incoming lepton charge
  int hid;//final hadron id
  int hic;//final hadron charge
  LHAPDF::PDF * pdfs;//pointer to collinear unpolarized PDFs
  LHAPDF::PDF * ffs;//pointer to collinear unpolarized FFs
  TLorentzVector Pl;//incoming lepton 
  TLorentzVector Plp;//outgoing lepton
  TLorentzVector Pq;//virtual boson
  TLorentzVector PP;//incoming nucleus
  TLorentzVector PPh;//outgoing hadron
  TLorentzVector PX;//outgoing undetected particles
  double Slepton;//incoming lepton polarization rate
  double SNL;//incoming nucleus longitudinal polarization rate
  double SNT;//incoming nucleus transverse polarization rate
  double x;//Bjorken scaling variable x
  double xn;//Nachtmann scaling variable xn
  double y;//lepton energy transferred fraction scalar variable
  double z;//fragmentation scalar variable
  double zn;//Nachtmann like light-front fragmentation fraction zn
  double Q2;//transferred momentum square, also used as factorization scale
  double Pt;//transverse momentum of final hadron in Trento convention
  double phih;//azimuthal angle of final hadron in Trento convention
  double phiS;//azimuthal angle of transverse polarization in Trento convention
  double W;//invariant mass of final state except scattered lepton
  double Wp;//invariant mass of final state except scattered lepton and the hadron
  double gamma;//
  double epsilon;//longi-trans-photon flux ratio
  double Rfactor;//a criterism for current fragmentation judgement
  double jacobian;//Jacobian from simulation space to cross section defined space
  bool physics_control;//
  double f1[6];//PDFs of u, d, s, ubar, dbar, sbar
  double D1[6];//FFs of u, d, s, ubar, dbar, sbar
  double g1[6];//Polarized PDF of u, d, s, ubar, dbar, sbar
  double h1[6];//Transversity PDF of u, d, s, ubar, dbar, sbar
//  double TMDpars[2];//Model parameters of gaussian type TMDs
  double Xmin[6];//simulation variables lower limits
  double Xmax[6];//simulation variables upper limits
  double Xseed[6];//start point for MCMC
  double volume;//simulation variables volume
  double sigmatotal;//total cross section within the simulation range [Xmin, Xmax]
  int FUUTmode;//mode for structure function
  TH1D * Xhisto[6];//for Gibbs sampler
  TF1 * TF_Maxwell;//for R1 sampling
  TF1 * TF_xi;
  TF1 * TF_zeta;
  TF1 * TF_kt;
 public:
  Lsidis();
  Lsidis(const TLorentzVector l, const TLorentzVector P);
  static double DEG(const double rad);//convert from rad to deg
  static double RAD(const double deg);//convert form deg to rad
  int SetNucleus(const double np, const double nn, const double pp, const double pn);//Set incoming nucleus
  int CheckNucleus();
  int SetHadron(const char * hadron);//Set detected hadron
  int GetHadronID();//Get pid of the detected hadron
  int GetHadronCharge();//Get charge of the detected hadron
  int CheckHadron();
  int SetFUUTmode(const int _FUUTmode);//Set mode for structure function FUU,T
  int SetInitialState(const TLorentzVector l, const TLorentzVector P);//Set incoming lepton and nucleon 4-momenta
  int SetLeptonBeam(const TLorentzVector l);//Set incoming lepton 4-momentum
  int SetLeptonBeam(const double Ebeam, const double m);//Set incoming lepton along z-direction
  int SetIonBeam(const TLorentzVector P);//Set incoming nucleon 4-momentum
  int SetIonBeam(const double Ebeam, const double M);//Set incoming nucleon along negative z-direction
  int SetTarget(const double np, const double nn, const double pp, const double pn);//Set fixed target
  int SetTarget();//Set fixed target
  int SetFinalState(const TLorentzVector lp, const TLorentzVector Ph);//Set final lepton and detected hadron 4-momenta
  int SetFinalLepton(const TLorentzVector lp);//Set final scattered lepton 4-momenta
  int SetFinalHadron(const TLorentzVector Ph);//Set final detected hadron 4-momenta
  int CalculateVariables();//Calculate Lorentz scalar variables
  int CalculateFinalState();//Calculate scattered lepton and detected hadron from x, y, z, Pt, phih, phiS
  int CalculateRfactor(const double kT2, const double MiT2, const double MfT2);//Calculate the Rfactor for current fragmentation criteria
  int R1SamplerStarter();//Initialize functions for R1 sampling
  double R1Sampler();//Random sampling R1 following Ted Rogers' method
  int SetVariables(const double x0, const double y0, const double z0, const double Pt0, const double phih0, const double phiS0);//Set variables x, y, z, Pt2, phih, phiS
  double GetVariable(const char * var);//Get particular variable of current event
  TLorentzVector GetLorentzVector(const char * part);//Get 4-momentum of a particle of current event
  int GetPDFs();//Get PDFs at current x, Q2
  int GetFFs();//Get FFs at current z, Q2
  double FUUT();//Get structure function F_UU,T at current x, z, Q2, Pt
  double getFUUT_p();//Get stucture function of the proton from the chosen mode
  double getFUUT_n();//Get stucture function of the neutron from the chosen mode	
  double dsigma(const int mode);//Differential cross section in dsigma/dx dy dz dPt2 dphih dphiS
  int SetRange(const double * mins, const double * maxs);//Uniform sampling range
  double GenerateEvent(const int mode, const int method);//Generate a sidis event
  double GetEventWeight(const int mode, const int method);//Get event weight
  int CalculateSigmaTotal(const int mode, const int method, const int precision);//Calculate total cross section within the simulation range [Xmin, Xmax]
  int GibbsStarter(const int mode, const int method);
  double GibbsSampler(const int mode, const int method);//Generate an event with Gibbs sampler
  int Test();
};

Lsidis::Lsidis(){//constructor
  st_nucleus = false;
  st_hadron = false;
  st_pdf = false;
  st_ff = false;
  st_l = false;
  st_P = false;
  st_lp = false;
  st_Ph = false;
  st_gibbs = false;
  physics_control = false;
//  TMDpars[0] = 0.57;//kt2
//  TMDpars[1] = 0.12;//pt2
  Slepton = 0;
  SNL = 0;
  SNT = 0;
  FUUTmode = 0; //default mode for structure function
}

Lsidis::Lsidis(const TLorentzVector l, const TLorentzVector P){//constructor with initial state
  Lsidis();
  SetInitialState(l, P);
}

double Lsidis::DEG(const double rad){//convert rad to deg
  return rad / PI * 180.0;
}

double Lsidis::RAD(const double deg){//convert deg to rad
  return deg / 180.0 * PI;
}

int Lsidis::SetNucleus(const double np, const double nn, const double pp = 0, const double pn = 0){//set number of protons and neutrons in the incoming nucleus
  Np = np;
  Nn = nn;
  Polp = pp;
  Poln = pn;
  st_nucleus = true;
  return 0;
}

int Lsidis::CheckNucleus(){//check the status of nucleus
  if (st_nucleus){
    printf("Number of protons:  %.1f\n", Np);
    printf("Number of neutrons: %.1f\n", Nn);
    return 0;
  }
  else {
    printf("No nucleus initialized!\n");
    return -1;
  }
}

int Lsidis::SetHadron(const char * hadron){//set final hadron
  if (strcmp(hadron, "pi+") == 0){
    hic = 1; hid = 211; 
    Mh = Mpion;
    MXminp = Mp;
    MXminn = Mp + Mpion;
  }
  else if (strcmp(hadron, "pi-") == 0){
    hic = -1; hid = -211;
    Mh = Mpion;
    MXminp = Mp + Mpion;
    MXminn = Mp;
  }
  else if (strcmp(hadron, "pi0") == 0){
    hic = 0; hid = 111;
    Mh = Mpi0;
    MXminp = Mp;
    MXminn = Mp;
  }
  else if (strcmp(hadron, "K+") == 0){
    hic = 1; hid = 321;
    Mh = Mkaon;
    MXminp = MLambda;
    MXminn = MSigmam;
  }
  else if (strcmp(hadron, "K-") == 0){
    hic = -1; hid = -321;
    Mh = Mkaon;
    MXminp = Mp + Mkaon;
    MXminn = Mp + Mkaon;
  }
  else if (strcmp(hadron, "K0") == 0){
    hic = 0; hid = 310;//Ks
    Mh = MK0;
    MXminp = MSigmap;
    MXminn = MLambda;
  }
  else if (strcmp(hadron, "rho+") == 0){
    hic = 1; hid = 213;
    Mh = Mrho;
    MXminp = Mp;
    MXminn = Mp + Mpion;
  }
  else if (strcmp(hadron, "rho-") == 0){
    hic = -1; hid = -213;
    Mh = Mrho;
    MXminp = Mp + Mpion;
    MXminn = Mp;
  }
  else if (strcmp(hadron, "rho0") == 0){
    hic = 0; hid = 113;
    Mh = Mrho;
    MXminp = Mp;
    MXminn = Mp;
  }
  else if (strcmp(hadron, "p") == 0){
    hic = 1; hid = 2212;
    Mh = Mp;
    MXminp = Mpi0;
    MXminn = Mpion;
  }
  else {
    perror("Invalid hadron option in SetHadron!");
    return -1;
  }
  strcpy(hname, hadron);
  StructureFunction::hid = hid;
  st_hadron = true;
  return 0;
}

int Lsidis::GetHadronID(){//get pid of the detected hadron
  return hid;
}

int Lsidis::GetHadronCharge(){//get charge of the detected hadron
  return hic;
}

int Lsidis::CheckHadron(){//check the status of final hadron
  if (st_hadron){
    printf("final hadron: %s\n", hname);
    printf("charge: %d,\t pid: %d\n", hic, hid);
    return 0;
  }
  else {
    perror("No hadron initialized!");
    return -1;
  }
}

int Lsidis::SetFUUTmode(const int _FUUTmode){
  FUUTmode = _FUUTmode;
  return 0;
}

int Lsidis::SetInitialState(const TLorentzVector l, const TLorentzVector P){//set incoming lepton and nucleon 4-momentums
  Pl = l;
  PP = P;
  st_l = true;
  st_P = true;
  return 0;
}

int Lsidis::SetLeptonBeam(const TLorentzVector l){//set incoming lepton 4-momentum
  Pl = l;
  st_l = true;
  return 0;
}

int Lsidis::SetLeptonBeam(const double Ebeam, const double m = Me){//set incoming lepton with beam energy Ebeam along z-direction
  Pl.SetXYZM(0, 0, sqrt(Ebeam * Ebeam - m * m), m);
  st_l = true;
  return 0;
}

int Lsidis::SetIonBeam(const TLorentzVector P){//set incoming nucleon 4-momentum
  PP = P;
  st_P = true;
  return 0;
}

int Lsidis::SetIonBeam(const double Ebeam, const double M = Mp){//set incoming nucleon with beam energy Ebeam along negative z-direction
  PP.SetXYZM(0, 0, -sqrt(Ebeam * Ebeam - M * M), Ebeam);
  st_P = true;
  return 0;
}

int Lsidis::SetTarget(const double np, const double nn, const double pp, const double pn){//set fixed target
  SetNucleus(np, nn, pp, pn);
  PP.SetXYZM(0, 0, 0, Mp);
  return 0;
}

int Lsidis::SetTarget(){//set fixed target
  if (st_nucleus){
    PP.SetXYZM(0, 0, 0, Mp);
    return 0;
  }
  else {
    perror("No nucleus initialized!");
    return -1;
  }
}

int Lsidis::SetFinalState(const TLorentzVector lp, const TLorentzVector Ph){//set final state
  SetFinalLepton(lp);
  SetFinalHadron(Ph);
  return 0;
}

int Lsidis::SetFinalLepton(const TLorentzVector lp){//set scattered lepton
  Plp = lp;
  st_lp = true;
  return 0;
}

int Lsidis::SetFinalHadron(const TLorentzVector Ph){//set final hadron
  PPh = Ph;
  st_Ph = true;
  return 0;
}

int Lsidis::CalculateVariables(){//calculate lorentz scalar variables
  if (st_l && st_P && st_lp && st_Ph && st_nucleus && st_hadron){
    Pq = Pl - Plp;//calculate virtual photon 4-momentum
    PX = Pq + PP - PPh;//calculate final X 4-momentum
    Q2 = - (Pq * Pq);
    if (Q2 < Mpion * Mpion){//Q2 below a cut
      physics_control = false;
      return 0;
    }
    double W2 = Mp * Mp + 2.0 * PP * Pq - Q2;
    if (W2 < pow(Mp + Mpi0, 2)){//below the lowest threshold
      physics_control = false;
      return 0;
    }
    W = sqrt(W2);
    double Wp2 = PX * PX;
    if (Wp2 < Mpi0 * Mpi0){//below the lowest threshold
      physics_control = false;
      return 0;
    }
    Wp = sqrt(Wp2);
    x = Q2 / (PP * Pq) / 2.0;
    if (x > 1.0 || x < 1e-6){//x below a cut
      physics_control = false;
      return 0;
    }
    y = (PP * Pq) / (PP * Pl);
    if (y > 1.0 || y < 0.0){
      physics_control = false;
      return 0;
    }
    z = (PP * PPh) / (PP * Pq);
    if (z > 1.0 || z < 0.05){//z below a cut
      physics_control = false;
      return 0;
    }
    gamma = 2.0 * x * Mp / sqrt(Q2);
    epsilon = (1.0 - y - 0.25 * gamma * gamma * y * y) / (1.0 - y + 0.5 * y * y + 0.25 * gamma * gamma * y * y);
    xn = 2.0 * x / (1.0 + sqrt(1.0 + gamma * gamma));
    double Pt2 = - (PPh * PPh) + 2.0 * (PP * PPh) * (Pq * PPh) / (PP * Pq) / (1.0 + gamma * gamma) - gamma * gamma / (1.0 + gamma * gamma) * (pow(Pq * PPh, 2) / Q2 - pow(PP * PPh, 2) / (Mp * Mp));
    if (Pt2 < 0) Pt2 = 0.0;
    Pt = sqrt(Pt2);
    zn = xn * z / (2.0 * x) * (1.0 + sqrt(1.0 - 4.0 * Mp * Mp * (Mh * Mh + Pt * Pt) * x * x / (z * z * Q2 * Q2)));//add since v3.0
    double lt2 = - (Pl * Pl) + 2.0 * (PP * Pl) * (Pq * Pl) / (PP * Pq) / (1.0 + gamma * gamma) - gamma * gamma / (1.0 + gamma * gamma) * (pow(Pq * Pl, 2) / Q2 - pow(PP * Pl, 2) / (Mp * Mp));
    double ch = - 1.0 / sqrt(lt2 * Pt * Pt) * ( (Pl * PPh) - ((Pq * Pl) * (PP * PPh) + (PP * Pl) * (Pq * PPh)) / (1.0 + gamma * gamma) / (PP * Pq) + gamma * gamma / (1.0 + gamma * gamma) * ( (Pq * Pl) * (Pq * PPh) / Q2 - (PP * Pl) * (PP * PPh) / (Mp * Mp)));
    TMatrixD eg5(4,4);
    eg5(0,0) = Pl.E(); eg5(0,1) = -Pl.X(); eg5(0,2) = -Pl.Y(); eg5(0,3) = -Pl.Z();
    eg5(1,0) = PPh.E(); eg5(1,1) = -PPh.X(); eg5(1,2) = -PPh.Y(); eg5(1,3) = -PPh.Z();
    eg5(2,0) = PP.E(); eg5(2,1) = -PP.X(); eg5(2,2) = -PP.Y(); eg5(2,3) = -PP.Z();
    eg5(3,0) = Pq.E(); eg5(3,1) = -Pq.X(); eg5(3,2) = -Pq.Y(); eg5(3,3) = -Pq.Z();
    double sh = -1.0 / sqrt(lt2 * Pt * Pt) * eg5.Determinant() / ((PP * Pq) * sqrt(1.0 + gamma * gamma));
    if (sh > 0) phih = acos(ch);
    else phih = -acos(ch);
    physics_control = true;
//    GetPDFs();
//    GetFFs();
    return 0;
  }
  else {
    physics_control = false;
    perror("Initial or final state missing!");
    return -1;
  }
}

int Lsidis::CalculateFinalState(){//Calculate scattered electron and detected hadron from x, y, z, Pt, phih, phiS
  if (st_l && st_P && st_nucleus && st_hadron){
    if (x > 1.0 || x < 1e-6){//x below a cut
      physics_control = false;
      return 0;
    }
    if (y > 1.0 || y < 0.0){//y
      physics_control = false;
      return 0;
    }
    if (z > 1.0 || z < 0.05){//z below a cut
      physics_control = false;
      return 0;
    }
    if (Pt < 0.0){//unphysical Pt
      physics_control = false;
      return 0;
    }
    Q2 = 2.0 * x * y * (PP * Pl);
    if (Q2 < Mpion * Mpion){//Q2 below a cut
      physics_control = false;
      return 0;
    }
    TLorentzVector Pl_2 = Pl;
    Pl_2.Boost(-PP.BoostVector());
    TLorentzVector Pl_1 = Pl_2;
    Pl_1.RotateZ(-Pl_2.Phi());
    Pl_1.RotateY(-Pl_2.Theta());
    double Elp_1 = Pl_1.E() - Q2 / (2.0 * x * Mp);
    double theta_lp_1 = acos(1.0 - Q2 / (2.0 * Pl_1.E() * Elp_1));
    TLorentzVector Plp_1(Elp_1 * sin(theta_lp_1), 0, Elp_1 * cos(theta_lp_1), Elp_1);
    TLorentzVector Pq_1 = Pl_1 - Plp_1;
    double shl = -cos(Pq_1.Theta()) * sin(phiS) / sqrt(1.0 - pow(sin(Pq_1.Theta()) * sin (phiS), 2));
    double chl = cos(phiS) / sqrt(1.0 - pow(sin(Pq_1.Theta()) * sin (phiS), 2));
    double phil;
    if (shl > 0) phil = acos(chl);
    else phil = -acos(chl);
    Plp_1.RotateZ(phil);
    Pq_1.RotateZ(phil);
    Plp = Plp_1;
    Plp.RotateY(Pl_2.Theta());
    Plp.RotateZ(Pl_2.Phi());
    Plp.Boost(PP.BoostVector());
    Pq = Pl - Plp;
    if (z < sqrt(Mh * Mh + Pt * Pt) / Pq_1.E()){//below hadron threshold
      physics_control = false;
      return 0;
    }
    TLorentzVector PPh_0;
    TLorentzVector Plp_0 = Plp_1;
    Plp_0.RotateZ(-Pq_1.Phi());
    Plp_0.RotateY(-Pq_1.Theta());
    PPh_0.SetXYZT(Pt * cos(phih + Plp_0.Phi()), Pt * sin(phih + Plp_0.Phi()), sqrt(pow(z * Pq_1.E(), 2) - Mh * Mh - Pt * Pt), z * Pq_1.E());
    TLorentzVector PPh_1 = PPh_0;
    PPh_1.RotateY(Pq_1.Theta());
    PPh_1.RotateZ(Pq_1.Phi());
    PPh = PPh_1;
    PPh.RotateY(Pl_2.Theta());
    PPh.RotateZ(Pl_2.Phi());
    PPh.Boost(PP.BoostVector());
    //CalculateVariables();
    double W2 = Mp * Mp + 2.0 * (PP * Pq) - Q2;
    if (W2 < pow(Mp + Mpi0, 2)){//below the lowest threshold
      physics_control = false;
      return 0;
    }
    W = sqrt(W2);
    PX = Pq + PP - PPh;
    double Wp2 = PX * PX;
    if (Wp2 < Mpi0 * Mpi0){//below the lowest threshold
      physics_control = false;
      return 0;
    }
    Wp = sqrt(Wp2);
    gamma = 2.0 * x * Mp / sqrt(Q2);
    epsilon = (1.0 - y - 0.25 * gamma * gamma * y * y) / (1.0 - y + 0.5 * y * y + 0.25 * gamma * gamma * y * y);
    xn = 2.0 * x / (1.0 + sqrt(1.0 + gamma * gamma));
    zn = xn * z / (2.0 * x) * (1.0 + sqrt(1.0 - 4.0 * Mp * Mp * (Mh * Mh + Pt * Pt) * x * x / (z * z * Q2 * Q2)));//add since v3.0
    physics_control = true;
//    GetPDFs();
//    GetFFs();
    return 0;
  }
  else {
    physics_control = false;
    perror("Initial state missing for final state calculation!");
    return -1;
  }
}

int Lsidis::CalculateRfactor(const double kT2 = 0.5, const double MiT2 = 0.5, const double MfT2 = 0.5){
  if (physics_control){
    double yi = 0.5 * log(Q2 / MiT2);
    double yf = -0.5 * log(Q2 / MfT2);
    double yh = log( (sqrt(Q2) * z * (Q2 - xn * xn * Mp * Mp)) / ( 2.0 * xn * xn * Mp * Mp * sqrt(Mh * Mh + Pt * Pt)) - sqrt(Q2) / (xn * Mp) * sqrt(pow(z * (Q2 - xn * xn * Mp * Mp), 2) / (4.0 * xn * xn  * Mp * Mp * (Mh * Mh + Pt * Pt)) - 1.0));
    double Rf = 0.5 * sqrt(Pt * Pt + Mh * Mh) * sqrt(MfT2) * (exp(yf - yh) + exp(yh - yf)) - sqrt(kT2) * Pt;
    double Ri = 0.5 * sqrt(Pt * Pt + Mh * Mh) * sqrt(MiT2) * (exp(yi - yh) - exp(yh - yi)) - sqrt(kT2) * Pt;
    Rfactor = std::abs(Rf / Ri);
    return 0;
  }
  else {
    Rfactor = 1. / 0.;//NaN
    return -1;
  }
}

int Lsidis::R1SamplerStarter(){//initialize the R1 sampler
  TF_Maxwell = new TF1("TF_Maxwell", "x*x*exp(-x*x/0.02)", 0.0, 1.0);
  TF_xi = new TF1("TF_xi", "x*pow(1.0-x,3)", 0.0, 1.0);
  TF_zeta = new TF1("TF_zeta", "x*pow(1.0-x,2)", 0.0, 1.0);
  TF_kt = new TF1("TF_kt", "exp(-x/0.4)", 0.0, 1.5); 
  return 0;
}

double Lsidis::R1Sampler(){//random sampling R1 following Ted Rogers' method
  double xi = TF_xi->GetRandom(xn, 1.0);
  double zeta = TF_zeta->GetRandom(zn, 1.0);
  double dkt = TF_kt->GetRandom();
  double ki = TF_Maxwell->GetRandom();
  double kf = TF_Maxwell->GetRandom();
  double angle = gRandom->Uniform(-M_PI, M_PI);
  double xnhat = xn / xi;
  double denominator = (zn * Q2) / (2.0 * xnhat) - pow(ki, 2) * xn * (Mh * Mh + Pt * Pt) / (2.0 * zn * Q2);
  double kft2 = pow(Pt / zeta, 2) + dkt * dkt - 2.0 * dkt * Pt / zeta * cos(angle);
  double numerator = (kft2 + kf * kf) * zeta / 2.0 + (Mh * Mh + Pt * Pt) / (2.0 * zeta) + Pt * Pt / zeta - dkt * Pt * cos(angle);
  return abs(numerator / denominator);
}

double Lsidis::GetVariable(const char * var){//get current variable
  if (!physics_control){
    std::cout << "Warning: unphysical event!" << std::endl;
    return -1;
  }
  if (strcmp(var, "x") == 0) return x;
  else if (strcmp(var, "y") == 0) return y;
  else if (strcmp(var, "z") == 0) return z;
  else if (strcmp(var, "Q2") == 0) return Q2;
  else if (strcmp(var, "Pt") == 0) return Pt;
  else if (strcmp(var, "phih") == 0) return phih;
  else if (strcmp(var, "phiS") == 0) return phiS;
  else if (strcmp(var, "W") == 0) return W;
  else if (strcmp(var, "Wp") == 0) return Wp;
  else if (strcmp(var, "xn") == 0) return xn;
  else if (strcmp(var, "gamma") == 0) return gamma;
  else if (strcmp(var, "epsilon") == 0) return epsilon;
  else if (strcmp(var, "Rfactor") == 0) return Rfactor;
  else if (strcmp(var, "sigmatotal") == 0) return sigmatotal;
  else if (strcmp(var, "Mh") == 0) return Mh;
  else {
    perror("No such variable!");
    return -1;
  }
}

int Lsidis::SetVariables(const double x0, const double y0, const double z0, const double Pt0, const double phih0, const double phiS0){//Set variables x, y, z, Pt2, phih, phiS
  x = x0;
  y = y0;
  z = z0;
  Pt = Pt0;
  phih = phih0;
  phiS = phiS0;
  st_lp = true;
  st_Ph = true;
  return 0;
}

TLorentzVector Lsidis::GetLorentzVector(const char * part){//Get 4-momentum of a particle in current event
  if (!physics_control){
    std::cout << "Warning: unphysical event!" << std::endl;
    TLorentzVector tt(0,0,0,0);
    return tt;
  }
  if (strcmp(part, "P") == 0) return PP;
  else if (strcmp(part, "l") == 0) return Pl;
  else if (strcmp(part, "lp") == 0) return Plp;
  else if (strcmp(part, "Ph") == 0) return PPh;
  else {
    perror("No such variable!");
    TLorentzVector tt(0,0,0,0);
    return tt;
  }
}

double Lsidis::FUUT(){//calculate the structure function F_UU,T at current x, z, Q2, Pt
  if (!(st_nucleus && st_hadron)){
    perror("Missing initialization of nucleus and hadron!");
    return -1;
  }
  double FFp = 0.0, FFn = 0.0;
  if (physics_control){
    if (Np > 0 && Wp > MXminp)
      FFp = Np * getFUUT_p();
    if (Nn > 0 && Wp > MXminn)
      FFn = Nn * getFUUT_n();
  }
  return FFp + FFn;
}

double Lsidis::getFUUT_p(){
  if (FUUTmode == 0){
    return StructureFunction::getFUUT_p_Gauss(x, z, Q2, Pt);
  }
  else if (FUUTmode == 1)
    return StructureFunction::getFUUT_n_grid(x, z, Q2, Pt);
  else {
    perror("No such struncture funtion mode");
    return -1;
  }
}

double Lsidis::getFUUT_n(){
  if (FUUTmode == 0){	
    return StructureFunction::getFUUT_n_Gauss(x, z, Q2, Pt);
  }
  else if (FUUTmode == 1)
    return StructureFunction::getFUUT_n_grid(x, z, Q2, Pt);
  else {
    perror("No such struncture funtion mode");
    return -1;
  }
}

double Lsidis::dsigma(const int mode = 0){//calculate the differential cross section at current x, y, z, Pt2, phih, phiS
  double ds = 0.0;
  if (physics_control){
    if (mode == 0){//No azimuthal modulations
      ds = alpha_em * alpha_em * y / (2.0 * x * Q2 * (1.0 - epsilon)) * (1.0 + gamma * gamma / (2.0 * x)) * FUUT();
    }
    else {
      perror("No such mode option!");
      return -1;
    }
  }
  return ds;
}

int Lsidis::SetRange(const double * mins, const double * maxs){//Set uniform sampling range
  volume = 1.0;
  for (int i = 0; i < 6; i++){
    Xmin[i] = mins[i];
    Xmax[i] = maxs[i];
    volume = volume * (maxs[i] - mins[i]);
  }
  return 0;
}

double Lsidis::GenerateEvent(const int mode = 0, const int method = 0){//Generate an event and return the cross section weight
  double var[6];
  double weight = 0;
  for (int i = 0; i < 6; i++){
      if(i==0 || i==1)
          var[i] = pow(10, gRandom->Uniform(log10(Xmin[i]), log10(Xmax[i])));
      else
          var[i] = gRandom->Uniform(Xmin[i], Xmax[i]);
  }
  if (method == 0){//generate in x, y, z, Pt, phih, phiS
    jacobian = 2.0 * var[3];
    SetVariables(var[0], var[1], var[2], var[3], var[4], var[5]);
    CalculateFinalState();
    weight = dsigma(mode);
  }
  else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
    double y0 = var[1] / (2.0 * var[0] * (PP * Pl));
    jacobian = 2.0 * var[3] * y0 / var[1];
    SetVariables(var[0], y0, var[2], var[3], var[4], var[5]);
    CalculateFinalState();
    weight = dsigma(mode);
  }
  else if (method == 2){//generate in log(x), Q2, z, Pt, phih, phiS
    double x0 = exp(var[0]);
    double y0 = var[1] / (2.0 * x0 * (PP * Pl));
    jacobian = 2.0 * var[3] * y0 / var[1] * x0;
    SetVariables(x0, y0, var[2], var[3], var[4], var[5]);
    CalculateFinalState();
    weight = dsigma(mode);
  }
  else {
    perror("No such method choice!");
    return -1;
  }
  return weight * volume * jacobian;
}

double Lsidis::GetEventWeight(const int mode = 0, const int method = 0){//Generate an event and return the cross section weight
  double weight = 0;
  if (method == 0){//generate in x, y, z, Pt, phih, phiS
    jacobian = 2.0 * Pt;
    weight = dsigma(mode);
  }
  else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
    double y0 = Q2 / (2.0 * x * (PP * Pl));
    jacobian = 2.0 * Pt * y0 / Q2;
    weight = dsigma(mode);
  }
  else {
    perror("No such method choice!");
    return -1;
  }
  return weight * volume * jacobian;
}

int Lsidis::CalculateSigmaTotal(const int mode = 0, const int method = 0, const int precision = 3){//Calculate the total cross section within the simulation range defined by [Xmin, Xmax] in unit of GeV^-2
  double Nsim = pow(10, precision * 2);
  double sum = 0.0;
  for (Long64_t i = 0; i < Nsim; i++){
    sum += GenerateEvent(mode, method);
  }
  sigmatotal = sum / Nsim;
  return 0;
}

int Lsidis::GibbsStarter(const int mode = 0, const int method = 0){//prepare for Gibbs sampler
  CalculateSigmaTotal(mode, method);//using default precision
  double var[6];
  double weight = 0;
  double judge = 0;
  const int itermax = 10;
  int iter = 0;
  while (judge < 1.0e-9 && iter < itermax){
    iter++;
    for (int ib = 0; ib < 1000; ib++){//search initial Xseed
      for (int i = 0; i < 6; i++){
	var[i] = gRandom->Uniform(Xmin[i], Xmax[i]);
      }
      if (method == 0){//generate in x, y, z, Pt, phih, phiS
	jacobian = 2.0 * var[3];
	SetVariables(var[0], var[1], var[2], var[3], var[4], var[5]);
	CalculateFinalState();
	weight = dsigma(mode) * jacobian;
      }
      else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
	double y0 = var[1] / (2.0 * var[0] * (PP * Pl));
	jacobian = 2.0 * var[3] * y0 / var[1];
	SetVariables(var[0], y0, var[2], var[3], var[4], var[5]);
	CalculateFinalState();
	weight = dsigma(mode) * jacobian;
      }
      else {
	perror("No such method choice!");
	return -1;
      }
      if (weight > judge){
	judge = weight;
	for (int ii = 0; ii < 6; ii++)
	  Xseed[ii] = var[ii];
      }
    }
  }
  if (judge == 0.0){
    perror("Gibbs starter failed!");
    return -1;
  }
  if (iter >= itermax){
    printf("Warning: max iteration reached by Gibbs starter!  %.3E\n", judge);
  }
  st_gibbs = true;
  for (int ih = 0; ih < 6; ih++){
    Xhisto[ih] = new TH1D(Form("X%d", ih), "", 1000, Xmin[ih], Xmax[ih]);
  }
  return 0;
}

double Lsidis::GibbsSampler(const int mode = 0, const int method = 0){//Generate an event with Gibbs sampler
  if (!st_gibbs){
    GibbsStarter(mode, method);
  }
  if (!st_gibbs) return -1;
  int dim = 6;
  if (mode == 0){
    dim = 4;
    Xseed[4] = gRandom->Uniform(-M_PI, M_PI);
    Xseed[5] = gRandom->Uniform(-M_PI, M_PI);
    for (int ih = 0; ih < dim; ih++){
      for (int ib = 1; ib <= Xhisto[ih]->GetNbinsX(); ib++){//set histogram
	Xseed[ih] = Xhisto[ih]->GetBinCenter(ib);
	if (method == 0){//generate in x, y, z, Pt, phih, phiS
	  jacobian = 2.0 * Xseed[3];
	  SetVariables(Xseed[0], Xseed[1], Xseed[2], Xseed[3], Xseed[4], Xseed[5]);
	}
	else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
	  double y0 = Xseed[1] / (2.0 * Xseed[0] * (PP * Pl));
	  jacobian = 2.0 * Xseed[3] * y0 / Xseed[1];
	  SetVariables(Xseed[0], y0, Xseed[2], Xseed[3], Xseed[4], Xseed[5]);
	}
	CalculateFinalState();
	Xhisto[ih]->SetBinContent(ib, dsigma(mode) * jacobian);
      }
      do {//Sample and double check
	if (Xhisto[ih]->Integral(1,-1) > 0){
	  Xseed[ih] = Xhisto[ih]->GetRandom();
	}
	else {
	  Xseed[ih] = gRandom->Uniform(Xmin[ih], Xmax[ih]);
	  //printf("Warning: GibbsSampler to zero point!\n");
	}
	if (method == 0){//generate in x, y, z, Pt, phih, phiS
	  jacobian = 2.0 * Xseed[3];
	  SetVariables(Xseed[0], Xseed[1], Xseed[2], Xseed[3], Xseed[4], Xseed[5]);
	}
	else if (method == 1){//generate in x, Q2, z, Pt, phih, phiS
	  double y0 = Xseed[1] / (2.0 * Xseed[0] * (PP * Pl));
	  jacobian = 2.0 * Xseed[3] * y0 / Xseed[1];
	  SetVariables(Xseed[0], y0, Xseed[2], Xseed[3], Xseed[4], Xseed[5]);
	}
	CalculateFinalState();
      } while(dsigma(mode) < 1.e-20);
    }
  }
  return sigmatotal;
}

int Lsidis::Test(){//Test code
  std::cout << x << "  " << y << "  " << z << "  " << Pt << "  " << phih << "  " << phiS << std::endl;
  return 0;
}

#endif

/* update log from v1.1 to 2.0
   implement Gibbs sampler for sidis event generator
   
   update log from v2.0 to 2.1
   fix pid bug

   update log from v2.2 to 2.3
   move static const double variables outside class definition

   update log from v2.3 to 3.0
   implement R1 Monte Carlo sampling follow Ted Rogers' method

*/


// Local Variables:
// mode: c++
// End:
