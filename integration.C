

{
TFile *file = new TFile("hists_f1_D1.root");
TH3D *f1u = (TH3D*)file->Get("f1u");
TH3D *D1u = (TH3D*)file->Get("D1u");

double x,z,Pt,Q2;


double val_u = 0;
double val_d = 0;

double k0 = 0.005;
double step = 0.01;
for(int i=0; i<100; i++){
double val_f1 = f1u->Intepolate(x, (k0 + i * step) * (k0 + i * step), Q2);
double val_D1 = D1u->Intepolate(z, Pt * Pt + z * z * (k0 + i * step) * (k0 + i * step), Q2);
val_u += val_f1 + val_D1;

double val_f1 = f1d->Intepolate(x,(k0+i*step)*(k0+i*step),Q2);
double val_D1 = D1d->Intepolate(z,Pt*Pt+z*z*(k0+i*step)*(k0+i*step),Q2);
val_d += val_f1 + val_D1;
}
val_u = val*x*eu2;
val_d = val*x*ed2;


return val_u + val_d;
}

double Tint(const double * kk, const double * par){//kk: theta, phi
  //par (every 3): p, pd
  TLorentzVector p, pd;
  p.SetXYZM(par[0] * sin(par[1]) * cos(par[2]), par[0] * sin(par[1]) * sin(par[2]), par[0] * cos(par[1]), MN);
  pd.SetXYZM(par[3] * sin(par[4]) * cos(par[5]), par[3] * sin(par[4]) * sin(par[5]), par[3] * cos(par[4]), Md);
  int L1 = (int) par[6];
  int L2 = (int) par[7];
  double result = Tk(kk, p, pd, L1, L2) * sin(kk[0]);
  return result;
}

double Tfi(TLorentzVector p, TLorentzVector pd, double L1, double L2){
  double par[8] = {p.P(), p.Theta(), p.Phi(), pd.P(), pd.Theta(), pd.Phi(), L1, L2};
  double xl[3] = {0, -M_PI};//integration lower bound
  double xu[3] = {M_PI, M_PI};//integration upper bound
  ROOT::Math::WrappedParamFunction<> wf(&Tint, 2, 8);
  wf.SetParameters(par);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0.0, 0.001, 100000);
  ig.SetFunction(wf);
  double result = ig.Integral(xl, xu);
  return result;
}

