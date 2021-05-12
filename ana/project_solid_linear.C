#include "include/all.h"

int main(){
  static const int kNPt = 5;
  static const int kNz = 4;
  const double kz_min[kNz] = {0.30, 0.40, 0.50, 0.60};
  const double kz_max[kNz] = {0.40, 0.50, 0.60, 0.70};
  const double kPt_min[kNPt] = {0.0, 0.2, 0.4, 0.6, 0.8}; 
  const double kPt_max[kNPt] = {0.2, 0.4, 0.6, 0.8, 1.0}; 
  const double scale = 10;

  TFile *file = new TFile("SoLID_projection/binsP11p.root");
  TTree *tree = (TTree*)file->Get("data");
  double x, y, z, Q2, Pt, stat;
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("Q2", &Q2);  
  tree->SetBranchAddress("Pt", &Pt);
  tree->SetBranchAddress("E0stat", &stat);

  int nevents = tree->GetEntries();
  TGraphErrors *gre[kNPt][kNz];

  for(int j=0; j<kNPt; j++){
    for(int k=0; k<kNz; k++){
      int count = 0;
      gre[j][k] = new TGraphErrors();
      for(int i=0; i<nevents; i++){
	tree->GetEntry(i);
	if(z>kz_min[k] && z<kz_max[k] && Pt>kPt_min[j] && Pt<kPt_max[j]){
	  gre[j][k]->SetPoint(count, x, log10(Q2));
	  gre[j][k]->SetPointError(count, 0, stat*scale);
	  count++;
	}
      }
    }	   
  }

  TFile *fout = new TFile("SoLID_projection_11GeVpiplus_P_full_linear.root","recreate");
  for(int j=0; j<kNPt; j++){
    for(int k=0; k<kNz; k++){
      gre[j][k]->SetNameTitle(Form("gre_%d_%d",j,k),Form("11 GeV, %.1f<Pt<%.1f, %.2f<z<%.2f",kPt_min[j],kPt_max[j],kz_min[k],kz_max[k]));
      gre[j][k]->Write();
    }
  }
  fout->Close();
}
