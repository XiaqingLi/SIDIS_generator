#include <iostream>
#include <fstream>
#include <cmath>
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"
#include "TGaxis.h"
#include "TText.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLegend.h"

using namespace std;

void DrawGraph(double logx_l, double logx_h, double logQ2_l, double logQ2_h, TGraphErrors *gre, int markerstyle, Color_t color, double size);


int main(){
  gStyle->SetOptStat(0);
  static const int kNconfig = 3;
  static const int kNPt = 2;
  static const int kNz = 4;
  const int kEe[kNconfig] = {/*10,*/ 5, 10, 5};
  const int kEp[kNconfig] = {/*275,*/ 100, 100, 41};
  const Color_t kColor[kNconfig] = {kRed+1, kGreen+2, kBlack};
  const double kz_min[kNz] = {0.30, 0.40, 0.50, 0.60};
  const double kz_max[kNz] = {0.40, 0.50, 0.60, 0.70};
  const double kPt_min[kNPt] = {0.2, 0.4}; 
  const double kPt_max[kNPt] = {0.4, 0.6}; 
  static const int kNQ2 = 5;
  double lgx_l = 0.002, lgx_h = 0.7;
  double lgQ2_l = -0.2, lgQ2_h = 2.15;

  double size = 0.08;
  double offset = 1.1;

  double scale = 10;

  TLegend *ll;
  TGraphErrors *gre[kNconfig][kNQ2], *gre_s11_NH3[kNconfig][kNQ2], *gre_s8_NH3[kNconfig][kNQ2], *gre_s11_He3[kNconfig][kNQ2], *gre_s8_He3[kNconfig][kNQ2];
  TFile *file;
  TFile *file_s11_NH3 = new TFile("SoLID_projection/gre_enhanced_11GeV_P_pi+_linearX.root");
  TFile *file_s8_NH3 = new TFile("SoLID_projection/gre_enhanced_8.8GeV_P_pi+_linearX.root");
  TFile *file_s11_He3 = new TFile("SoLID_projection/gre_enhanced_11GeV_N_pi+_linearX.root");
  TFile *file_s8_He3 = new TFile("SoLID_projection/gre_enhanced_8.8GeV_N_pi+_linearX.root");

  TCanvas *c = new TCanvas("c","c",2400,1000);
  c->Divide(kNz,kNPt,0,0);
  TH2F *hh = new TH2F("h","", 1,lgx_l,lgx_h, 1,lgQ2_l,lgQ2_h);
  hh->SetStats(0);
  hh->GetXaxis()->SetAxisColor(0);
  hh->GetXaxis()->SetLabelColor(0);
  hh->GetYaxis()->SetAxisColor(0);
  hh->GetYaxis()->SetLabelColor(0);
  TGaxis *axis_Q2 = new TGaxis(lgx_l,lgQ2_l,lgx_l,lgQ2_h,pow(10,lgQ2_l),pow(10,lgQ2_h),510,"G");
  axis_Q2->SetTitle("Q^{2} (GeV^{2})");
  axis_Q2->SetTitleSize(size);
  axis_Q2->SetLabelSize(size);
  axis_Q2->SetTitleOffset(offset);
  axis_Q2->SetTitleFont(62);
  TGaxis *axis_x = new TGaxis(lgx_l,lgQ2_l,lgx_h,lgQ2_l,lgx_l, lgx_h,505,"L");
  axis_x->SetTitle("x");
  axis_x->SetTitleSize(size);
  axis_x->SetLabelSize(size);
  axis_x->SetTitleOffset(offset-0.1);
  TGaxis *axis_asy = new TGaxis(lgx_h,lgQ2_l,lgx_h,lgQ2_h,(lgQ2_l/*-1*/)/scale,(lgQ2_h/*-1*/)/scale,505,"+L");
  axis_asy->SetTitle("SSA #pi^{+}");
  axis_asy->SetTitleSize(size);
  axis_asy->SetLabelSize(size);
  axis_asy->SetTitleOffset(offset+0.2);
  TLatex *tt1, *tt1p, *tt2;

  
  for(int j=0; j<kNPt; j++){
    for(int k=0; k<kNz; k++){
      c->cd(j*kNz+k+1);
      gPad->SetTopMargin(0);
      gPad->SetBottomMargin(0);
      gPad->SetLeftMargin(0);
      gPad->SetRightMargin(0);
      if(j==kNPt-1) gPad->SetBottomMargin(0.17);
      if(k==0) gPad->SetLeftMargin(0.2);
      if(k==kNz-1) gPad->SetRightMargin(0.2);

      hh->Draw("axis");
      
      for(int ii=0; ii<kNconfig; ii++){
	if(ii!=2) continue;//choose the EIC e-p energies
	file = new TFile(Form("EIC_e%d_p%d_10fb-1.root",kEe[ii],kEp[ii])); //contains EIC projections in TGraphErrors
	for(int i=0; i<kNQ2; i++){ 
	  gre[ii][i] = (TGraphErrors*)file->Get(Form("gre_Pt%.1f_%.1f_z%.2f_%.2f_Q2_%d",kPt_min[j],kPt_max[j],kz_min[k],kz_max[k],i));
	  int Ngre = gre[ii][i]->GetN();
	  for(int n_tmp=0; n_tmp<Ngre; n_tmp++){
	    if(gre[ii][i]->GetErrorY(n_tmp)/scale > 0.1){
	      gre[ii][i]->SetPoint(n_tmp, 0, -10);
	      gre[ii][i]->SetPointError(n_tmp, 0, 0);
	    }
	  }
	  DrawGraph(lgx_l, lgx_h, lgQ2_l, lgQ2_h, gre[ii][i], 25, kColor[ii], 1.5);
	} 
      }

      gre_s11_NH3[j][k] = (TGraphErrors*)file_s11_NH3->Get(Form("gre_%d_%d",j+1,k));
      for(int n_tmp=0; n_tmp<gre_s11_NH3[j][k]->GetN(); n_tmp++){
      	if(gre_s11_NH3[j][k]->GetErrorY(n_tmp)/scale > 0.02){
      	  gre_s11_NH3[j][k]->SetPoint(n_tmp, 0, -10);
      	  gre_s11_NH3[j][k]->SetPointError(n_tmp, 0, 0);
      	}
      }
      gre_s8_NH3[j][k] = (TGraphErrors*)file_s8_NH3->Get(Form("gre_%d_%d",j+1,k));
      for(int n_tmp=0; n_tmp<gre_s8_NH3[j][k]->GetN(); n_tmp++){
      	if(gre_s8_NH3[j][k]->GetErrorY(n_tmp)/scale > 0.02){
      	  gre_s8_NH3[j][k]->SetPoint(n_tmp, 0, -10);
      	  gre_s8_NH3[j][k]->SetPointError(n_tmp, 0, 0);
      	}
      }

      gre_s11_He3[j][k] = (TGraphErrors*)file_s11_He3->Get(Form("gre_%d_%d",j+1,k));
      for(int n_tmp=0; n_tmp<gre_s11_He3[j][k]->GetN(); n_tmp++){
      	if(gre_s11_He3[j][k]->GetErrorY(n_tmp)/scale > 0.02){
      	  gre_s11_He3[j][k]->SetPoint(n_tmp, 0, -10);
      	  gre_s11_He3[j][k]->SetPointError(n_tmp, 0, 0);
      	}
      }
      gre_s8_He3[j][k] = (TGraphErrors*)file_s8_He3->Get(Form("gre_%d_%d",j+1,k));
      for(int n_tmp=0; n_tmp<gre_s8_He3[j][k]->GetN(); n_tmp++){
      	if(gre_s8_He3[j][k]->GetErrorY(n_tmp)/scale > 0.02){
      	  gre_s8_He3[j][k]->SetPoint(n_tmp, 0, -10);
      	  gre_s8_He3[j][k]->SetPointError(n_tmp, 0, 0);
      	}
      }


      /////////// choose which SoLID target to draw
      //DrawGraph(lgx_l, lgx_h, lgQ2_l, lgQ2_h, gre_s11_NH3[j][k], 22, kBlue, 2);
      //DrawGraph(lgx_l, lgx_h, lgQ2_l, lgQ2_h, gre_s8_NH3[j][k], 22, kRed+1, 2);
      DrawGraph(lgx_l, lgx_h, lgQ2_l, lgQ2_h, gre_s11_He3[j][k], kFullCircle, kBlue, 1.7);
      DrawGraph(lgx_l, lgx_h, lgQ2_l, lgQ2_h, gre_s8_He3[j][k], kFullCircle, kRed+1, 1.7);

      if(j==1 && k==2){
	ll = new TLegend(.35, .7, 1, 1);
	//ll->AddEntry(gre_s8_NH3[j][k],"SoLID 8.8 GeV (NH_{3})  ","p");
	//ll->AddEntry(gre_s11_NH3[j][k],"SoLID 11 GeV (NH_{3})","p");
	ll->AddEntry(gre_s8_He3[j][k],"SoLID 8.8 GeV (He^{3})  ","p");
	ll->AddEntry(gre_s11_He3[j][k],"SoLID 11 GeV (He^{3})","p");
 	ll->AddEntry(gre[2][0],"EIC e-p #sqrt{s} = 29 GeV","p");
	ll->SetBorderSize(0);
	ll->SetFillStyle(0);
	ll->Draw();
      }

      axis_Q2->Draw();  
      axis_x->Draw();  
      if(k==kNz-1) 
      axis_asy->Draw();  

      tt1 = new TLatex(0.35,1.7,Form("%.1f GeV < P_{T} < %.1f GeV", kPt_min[j], kPt_max[j]));
      tt1->SetTextAlign(22);
      tt1->SetTextSize(0.09);
      if(k==0 && j==0) tt1->Draw();
      tt1p = new TLatex(0.35,1.9,Form("%.1f GeV < P_{T} < %.1f GeV", kPt_min[j], kPt_max[j]));
      tt1p->SetTextAlign(22);
      tt1p->SetTextSize(0.08);
      if(k==0 && j!=0) tt1p->Draw();
      tt2 = new TLatex(0.35,1.9,Form("%.2f < z < %.2f", kz_min[k], kz_max[k]));
      tt2->SetTextAlign(22);
      tt2->SetTextSize(0.09);
      if(j==0) tt2->Draw();
    }
  }

  c->SaveAs("SoLID_3He_20210217.pdf");
  return 0;
}

void DrawGraph(double logx_l, double logx_h, double logQ2_l, double logQ2_h, TGraphErrors *gre, int markerstyle, Color_t color, double size){
  gre->GetXaxis()->SetRangeUser(logx_l, logx_h);
  gre->GetYaxis()->SetRangeUser(logQ2_l, logQ2_h);
  gre->GetXaxis()->SetAxisColor(kRed);
  gre->GetXaxis()->SetLabelColor(kRed);
  gre->GetYaxis()->SetAxisColor(kRed);
  gre->GetYaxis()->SetLabelColor(kRed);
  gre->SetMarkerColor(color);
  gre->SetLineColor(color);
  gre->SetMarkerStyle(markerstyle);
  gre->SetMarkerSize(size);
  gre->SetLineWidth(1);
  gre->Draw("psame");
}




