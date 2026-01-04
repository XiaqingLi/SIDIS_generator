{
    double size = 0.05;
    gStyle->SetTitleSize(size,"X");
    gStyle->SetLabelSize(size,"X");
    gStyle->SetTitleSize(size,"Y");
    gStyle->SetLabelSize(size,"Y");
    gStyle->SetTitleSize(size,"Z");
    gStyle->SetLabelSize(size,"Z");
    gStyle->SetOptStat(0);
    
    
    const Int_t nbins = 100;
    double binwidth;
    double xmin = 1e-5;
    double xmax = 1;
    double logxmin = log10(xmin);
    double logxmax = log10(xmax);
    binwidth = (logxmax-logxmin)/nbins;
    double xbins[nbins+1];
    xbins[0] = xmin;
    for (Int_t i=1;i<=nbins;i++)
    xbins[i] = pow(10,logxmin+i*binwidth);
    double Q2min = 1;
    double Q2max = 1e4;
    double logQ2min = log10(Q2min);
    double logQ2max = log10(Q2max);
    binwidth = (logQ2max-logQ2min)/nbins;
    double Q2bins[nbins+1];
    Q2bins[0] = log10(Q2min);
    for (Int_t i=1;i<=nbins;i++)
    Q2bins[i] = pow(10,logQ2min+i*binwidth);
    
    
    TH2D *hQ2x = new TH2D("h_Q2_x", ";x;Q^{2} (GeV^{2})", nbins, xbins, nbins, Q2bins);
    
    double weight, Q2, x, y, z, W, Wp, phih, phis, Pt, lpE, lp, lptheta, lpphi, PhE, Ph, Phtheta, Phphi, rapid_e, rapid_h;
    
    TCanvas *c = new TCanvas();
    c->Divide(2,2);

    int Ee[4] = {5, 5, 10, 18};
    int Ep[4] = {41, 100, 100, 275};
    for(int n=0; n<4; n++){
        //    ifstream data("../output_mode0_CJ15lo/piplus_e5_p41_Pt0.4_0.6_z0.6_0.7_Q2_1_100_1e4_0009.dat");
        ifstream data(Form("../output_mode0_CJ15lo/piplus_e%d_p%d_1e7_log.dat", Ee[n], Ep[n]));
        string line;
        getline(data, line);
        getline(data, line);
        while(data>>weight>>Q2>>x>>y>>z>>W>>Wp>>phih>>phis>>Pt>>lpE>>lp>>lptheta>>lpphi>>PhE>>Ph>>Phtheta>>Phphi){
            double rapid_e = -log(tan(lptheta/2));
            double rapid_h = -log(tan(Phtheta/2));
            if(weight<=0) continue;
            if(Q2<1) continue;
            if(W*W<10) continue;
            if(rapid_e<-3.5 || rapid_e>3.5 || rapid_h<-3.5 || rapid_h>3.5) continue;
            if(y<=0.05 || y>=0.95) continue;
            if(z<0.3 || z>0.7) continue;
            if(Pt>1) continue;
            
            hQ2x->Fill(x,Q2,weight);
            
        }
        hQ2x->GetXaxis()->SetLimits(xmin, xmax);
        hQ2x->GetYaxis()->SetLimits(Q2min, Q2max);
        hQ2x->GetZaxis()->SetRangeUser(0.001, 1e5);
        
        c->cd(n+1);
        hQ2x->DrawCopy("colz");
        gPad->SetRightMargin(0.13);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        
        hQ2x->Reset();
    }
}
