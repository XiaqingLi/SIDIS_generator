#include "../Lsidis.h"

using namespace std;


int main(int argc, char * argv[]){

  if (argc < 3){
    cout << "./generator <Nsim> <Nfile>" << endl;
    return 0;
  }

  int Nsim = atoi(argv[1]);
  int Nfile = atoi(argv[2]);

  gRandom->SetSeed(0);
  
  TLorentzVector P(0, 0, -sqrt(60.0*60.0 - 0.938272*0.938272), 60.0);//Nucleon
  TLorentzVector l(0, 0, 11.0, 11.0);//electron
  Lsidis mysidis;
  mysidis.SetNucleus(1, 0);//proton number, neutron number
  mysidis.SetInitialState(l, P);//set initial state
  mysidis.SetHadron("pi+");//set detected hadron: pi+, pi-, pi0, K+, K-, K0, p
  mysidis.SetPDFset("CJ15lo");//choose a PDF set
  mysidis.SetFFset("DSSFFlo");//choose an FF set
  mysidis.SetSFmode(1); //chose a mode for structure function

  
  TLorentzVector lp, Ph;//scattered electron and detected hadron
  double weight;//weight of event
  
  double Xmin[6] = {0.0001, 0.05, 0.2, 0.0, -M_PI, -M_PI};//generator kinematic range x, y, z, Pt, phih, phiS
  double Xmax[6] = {1.0, 0.8, 0.8, 5.0, M_PI, M_PI};
  
  mysidis.SetRange(Xmin, Xmax);
  
  double Q2, x, y, z, W, Wp, phih, phis, Pt;
  FILE * fout;
  for (int j = 0; j < Nfile; j++){
    cout << "File #" << j << endl;
    fout = fopen(Form("../output/test_%.4d.dat", j),"w");
    fprintf(fout, "Number of events = %d, x = [%f,%f], Q2 = [%f,%f]\n", Nsim, Xmin[0],Xmax[0],Xmin[1],Xmax[1]);
    fprintf(fout, "weight\tQ2\tx\ty\tz\tW\tWp\tphih\tphis\tPt\tlpE\tlp\tlptheta\tlpphi\tPhE\tPh\tPhthete\tPhphi\n");

    for (int i = 0; i < Nsim/Nfile; i++){
      // if(i%(Nsim/100)==0) {
      // 	cout << i/(Nsim/100) << "%" <<endl;
      weight = mysidis.GenerateEvent(0, 0);//generate an event

      if (weight > 0){
	Q2 = mysidis.GetVariable("Q2");
	x = mysidis.GetVariable("x");
	y = mysidis.GetVariable("y");
	z = mysidis.GetVariable("z");
	W = mysidis.GetVariable("W");
	Wp = mysidis.GetVariable("Wp");
	phih = mysidis.GetVariable("phih");
	phis = mysidis.GetVariable("phiS");
	Pt = mysidis.GetVariable("Pt");
	lp = mysidis.GetLorentzVector("lp");
	Ph = mysidis.GetLorentzVector("Ph");
	  
	fprintf(fout, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		weight, Q2, x, y, z, W, Wp, phih, phis, Pt,
		lp.E(), lp.P(), lp.Theta(), lp.Phi(),
		Ph.E(), Ph.P(), Ph.Theta(), Ph.Phi());
	  

      }
      
    }
    fclose(fout);
  }
  
  return 0;
}
