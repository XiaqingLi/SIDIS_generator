#include "../include/Lsidis.h"
#include "../include/stru_func.h"


using namespace std;


int main(int argc, char * argv[]){

  if (argc < 12){
    cout << "./generator <Nsim> <Nfile> <electron beam energy in GeV> <ion beam energy in GeV> <xmin> <xmax> <z_min> <z_max> <Pt_min> <Pt_max> <Q2_min> <Q2_max> <path>" << endl;
    return 0;
  }

  int Nsim = atoi(argv[1]);
  int Nfile = atoi(argv[2]);
  double E_e = atof(argv[3]);
  double E_p = atof(argv[4]);
  double kx_min = atof(argv[5]);
  double kx_max = atof(argv[6]);
  double kz_min = atof(argv[7]);
  double kz_max = atof(argv[8]);
  double Pt_min = atof(argv[9]);
  double Pt_max = atof(argv[10]);
  double kQ2_min = atof(argv[11]);
  double kQ2_max = atof(argv[12]);
  TString path = argv[13];

  gRandom->SetSeed(0);

  TLorentzVector P(0, 0, sqrt(E_p*E_p - 0.938272*0.938272), E_p);//Nucleon, along z direction
  //TLorentzVector P(0, 0, 0, 0.938272);//Nucleon, fixed target
  TLorentzVector l(0, 0, -E_e, E_e);//electron, along -z direction
  Lsidis mysidis;
  mysidis.SetNucleus(1, 0);//proton number, neutron number
  //mysidis.SetNucleus(0.334*10+0.593*2, 0.334*7+0.593*2);//old SoLID NH3, use pol lumi

  /////////new SoLID NH3 (NH3 lumi = 4.785e35, 4He lumi = (0.69+0.47)e35), use unpol lumi
  //mysidis.SetNucleus(4.785*10/(17*5.945)+1.16*2/(4*5.945), 4.785*7/(17*5.945)+1.16*2/(4*5.945));

  ////////CLAS12 HD-ice (C12-11-111, 78% HD + 15% Al + 7% C2ClF3)
  //mysidis.SetNucleus(2+0.021*13, 1+0.021*14);//pol limi
  //mysidis.SetNucleus(2/3.567+0.021*13/3.567, 1/3.567+0.021*14/3.567);// unpol lumi

  mysidis.SetInitialState(l, P);//set initial state
  mysidis.SetHadron("pi+");//set detected hadron: pi+, pi-, pi0, K+, K-, K0, rho+, rho-, rho0, p
  mysidis.SetFUUTmode(0);//choose a mode for structure function, 0=Gauss, 1=grid

 
  TLorentzVector lp, Ph;//scattered electron and detected hadron
  double weight;//weight of event
  
  double Xmin[6] = {kx_min, kQ2_min, kz_min, Pt_min, -M_PI, -M_PI};//generator kinematic range x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {kx_max, kQ2_max, kz_max, Pt_max, M_PI, M_PI};
  // double Xmin[6] = {0.001, 1, 0.3, 0.0, -M_PI, -M_PI};//for SoLID
  // double Xmax[6] = {1.000, 50, 0.7, 2.0, M_PI, M_PI};//for SoLID
  // Xmin[0] = StructureFunction::xmin;
  // Xmin[1] = StructureFunction::Q2min;
  // Xmin[2] = StructureFunction::zmin;
  // Xmax[0] = StructureFunction::xmax;
  // Xmax[1] = StructureFunction::Q2max;
  // Xmax[2] = StructureFunction::zmax;

  
  mysidis.SetRange(Xmin, Xmax);
  
  double Q2, x, y, z, W, Wp, phih, phis, Pt;
  FILE * fout;
  for (int j = 0; j < Nfile; j++){
    cout << "File #" << j << endl;
    TString str_tmp = Form("1e%d_%.4d.dat",(int)log10(Nsim),j);
    str_tmp = path+str_tmp;
    fout = fopen(str_tmp.Data(),"w");

    fprintf(fout, "Number of events = %d, %d+%d GeV, x = [%.6f,%.2f], Q2 = [%.3f,%.2f], z = [%.2f,%.2f], Pt = [%.2f,%.2f]\n", Nsim, (int)E_e, (int)E_p, Xmin[0],Xmax[0],Xmin[1],Xmax[1],Xmin[2],Xmax[2],Xmin[3],Xmax[3]);
    fprintf(fout, "weight\tQ2\tx\ty\tz\tW\tWp\tphih\tphis\tPt\tlpE\tlp\tlptheta\tlpphi\tPhE\tPh\tPhthete\tPhphi\n");

    for (int i = 0; i < Nsim; i++){
      // if(i%1000000==0) 
      //   printf("%d events completed\n", i);      
      
      //cout << "\r"<<i/(Nsim/Nfile/100) << "%";
      weight = mysidis.GenerateEvent(0, 1);//generate an event (0,0)->in y, (0,1)->in Q2

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
