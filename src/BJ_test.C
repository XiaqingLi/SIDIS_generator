#include "../include/Lsidis.h"
#include "../include/stru_func.h"

using namespace std;

int main(int argc, char * argv[]){
    
  TLorentzVector P(0, 0, 0, 0.938272);//Nucleon 
  TLorentzVector l(0, 0, 10.6, 10.6);//electron 
  
  Lsidis mysidis;
  mysidis.SetNucleus(1, 0);//proton number, neutron number
  mysidis.SetInitialState(l, P);//set initial state
  mysidis.SetHadron("p");//set detected hadron: pi+, pi-, pi0, K+, K-, K0, p
  mysidis.SetFUUTmode(0);//choose a mode for structure function, 0=Gauss, 1=grid


  TLorentzVector lp, Ph;//scattered electron and detected hadron
  double weight;//weight of event
  
  double Xmin[6] = {0.01, 0.9, 0.05, 0.0, -M_PI, -M_PI};//generator kinematic range x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {0.999, 2.0, 0.999, 2.0, M_PI, M_PI};
  
  mysidis.SetRange(Xmin, Xmax);
  
  double x, Q2, z, Pt;
  for (int i = 0; i < 10; i++){
    weight = mysidis.GenerateEvent(0, 1);//generate an event
    if (weight > 0){
      lp = mysidis.GetLorentzVector("lp");
      Ph = mysidis.GetLorentzVector("Ph");
      x = mysidis.GetVariable("x");
      Q2 = mysidis.GetVariable("Q2");
      if (Q2 < 1.0) continue;
      z = mysidis.GetVariable("z");
      Pt = mysidis.GetVariable("Pt");

      cout << i << "\t" << x << "\t" << Q2 << "\t" << z << "\t" << Pt << endl;
    }
  }

  return 0;
}
