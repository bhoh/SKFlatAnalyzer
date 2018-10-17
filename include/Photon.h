#ifndef Photon_h
#define Photon_h

#include "Particle.h"

class Photon: public Particle {
public:

  Photon();
  ~Photon();
  
  void SetSC(double sceta, double scphi);
  inline double scEta() const { return j_scEta; }
  inline double scPhi() const { return j_scPhi; }
    
  
  void SetCutBasedIDVariables(
    double Full5x5_sigmaIetaIeta,
    double HoverE,
    double ChIsoWithEA,
    double NhIsoWithEA,
    double PhIsoWithEA
  ){
    j_Full5x5_sigmaIetaIeta = Full5x5_sigmaIetaIeta;
    j_HoverE = HoverE;
    j_ChIsoWithEA = ChIsoWithEA;
    j_NhIsoWithEA = NhIsoWithEA;
    j_PhIsoWithEA = PhIsoWithEA;
  }
  double Full5x5_sigmaIetaIeta() const { return j_Full5x5_sigmaIetaIeta; }
  double HoverE() const { return j_HoverE; }
  double ChIsoWithEA() const { return j_ChIsoWithEA; }
  double NhIsoWithEA() const { return j_NhIsoWithEA; }
  double PhIsoWithEA() const { return j_PhIsoWithEA; }
  

  void SetPOGIDs(std::vector<bool> bs);
  inline bool passLooseID() const { return j_passLooseID; }
  inline bool passMediumID() const { return j_passMediumID; }
  inline bool passTightID() const { return j_passTightID; }
  inline bool passMVAID_WP80() const { return j_passMVAID_WP80; }
  inline bool passMVAID_WP90() const { return j_passMVAID_WP90; }
  
  
  // === ID
  bool PassID(TString id);
  bool Pass_CutBasedLoose();
  bool Pass_CutBasedMedium();
  bool Pass_CutBasedTight();
  
  bool Pass_Prefire_ID();



  void SetRho(double r){ j_Rho = r; }
  double Rho() const { return j_Rho; }
  

private:
  double j_scEta;
  double j_scPhi;
  double j_scE;
  double j_Full5x5_sigmaIetaIeta;
  double j_HoverE;
  double j_ChIsoWithEA;
  double j_NhIsoWithEA;
  double j_PhIsoWithEA;
  bool j_passLooseID;
  bool j_passMediumID;
  bool j_passTightID;
  bool j_passMVAID_WP80;
  bool j_passMVAID_WP90;

  double j_Rho;

};

Photon::Photon() : Particle() {
  j_scEta = -999.;
  j_scPhi = -999.;
  j_scE = -999.;
  j_passLooseID = false;
  j_passMediumID = false;
  j_passTightID = false;
  j_passMVAID_WP80 = false;
  j_passMVAID_WP90 = false;
  j_Rho = -999.;
}

Photon::~Photon(){

}

void Photon::SetSC(double sceta, double scphi){
  j_scEta = sceta;
  j_scPhi = scphi;
}

void Photon::SetPOGIDs(std::vector<bool> bs){
  j_passLooseID = bs.at(0);
  j_passMediumID = bs.at(1);
  j_passTightID = bs.at(2);
  j_passMVAID_WP80 = bs.at(3);
  j_passMVAID_WP90 = bs.at(4);
  
}

bool Photon::PassID(TString id){
  
  if(id=="passLooseID") return passLooseID();
  if(id=="passMediumID") return passMediumID();
  if(id=="passTightID") return passTightID();
  if(id=="passMVAID_WP80") return passMVAID_WP80();
  if(id=="passMVAID_WP90") return passMVAID_WP90();
  return false;

}


#endif
