#ifndef TSCorr_input_h
#define TSCorr_input_h

#include "AnalyzerCore.h"

class TSCorr_input : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<TString> MuonIDs, MuonIDSFKeys;
  vector<Muon> AllMuons;
  vector<Electron> AllElectrons;
  vector<Jet> AllJets;
  vector<Gen> AllGens;

  double weight_Prefire;

  struct matchedPartonJet{
    Int_t truth_index;
    TLorentzVector matched_parton;
    TLorentzVector matched_jet;
  };

  TSCorr_input();
  ~TSCorr_input();

};



#endif

