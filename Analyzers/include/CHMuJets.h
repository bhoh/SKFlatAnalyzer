#ifndef CHMuJets_h
#define CHMuJets_h

#include "AnalyzerCore.h"

class CHMuJets : public AnalyzerCore {

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

  CHMuJets();
  ~CHMuJets();

};



#endif

