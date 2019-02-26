#ifndef FitterInput_h
#define FitterInput_h

#include <algorithm>
#include <numeric>
#include "AnalyzerCore.h"

class FitterInput : public AnalyzerCore {

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

  FitterInput();
  ~FitterInput();

};



#endif

