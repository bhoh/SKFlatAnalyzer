#ifndef AFBAnalyzer_h
#define AFBAnalyzer_h

#include "SMPAnalyzerCore.h"

class AFBAnalyzer : public SMPAnalyzerCore {

public:

  void executeEvent();
  void executeEventFromParameter(TString channelname,Event* ev);

  AFBAnalyzer();
  ~AFBAnalyzer();

  double GetCosThetaCS(const vector<Lepton*>& leps);
  //double GetCosThetaA(const vector<Lepton*>& leps);
  //double GetCosThetaB(const vector<Lepton*>& leps,const vector<Jet>& jets);
  //double GetCosThetaC(const vector<Lepton*>& leps,const vector<Jet>& jets);
  double GetCosTheta(const vector<Lepton*>& leps,const vector<Jet>& jets,TString option,double fcut);
  void FillAFBHists(TString pre,TString suf,const vector<Lepton*>& leps,const vector<Jet>& jets,double w);
  void FillAFBSystematicHists(TString pre,TString suf,const vector<Lepton*>& leps,const vector<Jet>& jets,map<TString,double> map_systematic);
  //static const int massbinnum=12;
  //const double massrange[massbinnum+1]={60,70,78,84,87,89,91,93,95,98,104,112,120};
  static const int massbinnum=14;
  const double massrange[massbinnum+1]={60,70,78,84,87,89,91,93,95,98,104,112,120,200,350};
  void FillHardHists(TString pre,TString suf,const Gen& genparton0,const Gen& genparton1,const Gen& genhardl0,const Gen& genhardl1,const Gen& genhardj0,double w);
  void FillGenAFBHists(TString pre,TString suf,const Gen& genl0,const Gen& genl1,const Gen& genphotons,double w);
  TRandom3* random;
  TString hardprefix;
};



#endif

