#include "CHMuJets.h"

CHMuJets::CHMuJets(){

}

void CHMuJets::initializeAnalyzer(){

  //==== if you've done "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");

  cout << "[CHMuJets::initializeAnalyzer] RunSyst = " << RunSyst << endl;

  //==== In the following examples, I will obtain dimuon Z-peak with two muon IDs
  //==== I defined "vector<TString> MuonIDs;" in Analyzers/include/CHMuJets.h
  MuonIDs = {
    "POGTight"
  };
  //==== corresponding Muon ID SF Keys for mcCorr.MuonID_SF()
  MuonIDSFKeys = {
    "NUM_TightID_DEN_genTracks",
  };

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/CHMuJets.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro
  if(DataYear==2016){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    TriggerSafePtCut = 26.;
  }
  else if(DataYear==2017){
    IsoMuTriggerName = "HLT_IsoMu27_v";
    TriggerSafePtCut = 29.;
  }

  cout << "[CHMuJets::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  cout << "[CHMuJets::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== add taggers and WP that you want to use in analysis                                                                                                                                                
  std::vector<Jet::Tagger> vtaggers;
  vtaggers.push_back(Jet::DeepCSV);

  std::vector<Jet::WP> v_wps;
  v_wps.push_back(Jet::Medium);

  //=== list of taggers, WP, setup systematics, use period SFs                                                                                                                                              
  SetupBTagger(vtaggers,v_wps, true, true);


}

CHMuJets::~CHMuJets(){

  //==== Destructor of this Analyzer

}

void CHMuJets::executeEvent(){

  //================================================================
  //====  (Example)
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/CHMuJets.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  //=== Jets too
  AllJets = GetAllJets();
  AllGens = GetGens();

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/CHMuJets.h
  weight_Prefire = GetPrefireWeight(0);

  //==== Declare AnalyzerParameter

  AnalyzerParameter param;

  //==== Loop over muon IDs

  for(unsigned int it_MuonID=0; it_MuonID<MuonIDs.size(); it_MuonID++){

    TString MuonID = MuonIDs.at(it_MuonID);
    TString MuonIDSFKey = MuonIDSFKeys.at(it_MuonID);

    //==== 1) First, let's run Central values of the systematics

    //==== clear parameter set
    param.Clear();

    //==== set which systematic sources you want to run this time
    //==== default syst_ is AnalyzerParameter::Central
    param.syst_ = AnalyzerParameter::Central;

    //==== set name of the parameter set
    //==== this will be used for the directory name of histograms
    param.Name = MuonID+"_"+"Central";

    //==== You can define lepton ID string here
    param.Muon_Tight_ID = MuonID;
    param.Muon_ID_SF_Key = MuonIDSFKey;
    param.Electron_Tight_ID = "passTightID";
    param.Electron_Veto_ID = "passVetoID";

    //==== And, Jet ID
    param.Jet_ID = "tightLepVeto"; 

    //==== Now, all parameters are set. Run executeEventFromParameter() with this parameter set
    executeEventFromParameter(param);

    //==== 2) Now, loop over systematic sources
    //==== without --userflag RunSyst, this will not be ran

    if(RunSyst){
      for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){
        //==== Everything else remains same, but only change syst_ and parameter name
        param.syst_ = AnalyzerParameter::Syst(it_syst);
        param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
        executeEventFromParameter(param);
      }
    }

  }

}

void CHMuJets::executeEventFromParameter(AnalyzerParameter param){

  //=============
  //==== No Cut
  //=============

  JSFillHist("CutFlow", "NoCut", 0., 1., 1, 0., 1.);

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  JSFillHist("CutFlow", "Event Cut", 0., 1., 1, 0., 1.);

  //==== Get Event object
  Event ev = GetEvent();
  int nVer = ev.nPV();
  Particle METv = ev.GetMETVector();

  //==============
  //==== Trigger
  //==============
  if(! (ev.PassTrigger(IsoMuTriggerName) )) return;

  JSFillHist("CutFlow", "Trigger Cut", 0., 1., 1, 0., 1.);
  //======================
  //==== Copy AllObjects
  //======================

  vector<Muon> this_AllMuons = AllMuons;
  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Jet> this_AllJets = AllJets;

  //==== Then, for each systematic sources
  //==== 1) Smear or scale them
  //==== 2) Then apply ID selections
  //==== This order should be explicitly followed
  //==== Below are all variables for available systematic sources

  if(param.syst_ == AnalyzerParameter::Central){

  }
  else if(param.syst_ == AnalyzerParameter::JetResUp){
    this_AllJets = SmearJets( this_AllJets, +1 );
    //this_AllFatJets = SmearFatJets( this_AllFatJets, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetResDown){
    this_AllJets = SmearJets( this_AllJets, -1 );
    //this_AllFatJets = SmearFatJets( this_AllFatJets, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetEnUp){
    this_AllJets = ScaleJets( this_AllJets, +1 );
    //this_AllFatJets = ScaleFatJets( this_AllFatJets, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::JetEnDown){
    this_AllJets = ScaleJets( this_AllJets, -1 );
    //this_AllFatJets = ScaleFatJets( this_AllFatJets, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::MuonEnUp){
    this_AllMuons = ScaleMuons( this_AllMuons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::MuonEnDown){
    this_AllMuons = ScaleMuons( this_AllMuons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronResUp){
    //this_AllElectrons = SmearElectrons( this_AllElectrons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronResDown){
    //this_AllElectrons = SmearElectrons( this_AllElectrons, -1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronEnUp){
    //this_AllElectrons = ScaleElectrons( this_AllElectrons, +1 );
  }
  else if(param.syst_ == AnalyzerParameter::ElectronEnDown){
    //this_AllElectrons = ScaleElectrons( this_AllElectrons, -1 );
  }
  else{
    cout << "[CHMuJets::executeEventFromParameter] Wrong syst" << endl;
    exit(EXIT_FAILURE);
  }

  //==================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  vector<Muon> muons = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 26., 2.4);
  vector<Electron> electrons_tight = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 30.,2.5);
  vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 15.,2.5);
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);

  //=======================
  //==== Sort in pt-order
  //=======================

  //==== 1) leptons : after scaling/smearing, pt ordring can differ from MINIAOD
  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(electrons_tight.begin(), electrons_tight.end(), PtComparing);
  std::sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
  //==== 2) jets : similar, but also when applying new JEC, ordering is changes. This is important if you use leading jets
  std::sort(jets.begin(), jets.end(), PtComparing);

  //=========================
  //==== Event selections..
  //=========================

  //==== single muon
  if(muons.size() != 1) return;
  //==== leading muon has trigger-safe pt
  if( muons.at(0).Pt() <= TriggerSafePtCut ) return;
  JSFillHist("CutFlow", "Single Muon", 0., 1., 1, 0., 1.);
  //==== extra lepton veto
  if(electrons_veto.size() != 0) return;
  JSFillHist("CutFlow", "Extra Lepton Veto", 0., 1., 1, 0., 1.);
  //==== MET cut
  if(METv.E()<20) return;
  JSFillHist("CutFlow", "MET > 20", 0., 1., 1, 0., 1.);
  //==== 4jet cut
  int njets = jets.size();
  if(njets<4) return;
  JSFillHist("CutFlow", "more then 4jets", 0., 1., 1, 0., 1.);


  //==== btag cut
  Int_t nbtags=GetNBTags(jets,Jet::DeepCSV, Jet::Medium);

  if(nbtags<2) return;
  JSFillHist("CutFlow", "more then 2btaged", 0., 1., 1, 0., 1.);
  if(nbtags>=3){
    JSFillHist("CutFlow", "more then 3btaged", 0., 1., 1, 0., 1.);
  }

  //===================
  //==== Event weight
  //===================

  double weight = 1.;

  JSFillHist("no_weight", "leading_jet_pT", jets.at(0).Pt(), weight,60,0.,300.);
  JSFillHist("no_weight", "leading_jet_eta", jets.at(0).Eta(), weight,60,-2.4,2.4);
  JSFillHist("no_weight", "muon_pT", muons.at(0).Pt(), weight,60,0.,300.);
  JSFillHist("no_weight", "muon_eta", muons.at(0).Eta(), weight,60,-2.4,2.4);
  JSFillHist("no_weight", "nbtags", (double)nbtags, weight,3,2.,5.);
  JSFillHist("no_weight", "njets", (double)njets, weight,10,4.,14.);
  JSFillHist("no_weight", "nVertices", (double)nVer, weight, 50,0.,50.);

  //==== If MC
  if(!IsDATA){

    //==== weight_norm_1invpb is set to be event weight normalized to 1 pb-1
    //==== So, you have to multiply trigger luminosity
    //==== you can pass trigger names to ev.GetTriggerLumi(), but if you are using unprescaled trigger, simply pass "Full"
    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");

    //==== MCweight is +1 or -1. Should be multiplied if you are using e.g., aMC@NLO NLO samples
    weight *= ev.MCweight();

    JSFillHist("normalization", "w1_leading_jet_pT", jets.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("normalization", "w1_leading_jet_eta", jets.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("normalization", "w1_muon_pT", muons.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("normalization", "w1_muon_eta", muons.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("normalization", "w1_nbtags", (double)nbtags, weight,3,2.,5.);
    JSFillHist("normalization", "w1_njets", (double)njets, weight,10,4.,14.);
    JSFillHist("normalization", "w1_nVertices", (double)nVer, weight, 50,0.,50.);

    //==== L1Prefire reweight
    weight *= weight_Prefire;

    JSFillHist("prefire", "w2_leading_jet_pT", jets.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("prefire", "w2_leading_jet_eta", jets.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("prefire", "w2_muon_pT", muons.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("prefire", "w2_muon_eta", muons.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("prefire", "w2_nbtags", (double)nbtags, weight,3,2.,5.);
    JSFillHist("prefire", "w2_njets", (double)njets, weight,10,4.,14.);
    JSFillHist("prefire", "w2_nVertices", (double)nVer, weight, 50,0.,50.);


    //==== trigger SF
    double this_trigsf = mcCorr->MuonTrigger_SF(param.Muon_ID_SF_Key, IsoMuTriggerName, muons);
    weight *= this_trigsf;

    JSFillHist("trigger_SF", "w3_leading_jet_pT", jets.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("trigger_SF", "w3_leading_jet_eta", jets.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("trigger_SF", "w3_muon_pT", muons.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("trigger_SF", "w3_muon_eta", muons.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("trigger_SF", "w3_nbtags", (double)nbtags, weight,3,2.,5.);
    JSFillHist("trigger_SF", "w3_njets", (double)njets, weight,10,4.,14.);
    JSFillHist("trigger_SF", "w3_nVertices", (double)nVer, weight, 50,0.,50.);


    //==== Example of applying Muon scale factors
    for(unsigned int i=0; i<muons.size(); i++){

      double this_idsf  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt());

      //==== If you have iso SF, do below. Here we don't.
      //double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt());
      double this_isosf = 1.;

      weight *= this_idsf*this_isosf;

    }

    JSFillHist("muon_SF", "w4_leading_jet_pT", jets.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("muon_SF", "w4_leading_jet_eta", jets.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("muon_SF", "w4_muon_pT", muons.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("muon_SF", "w4_muon_eta", muons.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("muon_SF", "w4_nbtags", (double)nbtags, weight,3,2.,5.);
    JSFillHist("muon_SF", "w4_njets", (double)njets, weight,10,4.,14.);
    JSFillHist("muon_SF", "w4_nVertices", (double)nVer, weight, 50,0.,50.);


    //==== PileUpWeight
    double pileup_weight = GetPileUpWeight(nVer, 0);
    weight *= pileup_weight;

    JSFillHist("pileup_reweight", "w5_leading_jet_pT", jets.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("pileup_reweight", "w5_leading_jet_eta", jets.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("pileup_reweight", "w5_muon_pT", muons.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("pileup_reweight", "w5_muon_eta", muons.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("pileup_reweight", "w5_nbtags", (double)nbtags, weight,3,2.,5.);
    JSFillHist("pileup_reweight", "w5_njets", (double)njets, weight,10,4.,14.);
    JSFillHist("pileup_reweight", "w5_nVertices", (double)nVer, weight, 50,0.,50.);

    //==== top-pt-reweight
    double top_pt_reweight = mcCorr->GetTopPtReweight(AllGens);
    weight *= top_pt_reweight;

    JSFillHist("top_pt_reweight", "w6_leading_jet_pT", jets.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("top_pt_reweight", "w6_leading_jet_eta", jets.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("top_pt_reweight", "w6_muon_pT", muons.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("top_pt_reweight", "w6_muon_eta", muons.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("top_pt_reweight", "w6_nbtags", (double)nbtags, weight,3,2.,5.);
    JSFillHist("top_pt_reweight", "w6_njets", (double)njets, weight,10,4.,14.);
    JSFillHist("top_pt_reweight", "w6_nVertices", (double)nVer, weight, 50,0.,50.);

    //==== btag SF
    float btag_sf=1., mistag_sf=1.;
    BtaggingSFEvtbyEvt(jets, Jet::DeepCSV, Jet::Medium, 0, btag_sf, mistag_sf);
    weight *= btag_sf*mistag_sf;

    JSFillHist("btag_sf", "w7_leading_jet_pT", jets.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("btag_sf", "w7_leading_jet_eta", jets.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("btag_sf", "w7_muon_pT", muons.at(0).Pt(), weight,60,0.,300.);
    JSFillHist("btag_sf", "w7_muon_eta", muons.at(0).Eta(), weight,60,-2.4,2.4);
    JSFillHist("btag_sf", "w7_nbtags", (double)nbtags, weight,3,2.,5.);
    JSFillHist("btag_sf", "w7_njets", (double)njets, weight,10,4.,14.);
    JSFillHist("btag_sf", "w7_nVertices", (double)nVer, weight, 50,0.,50.);

  }
}

