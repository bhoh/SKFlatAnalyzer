#include "FitterInput.h"

FitterInput::FitterInput(){

}

void FitterInput::initializeAnalyzer(){

  //==== if you've done "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");

  cout << "[FitterInput::initializeAnalyzer] RunSyst = " << RunSyst << endl;

  //==== In the following examples, I will obtain dimuon Z-peak with two muon IDs
  //==== I defined "vector<TString> MuonIDs;" in Analyzers/include/FitterInput.h
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
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/FitterInput.h 
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

  cout << "[FitterInput::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  cout << "[FitterInput::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== add taggers and WP that you want to use in analysis                                                                                                                                                
  std::vector<Jet::Tagger> vtaggers;
  vtaggers.push_back(Jet::DeepCSV);

  std::vector<Jet::WP> v_wps;
  v_wps.push_back(Jet::Medium);

  //=== list of taggers, WP, setup systematics, use period SFs                                                                                                                                              
  SetupBTagger(vtaggers,v_wps, true, true);


}

FitterInput::~FitterInput(){

  //==== Destructor of this Analyzer

}

void FitterInput::executeEvent(){

  //================================================================
  //====  (Example)
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/FitterInput.h,
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
  //==== I defined "double weight_Prefire;" in Analyzers/include/FitterInput.h
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

void FitterInput::executeEventFromParameter(AnalyzerParameter param){

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
    cout << "[FitterInput::executeEventFromParameter] Wrong syst" << endl;
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
  if(METv.E()<20) return; //FIXME: anti-isolation method
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
  //==== get btag vector
  std::vector<bool> btag_vector = GetBTagVector(jets,Jet::DeepCSV, Jet::Medium);
  //===================
  //==== Event weight
  //===================

  double weight = 1.;
  double weight_norm_lumi=-1.;
  double trigger_sf=-1.;
  double lepton_sf=-1.;
  double pileup_reweight=-1.;
  double toppt_reweight=-1.;
  double btag_sf=-1.;
  double mistag_sf=-1.;

  //==== If MC
  if(!IsDATA){

    //==== weight_norm_1invpb is set to be event weight normalized to 1 pb-1
    //==== So, you have to multiply trigger luminosity
    //==== you can pass trigger names to ev.GetTriggerLumi(), but if you are using unprescaled trigger, simply pass "Full"
    weight_norm_lumi = weight_norm_1invpb*ev.GetTriggerLumi("Full");
    //==== MCweight is +1 or -1. Should be multiplied if you are using e.g., aMC@NLO NLO samples
    weight_norm_lumi *= ev.MCweight();
    weight *= weight_norm_lumi;
  
    //==== L1Prefire reweight
    weight *= weight_Prefire;


    //==== trigger SF
    trigger_sf = mcCorr->MuonTrigger_SF(param.Muon_ID_SF_Key, IsoMuTriggerName, muons);
    weight *= trigger_sf;

    //==== Example of applying Muon scale factors
    for(unsigned int i=0; i<muons.size(); i++){

      double this_idsf  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt());

      //==== If you have iso SF, do below. Here we don't.
      //double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt());
      double this_isosf = 1.;

      lepton_sf = this_idsf*this_isosf;

    }
    weight *= lepton_sf;

    //==== PileUpWeight
    pileup_reweight = GetPileUpWeight(nVer, 0);
    weight *= pileup_reweight;

    //==== top-pt-reweight
    toppt_reweight = mcCorr->GetTopPtReweight(AllGens);
    weight *= toppt_reweight;

    //==== btag SF
    float btag_sf_=1., mistag_sf_=1.;
    BtaggingSFEvtbyEvt(jets, Jet::DeepCSV, Jet::Medium, 0, btag_sf_, mistag_sf_);
    btag_sf = btag_sf_;
    mistag_sf = mistag_sf_;
    weight *= btag_sf*mistag_sf;
  }

  //=== make jet vector, mu/el and MET vector in TLorentzVector
  //=== jet vector
  std::vector<TLorentzVector> jet_vector;
  std::copy(jets.begin(),jets.end(),jet_vector.begin());
  //=== muon vector
  TLorentzVector muon_vector = muons.at(0);
  //=== met vector
  TLorentzVector met_vector = METv;
  //=== make tree and fill tree
  TString tree_name = "tree_var";
  //=== if key not exist, initialize tree
  if(mapTree.find(tree_name) == mapTree.end()){
    mapTree[tree_name] = new TTree(tree_name,tree_name);

    outfile->SetCompressionLevel(1);
    mapTree[tree_name]->Branch("jet_vector","std::vector<TLorentzVector>", &jet_vector);
    mapTree[tree_name]->Branch("btag_vector","std::vector<bool>", &btag_vector);
    mapTree[tree_name]->Branch("muon_vector","TLorentzVector", &muon_vector);
    mapTree[tree_name]->Branch("met_vector","TLorentzVector", &met_vector);
    mapTree[tree_name]->Branch("weight_total", &weight);
    mapTree[tree_name]->Branch("weight_norm_lumi", &weight_norm_lumi);
    mapTree[tree_name]->Branch("weight_Prefire", &weight_Prefire);
    mapTree[tree_name]->Branch("trigger_sf", &trigger_sf);
    mapTree[tree_name]->Branch("lepton_sf", &lepton_sf);
    mapTree[tree_name]->Branch("pileup_reweight", &pileup_reweight);
    mapTree[tree_name]->Branch("toppt_reweight", &toppt_reweight);
    mapTree[tree_name]->Branch("btag_sf", &btag_sf);
    mapTree[tree_name]->Branch("mistag_sf", &mistag_sf);
  }
  else{
    mapTree[tree_name]->Fill();
  }

}

