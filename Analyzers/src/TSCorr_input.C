#include "TSCorr_input.h"

TSCorr_input::TSCorr_input(){

}

void TSCorr_input::initializeAnalyzer(){

  //==== if you've done "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
  RunSyst = HasFlag("RunSyst");

  cout << "[TSCorr_input::initializeAnalyzer] RunSyst = " << RunSyst << endl;

  //==== In the following examples, I will obtain dimuon Z-peak with two muon IDs
  //==== I defined "vector<TString> MuonIDs;" in Analyzers/include/TSCorr_input.h
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
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/TSCorr_input.h 
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

  cout << "[TSCorr_input::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  cout << "[TSCorr_input::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== add taggers and WP that you want to use in analysis                                                                                                                                                
  std::vector<Jet::Tagger> vtaggers;
  vtaggers.push_back(Jet::DeepCSV);

  std::vector<Jet::WP> v_wps;
  v_wps.push_back(Jet::Medium);

  //=== list of taggers, WP, setup systematics, use period SFs                                                                                                                                              
  SetupBTagger(vtaggers,v_wps, true, true);


}

TSCorr_input::~TSCorr_input(){

  //==== Destructor of this Analyzer

}

void TSCorr_input::executeEvent(){

  //================================================================
  //====  (Example)
  //====  Dimuon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/TSCorr_input.h,
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
  //==== I defined "double weight_Prefire;" in Analyzers/include/TSCorr_input.h
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

void TSCorr_input::executeEventFromParameter(AnalyzerParameter param){

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
  //int nVer = ev.nPV();
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
    cout << "[TSCorr_input::executeEventFromParameter] Wrong syst" << endl;
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


  //===================
  //==== parton-jet matching
  //===================

  parton_jet_matching_CHToCB->SetGens(AllGens);
  parton_jet_matching_CHToCB->SetJets(jets);
  if(!parton_jet_matching_CHToCB->FindHardProcessParton()){
    return;
  }
  if(!parton_jet_matching_CHToCB->MatchJets()){
    return;
  }

  //=== get TLorentzVectors
  TLorentzVector hadronic_top_b, leptonic_top_b, up_type_jet, down_type_jet, neutrino;
  hadronic_top_b = parton_jet_matching_CHToCB->Get_hadronic_top_b_jet();
  leptonic_top_b = parton_jet_matching_CHToCB->Get_leptonic_top_b_jet();
  up_type_jet = parton_jet_matching_CHToCB->Get_up_type_jet();
  down_type_jet = parton_jet_matching_CHToCB->Get_down_type_jet();
  neutrino = parton_jet_matching_CHToCB->Get_neutrino();

  TLorentzVector Mu_vec;
  Muon Mu = muons.at(0);
  Mu_vec.SetPtEtaPhiE(Mu.Pt(),Mu.Eta(),Mu.Phi(),Mu.E());

  // Pz of neutrino.
  double PzSol[2] = {0,0};
  double d = 0.5*(80.385*80.385-Mu_vec.M()*Mu_vec.M())+Mu_vec.Px()*METv.Px()+METv.Py()*METv.Py();
  double a = Mu_vec.E()*Mu_vec.E() - Mu_vec.Pz()*Mu_vec.Pz();
  //Double_t b = 2*d*Mu_vec.Pz();
  double b = 2*d*Mu_vec.Pz();
  double c = Mu_vec.E()*Mu_vec.E()*METv.Pt()*METv.Pt()-d*d;
  //
  if (b*b-4*a*c>=0){
    PzSol[0] = (b+sqrt(b*b-4*a*c))/(2*a);
    PzSol[1] = (b-sqrt(b*b-4*a*c))/(2*a);
  }
  else { //	taking real part for complex solution.
    return;
    PzSol[0] = (b)/(2*a);
    PzSol[1] = (b)/(2*a);
  }
  JSFillHist("CutFlow", "real Neu Pz", 0., 1., 1, 0., 1.);
  TLorentzVector Neu_vec;

  double neupz =0;
  //double neupz_false=0;
  double truth_neupz = neutrino.Pz();
  double d_pz[2];
  d_pz[0] = abs(PzSol[0]-truth_neupz);
  d_pz[1] = abs(PzSol[1]-truth_neupz);
	
  if (d_pz[0] < d_pz[1]){
	neupz = PzSol[0];
	//neupz_false = PzSol[1];
  }
  else{
	neupz = PzSol[1];
	//neupz_false = PzSol[0];
  }

  Neu_vec.SetPxPyPzE(METv.Px(),METv.Py(),neupz,sqrt(METv.E()*METv.E()+neupz*neupz));
  TLorentzVector leptonic_w = Mu_vec+ Neu_vec;

  TString tree_name = "tree_var";
  //=== if key not exist, initialize tree
  if(mapTree.find(tree_name) == mapTree.end()){
    mapTree[tree_name] = new TTree(tree_name,tree_name);

    outfile->SetCompressionLevel(1);
    mapTree[tree_name]->Branch("hadronic_top_b_jet","TLorentzVector", &hadronic_top_b);
    mapTree[tree_name]->Branch("leptonic_top_b_jet","TLorentzVector", &leptonic_top_b);
    mapTree[tree_name]->Branch("up_type_jet","TLorentzVector", &up_type_jet);
    mapTree[tree_name]->Branch("down_type_jet","TLorentzVector", &down_type_jet);
    mapTree[tree_name]->Branch("muon","TLorentzVector", &Mu_vec);
    mapTree[tree_name]->Branch("neutrino","TLorentzVector", &Neu_vec);
  }
  else{
    mapTree[tree_name]->Fill();
  }

}

