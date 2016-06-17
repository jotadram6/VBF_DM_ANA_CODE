#include "Analyzer.h"
//#include "Histo.h"
#define ival(x) static_cast<int>(x)
//#define at(x) operator[](x)
typedef vector<int>::iterator vec_iter;
const string FILESPACE = "PartDet/";



void Analyzer::clear_values() {
  for(int i=0; i < (int)goodParts.size(); i++) {
    goodParts[i].clear();
  }
  deltaMEx=0;
  deltaMEy=0;
  sumpxForMht = 0.0;
  sumpyForMht = 0.0;
  sumptForHt  = 0.0;
  leadIndex=-1;
}


void Analyzer::preprocess(int event) {
  BOOM->GetEntry(event);
  
  theMETVector.SetPxPyPzE(Met_px, Met_py, Met_pz, sqrt(pow(Met_px,2) + pow(Met_py,2)));
  pu_weight = (!isData && CalculatePUSystematics) ? getPileupWeight(nTruePU) : 1.0;
  // SET NUMBER OF GEN PARTICLES
  // TODOGeneralize to remove magic numbers

  if(!isData){
    PartStats genStat = _Gen->pstats["Gen"];
    ExtractNumberOfGoodGen(15, 2, CUTS::eGTau, &genStat);
    ExtractNumberOfGoodGen(6, 2, CUTS::eGTop, &genStat);
    ExtractNumberOfGoodGen(11, 1, CUTS::eGElec, &genStat);
    ExtractNumberOfGoodGen(13, 1, CUTS::eGMuon, &genStat);
    ExtractNumberOfGoodGen(23, 2, CUTS::eGZ, &genStat);
    ExtractNumberOfGoodGen(24, 2, CUTS::eGW, &genStat);
    ExtractNumberOfGoodGen(25, 2, CUTS::eGHiggs, &genStat);
    ExtractNumberOfTauNu();
  }

  //////Smearing  
  SmearLepton(_Electron, CUTS::eGElec, &_Electron->pstats["Smear"]);
  SmearLepton(_Muon, CUTS::eGMuon, &_Muon->pstats["Smear"]);
  SmearLepton(_Tau, CUTS::eGTau, &_Tau->pstats["Smear"]);
  SmearJet();


  //////Triggers and Vertices
  goodParts[ival(CUTS::eRVertex)].resize(bestVertices);
  if(passTriggerCuts("Trigger1")) goodParts[ival(CUTS::eRTrig1)].resize(1);
  if(passTriggerCuts("Trigger2")) goodParts[ival(CUTS::eRTrig2)].resize(1);

  // // SET NUMBER OF RECO PARTICLES
  // // MUST BE IN ORDER: Muon/Electron, Tau, Jet
  ExtractNumberOfGoodReco(_Electron, CUTS::eRElec1, CUTS::eGElec, &_Electron->pstats["Elec1"]);
  ExtractNumberOfGoodReco(_Electron, CUTS::eRElec2, CUTS::eGElec, &_Electron->pstats["Elec2"]);
  ExtractNumberOfGoodReco(_Muon, CUTS::eRMuon1, CUTS::eGMuon, &_Muon->pstats["Muon1"]);
  ExtractNumberOfGoodReco(_Muon, CUTS::eRMuon2, CUTS::eGMuon, &_Muon->pstats["Muon2"]);
  ExtractNumberOfGoodReco(_Tau, CUTS::eRTau1, CUTS::eGTau, &_Tau->pstats["Tau1"]);
  ExtractNumberOfGoodReco(_Tau, CUTS::eRTau2, CUTS::eGTau, &_Tau->pstats["Tau2"]);

  ExtractNumberOfGoodRecoJets(CUTS::eRJet1, &_Jet->pstats["Jet1"]);
  ExtractNumberOfGoodRecoJets(CUTS::eRJet2, &_Jet->pstats["Jet2"]);
  ExtractNumberOfGoodRecoJets(CUTS::eRCenJet, &_Jet->pstats["CentralJet"]);
  ExtractNumberOfGoodRecoJets(CUTS::eRBJet, &_Jet->pstats["BJet"]);

  ExtractNumberOfGoodRecoJets(CUTS::eR1stJet, &_Jet->pstats["FirstLeadingJet"]);
  leadIndex = (goodParts[ival(CUTS::eR1stJet)].size() != 0) ? goodParts[ival(CUTS::eR1stJet)].at(0) : -1;
  ExtractNumberOfGoodRecoJets(CUTS::eR2ndJet, &_Jet->pstats["SecondLeadingJet"]);

  updateMet();

  /////  SET NUMBER OF RECO MET TOPOLOGY PARTICLES
  passRecoLeptonMetTopologyCuts(_Electron, CUTS::eRElec1, CUTS::eTElec1, &_Electron->pstats["Elec1"]);
  passRecoLeptonMetTopologyCuts(_Electron, CUTS::eRElec2, CUTS::eTElec2, &_Electron->pstats["Elec2"]);
  passRecoLeptonMetTopologyCuts(_Muon, CUTS::eRMuon1, CUTS::eTMuon1, &_Muon->pstats["Muon1"]);
  passRecoLeptonMetTopologyCuts(_Muon, CUTS::eRMuon2, CUTS::eTMuon2, &_Muon->pstats["Muon2"]);
  passRecoLeptonMetTopologyCuts(_Tau, CUTS::eRTau1, CUTS::eTTau1, &_Tau->pstats["Tau1"]);
  passRecoLeptonMetTopologyCuts(_Tau, CUTS::eRTau2, CUTS::eTTau2, &_Tau->pstats["Tau2"]);

  /////fix up susy cut
  /////need leeding jets and then maybe change from pass to extract
  if(goodParts[ival(CUTS::eR1stJet)].size() != 0 && goodParts[ival(CUTS::eR2ndJet)].size() != 0) SusyTopologyCuts();

  /////lepton lepton topology cuts
  passLeptonComboTopologyCut(_Tau, _Electron, CUTS::eRTau1, CUTS::eRElec1, CUTS::eElec1Tau1, &distats["Electron1Tau1"]);
  passLeptonComboTopologyCut(_Tau, _Electron, CUTS::eRTau1, CUTS::eRElec2, CUTS::eElec2Tau1, &distats["Electron2Tau1"]);
  passLeptonComboTopologyCut(_Tau, _Muon, CUTS::eRTau1, CUTS::eRMuon1, CUTS::eMuon1Tau1, &distats["Muon1Tau1"]);
  passLeptonComboTopologyCut(_Tau, _Muon, CUTS::eRTau1, CUTS::eRMuon2, CUTS::eMuon2Tau1, &distats["Muon2Tau1"]);
  passLeptonComboTopologyCut(_Tau, _Electron, CUTS::eRTau2, CUTS::eRElec1, CUTS::eElec1Tau2, &distats["Electron1Tau2"]);
  passLeptonComboTopologyCut(_Tau, _Electron, CUTS::eRTau2, CUTS::eRElec2, CUTS::eElec2Tau2, &distats["Electron2Tau2"]);
  passLeptonComboTopologyCut(_Tau, _Muon, CUTS::eRTau2, CUTS::eRMuon1, CUTS::eMuon1Tau2, &distats["Muon1Tau2"]);
  passLeptonComboTopologyCut(_Tau, _Muon, CUTS::eRTau2, CUTS::eRMuon2, CUTS::eMuon2Tau2, &distats["Muon2Tau2"]);

  ////DIlepton topology cuts
  passLeptonComboTopologyCut(_Tau, _Tau, CUTS::eRTau1, CUTS::eRTau2, CUTS::eDiTau, &distats["DiTau"]);
  passLeptonComboTopologyCut(_Electron, _Electron, CUTS::eRElec1, CUTS::eRElec2, CUTS::eDiElec, &distats["DiElectron"]);
  passLeptonComboTopologyCut(_Muon, _Muon, CUTS::eRMuon1, CUTS::eRMuon2, CUTS::eDiMuon, &distats["DiMuon"]);

  ///////////////fix up dijet stuff
  passDiJetTopologyCuts(&distats["DiJet"]);

  // if(event % 50000 == 0) {
  //   cout << "Event #" << event << endl;
  // }
}

int Analyzer::fillCuts() {
  unordered_map<string,pair<int,int> >* cut_info = histo->get_cuts();
  vector<string>* cut_order = histo->get_order();

  string cut;
  int min, max;
  bool prevTrue = true;
  int nparticles, i=0;
  int maxCut=0;

  for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
    if(isData && it->find("Gen") != string::npos) continue;
    cut = *it;
    min= cut_info->at(cut).first;
    max= cut_info->at(cut).second;
    nparticles = goodParts[ival(cut_num[cut])].size();
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      cuts_per[i]++;
      cuts_cumul[i] += (prevTrue) ? 1 : 0;
      maxCut += (prevTrue) ? 1 : 0;
    } else prevTrue = false;
  }
  return maxCut;
}

void Analyzer::printCuts() {
  vector<string>* cut_order = histo->get_order();
  int i =0;

  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(3);
  cout << "\n";
  cout << "Selection Efficiency " << "\n";
  cout << "Total events: " << nentries << "\n";
  cout << "         Name                     Indiv.      Cumulative\n";
  cout << "---------------------------------------------------------------------------\n";
  for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
    cout << setw(28) << *it << " ";
    if(isData && it->find("Gen") != string::npos) cout << "Skipped" << endl;
    else cout << setw(8) << cuts_per.at(i) << " (" << setw(8) << ((float)cuts_per.at(i)) / nentries << ") "
	      << setw(8) << cuts_cumul.at(i) << "( " << setw(8) << ((float)cuts_cumul.at(i)) / nentries << ") " << endl;
  }
  cout << "---------------------------------------------------------------------------\n";  
}




void Analyzer::updateMet() {
  ////// Neutrino update before calculation
  if(distats["Run"].bmap["TreatMuonsAsNeutrinos"]) {
    for(vec_iter it=goodParts[ival(CUTS::eRMuon1)].begin(); it!=goodParts[ival(CUTS::eRMuon1)].end(); it++) {
      if(find(goodParts[ival(CUTS::eRMuon2)].begin(), goodParts[ival(CUTS::eRMuon2)].end(), (*it)) != goodParts[ival(CUTS::eRMuon2)].end() ) continue;
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }    
    for(vec_iter it=goodParts[ival(CUTS::eRMuon2)].begin(); it!=goodParts[ival(CUTS::eRMuon2)].end(); it++) {
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }
  }
  ///---MHT and HT calculations----////
  int i=0;
  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it!=_Jet->smearP.end(); it++, i++) {
    if( (it->Pt() > distats["Run"].dmap["JetPtForMhtAndHt"]) && (fabs(it->Eta()) < distats["Run"].dmap["JetEtaForMhtAndHt"]) ) {
      if(distats["Run"].bmap["ApplyJetLooseIDforMhtAndHt"] ) continue;//&& !passLooseJetID() continue;
      
      sumpxForMht -= it->Px();
      sumpyForMht -= it->Py();
      sumptForHt  += it->Pt();
    }
  }
  phiForMht = atan(sumpyForMht/sumpxForMht);
  if(sumpxForMht < 0) phiForMht += (sumpyForMht >= 0) ? TMath::Pi() : -TMath::Pi();

  theMETVector.SetPxPyPzE(theMETVector.Px()+deltaMEx, theMETVector.Py()+deltaMEy, theMETVector.Pz(), 
  			  TMath::Sqrt(pow(theMETVector.Px()+deltaMEx,2) + pow(theMETVector.Py()+deltaMEy,2)));
}


Analyzer::Analyzer(string infile, string outfile) {

  setupJob(infile);
  
  setupGeneral(BOOM,infile);

  isData = distats["Run"].bmap["isData"];
  CalculatePUSystematics = distats["Run"].bmap["CalculatePUSystematics"];
  prevTrig["Trigger1"] = make_pair(0,0);
  prevTrig["Trigger2"] = make_pair(0,0);

  initializePileupInfo(distats["Run"].smap["MCHistos"], distats["Run"].smap["DataHistos"]);

  //////need to initialize histo and get values for cut arrays

  histo = new Histogramer(1, FILESPACE+"Hist_entries.in", FILESPACE+"Cuts.in", outfile);
  cuts_per.resize(histo->get_cuts()->size());
  cuts_cumul.resize(histo->get_cuts()->size());

  _Gen = new Generated(BOOM, FILESPACE + "Gen_info.in");
  _Electron = new Electron(BOOM, FILESPACE + "Electron_info.in");
  _Muon = new Muon(BOOM, FILESPACE + "Muon_info.in");
  _Tau = new Taus(BOOM, FILESPACE + "Tau_info.in");
  _Jet = new Jet(BOOM, FILESPACE + "Jet_info.in");
  
  cout << "finished setup" << endl;

}

void Analyzer::writeout() {
  delete histo;
}

Analyzer::~Analyzer() {
  delete _Gen;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
}

void Analyzer::setupJob(string filename) {
  f = TFile::Open(filename.c_str());
  f->cd("TNT");
  BOOM = (TTree*)f->Get("TNT/BOOM");
  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "setup" << std::endl;
  std::cout << nentries << std::endl;
}

void Analyzer::setupGeneral(TTree* BOOM, string infile) {
  BOOM->SetBranchStatus("Trigger_decision", 1);
  BOOM->SetBranchStatus("Trigger_names", 1);
  BOOM->SetBranchStatus("nTruePUInteractions", 1);
  BOOM->SetBranchStatus("bestVertices", 1);
  BOOM->SetBranchStatus("Met_type1PF_px", 1);
  BOOM->SetBranchStatus("Met_type1PF_py", 1);
  BOOM->SetBranchStatus("Met_type1PF_pz", 1);

  BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision);
  BOOM->SetBranchAddress("Trigger_names", &Trigger_names);
  BOOM->SetBranchAddress("nTruePUInteractions", &nTruePU);
  BOOM->SetBranchAddress("bestVertices", &bestVertices);
  BOOM->SetBranchAddress("Met_type1PF_px", &Met_px);
  BOOM->SetBranchAddress("Met_type1PF_py", &Met_py);
  BOOM->SetBranchAddress("Met_type1PF_pz", &Met_pz);

  read_info(FILESPACE + "ElectronTau_info.in");
  read_info(FILESPACE + "MuonTau_info.in");
  read_info(FILESPACE + "DiParticle_info.in");
  read_info(FILESPACE + "VBFCuts_info.in");
  read_info(FILESPACE + "Run_info.in");
}

void Analyzer::read_info(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }

  vector<string> stemp;
  string group, line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      cout << "error in " << filename << "; no groups specified for data" << endl;
      exit(1);
    } else if(stemp.size() == 2) {
      if(stemp[1].find(".") != string::npos && stemp[1].find("root") == string::npos) distats[group].dmap[stemp[0]]=stod(stemp[1]);
      else if(stemp[1] == "1" || stemp[1] == "true") distats[group].bmap[stemp[0]]=true;
      else if(stemp[1] == "0" || stemp[1] == "false") distats[group].bmap[stemp[0]]=false; 
      else distats[group].smap[stemp[0]] = stemp[1];

    } else  distats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
  }
  info_file.close();
}


void Analyzer::SmearLepton(Lepton* lepton, CUTS eGenPos, PartStats* stats) {
  lepton->smearP.clear();

  double smearedPt;
  double smearedEta;
  double smearedPhi;
  double smearedEnergy;


  for(int i = 0; i < (int)lepton->pt->size(); i++) {
    TLorentzVector tmpSmear;
    tmpSmear.SetPtEtaPhiE(lepton->pt->at(i), lepton->eta->at(i), lepton->phi->at(i), lepton->energy->at(i));

    if(isData || !stats->bmap["SmearTheParticle"]) {
      lepton->smearP.push_back(tmpSmear);
      continue;
    }

    TLorentzVector* genVec =  matchLeptonToGen(&tmpSmear, lepton,eGenPos);
    if(genVec == nullptr) {      
      lepton->smearP.push_back(tmpSmear);
      continue;
    }

    smearedPt = (genVec->Pt()*stats->dmap["PtScaleOffset"]) + (tmpSmear.Pt() - genVec->Pt())*stats->dmap["PtSigmaOffset"];
    smearedEta =(genVec->Eta()*stats->dmap["EtaScaleOffset"]) + (tmpSmear.Eta() - genVec->Eta())*stats->dmap["EtaSigmaOffset"];
    smearedPhi = (genVec->Phi() * stats->dmap["PhiScaleOffset"]) + (tmpSmear.Phi() - genVec->Phi())*stats->dmap["PhiSigmaOffset"];
    smearedEnergy = (genVec->Energy()*stats->dmap["EnergyScaleOffset"]) + (tmpSmear.Energy() - genVec->Energy())*stats->dmap["EnergySigmaOffset"];
    
    TLorentzVector final;
    final.SetPtEtaPhiE(smearedPt, smearedEta, smearedPhi, smearedEnergy);
    lepton->smearP.push_back(final);
    deltaMEx += tmpSmear.Px() - final.Px();
    deltaMEy += tmpSmear.Py() - final.Py();
    delete genVec;
  }
}


void Analyzer::SmearJet() {
  _Jet->smearP.clear();
  TLorentzVector jetV;
  PartStats* stats = &_Jet->pstats["Smear"];

  for(int i=0; i< (int)_Jet->pt->size(); i++) {
    jetV.SetPtEtaPhiE(_Jet->pt->at(i), _Jet->eta->at(i), _Jet->phi->at(i), _Jet->energy->at(i));

    if(isData || !stats->bmap["SmearTheParticle"]) {
      _Jet->smearP.push_back(jetV);
      continue;
    }
    
    if(JetMatchesLepton(_Muon, jetV, stats->dmap["MuonMatchingDeltaR"], CUTS::eGMuon) ||
       JetMatchesLepton(_Tau, jetV, stats->dmap["TauMatchingDeltaR"], CUTS::eGTau) ||       
       JetMatchesLepton(_Electron, jetV,stats->dmap["ElectronMatchingDeltaR"], CUTS::eGElec)){

      _Jet->smearP.push_back(jetV);
      continue;
    }
    
    _Jet->smearP.push_back(stats->dmap["JetEnergyScaleOffset"] * jetV);

  }
}


bool Analyzer::JetMatchesLepton(Lepton* lepton, TLorentzVector jetV, double partDeltaR, CUTS eGenPos) {
  TLorentzVector tempV;
  for(int j = 0; j < (int)lepton->pt->size(); j++) {
    tempV.SetPtEtaPhiE(lepton->pt->at(j), lepton->eta->at(j), lepton->phi->at(j), lepton->energy->at(j));
    if(jetV.DeltaR(tempV) < partDeltaR && matchLeptonToGen(&tempV, lepton, eGenPos) != nullptr) return true;
  }
  return false;
}

TLorentzVector* Analyzer::matchLeptonToGen(TLorentzVector* lvec, Lepton* lep, CUTS ePos) {
  PartStats* stats = &lep->pstats["Smear"];

  if(ePos == CUTS::eGTau) {
    return matchTauToGen(lvec, stats->dmap["matchingDeltaR"]);
  }
  TLorentzVector* genVec = new TLorentzVector(0,0,0,0);
  
  for(vec_iter it=goodParts[ival(ePos)].begin(); it !=goodParts[ival(ePos)].end();it++) {
    genVec->SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
    if(lvec->DeltaR(*genVec) <= stats->dmap["matchingDeltaR"]) {
      if(stats->bmap["UseMotherID"] && abs(_Gen->motherpdg_id->at(*it)) != stats->dmap["MotherID"]) continue; //MOTHER ID OPTION
      return genVec;
    }
  }
  delete genVec;
  return nullptr;
}


TLorentzVector* Analyzer::matchTauToGen(TLorentzVector* lvec, double lDeltaR) {
  TLorentzVector* genVec = new TLorentzVector(0,0,0,0);
  bool leptonicDecay = false;
  
  for(vec_iter it=goodParts[ival(CUTS::eGTau)].begin(); it !=goodParts[ival(CUTS::eGTau)].end();it++) {
    leptonicDecay = false;
    for(int j = 0; j < (int)_Gen->pt->size(); j++) {
      if( ((abs(_Gen->pdg_id->at(j)) == 12) || (abs(_Gen->pdg_id->at(j)) == 14)) && (_Gen->BmotherIndex->at(j) == (*it)) ) {
	leptonicDecay = true; 
	break;
      }
    }
    if(leptonicDecay) continue;

    for(vec_iter inu=goodParts[ival(CUTS::eNuTau)].begin(); inu !=goodParts[ival(CUTS::eNuTau)].end();inu++) {
      if(_Gen->BmotherIndex->at(*inu) != (*it)) continue;
      
      genVec->SetPtEtaPhiE(_Gen->pt->at(*it)-_Gen->pt->at(*inu), _Gen->eta->at(*it)-_Gen->eta->at(*inu), 
			   _Gen->phi->at(*it)-_Gen->phi->at(*inu), _Gen->energy->at(*it)-_Gen->energy->at(*inu));
      if(lvec->DeltaR(*genVec) <= lDeltaR) {
	return genVec;
      }
    }
  }
  delete genVec;
  return nullptr;
}


void Analyzer::ExtractNumberOfGoodGen(int particle_id, int particle_status, CUTS ePos, PartStats* stats) {
  for(int j = 0; j < (int)_Gen->pt->size(); j++) {
    if(particle_id == 15 && (_Gen->pt->at(j) < stats->dmap["TauPtMinCut"]) && (abs(_Gen->eta->at(j)) > stats->dmap["TauEtaMaxCut"])) continue;
    
    if((abs(_Gen->pdg_id->at(j)) == particle_id) && (_Gen->status->at(j) == particle_status)) {
      goodParts[ival(ePos)].push_back(j);
    }
  }
}

void Analyzer::ExtractNumberOfTauNu() {
  for(int j = 0; j < (int)_Gen->pt->size(); j++) {
    
    if( (abs(_Gen->pdg_id->at(j)) == 16) && (abs(_Gen->pdg_id->at(_Gen->BmotherIndex->at(j))) == 15) && (_Gen->status->at(_Gen->BmotherIndex->at(j)) == 2) ) {
      goodParts[ival(CUTS::eNuTau)].push_back(j);

    }
  }
}


void Analyzer::ExtractNumberOfGoodReco(Lepton* lep, CUTS ePos, CUTS eGenPos, PartStats* stats) {
  Electron* partE = dynamic_cast<Electron*>(lep);
  Muon* partM = dynamic_cast<Muon*>(lep);
  Taus* partT = dynamic_cast<Taus*>(lep);
  int i = 0;
  for(vector<TLorentzVector>::iterator it=lep->smearP.begin(); it != lep->smearP.end(); it++, i++) {
    TLorentzVector lvec = (*it);

    if (fabs(lvec.Eta()) > stats->dmap["EtaCut"]) continue;
    if (lvec.Pt() < stats->pmap["PtCut"].first || lvec.Pt() > stats->pmap["PtCut"].second) continue;

    if((stats->bmap["MatchToGen"]) && (!isData)) {   /////check
      if(matchLeptonToGen(&lvec, lep ,eGenPos) == nullptr) continue;
    }

    if(partM != nullptr) {      ////////////////MUON CUTS/////////////
      if(stats->bmap["DoDiscrByTightID"] && (partM->tight->at(i) == 0)) continue;
      if(stats->bmap["DoDiscrBySoftID"] && (partM->soft->at(i) == 0)) continue;


      if (stats->bmap["DoDiscrByIsolation"]) {
	double maxIsoval = max(0.0, partM->isoNeutralHadron->at(i) + partM->isoPhoton->at(i) - 0.5 * partM->isoPU->at(i));
	double isoSum = (partM->isoCharged->at(i) + maxIsoval) / lvec.Pt();
	if(isoSum < stats->pmap["IsoSumPtCutValue"].first || isoSum >= stats->pmap["IsoSumPtCutValue"].second) continue;
      }
    } else if(partE != nullptr) {    ///////////////ELECTRON CUT///////////
      if (stats->bmap["DoDiscrByIsolation"]) {
	double maxIsoval = std::max(0.0 , partE->isoNeutralHadrons->at(i) + partE->isoPhotons->at(i) - 0.5 * partE->isoPU->at(i) );
	double isoSum = (partE->isoChargedHadrons->at(i) + maxIsoval) / lvec.Pt();
	if(isoSum < stats->pmap["IsoSumPtCutValue"].first || isoSum >= stats->pmap["IsoSumPtCutValue"].second) continue;
      }

      //----Require electron to pass ID discriminators
      if(stats->bmap["DoDiscrByVetoID"] && (partE->isPassVeto->at(i) == 0)) continue;
      if(stats->bmap["DoDiscrByLooseID"] && (partE->isPassLoose->at(i) == 0)) continue;
      if(stats->bmap["DoDiscrByMediumID"] && (partE->isPassMedium->at(i) == 0)) continue;
      if(stats->bmap["DoDiscrByTightID"] && (partE->isPassTight->at(i) == 0)) continue;
      if(stats->bmap["DoDiscrByHEEPID"] && (partE->isPassHEEPId->at(i) == 0)) continue;

    } else if(partT != nullptr) {   /////////////TAU CUT/////////////////
      if (stats->bmap["DoDiscrByLeadTrack"]) {
	if(partT->leadChargedCandPt->at(i) < stats->dmap["LeadTrackThreshold"]) continue;
      }
      // ----Isolation requirement
      if (stats->bmap["DoDiscrByIsolation"]) {
	//--- max isolation requirement
	maxIso = (ePos == CUTS::eRTau1) ? partT->maxIso.first->at(i) : partT->maxIso.second->at(i);
	if(maxIso < 0.5) continue;
	if(stats->smap["DiscrByMinIsolation"] != "ZERO") {
	  minIso = (ePos == CUTS::eRTau1) ? partT->minIso.first->at(i) : partT->minIso.second->at(i);
	  if (minIso > 0.5) continue;
	}
      }

      // ----Require 1 or 3 prongs
      if(stats->smap["DiscrByProngType"].find("hps") != string::npos && partT->decayModeFindingNewDMs->at(i) < 0.5) continue;
      if(!passProng(stats->smap["DiscrByProngType"], partT->nProngs->at(i))) continue;

      // ----Electron and Muon vetos
      if (stats->bmap["DoDiscrAgainstElectron"] && !stats->bmap["SelectTausThatAreElectrons"]) {
	againstElectron = (ePos == CUTS::eRTau1) ? partT->againstElectron.first->at(i) : partT->againstElectron.second->at(i);
	if (againstElectron < 0.5) continue;
      } else if (stats->dmap["SelectTausThatAreElectrons"]) {
	if (againstElectron > 0.5) continue;
      }
      
      if (stats->bmap["DoDiscrAgainstMuon"] && !stats->bmap["SelectTausThatAreMuons"]) {
	againstMuon = (ePos == CUTS::eRTau1) ? partT->againstMuon.first->at(i) : partT->againstMuon.second->at(i);
	if (againstMuon < 0.5) continue;
      } else if (stats->bmap["SelectTausThatAreMuons"]) {
	againstMuon = (ePos == CUTS::eRTau1) ? partT->againstMuon.first->at(i) : partT->againstMuon.second->at(i);
	if (againstMuon > 0.5) continue;
      }

      if (stats->bmap["DoDiscrByCrackCut"] && isInTheCracks(lvec.Eta())) continue;

      // ----anti-overlap requirements
      if (stats->bmap["RemoveOverlapWithMuon1s"] && isOverlaping(lvec, _Muon, CUTS::eRMuon1, stats->dmap["Muon1MatchingDeltaR"])) continue;
      if (stats->bmap["RemoveOverlapWithMuon2s"] && isOverlaping(lvec, _Muon, CUTS::eRMuon2, stats->dmap["Muon2MatchingDeltaR"])) continue;
      if (stats->bmap["RemoveOverlapWithElectron1s"] && isOverlaping(lvec, _Electron, CUTS::eRElec1, stats->dmap["Electron1MatchingDeltaR"])) continue;
      if (stats->bmap["RemoveOverlapWithElectron2s"] && isOverlaping(lvec, _Electron, CUTS::eRElec2, stats->dmap["Electron2MatchingDeltaR"])) continue;
    }
    goodParts[ival(ePos)].push_back(i);    
  }
}


void Analyzer::ExtractNumberOfGoodRecoJets(CUTS ePos, PartStats* stats) {
  int i=0, prev=-1;
  double prevPt=0;
  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it != _Jet->smearP.end(); it++, i++) {
    TLorentzVector lvec = (*it);
    ///if else loop for central jet requirements
    if( ePos == CUTS::eRCenJet) {
      if(lvec.Pt() > 2.5) continue;
    } else if (fabs(lvec.Eta()) > stats->pmap["RecoEtaCut"].first || fabs(lvec.Eta()) > stats->pmap["RecoEtaCut"].second) continue;

    if (lvec.Pt() < stats->dmap["RecoPtCut"]) continue;

    if (stats->bmap["ApplyLeadingJetsLooseID"] && !passedLooseJetID(i)) continue;
    
  // ----anti-overlap requirements
    if(stats->bmap["RemoveOverlapWithMuon1s"] && isOverlaping(lvec, _Muon, CUTS::eRMuon1, stats->dmap["Muon1MatchingDeltaR"])) continue;
    if(stats->bmap["RemoveOverlapWithMuon2s"] && isOverlaping(lvec, _Muon, CUTS::eRMuon2, stats->dmap["Muon2MatchingDeltaR"])) continue;
    if(stats->bmap["RemoveOverlapWithElectron1s"] && isOverlaping(lvec, _Electron, CUTS::eRElec1, stats->dmap["Electron1MatchingDeltaR"])) continue;
    if(stats->bmap["RemoveOverlapWithElectron2s"] && isOverlaping(lvec, _Electron, CUTS::eRElec2, stats->dmap["Electron2MatchingDeltaR"])) continue;
    if(stats->bmap["RemoveOverlapWithTau1s"] && isOverlaping(lvec, _Tau, CUTS::eRTau1, stats->dmap["Tau1MatchingDeltaR"])) continue;
    if(stats->bmap["RemoveOverlapWithTau2s"] && isOverlaping(lvec, _Tau, CUTS::eRTau2, stats->dmap["Tau2MatchingDeltaR"])) continue;

    /// BJet specific
    if(stats->bmap["ApplyJetBTagging"] && _Jet->bDiscriminator->at(i) <= stats->dmap["JetBTaggingCut"]) continue;
    if((stats->bmap["MatchBToGen"]) && !isData && abs(_Jet->partonFlavour->at(i)) != 5) continue;

    /////fill up array
    if(ePos == CUTS::eR1stJet && it->Pt() > prevPt) prev = i;
    else if(ePos == CUTS::eR2ndJet && i != leadIndex && it->Pt() > prevPt) prev = i;
    else goodParts[ival(ePos)].push_back(i);    
  }
  if(prev != -1) goodParts[ival(ePos)].push_back(prev);    
} 


bool Analyzer::isOverlaping(TLorentzVector lvec, Lepton* overlapper, CUTS ePos, double MatchingDeltaR) {
  for(vec_iter it=goodParts[ival(ePos)].begin(); it < goodParts[ival(ePos)].end(); it++) {
    if(lvec.DeltaR(overlapper->smearP.at(*it)) < MatchingDeltaR) return true;
  }
  return false;
}

bool Analyzer::passProng(string prong, int value) {
  return ( (prong.find("1") != string::npos && value == 1) ||
	   (prong.find("2") != string::npos && value == 2) ||
	   (prong.find("3") != string::npos && value == 3) );
}

bool Analyzer::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 ||
          (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
          (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
          (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
          (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}
 
 
bool Analyzer::passedLooseJetID(int nobj) {
  if (_Jet->neutralHadEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->neutralEmEmEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->numberOfConstituents->at(nobj) <= 1) return false;
  if (_Jet->muonEnergyFraction->at(nobj) >= 0.80) return false;
  if ( (fabs(_Jet->smearP.at(nobj).Eta()) < 2.4) && 
           ((_Jet->chargedHadronEnergyFraction->at(nobj) <= 0.0) || 
	    (_Jet->chargedMultiplicity->at(nobj) <= 0.0) || 
	    (_Jet->chargedEmEnergyFraction->at(nobj) >= 0.99) )) return false;
  return true;
}

bool Analyzer::passTriggerCuts(string TriggerN) {
  if(prevTrig[TriggerN].first >= (int)Trigger_names->size() ||
     distats["Run"].smap[TriggerN+"FirstRequirement"] != Trigger_names->at(prevTrig[TriggerN].first)) {

    int temp = 0;
    vector<string>::iterator it = find(Trigger_names->begin(), Trigger_names->end(), distats["Run"].smap[TriggerN+"FirstRequirement"]);
    while(it != Trigger_names->begin()) {
      temp++; it--;
    }
    prevTrig[TriggerN].first = temp;
  }
  if(prevTrig[TriggerN].second >= (int)Trigger_names->size() ||
     distats["Run"].smap[TriggerN+"SecondRequirement"] != Trigger_names->at(prevTrig[TriggerN].second)) {
    int temp = 0;
    vector<string>::iterator it = find(Trigger_names->begin(), Trigger_names->end(), distats["Run"].smap[TriggerN+"SecondRequirement"]);
    while(it != Trigger_names->begin()) {
      temp++; 
      it--;
    }
    prevTrig[TriggerN].second = temp;
  }     
  if( (prevTrig[TriggerN].first < (int)Trigger_names->size() && Trigger_decision->at(prevTrig[TriggerN].first) == 1) || 
	(prevTrig[TriggerN].second < (int)Trigger_names->size() && Trigger_decision->at(prevTrig[TriggerN].second) == 1) ) return true;

  return false;
}


void Analyzer::SusyTopologyCuts() {
  PartStats* stats = &distats["VBFSUSY"];
  TLorentzVector ljet1 = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0));
  TLorentzVector ljet2 = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));
  
  if(stats->bmap["DiscrByLeadDiJetMass"]) {
    if((ljet1 + ljet2).M() < stats->pmap["LeadDiJetMassCut"].first) return;
    if((ljet1 + ljet2).M() > stats->pmap["LeadDiJetMassCut"].second) return;
  }

  if(stats->bmap["DiscrByLeadDiJetPt"]) {
    if((ljet1 + ljet2).Pt() < stats->pmap["LeadDiJetPtCut"].first) return;
    if((ljet1 + ljet2).Pt() > stats->pmap["LeadDiJetPtCut"].second) return;
  }
  
  if(stats->bmap["DiscrByLeadDiJetDeltaEta"]) {
    if(fabs(ljet1.Eta() - ljet2.Eta()) < stats->pmap["LeadDiJeTDeltaEtaCut"].first) return;
    if(fabs(ljet1.Eta() - ljet2.Eta()) > stats->pmap["LeadDiJetDeltaEtaCut"].second) return;
  }

  if(stats->bmap["DiscrByLeadDiJetDeltaPhi"]) {
    if(absnormPhi(ljet1.Phi() - ljet2.Phi()) < stats->pmap["LeadDiJetDeltaPhiCut"].first) return;
    if(absnormPhi(ljet1.Phi() - ljet2.Phi()) > stats->pmap["LeadDiJetDeltaPhiCut"].second) return;
  }
  
  if(stats->bmap["DiscrByLeadDiJetOSEta"]) {
    if((ljet1.Eta() * ljet2.Eta()) >= 0) return;
  }
  if(stats->bmap["DiscrByMet"]) {
    if(theMETVector.Pt() < stats->pmap["RecoMetCut"].first) return;
    if(theMETVector.Pt() > stats->pmap["RecoMetCut"].second) return;
  }
  if(stats->bmap["DiscrByMHT"] && sqrt(pow(sumpxForMht,2) + pow(sumpyForMht,2)) < stats->dmap["MhtCut"]) return;

  if(stats->bmap["DiscrByHT"] && sumptForHt < stats->dmap["HtCut"]) return; 

  double dphi1 = normPhi(ljet1.Phi() - theMETVector.Phi());
  double dphi2 = normPhi(ljet2.Phi() - theMETVector.Phi());
  double r1, r2, alpha;
  
  if(stats->bmap["DiscrByR1"]) {
    r1 = sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) );
    if(r1 < stats->pmap["R1Cut"].first || r1 > stats->pmap["R1Cut"].second) return;

  }
  if(stats->bmap["DiscrByR2"]) {
    r2 = sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0) );
    if(r2 < stats->pmap["R2Cut"].first || r2 > stats->pmap["R2Cut"].second) return;
  }
  if(stats->bmap["DiscrByAlpha"]) {
    TLorentzVector addVec = ljet1 + ljet2;
    alpha = (addVec.M() > 0) ? ljet2.Pt() / addVec.M() : -1;
    if(alpha < stats->pmap["AlphaCut"].first || alpha > stats->pmap["AlphaCut"].second) return;
  }
  if(stats->bmap["DiscrByDphi1"]) {
    if(abs(dphi1) < stats->pmap["Dphi1Cut"].first || abs(dphi1) > stats->pmap["Dphi1Cut"].second) return;
  }
  if(stats->bmap["DiscrByDphi2"]) {
    if(abs(dphi2) < stats->pmap["Dphi2Cut"].first || abs(dphi2) > stats->pmap["Dphi2Cut"].second) return;
  }

  goodParts[ival(CUTS::eSusyCom)].push_back(0);
}


void Analyzer::passRecoLeptonMetTopologyCuts(Lepton* lep, CUTS eReco, CUTS ePos, PartStats* stats) {
  Taus* tauTest = dynamic_cast<Taus*>(lep);

  for(vec_iter it=goodParts[ival(eReco)].begin(); it != goodParts[ival(eReco)].end(); it++) {
    if (tauTest == nullptr && stats->bmap["DiscrByIsZllCut"]) { 
      if(isZdecay(lep->smearP.at(*it), lep)) continue;
    }
    if (stats->bmap["DiscrByMetDphi"]) {
      double Dphi = absnormPhi(lep->smearP.at(*it).Phi() - theMETVector.Phi());
      if(Dphi < stats->pmap["MetDphiCut"].first || Dphi > stats->pmap["MetDphiCut"].second) continue;
    }
    if (stats->bmap["DoDiscrByMetMt"]) {
      double LeptonMetMt = CalculateLeptonMetMt(lep->smearP.at(*it));
      if( LeptonMetMt < stats->pmap["MetMtCut"].first || LeptonMetMt > stats->pmap["MetMtCut"].second) continue;
    }
    goodParts[ival(ePos)].push_back(*it);
  }
}


//-----Calculate lepton+met transverse mass
double Analyzer::CalculateLeptonMetMt(TLorentzVector Tobj) {
  double px = Tobj.Px() + theMETVector.Px();
  double py = Tobj.Py() + theMETVector.Py();
  double et = Tobj.Et() + TMath::Sqrt((theMETVector.Px() * theMETVector.Px()) + (theMETVector.Py() * theMETVector.Py()));
  double mt2 = et*et - (px*px + py*py);
  return (mt2 >= 0) ? sqrt(mt2) : -1;
}

/////all it does is add lorentz vectors. 
/////keep in case needed later
// TLorentzVector Analyzer::CalculateTheDiJet4Momentum(TLorentzVector* Tobj1, TLorentzVector* Tobj2) {
//   return (*Tobj1) + (*Tobj2);
// }



double Analyzer::DiParticleMass(TLorentzVector Tobj1, TLorentzVector Tobj2, string howCalc) {
  bool ratioNotInRange = false;
  TLorentzVector The_LorentzVect;

  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + theMETVector.Px())) - (Tobj2.Px() * (Tobj1.Py() + theMETVector.Py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + theMETVector.Py())) - (Tobj1.Py() * (Tobj2.Px() + theMETVector.Px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    ratioNotInRange=!((x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px() / x1) + (Tobj2.Px() / x2), (Tobj1.Py() / x1) + (Tobj2.Py() / x2), (Tobj1.Pz() / x1) + (Tobj2.Pz() / x2), (Tobj1.Energy() / x1) + (Tobj2.Energy() / x2) );
      return The_LorentzVect.M();
    }
  } 

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    double px = Tobj1.Px() + Tobj2.Px() + theMETVector.Px();
    double py = Tobj1.Py() + Tobj2.Py() + theMETVector.Py();
    double pz = Tobj1.Pz() + Tobj2.Pz();
    double e  = Tobj1.Energy() + Tobj2.Energy() + TMath::Sqrt(pow(theMETVector.Px(),2) + pow(theMETVector.Py(),2));
    The_LorentzVect.SetPxPyPzE(px, py, pz, e);
    return The_LorentzVect.M();
  }

  return (Tobj1 + Tobj2).M();
}

//////////////////////rename
bool Analyzer::passDiParticleApprox(TLorentzVector Tobj1, TLorentzVector Tobj2, string howCalc) {
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + theMETVector.Px())) - (Tobj2.Px() * (Tobj1.Py() + theMETVector.Py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + theMETVector.Py())) - (Tobj1.Py() * (Tobj2.Px() + theMETVector.Px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    return (x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.);
  } else {
    return true;
  }
}




/////abs for values
void Analyzer::passLeptonComboTopologyCut(Lepton* lep1, Lepton* lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, PartStats* stats) {
  for(vec_iter i1=goodParts[ival(ePos1)].begin(); i1 != goodParts[ival(ePos1)].end(); i1++) {
    for(vec_iter i2=goodParts[ival(ePos2)].begin(); i2 != goodParts[ival(ePos2)].end(); i2++) {

      if(stats->bmap["DiscrByDeltaR"] && (lep1->smearP.at(*i1).DeltaR(lep2->smearP.at(*i2))) < stats->dmap["DeltaRCut"]) continue;
      /////think about type for osls
      if(stats->dmap["DiscrByOSLSType"] == -1 && (lep1->charge->at(*i1) * lep2->charge->at(*i2) >= 0)) continue;
      else if(stats->dmap["DiscrByOSLSType"] == 1 && (lep1->charge->at(*i1) * lep2->charge->at(*i2) <= 0)) continue;

      if(stats->bmap["DiscrByCosDphi"]) {
	double Dphi = absnormPhi( lep1->smearP.at(*i1).Phi() - lep2->smearP.at(*i2).Phi());
	if(cos(Dphi) > stats->pmap["CosDphiCut"].first || cos(Dphi) < stats->pmap["CosDphiCut"].second) continue;
      }
  // ----Mass window requirement
      if (stats->bmap["DiscrByMassReco"]) {
      	double diMass = DiParticleMass(lep1->smearP.at(*i1),lep2->smearP.at(*i2), stats->smap["HowCalculateMassReco"]);
      	if( diMass < stats->pmap["MassCut"].first || diMass > stats->pmap["MassMaxCut"].second) continue;
      }

      if (stats->bmap["DiscrByCDFzeta2D"]) {
      	double CDFzeta = stats->dmap["PZetaCutCoefficient"] * CalculatePZeta(lep1->smearP.at(*i1), lep2->smearP.at(*i2)) 
	  + stats->dmap["PZetaVisCutCoeffiecient"] * CalculatePZetaVis(lep1->smearP.at(*i1), lep2->smearP.at(*i2));
      	if( CDFzeta < stats->pmap["CDFzeta2DCutValue"].first || CDFzeta > stats->pmap["CDFzeta2DCutValue"].second ) continue;
      }
      //////////abs on the difference????
      ///////////////////
      if (stats->bmap["DiscrByDeltaPtDivSumPt"]) {
	double ptDiv = (lep1->smearP.at(*i1).Pt() - lep2->smearP.at(*i2).Pt()) / (lep1->smearP.at(*i1).Pt() + lep2->smearP.at(*i2).Pt());
	if( ptDiv < stats->pmap["DeltaPtDivSumPtCut"].first || ptDiv > stats->pmap["DeltaPtDivSumPtCut"].second) continue;
      }

      if (stats->bmap["DiscrByDeltaPt"]) {
	double deltaPt = lep1->smearP.at(*i1).Pt() - lep2->smearP.at(*i2).Pt(); 
	if(deltaPt < stats->pmap["DeltaPtCutValue"].first || deltaPt > stats->pmap["DeltaPtCutValue"].second) continue;
      }
      ///Particlesp that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2 
      goodParts[ival(ePosFin)].push_back((*i1)*_Gen->pt->size() + (*i2));
    }
  }
}

/////still work out values stuff
void Analyzer::passDiJetTopologyCuts(PartStats* stats) {
  // ----Separation cut between jets (remove overlaps)
  for(vec_iter ij1=goodParts[ival(CUTS::eRJet1)].begin(); ij1 != goodParts[ival(CUTS::eRJet1)].end(); ij1++) {
    for(vec_iter ij2=goodParts[ival(CUTS::eRJet2)].begin(); ij2 != goodParts[ival(CUTS::eRJet2)].end(); ij2++) {
      if (stats->bmap["DoDiscrByDeltaR"]) {
	if(_Jet->smearP.at(*ij1).DeltaR(_Jet->smearP.at(*ij2)) < stats->dmap["DeltaRCut"]) continue;
      }
      if (stats->bmap["DoDiscrByDeltaEta"]) {
	if(fabs(_Jet->smearP.at(*ij1).Eta() - _Jet->smearP.at(*ij2).Eta()) < stats->pmap["DeltaEtaCut"].first) continue;
	if(fabs(_Jet->smearP.at(*ij1).Eta() - _Jet->smearP.at(*ij2).Eta()) > stats->pmap["DeltaEtaCut"].second) continue;
      }
      if (stats->bmap["DoDiscrByDeltaPhi"]) {
	if(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi()) < stats->pmap["DeltaPhiCut"].first) continue;
	if(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi()) > stats->pmap["DeltaPhiCut"].second) continue;
      } 
      if (stats->bmap["DoDiscrByOSEta"]) {
	if((_Jet->smearP.at(*ij1).Eta() * _Jet->smearP.at(*ij2).Eta()) >= 0) continue;
      }
      // ----Require both legs to be almost back-to-back in phi
      if (stats->bmap["DoDiscrByCosDphi"]) {
	if(cos(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi())) > stats->pmap["CosDphiMaxCut"].first) continue;
	if(cos(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi())) < stats->pmap["CosDphiMinCut"].second) continue;
      }
      // ----Mass window requirement
      if (stats->bmap["DoDiscrByDiJetMassReco"]) {
	if( ((_Jet->smearP.at(*ij1) + _Jet->smearP.at(*ij2)).M() < stats->pmap["MassCut"].first) || ((_Jet->smearP.at(*ij1) + _Jet->smearP.at(*ij2)).M() > stats->pmap["MassMaxCut"].first) ) continue;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2 
      goodParts[ival(CUTS::eDiJet)].push_back((*ij1)*_Jet->smearP.size() + (*ij2));
    }
  }
}


///////Only tested for if is Zdecay, can include massptasymmpair later?
bool Analyzer::isZdecay(TLorentzVector theObject, Lepton* lep) {
  bool eventIsZdecay = false;
  const float zMass = 90.1876;
  const float zWidth = 2.4952;
  float zmmPtAsymmetry = -10.;

  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.
  for(vector<TLorentzVector>::iterator lepit= lep->smearP.begin(); lepit != lep->smearP.end(); lepit++) {
    if(theObject.DeltaR(*lepit) < 0.3) continue;
    if(theObject == (*lepit)) continue;

    TLorentzVector The_LorentzVect = theObject + (*lepit);
    zmmPtAsymmetry = (theObject.Pt() - lepit->Pt()) / (theObject.Pt() + lepit->Pt());
    //    theMassPtAsymmPair = std::make_pair<float, float>(The_LorentzVect.M(), float(zmmPtAsymmetry));

    if( (abs(The_LorentzVect.M() - zMass) < 3.*zWidth) || (fabs(zmmPtAsymmetry) < 0.20) ) {
      eventIsZdecay = true;
      break;
    }
  }

  return eventIsZdecay;
}



double Analyzer::CalculatePZeta(TLorentzVector Tobj1, TLorentzVector Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double px = visPx + theMETVector.Px();
  double py = visPy + theMETVector.Py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}

double Analyzer::CalculatePZetaVis(TLorentzVector Tobj1, TLorentzVector Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}

double Analyzer::normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}

double Analyzer::absnormPhi(double phi) {
  return abs(normPhi(phi));
}



// void BSM3GAnalyzer::fillHistograms(string group, vector<int> group_list) {
//   for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++){
//     double weight = isrgluon_weight * pdfWeightVector.at(NpdfID);
//       _hNPVertices[i][NpdfID]->Fill(bestVertices,weight);
//     }

    
    
//       int nGenTaus = 0;
//       for(int j = 0; j < Gen_pt->size(); j++) {
//         if(abs(Gen_pdg_id->at(j)) == 15) {

//           _hGenTauStatusCode[i][NpdfID]->Fill(Gen_status->at(j),weight);
//         }
//         if((abs(Gen_pdg_id->at(j)) == 15) && (Gen_status->at(j) == 2)) {
//           nGenTaus++;
// 	}
//       }
    
//       _hNGenTau[i][NpdfID]->Fill(nGenTaus,weight);

//       int nGenHadTaus = 0;
//       bool IsItAHadronicDecay; 
//       TLorentzVector theGenObject(0,0,0,0); 
//       TLorentzVector theNeutrinoObject(0,0,0,0);
//       vector<bool> IsItAHadronicDecayVector; 
//       IsItAHadronicDecayVector.clear();
//       vector<int> tempTauIndexVector; 
//       tempTauIndexVector.clear();
//       vector<TLorentzVector> tempNeutrinoMomentumVector; 
//       tempNeutrinoMomentumVector.clear();
//       for(int j = 0; j < Gen_pt->size(); j++) {
//         if( (abs(Gen_pdg_id->at(j)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(j))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(j)) == 2) ) {
//           tempTauIndexVector.push_back(Gen_BmotherIndex->at(j));
//           theNeutrinoObject.SetPtEtaPhiE(Gen_pt->at(j), Gen_eta->at(j), Gen_phi->at(j), Gen_energy->at(j));
//           tempNeutrinoMomentumVector.push_back(theNeutrinoObject);
//         }
//       }
//       if(tempTauIndexVector.size() > 0) {
//         for(int jj = 0; jj < tempTauIndexVector.size(); jj++) {
//           IsItAHadronicDecay = true;
//           for(int j = 0; j < Gen_pt->size(); j++) {
//             if( ((abs(Gen_pdg_id->at(j)) == 12) || (abs(Gen_pdg_id->at(j)) == 14)) && (Gen_BmotherIndex->at(j) == tempTauIndexVector.at(jj)) ) {
//               IsItAHadronicDecay = false; // it is not a hadronic tau decay since it decayed to a electron/muon neutrino
//             }
//           }
//           IsItAHadronicDecayVector.push_back(IsItAHadronicDecay);
//         }
//         for(int jj = 0; jj < tempTauIndexVector.size(); jj++) {
//           for(int j = 0; j < Gen_pt->size(); j++) {
//             if(j == tempTauIndexVector.at(jj)) {
//               theGenObject.SetPtEtaPhiE(Gen_pt->at(j), Gen_eta->at(j), Gen_phi->at(j), Gen_energy->at(j));
//               theGenObject = theGenObject - tempNeutrinoMomentumVector.at(jj);
//               if( (IsItAHadronicDecayVector.at(jj)) ) {
//                 nGenHadTaus++;
//                 _hGenHadTauPt[i][NpdfID]->Fill(theGenObject.Pt(),weight);
//                 _hGenHadTauEta[i][NpdfID]->Fill(theGenObject.Eta(),weight);
//               }
//             }
//           }
//         }
//       }
//       _hNGenHadTau[i][NpdfID]->Fill(nGenHadTaus,weight);

//       _hNGenMuon[i][NpdfID]->Fill(nGenMuons,weight);

//       for(int j = 0; j < Gen_pt->size(); j++) {
//         if(abs(Gen_pdg_id->at(j)) == 32) {
//             _hGenZprimeStatusCode[i][NpdfID]->Fill(Gen_status->at(j),weight);
//         }
//         if((abs(Gen_pdg_id->at(j)) == 32) && (Gen_status->at(j) != 22)) {
//           TLorentzVector genObjt1;
//           genObjt1.SetPtEtaPhiE(Gen_pt->at(j), Gen_eta->at(j), Gen_phi->at(j), Gen_energy->at(j));
//           _hGenZprimeMass[i][NpdfID]->Fill(genObjt1.M(),weight);
//         }
//       }


//       for(int j = 0; j < Gen_pt->size(); j++) {
//         for(int jj = 0; jj < Gen_pt->size(); jj++) {
//           if((abs(Gen_pdg_id->at(j)) == 15) && (Gen_status->at(j) == 2) && (abs(Gen_pdg_id->at(jj)) == 15) && (Gen_status->at(jj) == 2) && (jj > j)) {
//             TLorentzVector genObjt1;
//             genObjt1.SetPtEtaPhiE(Gen_pt->at(j), Gen_eta->at(j), Gen_phi->at(j), Gen_energy->at(j));
//             TLorentzVector genObjt2;
//             genObjt2.SetPtEtaPhiE(Gen_pt->at(jj), Gen_eta->at(jj), Gen_phi->at(jj), Gen_energy->at(jj));
//             _hGenDiTauMass[i][NpdfID]->Fill(CalculateTheDiJet4Momentum(genObjt1,genObjt2).second.M(),weight);
//           }
//         }
//       }

//     }



//       // ------Central Jet Histograms
//       int nCJets = 0;
//       for(int j = 0; j < Jet_pt->size(); j++) {
// 	if (!passRecoCentralJetCuts(j)) continue;
// 	_hCentralJetPt[i][NpdfID]->Fill(smearedJetMomentumVector.at(j).Pt(),weight);
// 	_hCentralJetEta[i][NpdfID]->Fill(smearedJetMomentumVector.at(j).Eta(),weight);
// 	nCJets++;

//       _hNCentralJet[i][NpdfID]->Fill(nCJets,weight);

  
  // if(group == "FillGenTau") {
  //   for(int i = 0; i < goodParts[ival(eGTau)].size(); i++) {
  //     histo.addVal(Gen->energy->at(i), group, );   //GenTauEnergy
  //     histo.addVal(Gen->pt->at(i), group, );   //GenTauPt
  //     histo.addVal(Gen->eta->at(i),group, );   //GenTauEta
  //     histo.addVal(Gen->phi->at(i),group, );   //GenTauPhi
  //   }

  // } else if(group == "FillGenMuon") {
  //   for(int i = 0; i < goodParts[ival(eGMuon)].size(); i++) {

  // 	switch(j) {

  // 	case 0: Fill(Gen->energy->at(i));   //GenMuonEnergy
  //       case 1: Fill(Gen_pt->at(i));   //GenMuonPt
  //       case 2: Fill(Gen_eta->at(i));   //GenMuonEta
  //       case 3: Fill(Gen_phi->at(i));   //GenMuonPhi

  // 	}
  //     }
  //   }

  //     if(Muon->smearP.at(j).Pt() >= leadingmuonpt) {
  // 	leadingmuonpt = Muon->smearP.at(j).Pt();
  // 	leadingmuoneta = Muon->smearP.at(j).Eta();
  //   if(nMuons > 0) {
  //     _hFirstLeadingMuon1Pt->Fill(leadingmuonpt);
  //     _hFirstLeadingMuon1Eta->Fill(leadingmuoneta);
  //   }
  //////////num with pt/eta/phi/e 


void Analyzer::fill_histogram() {
  int maxCut = fillCuts();
  vector<string>* groups = histo->get_groups();

  for(vector<string>::iterator it = groups->begin(); it!=groups->end(); it++) {
    fill_Folder(*it, maxCut);
  }
}


void Analyzer::fill_Folder(string group, int max) {

  if(group == "FillTauJet1" || group == "FillTauJet2" || group == "FillMuon1" || group == "FillMuon2" || group == "FillJet1" || group == "FillJet2" || group == "FillBJet") {
    Particle* part;
    if(group == "FillTauJet1" || group == "FillTauJet2") part=_Tau;
    else if(group == "FillMuon1" || group == "FillMuon2") part=_Muon;
    else part = _Jet;
    CUTS ePos = fill_num[group];

    for(vec_iter it=goodParts[ival(ePos)].begin(); it!=goodParts[ival(ePos)].end(); it++) {
      histo->addVal(part->smearP.at(*it).Energy(), group,max, "Energy");
      histo->addVal(part->smearP.at(*it).Pt(), group,max, "Pt");
      histo->addVal(part->smearP.at(*it).Eta(), group,max, "Eta");
      histo->addVal(part->smearP.at(*it).Phi(), group,max, "Phi");
      if(dynamic_cast<Taus*>(part) != NULL) {
	histo->addVal(_Tau->nProngs->at(*it), group,max, "NumSignalTracks");
  	histo->addVal(_Tau->charge->at(*it), group,max, "Charge");
	histo->addVal(_Tau->leadChargedCandPt->at(*it), group,max, "SeedTracks");
      } else if(dynamic_cast<Muon*>(part)) {
	histo->addVal(CalculateLeptonMetMt(_Muon->smearP.at(*it)), group,max, "MetMt");  
      }
    }

    histo->addVal(goodParts[ival(ePos)].size(), group,max, "N");

  } else if(group == "FillSusyCuts") {

    histo->addVal(sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)),group,max, "MHT");
    histo->addVal(sumptForHt,group,max, "HT");  
    histo->addVal(sumptForHt + sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)),group,max, "Meff");
    histo->addVal(theMETVector.Pt(), group,max, "Met");
    if(goodParts[ival(CUTS::eR1stJet)].size() * goodParts[ival(CUTS::eR2ndJet)].size() != 0) {
      TLorentzVector DiJet = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0)) + _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));
      histo->addVal(absnormPhi(theMETVector.Phi() - DiJet.Phi()), group,max, "MetDiJetDeltaPhi");
    }
  


  } else if(group == "FillLeadingJet" && goodParts[ival(CUTS::eSusyCom)].size() != 0) {
    TLorentzVector first = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0));
    TLorentzVector second = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));

    histo->addVal(first.Pt(),group,max, "FirstPt");
    histo->addVal(second.Pt(),group,max, "SecondPt");
    histo->addVal(first.Eta(),group,max, "FirstEta");
    histo->addVal(second.Eta(),group,max, "SecondEta");
    
    TLorentzVector LeadDiJet = first + second;
    
    histo->addVal(LeadDiJet.M(), group,max, "Mass"); 
    histo->addVal(LeadDiJet.Pt(), group,max, "Pt");  
    histo->addVal(fabs(first.Eta() - second.Eta()), group,max, "DeltaEta"); 
    histo->addVal(first.DeltaR(second),group,max, "DeltaR");  

    double dphiDijets = absnormPhi(second.Phi() - second.Phi());
    double dphi1 = normPhi(first.Phi() - theMETVector.Phi());
    double dphi2 = normPhi(second.Phi() - theMETVector.Phi());
    double alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histo->addVal(dphiDijets,group,max, "LeadSublDijetDphi"); 
    //////histo->addVal(theMETVector.Pt(),dphiDijets, group,max, "");   //MetVsDiJetDeltaPhiLeadSubl
    /////histo->addVal(fabs(first.Eta() - second.Eta()), dphiDijets, group,max, "");   //DeltaEtaVsDeltaPhiLeadSubl

    histo->addVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), group,max, "R1");
    histo->addVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), group,max, "R2");
    histo->addVal(normPhi(first.Phi() - phiForMht), group,max, "Dphi1MHT"); 
    histo->addVal(normPhi(second.Phi() - phiForMht), group,max, "Dphi2MHT");
    histo->addVal(dphi1,group,max, "Dphi1");
    histo->addVal(dphi2,group,max, "Dphi2");
    ///histo->addVal(dphi1,dphi2,group,max, "");   //Dphi1VsDphi2
    histo->addVal(alpha,group,max, "Alpha");


    //dijet info
  } else if(group == "FillDiJet") {
    double leaddijetmass = 0;
    double leaddijetpt = 0;
    double leaddijetdeltaR = 0;
    double leaddijetdeltaEta = 0;
    double etaproduct = -100;
    for(vec_iter it=goodParts[ival(CUTS::eDiJet)].begin(); it!=goodParts[ival(CUTS::eDiJet)].end(); it++) {
      int p1 = (*it) / _Jet->smearP.size();
      int p2 = (*it) % _Jet->smearP.size();
      TLorentzVector jet1 = _Jet->smearP.at(p1);
      TLorentzVector jet2 = _Jet->smearP.at(p2);
      TLorentzVector DiJet = jet1 + jet2;
	  
      if(DiJet.M() > leaddijetmass) {
	leaddijetmass = DiJet.M();
	etaproduct = (jet1.Eta() * jet2.Eta() > 0) ? 1 : -1;
      }
      leaddijetpt = (DiJet.Pt() > leaddijetpt) ? DiJet.Pt() : leaddijetpt;
      leaddijetdeltaEta = (fabs(jet1.Eta() - jet2.Eta()) > leaddijetdeltaEta) ? fabs(jet1.Eta() - jet2.Eta()) : leaddijetdeltaEta;
      leaddijetdeltaR = (jet1.DeltaR(jet2) > leaddijetdeltaR) ? jet1.DeltaR(jet2) : leaddijetdeltaR;

      histo->addVal(DiJet.M(), group,max, "Mass");
      histo->addVal(DiJet.Pt(), group,max, "Pt");
      histo->addVal(fabs(jet1.Eta() - jet2.Eta()), group,max, "DetlaEta");
      histo->addVal(absnormPhi(jet1.Phi() - jet2.Phi()), group,max, "DeltaPhi");
      histo->addVal(jet1.DeltaR(jet2), group,max, "DeltaR");
    }
  
    histo->addVal(leaddijetmass, group,max, "LeadMass");
    histo->addVal(leaddijetpt, group,max, "LeadPt");  
    histo->addVal(leaddijetdeltaEta, group,max, "LeadDeltaEta");
    histo->addVal(leaddijetdeltaR, group,max, "LeadDeltaR");
    histo->addVal(etaproduct, group,max, "LeadEtaProduct");


    ////diparticle stuff
  } else if(group == "FillDiMuon" || group == "FillDiTau" || group == "Muon1Tau1" || group == "Muon1Tau2" || group == "Muon2Tau1" || group == "Muon2Tau2") {  ///mumu/mutau/tautau
    Lepton* lep1 = NULL;
    Lepton* lep2 = NULL;
    CUTS ePos = fill_num[group];
    if(ePos == CUTS::eMuon1Tau1 || ePos == CUTS::eMuon1Tau2 || ePos == CUTS::eMuon2Tau1 || ePos == CUTS::eMuon2Tau2) {
      lep1 = _Muon; lep2 = _Tau;
    } else if(ePos == CUTS::eElec1Tau1 || ePos == CUTS::eElec1Tau2 || ePos == CUTS::eElec2Tau1 || ePos == CUTS::eElec2Tau2) {
      lep1 = _Electron; lep2 = _Tau;
    } else if(ePos == CUTS::eDiMuon) {
      lep1 = _Muon; lep2 = _Muon;
    } else if(ePos == CUTS::eDiTau) { lep1 = _Tau; lep2 = _Tau; 
    } else if (ePos == CUTS::eDiElec) { lep1 = _Electron; lep2 = _Electron; }

    
    for(vec_iter it=goodParts[ival(ePos)].begin(); it!=goodParts[ival(ePos)].end(); it++) {

      int p1= (*it) / _Gen->pt->size();
      int p2= (*it) % _Gen->pt->size();

      //    histo->addVal(Muon->smearP.at(mj).Pt(),Tau->smearP.at(tj).Pt());   //Muon1PtVsTau1Pt######################
      histo->addVal(lep1->smearP.at(p1).DeltaR(lep2->smearP.at(p2)), group,max, "DeltaR"); 
      histo->addVal((lep1->smearP.at(p1).Pt() - lep1->smearP.at(p2).Pt()) / (lep1->smearP.at(p1).Pt() + lep2->smearP.at(p2).Pt()), group,max, "DeltaPtDivSumPt");   //Muon1Tau1DeltaPtDivSumPt
      histo->addVal(lep1->smearP.at(p1).Pt() - lep2->smearP.at(p2).Pt(), group,max, "DeltaPt");
      histo->addVal(cos(absnormPhi(lep2->smearP.at(p2).Phi() - lep1->smearP.at(p1).Phi())), group,max, "CosDphi");
      histo->addVal(absnormPhi(lep2->smearP.at(p2).Phi() - theMETVector.Phi()), group,max, "Part1MetDeltaPhi");
      ///      histo->addVal(absnormPhi(lep2->smearP.at(p2).Phi() - theMETVector.Phi()), cos(absnormPhi(lep2->smearP.at(p2).Phi() - lep1->smearP.at(p1).Phi())));   //Muon1MetDeltaPhiVsMuon1Tau1CosDphi
      histo->addVal(absnormPhi(lep1->smearP.at(p1).Phi() - theMETVector.Phi()), group,max, "Part2MetDeltaPhi");
      // if(CalculateTheDiTau4Momentum(lep1->smearP.at(p1),lep2->smearP.at(p2)).first) {
      // 	histo->addVal(DiParticleMass(lep1->smearP.at(p1),lep2->smearP.at(p2)), group,max, "ReconstructableMass");
      // } else {
      // 	histo->addVal(DiParticleMass(lep1->smearP.at(p1),lep2->smearP.at(p2)), group,max, "NotReconstructableMass");
      // }
      double PZeta = CalculatePZeta(lep1->smearP.at(p1),lep2->smearP.at(p2));
      double PZetaVis = CalculatePZetaVis(lep1->smearP.at(p1),lep2->smearP.at(p2));
      histo->addVal(CalculateLeptonMetMt(lep2->smearP.at(p2)), group,max, "Part1MetMt");
      histo->addVal(CalculateLeptonMetMt(lep1->smearP.at(p1)), group,max, "Part2MetMt"); 
      histo->addVal(lep2->charge->at(p2) * lep1->charge->at(p1), group,max, "OSLS");  
      histo->addVal(PZeta, group,max, "PZeta"); 
      histo->addVal(PZetaVis, group,max, "PZetaVis");  
      // histo->addVal(PZetaVis,PZeta, group,max, "Zeta2D");   //Muon1Tau1Zeta2D
      ///histo->addVal((_Muon1Tau1PZetaCutCoefficient * PZeta) + (_Muon1Tau1PZetaVisCutCoefficient * PZetaVis));   //Muon

      if ((goodParts[ival(CUTS::eR1stJet)].size() != 0) && (goodParts[ival(CUTS::eR2ndJet)].size() != 0)) {
	TLorentzVector TheLeadDiJetVect = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0)) + _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].size());
	histo->addVal(absnormPhi(lep1->smearP.at(p1).Phi() - TheLeadDiJetVect.Phi()), group,max, "Part1DiJetDeltaPhi");
	histo->addVal(absnormPhi(lep2->smearP.at(p2).Phi() - TheLeadDiJetVect.Phi()), group,max, "Part2DiJetDeltaPhi");
	////////#####histo->addVal(DiParticleMass(TheLeadDiJetVect, Muon->smearP.at(mj)+Tau->smearP.at(tj) ));   //Muon1Tau1DiJetReconstructableMass
      }
      // Fill(isZmm(Muon->smearP.at(mj)).first);   //Muon1Tau1_Muon1IsZmm
    }
  }
}


void Analyzer::initializePileupInfo(string MCHisto, string DataHisto) {
  // Filenames must be c_strings below. Here is the conversion from strings to c_strings
  // As you can see above cstr1 corresponds to MC and cstr2 corresponds to data.

  TFile *file = TFile::Open(MCHisto.c_str());
  TH1* histmc = dynamic_cast<TH1*>(file->Get("analyzeHiMassTau/NVertices_0"));
  if(!histmc) {throw std::runtime_error("failed to extract histogram");}
  for(int bin=0; bin<=(histmc->GetXaxis()->GetNbins() + 1); bin++) {
    hPUmc->SetBinContent(bin,histmc->GetBinContent(bin));
  }
  file->Close();

  file = TFile::Open(DataHisto.c_str());
  TH1* histdata = dynamic_cast<TH1*>(file->Get("analyzeHiMassTau/NVertices_0"));
  if(!histdata) {throw std::runtime_error("failed to extract histogram");}
  for(int bin=0; bin<=(histdata->GetXaxis()->GetNbins() + 1); bin++) {
    hPUdata->SetBinContent(bin,histdata->GetBinContent(bin));
  }
  file->Close();
}

double Analyzer::getPileupWeight(float ntruePUInt) {
  int bin;
  double MCintegral;
  double MCvalue;
  double Dataintegral;
  double Datavalue;

  // The probability that data (or MC) has N pileup interactions is value / integral
  // The ratio of the data and MC probability density functions gives us our pileup weights

  //std::cout << "Grabbing pileup info. " << std::endl;
  bin = hPUmc->GetBin(ntruePUInt+1);
  MCvalue = hPUmc->GetBinContent(bin);
  MCintegral = hPUmc->Integral();
  Datavalue = hPUdata->GetBinContent(bin);
  Dataintegral = hPUdata->Integral();

  // printouts for debugging
  //std::cout << "Number of true pileup interactions = " << ntruePUInt << std::endl;
  //std::cout << "Histogram bin, given the number of true pileup interactions = " << bin << std::endl;
  //std::cout << "MC PU probability density, given the number of true pileup interactions = " << MCvalue << std::endl;
  //std::cout << "Data PU probability density, given the number of true pileup interactions = " << Datavalue << std::endl;

  //std::cout << "Grabbing pileup weight. " << std::endl;
  //Ratio of normalized histograms in given bin

  return ((MCvalue * Dataintegral) != 0) ? (Datavalue * MCintegral) / (MCvalue * Dataintegral) : 1.0;
}





