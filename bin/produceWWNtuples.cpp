#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"

#include "../interface/setInputTree.h"
//#include "../interface/initInputTree.h"
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/analysisUtils.h"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  bool isMC = argv[3];
  std::string leptonName = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }
    
  //--------input tree-----------
  //  TChain* chain = new TChain("TreeMaker2/PreSelection");
  //  TChain* chain = new TChain(inputTreeName.c_str());
  //  InitTree(chain);
  //  chain->Add((inputFolder+inputFile).c_str());

  TFile *MyFile = new TFile((inputFolder+inputFile).c_str(),"READ");
  setInputTree *WWTree = new setInputTree (MyFile, inputTreeName.c_str());
  if (WWTree->fChain == 0) return (-1);
  WWTree->Init();
  //  TTree *chain = (TTree *) MyFile->Get(inputTreeName.c_str());
  //  setInputTree(chain);

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_")+leptonName+std::string("/")+outputFile).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);
  SetOutTree(outTree);

  //---------start loop on events------------
  for (Long64_t jentry=0; jentry<WWTree->fChain->GetEntries();jentry++) {

    Long64_t iEntry = WWTree->LoadTree(jentry);
    if (iEntry < 0) break;
    int nb = WWTree->fChain->GetEntry(jentry);   
    // if (Cut(ientry) < 0) continue;                                                                                                                           

    if(iEntry % 1000 == 0)    
      cout << "read entry: " << iEntry << endl;

    init(); //initialize all variables

    //require exactly one lepton
    if ( strcmp(leptonName.c_str(),"el")==0 && WWTree->selectedIDIsoElectronsNum!=1) continue; 
    if ( strcmp(leptonName.c_str(),"mu")==0 && WWTree->selectedIDIsoMuonsNum!=1) continue;      
    
    if (WWTree->AK8JetsNum < 1) continue; //at least one jet

    if (WWTree->AK8JetsPt[0] < 150) continue; 
    if (WWTree->METPt < 50) continue;

    //lepton Pt selection
    if ( strcmp(leptonName.c_str(),"el")==0 && WWTree->selectedIDIsoElectronsPt[0]<35) continue; 
    if ( strcmp(leptonName.c_str(),"mu")==0 && WWTree->selectedIDIsoMuonsPt[0]<30) continue; 
    
    //save event variables
    event_runNo   = WWTree->RunNum;
    event = WWTree->EvtNum;
    njets = WWTree->NJets;
    nPV  = WWTree->NVtx;
    
    /////////////////LEPTON
    if (strcmp(leptonName.c_str(),"el")==0) {
	l_pt  = WWTree->selectedIDIsoElectronsPt[0];
	l_eta = WWTree->selectedIDIsoElectronsEta[0];
	l_phi = WWTree->selectedIDIsoElectronsPhi[0];	
	l_e = WWTree->selectedIDIsoElectronsE[0];	
      }
    else if (strcmp(leptonName.c_str(),"mu")==0) {
	l_pt  = WWTree->selectedIDIsoMuonsPt[0];
	l_eta = WWTree->selectedIDIsoMuonsEta[0];
	l_phi = WWTree->selectedIDIsoMuonsPhi[0];
	l_e = WWTree->selectedIDIsoMuonsE[0];
      }

    //////////////MET

    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*

    float Wmass = 80.385;

    TLorentzVector W_mu, W_Met;

    W_mu.SetPtEtaPhiE(l_pt,l_eta,l_phi,l_e);
    W_Met.SetPxPyPzE(WWTree->METPt * TMath::Cos(WWTree->METPhi), WWTree->METPt * TMath::Sin(WWTree->METPhi), 0., sqrt(WWTree->METPt*WWTree->METPt));

    if(W_mu.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }


    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_mu);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());

    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    double pz2_type0 = NeutrinoPz_type0.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    //    W_mass_type0_met = (W_neutrino_type0_met+W_mu).M();
    //    W_pz_type0_met = (W_neutrino_type0_met+W_mu).Pz();
    //    W_nu1_pz_type0_met = pz1_type0;
    //    W_nu2_pz_type0_met = pz2_type0;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(WWTree->METPhi),
			      nu_pt1 * TMath::Sin(WWTree->METPhi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(WWTree->METPhi),
			      nu_pt2 * TMath::Sin(WWTree->METPhi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }

    //    W_mass_type0 = (W_mu+W_neutrino_type0).M();
    //    W_pz_type0 = (W_mu+W_neutrino_type0).Pz();
    //    W_nu1_pz_type0 = pz1_type0;
    //    W_nu2_pz_type0 = pz2_type0;

    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(W_mu);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());

    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    double pz2_type2 = NeutrinoPz_type2.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    //    W_mass_type2_met = (W_neutrino_type2_met+W_mu).M();
    //    W_pz_type2_met = (W_neutrino_type2_met+W_mu).Pz();
    //    W_nu1_pz_type2_met = pz1_type2;
    //    W_nu2_pz_type2_met = pz2_type2;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2; 
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    if (NeutrinoPz_type2.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(WWTree->METPhi),
			      nu_pt1 * TMath::Sin(WWTree->METPhi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(WWTree->METPhi),
			      nu_pt2 * TMath::Sin(WWTree->METPhi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }

    //    W_mass_type2 = (W_mu+W_neutrino_type2).M();
    //    W_pz_type2 = (W_mu+W_neutrino_type2).Pz();
    //    W_nu1_pz_type2 = pz1_type2;
    //    W_nu2_pz_type2 = pz2_type2;

    pfMET   = sqrt(WWTree->METPt*WWTree->METPt);
    pfMET_Phi = WWTree->METPhi;
    nu_pz_type0 = pz1_type0;
    nu_pz_type2 = pz1_type2;


    /////////////////LEPTONIC W

    TLorentzVector *W = new TLorentzVector();
    TLorentzVector *LEP = new TLorentzVector();
    TLorentzVector *NU0  = new TLorentzVector();
    TLorentzVector *NU2  = new TLorentzVector();
    
    LEP->SetPtEtaPhiE(l_pt,l_eta,l_phi,l_e);
    NU0->SetPxPyPzE(WWTree->METPt*TMath::Cos(WWTree->METPhi),WWTree->METPt*TMath::Sin(WWTree->METPhi),nu_pz_type0,pfMET);
    NU2->SetPxPyPzE(WWTree->METPt*TMath::Cos(WWTree->METPhi),WWTree->METPt*TMath::Sin(WWTree->METPhi),nu_pz_type2,pfMET);
    *W = *LEP + *NU0;
    
    v_pt = W->Pt();
    v_eta = W->Eta();
    v_phi = W->Phi();
    v_mt = TMath::Sqrt(2*LEP->Et()*NU0->Et()*(1-TMath::Cos(LEP->DeltaPhi(*NU0))));
    //    W_mt = W->Mt();

    //////////////////ANGULAR VARIABLES

    TLorentzVector *JET = new TLorentzVector();
    JET->SetPtEtaPhiE(WWTree->AK8JetsPt[0],WWTree->AK8JetsEta[0],WWTree->AK8JetsPhi[0],WWTree->AK8JetsE[0]);
    deltaR_lak8jet = JET->DeltaR(*LEP);
    deltaphi_METak8jet = JET->DeltaPhi(*NU0);
    deltaphi_Vak8jet = JET->DeltaPhi(*W);

    //FOUR-BODY INVARIANT MASS
    mass_lvj_type0 = (*LEP + *NU0 + *JET).M();
    mass_lvj_type2 = (*LEP + *NU2 + *JET).M();

    //delete all the TLorentzVector before a new selection
    delete W;
    delete LEP;
    delete NU0;
    delete NU2;
    delete JET;

    if (v_pt < 150) continue;
    if (deltaR_lak8jet < (TMath::Pi()/2.0))   continue;

    ///////////JETS
    float tempPt=0.;
    for (unsigned int i=0; i<WWTree->AK8JetsNum; i++)
      {
	if (WWTree->AK8JetsPt[i]<30 || WWTree->AK8JetsEta[i]>4.7)  continue;
	if (WWTree->AK8JetsPt[i]<=tempPt) continue; //to save the jet with largest pt
	ungroomed_jet_pt  = WWTree->AK8JetsPt[i];
	ungroomed_jet_eta = WWTree->AK8JetsEta[i];
	ungroomed_jet_phi = WWTree->AK8JetsPhi[i];
	ungroomed_jet_e   = WWTree->AK8JetsE[i];
	jet_mass_pr   = WWTree->AK8Jets_prunedMass[i];
	jet_mass_tr   = WWTree->AK8Jets_trimmedMass[i];
	jet_mass_fi   = WWTree->AK8Jets_filteredMass[i];
	jet_tau2tau1   = WWTree->AK8Jets_tau2[i]/WWTree->AK8Jets_tau1[i];
	tempPt = ungroomed_jet_pt;
      }

    /////////VBF PART
    bool fillVBF = true;

    TLorentzVector *HADW = new TLorentzVector();
    HADW->SetPtEtaPhiE(ungroomed_jet_pt,ungroomed_jet_eta,ungroomed_jet_phi,ungroomed_jet_e); //AK8 fat jet (hadronic W)
    TLorentzVector *AK4 = new TLorentzVector();
    std::vector<int> indexGoodJets;
    indexGoodJets.clear();
    if (indexGoodJets.size()!=0)  fillVBF=false;

    for (unsigned int i=0; i<WWTree->JetsNum; i++) //loop on AK4 jet
      {
	if (WWTree->JetsPt[i]<30 || WWTree->JetsEta[i]>4.7)  continue;
	AK4->SetPtEtaPhiE(WWTree->JetsPt[i],WWTree->JetsEta[i],WWTree->JetsPhi[i],WWTree->JetsE[i]);
	float deltaR = HADW->DeltaR(*AK4);
	if (deltaR<0.8) continue; //the vbf jets must be outside the had W cone
	indexGoodJets.push_back(i); //save index of the "good" vbf jets candidate
      }
    if (indexGoodJets.size()<2)  fillVBF=false; //check if at least 2 jets are inside the collection
    delete HADW;
    delete AK4;

    if (fillVBF) 
      {
	TLorentzVector *VBF1 = new TLorentzVector();
	TLorentzVector *VBF2 = new TLorentzVector();
	TLorentzVector *TOT = new TLorentzVector();
	float tempPtMax=0.;
	int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
    
	for (unsigned int i=0; i<indexGoodJets.size()-1; i++) {
	  for (unsigned int ii=i+1; ii<indexGoodJets.size(); ii++) {
	    VBF1->SetPtEtaPhiE(WWTree->JetsPt[indexGoodJets.at(i)],WWTree->JetsEta[indexGoodJets.at(i)],WWTree->JetsPhi[indexGoodJets.at(i)],WWTree->JetsE[indexGoodJets.at(i)]);
	    VBF2->SetPtEtaPhiE(WWTree->JetsPt[indexGoodJets.at(ii)],WWTree->JetsEta[indexGoodJets.at(ii)],WWTree->JetsPhi[indexGoodJets.at(ii)],WWTree->JetsE[indexGoodJets.at(ii)]);
	    *TOT = *VBF1 + *VBF2;
	    if (TOT->Pt() < tempPtMax) continue;
	    tempPtMax = TOT->Pt(); //take the jet pair with largest Pt
	    nVBF1 = indexGoodJets.at(i); //save position of the 1st vbf jet
	    nVBF2 = indexGoodJets.at(ii); //save position of the 2nd vbf jet
	  }
	}
	
	if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
	  {
	    VBF1->SetPtEtaPhiE(WWTree->JetsPt[nVBF1],WWTree->JetsEta[nVBF1],WWTree->JetsPhi[nVBF1],WWTree->JetsE[nVBF1]);
	    VBF2->SetPtEtaPhiE(WWTree->JetsPt[nVBF2],WWTree->JetsEta[nVBF2],WWTree->JetsPhi[nVBF2],WWTree->JetsE[nVBF2]);
	    *TOT = *VBF1 + *VBF2;
	    
	    vbf_maxpt_j1_pt = WWTree->JetsPt[nVBF1];
	    vbf_maxpt_j1_eta = WWTree->JetsEta[nVBF1];
	    vbf_maxpt_j1_phi = WWTree->JetsPhi[nVBF1];
	    vbf_maxpt_j1_e = WWTree->JetsE[nVBF1];
	    vbf_maxpt_j1_bDiscriminatorCSV = WWTree->Jets_bDiscriminator[nVBF1];
	    vbf_maxpt_j2_pt = WWTree->JetsPt[nVBF2];
	    vbf_maxpt_j2_eta = WWTree->JetsEta[nVBF2];
	    vbf_maxpt_j2_phi = WWTree->JetsPhi[nVBF2];
	    vbf_maxpt_j2_e = WWTree->JetsE[nVBF2];
	    vbf_maxpt_j2_bDiscriminatorCSV = WWTree->Jets_bDiscriminator[nVBF2];
	    vbf_maxpt_jj_pt = TOT->Pt();
	    vbf_maxpt_jj_eta = TOT->Eta();
	    vbf_maxpt_jj_phi = TOT->Phi();
	    vbf_maxpt_jj_m = TOT->M();	
	  }
	
	delete VBF1;
	delete VBF2;
	delete TOT;
      }

    /////////////////MC Infos
    if (isMC)
      {
	TLorentzVector temp, temp2;
	//	std::cout<<"entry: "<<iEntry<<" "<<GenNuNum<<std::endl;
	double deltaPhiOld=100.;
	for (int i=0; i<WWTree->GenBosonNum; i++) {
	  double deltaPhi = getDeltaPhi(WWTree->GenBosonPhi[i],v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;
	  //	  std::cout<<"bosone: "<<i<<" "<<WWTree->GenBosonPhi[i]<<" "<<v_phi<<std::endl;
	  temp.SetPtEtaPhiE(WWTree->GenBosonPt[i],WWTree->GenBosonEta[i],WWTree->GenBosonPhi[i],WWTree->GenBosonE[i]);
	  W_pt_gen = WWTree->GenBosonPt[i];
	  W_pz_gen = temp.Pz();
	  deltaPhiOld = deltaPhi;
	}	
	if (WWTree->GenBosonNum==2) {
	  temp.SetPtEtaPhiE(WWTree->GenBosonPt[0],WWTree->GenBosonEta[0],WWTree->GenBosonPhi[0],WWTree->GenBosonE[0]);
	  temp2.SetPtEtaPhiE(WWTree->GenBosonPt[1],WWTree->GenBosonEta[1],WWTree->GenBosonPhi[1],WWTree->GenBosonE[1]);
	  gen_GravMass=(temp+temp2).M();	
	}

	deltaPhiOld=100.;
       	for (int i=0; i<WWTree->GenNuNum; i++) {
	  double deltaPhi = getDeltaPhi(WWTree->GenNuPhi[i],v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;	  
	  temp.SetPtEtaPhiE(WWTree->GenNuPt[i],WWTree->GenNuEta[i],WWTree->GenNuPhi[i],WWTree->GenNuE[i]);
	  nu_pz_gen=temp.Pz();	  
	  deltaPhiOld = deltaPhi;
	}		
      }
    
    //fill the tree
    outTree->Fill();
  }

  //--------close everything-------------
  WWTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
