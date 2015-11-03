#define TauTree_cxx
#include "TauTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;

static const int nTriggers =18;
TString triggerlist[nTriggers]= {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1",
	"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1",
	"HLT_IsoMu24_eta2p1_v1", 
	"HLT_IsoMu27_v1", 
	"HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1",
	"HLT_Ele32_eta2p1_WP75_Gsf_v1",
	"HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1", 
	"HLT_IsoMu24_eta2p1_v1", 
	"HLT_IsoMu27_v1",
	"HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1", 
	"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2",
	"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2", 
	"HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1", 
	"HLT_Ele32_eta2p1_WPTight_Gsf_v1",
	"HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2", 
	"HLT_IsoMu24_eta2p1_v2",
	"HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v2"
		//	"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1",
		//	"HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1",
		//	"HLT_IsoMu24_eta2p1_v1",
		//	"HLT_IsoMu27_v1",
		//	"HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1",
		//	"HLT_Ele32_eta2p1_WP75_Gsf_v1",
		//	"HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1",
		//	"HLT_IsoMu24_eta2p1_v2",
		//	"HLT_IsoMu24_eta2p1_IterTrk02_v1",
		//	"HLT_IsoMu24_eta2p1_IterTrk02_v2"
};

void TauTree::Loop(string isoTau)
{


	// intialization of gloab variables
	Nmutaudecay = NdelR = Nmujetoverlap = Ntaujetoverlap = Nmueleoverlap = Ntaueleoverlap = Nmucuts = Ntaucuts = Ngoodpair  = NextraElectron = NmtCut = NpzetA = Ndilepton = NpassDMF= Nextramu =Ntauiso= NchargeReq=0;
	Nbjets = evenT = NtrigFired = 0;
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		int GoodPair = -1;
		mutaudecay = 0;
		delR = 0;

		mutaudecay = delR = mujetoverlap = taujetoverlap = mueleoverlap = taueleoverlap = mucuts =  extraElectron = taucuts= mtCut = pzetA = dilepton = passDMF= extramu =tauiso= chargeReq=0;


		bool goodpair = false;;
		bool PassMuTrigger=false;
		evenT++;
			if (IsTriggerFired(triggerbit,FindTriggerNumber("HLT_IsoMu24_eta2p1_v1"))) PassMuTrigger=true;
//		if (IsTriggerFired(15,FindTriggerNumber("HLT_IsoMu24_eta2p1_v2"))) PassMuTrigger=true;

		if((!PassMuTrigger)) continue;
		cout<<"HELLO"<<endl;
		NtrigFired++;
		// looping over candidates

		for (unsigned int iMoth = 0; iMoth < mothers_px->size(); iMoth++){
			int dau1index = indexDau1->at(iMoth);
			int dau2index = indexDau2->at(iMoth);


			if (fabs(particleType->at(dau1index) - particleType->at(dau2index)) ==2) {

				mutaudecay++;

				TLorentzVector DauTau, DauMu ;
				DauMu.SetPxPyPzE(daughters_px->at(dau1index), daughters_py->at(dau1index), daughters_pz->at(dau1index),daughters_e->at(dau1index));
				DauTau.SetPxPyPzE(daughters_px->at(dau2index), daughters_py->at(dau2index), daughters_pz->at(dau2index),daughters_e->at(dau2index));
 bool matchMuTrigger = false;
                                if(IsTriggerFired(daughters_FilterFired->at(dau1index),FindTriggerNumber("HLT_IsoMu24_eta2p1_v1"))) matchMuTrigger= true;
                                if(IsTriggerFired(daughters_FilterFired->at(dau1index),FindTriggerNumber("HLT_IsoMu24_eta2p1_v2"))) matchMuTrigger= true;
 
                                if (!matchMuTrigger) continue;

				if(OverLap05(DauMu,DauTau,0.5)) continue;
				delR++;
				if(!(PassMuSelections(dau1index))) continue;
				mucuts++;
				if(!(PreselectionTauCuts(dau2index))) continue;
				taucuts++;
				if(OverlapWithElectrons(dau1index))continue;
				mueleoverlap++;
				if(OverlapWithJets(dau1index)) continue;
				mujetoverlap++;
				if(OverlapWithElectrons(dau2index))continue;
				taueleoverlap++;
				if(OverlapWithJets(dau2index)) continue;
				taujetoverlap++;
				if(!(daughters_decayModeFindingOldDMs->at(dau2index) == 1)) continue;                
				passDMF++;     
				if(isoTau == "LooseDB") if(!(daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(dau2index) == 1)) continue;
				if(isoTau == "MediumDB") if(!(daughters_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(dau2index) ==1)) continue;
				if(isoTau == "TightDB") if(!(daughters_byTightCombinedIsolationDeltaBetaCorr3Hits->at(dau2index) == 1)) continue;
				tauiso++;

				if(!((daughters_charge->at(dau1index) * daughters_charge->at(dau2index)) == -1)) continue;
				chargeReq++;

				if(DiLeptonVeto(dau1index)) continue; 
				dilepton++;

				mtvalue =mTCalculation(METx->at(iMoth), METy->at(iMoth), DauMu.Px(),DauMu.Py(),DauMu.Pt());
				pzeta = PZeta(dau1index, dau2index ,METx->at(iMoth), METy->at(iMoth));
				pzetavis = PZetaVis(dau1index, dau2index);
				HistoFiller(hmt_woCut, mtvalue);
				if(mtvalue > 40.) continue;

				mtCut++;
				HistoFiller(hpzeta_woCut,(pzeta - 1.85*pzetavis));

				if(!(pzeta - 1.85*pzetavis > -20)) continue;
				pzetA++;

				if(ExtraMuon(dau1index)) continue;
				extramu++;

				GoodPair = iMoth;
				goodpair = true;
				break;
			}
		} // daughter


		if(mutaudecay > 0) Nmutaudecay++;
		if(delR > 0) NdelR++;
		if(mucuts >0 ) Nmucuts++;
		if(taucuts >0 )  Ntaucuts++;
		if(mueleoverlap > 0) Nmueleoverlap++;
		if(mujetoverlap > 0) Nmujetoverlap++;
		if(taueleoverlap > 0) Ntaueleoverlap++;
		if(taujetoverlap > 0) Ntaujetoverlap++;
		if(passDMF > 0) NpassDMF++;
		if(tauiso > 0) Ntauiso++;
		if(chargeReq > 0) NchargeReq++;
		if(dilepton > 0) Ndilepton++;
		if(mtCut > 0) NmtCut++;
		if(pzetA > 0) NpzetA++;
		if(extramu > 0) Nextramu++;

		if(!(goodpair)) continue;
		Ngoodpair++;
		HistoFiller(h_nbjets, NumberBJets());
		HistoFiller(h_nextraElectron, NumberExtraElectron());



		if(BJets()) continue;
		Nbjets++;
		if(ExtraElectron()) continue;
		NextraElectron++;

		int tindex = indexDau1->at(GoodPair);
		int sindex=indexDau2->at(GoodPair);   


		TLorentzVector FirstObj, SecondObj, ThirdObj ;
		FirstObj.SetPxPyPzE(daughters_px->at(tindex), daughters_py->at(tindex),daughters_pz->at(tindex),daughters_e->at(tindex));
		SecondObj.SetPxPyPzE(daughters_px->at(sindex), daughters_py->at(sindex),daughters_pz->at(sindex),daughters_e->at(sindex));

                 ThirdObj=FirstObj+SecondObj;

                HistoFiller(hmutaumass, ThirdObj.M()); 
		// histograms
		HistoFiller(hmu_px,daughters_px->at(tindex));
		HistoFiller(hmu_py,daughters_py->at(tindex));
		HistoFiller(hmu_pz,daughters_pz->at(tindex));
		HistoFiller(hmu_en,daughters_e->at(tindex));

		HistoFiller(hmu_pt,FirstObj.Pt());
		HistoFiller(hmu_eta,FirstObj.Eta());
		HistoFiller(hmu_phi,FirstObj.Phi());
		HistoFiller(hmu_charge,daughters_charge->at(tindex));

		HistoFiller(htau_px,daughters_px->at(sindex));
		HistoFiller(htau_py,daughters_py->at(sindex));
		HistoFiller(htau_pz,daughters_pz->at(sindex));
		HistoFiller(htau_en,daughters_e->at(sindex));

		HistoFiller(htau_pt,SecondObj.Pt());
		HistoFiller(htau_eta,SecondObj.Eta());
		HistoFiller(htau_phi,SecondObj.Phi());
		HistoFiller(htau_charge,daughters_charge->at(sindex));

		HistoFiller(htau_oldDM,daughters_decayModeFindingOldDMs->at(sindex));
		HistoFiller(htau_newDM,daughters_decayModeFindingNewDMs->at(sindex));

		HistoFiller(htau_Loose3Hits,daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(sindex));

		HistoFiller(htau_Medium3Hits,daughters_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(sindex));

		HistoFiller(htau_Tight3Hits,daughters_byTightCombinedIsolationDeltaBetaCorr3Hits->at(sindex));

		HistoFiller(htau_RawIso,daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(sindex));

		HistoFiller(htau_oldDMwoLTraw, daughters_byIsolationMVA3oldDMwoLTraw->at(sindex));
		HistoFiller(htau_newDMwoLTraw, daughters_byIsolationMVA3newDMwoLTraw->at(sindex));
		HistoFiller(htau_oldDMwLT, daughters_byIsolationMVA3oldDMwLTraw->at(sindex));
		HistoFiller(htau_newDMwLT, daughters_byIsolationMVA3newDMwLTraw->at(sindex));

		HistoFiller(h_chargedIso, daughters_chargedIsoPtSum->at(sindex));
		HistoFiller(h_neutral, daughters_neutralIsoPtSum->at(sindex));

		HistoFiller(h_numIsoCone, daughters_numParticlesIsoCone->at(sindex));
		HistoFiller(h_numIsoPhotons, daughters_numPhotonsIsoCone->at(sindex));
		HistoFiller(h_numIsoNeutral, daughters_numNeutralHadronsIsoCone->at(sindex));
		HistoFiller(h_numCharged, daughters_numChargedParticlesIsoCone->at(sindex));
		HistoFiller(h_numPhotons, daughters_numPhotonsSignalCone->at(sindex));
		HistoFiller(h_numNeuHadr, daughters_numNeutralHadronsSignalCone->at(sindex));
		HistoFiller(h_numCharSignal, daughters_numChargedParticlesSignalCone->at(sindex));
		HistoFiller(h_puCorr , daughters_puCorrPtSum->at(sindex));
		HistoFiller(h_muonloose, daughters_againstMuonLoose3->at(sindex));
		HistoFiller(h_muontight, daughters_againstMuonTight3->at(sindex));
		HistoFiller(h_evloose, daughters_againstElectronVLooseMVA5->at(sindex));
		HistoFiller(h_eloose, daughters_againstElectronLooseMVA5->at(sindex));
		HistoFiller(h_emedium, daughters_againstElectronMediumMVA5->at(sindex));
		HistoFiller(h_etight, daughters_againstElectronTightMVA5->at(sindex));
		HistoFiller(h_evtight, daughters_againstElectronVTightMVA5->at(sindex));
		HistoFiller(hZmt,mtvalue);
		HistoFiller(hpzetacut, (pzeta - 1.85*pzetavis));

		// new histos


		//		int tau_id = -1; 
		//		if(!(GenSelection(tindex, sindex,tau_id))) continue;
		//		if(tau_id == -1) continue;		
		//		cout<<"PZETA"<<pzeta<<":pzetavis:"<<pzetavis<<endl; 
		//		cout<<"tau_id"<<tau_id<<endl;
		int NoCharged = daughters_numChargedParticlesSignalCone->at(sindex) + daughters_numChargedParticlesIsoCone->at(sindex);
		HistoFiller(NTracks,NoCharged);
		HistoFiller(NTracksSignal, daughters_numChargedParticlesSignalCone->at(sindex));
		HistoFiller(NTracksIsolation, daughters_numChargedParticlesIsoCone->at(sindex));




	}


	hEventCounter->SetBinContent(1, evenT);
	hEventCounter->SetBinContent(2, NtrigFired);
	hEventCounter->SetBinContent(3, Nmutaudecay);
	hEventCounter->SetBinContent(4, NdelR);
	hEventCounter->SetBinContent(5, Nmucuts);
	hEventCounter->SetBinContent(6, Ntaucuts);
	hEventCounter->SetBinContent(7, Nmueleoverlap);
	hEventCounter->SetBinContent(8, Nmujetoverlap);
	hEventCounter->SetBinContent(9, Ntaueleoverlap);
	hEventCounter->SetBinContent(10, Ntaujetoverlap);
	hEventCounter->SetBinContent(11, NpassDMF);
	hEventCounter->SetBinContent(12, Ntauiso);
	hEventCounter->SetBinContent(13, NchargeReq);
	hEventCounter->SetBinContent(14, Ndilepton);
	hEventCounter->SetBinContent(15, NmtCut);
	hEventCounter->SetBinContent(16, NpzetA);
	hEventCounter->SetBinContent(17,Nextramu);
	hEventCounter->SetBinContent(18, Ngoodpair);
	hEventCounter->SetBinContent(19, Nbjets);
	hEventCounter->SetBinContent(20, NextraElectron);

	output_file->cd();
	output_file->Write();
	output_file->Close();

}


bool TauTree::PreselectionTauCuts(int tauindex) {

	TLorentzVector taulep;
	taulep.SetPxPyPzE(daughters_px->at(tauindex) , daughters_py->at(tauindex),daughters_pz->at(tauindex),daughters_e->at(tauindex));

	if((taulep.Pt() <= 20)) return false;
	if(fabs(taulep.Eta()) >= 2.3) return false;
	if(daughters_trackRefPt->at(tauindex) < 5.) return false;
	//// dz of leading track of tau                             
	if(fabs(dz->at(tauindex)) >= 0.2) return false;

	// for synchronization with cecile 30 sept 
	if(daughters_againstElectronLooseMVA5->at(tauindex) != 1) return false;
	if(daughters_againstMuonTight3->at(tauindex) != 1) return false;
	//  for synchronization with cecile 30 sept                                
	// remove overlap with global muon                                         
	if(muonTauOverlap(tauindex)) return false;

	return true;
}

bool TauTree::PassMuSelections(int muindex){

	TLorentzVector mu;
	mu.SetPxPyPzE(daughters_px->at(muindex) , daughters_py->at(muindex),daughters_pz->at(muindex),daughters_e->at(muindex));

	if(mu.Pt() <= 25.) return false;

	if((fabs(mu.Eta()) >= 2.1)) return false;

	//if(!(IsTriggerFired(daughters_FilterFired->at(muindex),FindTriggerNumber("HLT_IsoMu24_eta2p1_v1")))) return false;

	if(!((daughters_typeOfMuon->at(muindex) >> 1 ) & 1)) return false; 
	// medium wp                                                                                                   
	if(!((daughters_muonID->at(muindex) >> 2) & 1)) return false;
	if((fabs(dxy->at(muindex)) >= 0.045)) return false;
	if((fabs(dz->at(muindex)) >= 0.2)) return false;
	if(!(combreliso->at(muindex) < 0.1)) return false;                              	return true;

}

bool TauTree::OverLap05(TLorentzVector l1 , TLorentzVector l2, float conesize) {
	if(dR(l1.Eta(), l1.Phi(), l2.Eta(), l2.Phi()) <= conesize) return true;
	else return false;
}

/// delta Phi                                                                                                          
//                                                                                                                     
float TauTree::deltaPhi( float a, float b) {
	float result = a-b;
	while (result > M_PI) result -= 2* M_PI;
	while (result <= -M_PI) result += 2* M_PI;
	return (fabs(result));

}

float TauTree::dR(float l1eta, float l1phi, float l2eta, float l2phi ) {
	float deta = l1eta - l2eta;
	float dphi = deltaPhi(l1phi,l2phi);
	return sqrt(deta*deta + dphi*dphi);
}

bool TauTree::muselection(int index1, float ptcut, float etacut, float isolation) {

	TLorentzVector lep;
	lep.SetPxPyPzE(daughters_px->at(index1) , daughters_py->at(index1),daughters_pz->at(index1),daughters_e->at(index1));
	// global muon                                                                                                 
	if(lep.Pt() > ptcut && fabs(lep.Eta()) < etacut && ((daughters_typeOfMuon->at(index1) >> 1 ) & 1) && combreliso->at(index1) < isolation && ((daughters_muonID->at(index1) >> 3) & 1)) {
		return true;
	}else return false;
}
// Third lepton vet0    

bool TauTree::tauselection(int index1, float ptcut, float etacut, float isolation, string isoWP,string muWP, string eleWP, string dmf) {
	TLorentzVector tau;
	tau.SetPxPyPzE(daughters_px->at(index1) , daughters_py->at(index1),daughters_pz->at(index1),daughters_e->at(index1));
	bool kinematic = false;
	bool dmfbool = false;
	bool elebool = false;
	bool mubool = false;
	bool isobool = false;
	if(tau.Pt() > ptcut && fabs(tau.Eta()) < etacut) kinematic=true;
	if(dmf == "Old"){
		if((daughters_decayModeFindingOldDMs->at(index1)) == 1) dmfbool = true;
	} else if(dmf == "New") {
		if((daughters_decayModeFindingNewDMs->at(index1)) == 1) dmfbool = true;
	}else {cout<<"give correct dmf wp string"<<endl;}


	if(eleWP == "vloose") { if(daughters_againstElectronVLooseMVA5->at(index1) == 1) elebool = true;}
	else if(eleWP == "loose") { if(daughters_againstElectronLooseMVA5->at(index1) == 1) elebool = true;}
	else if(eleWP == "medium") { if(daughters_againstElectronMediumMVA5->at(index1) == 1) elebool = true;}
	else if(eleWP == "tight") { if(daughters_againstElectronTightMVA5->at(index1) == 1) elebool = true;}
	else if(eleWP == "vtight") { if(daughters_againstElectronVTightMVA5->at(index1) == 1) elebool = true;}
	else {cout<<"give correct elec wp string"<<endl;}


	if(muWP == "loose") { if(daughters_againstMuonLoose3->at(index1) == 1) mubool = true;}
	else if (muWP == "tight") { if(daughters_againstMuonTight3->at(index1) == 1) mubool = true;}
	else {cout<<"give correct muon  wp string"<<endl;}

	if(isoWP == "isoDBloose") {if(daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(index1) == 1) isobool = true;}

	else if(isoWP == "isoDBmedium") {if(daughters_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(index1) ==1) isobool = true;}

	else if(isoWP == "isoDBtight") {if(daughters_byTightCombinedIsolationDeltaBetaCorr3Hits->at(index1) == 1) isobool = true;}

	else {cout<<"give correct iso wp string"<<endl;}
	if(kinematic&&dmfbool&&elebool&&isobool) return true; else return false;
}


bool TauTree::muonTauOverlap( int index1 ) {
	bool muon1;
	muon1= false;
	TLorentzVector lep, DauTau;

	DauTau.SetPxPyPzE(daughters_px->at(index1),daughters_py->at(index1),daughters_pz->at(index1),daughters_e->at(index1));
	for (unsigned int i = 0 ; i < daughters_px->size(); i++) {
		lep.SetPxPyPzE(daughters_px->at(i),daughters_py->at(i),daughters_pz->at(i),daughters_e->at(i));
		if(particleType->at(i) == 0 && lep.Pt() > 5 && fabs(lep.Eta()) < 2.4) {
			// type of muon
			if((daughters_typeOfMuon->at(i) >> 0) & 1) {
				if(OverLap05(lep,DauTau,0.5)) {
					muon1 = true;
					break;
				}

			}
		}
	}
	return muon1;
}



int TauTree::NumberExtraElectron(){
	TLorentzVector ele;
	int nelectron;
	nelectron = 0;
	for (unsigned int j = 0 ; j < daughters_px->size(); j++) {
		ele.SetPxPyPzE(daughters_px->at(j),daughters_py->at(j),daughters_pz->at(j),daughters_e->at(j));
		if(particleType->at(j) == 1 && ele.Pt() > 15 && fabs(ele.Eta()) < 2.4 && fabs(dxy->at(j)) < 0.045 && fabs(dz->at(j)) < 0.2 && combreliso->at(j) < 0.3 && (daughters_passConversionVeto->at(j)) && (daughters_iseleWP90->at(j)) ){ 
			nelectron++;
		} 
	}
	return nelectron;
}


bool TauTree::ExtraElectron(){
	TLorentzVector ele;
	bool iselectron;
	iselectron = false;
	for (unsigned int j = 0 ; j < daughters_px->size(); j++) {
		ele.SetPxPyPzE(daughters_px->at(j),daughters_py->at(j),daughters_pz->at(j),daughters_e->at(j));
		if(particleType->at(j) == 1 && ele.Pt() > 15 && fabs(ele.Eta()) < 2.4 && fabs(dxy->at(j)) < 0.045 && fabs(dz->at(j)) < 0.2 && combreliso->at(j) < 0.3 && (daughters_passConversionVeto->at(j)) && (daughters_iseleWP90->at(j)) ){ 
			iselectron = true;
			break;
		} 
	}
	return iselectron;
}


int TauTree::ExtraMuonNumber(int dau1index){
	TLorentzVector mu, mu1;
	mu1.SetPxPyPzE(daughters_px->at(dau1index),daughters_py->at(dau1index),daughters_pz->at(dau1index),daughters_e->at(dau1index));
	int nmuon;
	nmuon = 0;
	for (unsigned int j = 0 ; j < daughters_px->size(); j++) {
		mu.SetPxPyPzE(daughters_px->at(j),daughters_py->at(j),daughters_pz->at(j),daughters_e->at(j));

		if( (!OverLap05(mu,mu1,0.5)) && particleType->at(j) == 0 && mu.Pt() > 10 && fabs(mu.Eta()) < 2.4 && fabs(dxy->at(j)) < 0.045 && fabs(dz->at(j)) < 0.2 && combreliso->at(j) < 0.3 ) { 
			nmuon++;

		}
	}
	return nmuon;
}



bool TauTree::ExtraMuon(int dau1index){
	TLorentzVector mu, mu1;
	mu1.SetPxPyPzE(daughters_px->at(dau1index),daughters_py->at(dau1index),daughters_pz->at(dau1index),daughters_e->at(dau1index));
	bool ismuon;
	ismuon = false;
	for (unsigned int j = 0 ; j < daughters_px->size(); j++) {
		mu.SetPxPyPzE(daughters_px->at(j),daughters_py->at(j),daughters_pz->at(j),daughters_e->at(j));

		if( (!OverLap05(mu,mu1,0.5)) && particleType->at(j) == 0 && mu.Pt() > 10 && fabs(mu.Eta()) < 2.4 && fabs(dxy->at(j)) < 0.045 && fabs(dz->at(j)) < 0.2 && combreliso->at(j) < 0.3 ) { 
			ismuon =  true;

		}
	}
	return ismuon;
}

bool TauTree::BJets() {
	TLorentzVector bjet;
	bool isbjet;
	isbjet=false;
	for(unsigned int i = 0 ; i < jets_px->size(); i++) {
		bjet.SetPxPyPzE(jets_px->at(i), jets_py->at(i), jets_pz->at(i), jets_e->at(i));
		if(bjet.Pt() > 20 && (fabs(bjet.Eta()) < 2.4) && (PFjetID->at(i) == 2) && (bCSVscore->at(i) > 0.814) ) {
			isbjet = true;
		}
	}
	return isbjet;
}


int TauTree::NumberBJets() {
	TLorentzVector bjet;
	int bjetn;
	bjetn = 0;
	for(unsigned int i = 0 ; i < jets_px->size(); i++) {
		bjet.SetPxPyPzE(jets_px->at(i), jets_py->at(i), jets_pz->at(i), jets_e->at(i));
		if(bjet.Pt() > 20 && (fabs(bjet.Eta()) < 2.4) && (PFjetID->at(i) == 2) && (bCSVscore->at(i) > 0.814) ) {
			bjetn++;
		}
	}
	return bjetn;
}



bool TauTree::IsTriggerFired(int triggerbit, int triggernumber){ 
	if(triggernumber>=0 && triggernumber<nTriggers) { return triggerbit;}// & (1 << triggernumber);}
	else return false;
} 


int TauTree::FindTriggerNumber(TString triggername){              
	for(unsigned int it=0;it<nTriggers;it++){
		if(triggerlist[it].CompareTo(triggername.Data())==0) return it;
		else {              
			TString newName=triggername.Data();
			if(triggerlist[it].CompareTo(newName.Data())==0)return it;
		}                   
	}                         
	return -1;                
}

bool TauTree::OverlapWithMuons( int dau1index){
	TLorentzVector mu,mu1;
	mu1.SetPxPyPzE(daughters_px->at(dau1index),daughters_py->at(dau1index),daughters_pz->at(dau1index),daughters_e->at(dau1index));
	bool ismuon;
	ismuon = false;
	for (unsigned int j = 0 ; j < daughters_px->size(); j++) {
		mu.SetPxPyPzE(daughters_px->at(j),daughters_py->at(j),daughters_pz->at(j),daughters_e->at(j));

		if( (!(OverLap05(mu1,mu,0.5))) && particleType->at(j) == 0 && mu.Pt() > 5. && fabs(mu.Eta()) < 2.4 && fabs(dxy->at(j)) < 0.045 && fabs(dz->at(j)) < 0.2 && combreliso->at(j) < 0.3 ) { 
			ismuon =  true;

		}
	}
	return ismuon;
}      

bool TauTree::OverlapWithElectrons(int dau1index){
	TLorentzVector ele, ele1;
	ele1.SetPxPyPzE(daughters_px->at(dau1index),daughters_py->at(dau1index),daughters_pz->at(dau1index),daughters_e->at(dau1index));
	bool iselectron;
	iselectron = false;
	for (unsigned int j = 0 ; j < daughters_px->size(); j++) {
		ele.SetPxPyPzE(daughters_px->at(j),daughters_py->at(j),daughters_pz->at(j),daughters_e->at(j));
		if((!(OverLap05(ele1,ele,0.5))) && particleType->at(j) == 1 && ele.Pt() > 15 && fabs(ele.Eta()) < 2.4 && fabs(dxy->at(j)) < 0.045 && fabs(dz->at(j)) < 0.2 && combreliso->at(j) < 0.3 && (daughters_passConversionVeto->at(j)) && (daughters_iseleWP90->at(j)) ){ 
			iselectron = true;
		} 
	}
	return iselectron;
}  

bool TauTree::OverlapWithTaus(int dau1index){
	TLorentzVector tau, tau1;
	bool istau;
	istau=false;
	tau1.SetPxPyPzE(daughters_px->at(dau1index),daughters_py->at(dau1index),daughters_pz->at(dau1index),daughters_e->at(dau1index));
	for (unsigned int j = 0 ; j < daughters_px->size(); j++) {
		tau.SetPxPyPzE(daughters_px->at(j),daughters_py->at(j),daughters_pz->at(j),daughters_e->at(j));
		if((!(OverLap05(tau1,tau,0.5))) && (tau.Pt() > 20) && (fabs(tau.Eta()) < 2.3) &&  ((daughters_decayModeFindingNewDMs->at(dau1index)) == 1) && (daughters_againstElectronTightMVA5->at(dau1index) == 1) && (daughters_againstMuonTight3->at(dau1index) == 1) && (daughters_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(dau1index) == 1) ) {
			istau = true;
		}
	}
	return istau;
}       

bool TauTree::OverlapWithJets( int dau1index){
	TLorentzVector jet, object1;
	object1.SetPxPyPzE(daughters_px->at(dau1index),daughters_py->at(dau1index),daughters_pz->at(dau1index),daughters_e->at(dau1index));

	bool isjet;
	isjet=false;
	for(unsigned int i = 0 ; i < jets_px->size(); i++) {
		jet.SetPxPyPzE(jets_px->at(i), jets_py->at(i), jets_pz->at(i), jets_e->at(i));
		if(jet.Pt() > 30. && (fabs(jet.Eta()) < 5.0) && (PFjetID->at(i) == 2)  ) {
			isjet = true;
		}
	}
	return isjet;

}       



float TauTree::mTCalculation(float metx, float mety, float mupx, float mupy, float mupt){
	float mt = -1;
	float pX = mupx+metx;
	float pY = mupy+mety;
	float et = mupt + TMath::Sqrt(metx*metx + mety*mety);
	mt = TMath::Sqrt(et*et-(pX*pX + pY*pY));
	return mt;

}


float TauTree::PZetaVis( int muindex, int tauindex){
	float pzetavis;
	pzetavis = 999;
	TLorentzVector tau, mu;  
	tau.SetPxPyPzE(daughters_px->at(tauindex),daughters_py->at(tauindex),daughters_pz->at(tauindex),daughters_e->at(tauindex));
	mu.SetPxPyPzE(daughters_px->at(muindex),daughters_py->at(muindex),daughters_pz->at(muindex),daughters_e->at(muindex));
	float zetax = TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()) ;
	float zetay = TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()) ;
	float zetaR = TMath::Sqrt(pow(zetax,2) + pow(zetay,2));
	zetax = zetax/zetaR;
	zetay = zetay/zetaR;

	float visPx = mu.Px() + tau.Px() ;
	float visPy = mu.Py() + tau.Py() ;

	pzetavis = visPx*zetax + visPy*zetay;
	return pzetavis;

}

float TauTree::PZeta(int muindex, int tauindex , float metpx, float metpy){
	float pzeta;
	pzeta = 999;
	TLorentzVector tau, mu;                                         
	tau.SetPxPyPzE(daughters_px->at(tauindex),daughters_py->at(tauindex),daughters_pz->at(tauindex),daughters_e->at(tauindex));
	mu.SetPxPyPzE(daughters_px->at(muindex),daughters_py->at(muindex),daughters_pz->at(muindex),daughters_e->at(muindex));
	float zetax = TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()) ;
	float zetay = TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()) ;
	float zetaR = TMath::Sqrt(pow(zetax,2) + pow(zetay,2)); 
	zetax = zetax/zetaR;                                    
	zetay = zetay/zetaR;         

	float vPx = mu.Px() + tau.Px()+metpx ;
	float vPy = mu.Py() + tau.Py()+metpy ;

	pzeta = vPx*zetax + vPy*zetay;
	return pzeta;


}

bool TauTree::DiLeptonVeto(int muindex){
	TLorentzVector mu;
	bool dilepton=false;
	mu.SetPxPyPzE(daughters_px->at(muindex),daughters_py->at(muindex),daughters_pz->at(muindex),daughters_e->at(muindex));
	for (unsigned int j = 0 ; j < daughters_px->size(); j++) {
		TLorentzVector tmp ;
		tmp.SetPxPyPzE(daughters_px->at(j),daughters_py->at(j),daughters_pz->at(j),daughters_e->at(j));
		if(OverLap05(mu ,tmp, 0.5)) continue;
		if(!(daughters_charge->at(muindex) * daughters_charge->at(j) < 0)) continue;
		if(!(particleType->at(j) == 0)) continue;
		if(!(combreliso->at(j) < 0.3)) continue;
		if(!(fabs(dxy->at(j)) < 0.045)) continue;
		if(!(fabs(dz->at(j)) < 0.2)) continue;
		if(!(fabs(tmp.Eta()) < 2.4)) continue;
		if(!(tmp.Pt() > 15)) continue;
		if(((daughters_typeOfMuon->at(j) >> 0) & 1) && ((daughters_typeOfMuon->at(j) >> 1) & 1) && ((daughters_typeOfMuon->at(j) >> 2) & 1)) {
			dilepton= true;
		}  
	}
	return dilepton;
}
/*
   bool TauTree::GenSelection(int muindex, int tauindex, int &tau_id){
   TLorentzVector tau2;
   tau_id = -1;
   tau2.SetPxPyPzE(daughters_px->at(tauindex),daughters_py->at(tauindex),daughters_pz->at(tauindex),daughters_e->at(tauindex));
   bool isZmt=false;            
   double newdr = 0.5;
   int indice_gen=0;
   for (unsigned int iGen = 2; iGen < genpart_pdg->size(); iGen++){
   TLorentzVector g;
   g.SetPxPyPzE(genpart_px->at(iGen),genpart_py->at(iGen),genpart_pz->at(iGen),genpart_e->at(iGen));     
   if ((fabs(genpart_pdg->at(iGen))==11 or fabs(genpart_pdg->at(iGen))==13 or fabs(genpart_pdg->at(iGen))==15) && g.DeltaR(tau2)<0.5 && g.DeltaR(tau2)<newdr){
   newdr=g.DeltaR(tau2);
   indice_gen=iGen;
   }
   }
   if (indice_gen==0){
   for (unsigned int iGen = 2; iGen < genpart_pdg->size(); iGen++){
   TLorentzVector g;
   g.SetPxPyPzE(genpart_px->at(iGen),genpart_py->at(iGen),genpart_pz->at(iGen),genpart_e->at(iGen));
   if ((fabs(genpart_pdg->at(iGen))<30) && g.DeltaR(tau2)<0.5 && g.DeltaR(tau2)<newdr){
   newdr=g.DeltaR(tau2);
   indice_gen=iGen;
   }
   }
   }
   tau_id=genpart_pdg->at(indice_gen);

   for (unsigned int iGen = 0; iGen < genpart_px->size(); iGen++){
   if (genpart_pdg->at(iGen)==23){
   if (genpart_HZDecayMode->at(iGen)==0)
   isZmt=true;
   } 
   }
   return isZmt;

   }
   */
void TauTree::HistoFiller(TH1F *histo, double value){

	histo->Fill(value);

}


void TauTree::HistoDec(const char *fname){
	output_file = TFile::Open(fname,"RECREATE");


	string counE[] = {"Total","TriggerPassed","MuTauDecay","DeltaR","Mu Cuts","Tau Cuts", "MuEle overlap","MuJet overlap","TauEle overlap","TauJet overlap","tau DMF","tau Isolation", "Charge Req.", "NDilepton","mt cut", "pzeta cut", "extra muon", "GoodPair","BJet veto","extra electron"};

	int hhh = sizeof(counE)/sizeof(string);
	hEventCounter =  new TH1F("hEventCounter","hEventCounter",hhh,0,hhh);
	for (int i = 0; i <hhh ; i++){

		hEventCounter->GetXaxis()->SetBinLabel(i+1,counE[i].c_str());
	}


	hmutaumass = new TH1F("hmutaumass","hmutaumass",1000,0,1000);
	hmu_pt = new TH1F("hmu_pt","hmu_pt",1000,0,1000);
	hmu_px = new TH1F("hmu_px","hmu_px",1000,-500,500);
	hmu_py = new TH1F("hmu_py","hmu_py",1000,-500,500);
	hmu_pz = new TH1F("hmu_pz","hmu_pz",1000,-500,500);
	hmu_en = new TH1F("hmu_en","hmu_en",1000,0,1000);
	hmu_eta = new TH1F("hmu_eta","hmu_eta",100,-5,5);
	hmu_phi = new TH1F("hmu_phi","hmu_phi",100,-5,5);
	hmu_charge = new TH1F("hmu_charge","hmu_charge",10,-5,5);

	/////////////
	htau_pt = new TH1F("htau_pt","htau_pt",1000,0,1000);
	htau_px = new TH1F("htau_px","htau_px",1000,-500,500);
	htau_py = new TH1F("htau_py","htau_py",1000,-500,500);
	htau_pz = new TH1F("htau_pz","htau_pz",1000,-500,500);
	htau_en = new TH1F("htau_en","htau_en",1000,0,1000);
	htau_eta = new TH1F("htau_eta","htau_eta",100,-5,5);
	htau_phi = new TH1F("htau_phi","htau_phi",100,-5,5);
	htau_charge = new TH1F("htau_charge","htau_charge",10,-5,5);
	htau_oldDM = new TH1F("htau_oldDM","htau_oldDM",10,-5,5);
	htau_newDM = new TH1F("htau_newDM","htau_newDM",10,-5,5);

	hZmt  = new TH1F("hZmt","hZmt",1000,0,1000);

	hpzetacut = new TH1F("hpzetacut","hpzetacut",1000,-500,500);
	hpzeta_woCut  = new TH1F("hpzeta_woCut","hpzeta_woCut",1000,-500,500);
	hmt_woCut = new TH1F("hmt_woCut","hmt_woCut",1000,0,1000);
	h_dilepton  = new TH1F("h_dilepton","h_dilepton",2,0,2);
	h_chargedIso = new TH1F("h_chargedIso","h_chargedIso",200,0,20);
	h_neutral = new TH1F("h_neutral","h_neutral",200,0,20);
	h_numIsoCone = new TH1F("h_numIsoCone","h_numIsoCone",20,0,20);
	h_numIsoPhotons  = new TH1F("h_numIsoPhotons","h_numIsoPhotons",20,0,20);
	h_numIsoNeutral = new TH1F("h_numIsoNeutral","h_numIsoNeutral",20,0,20);
	h_numCharged  = new TH1F("h_numCharged","h_numCharged",20,0,20);
	h_numPhotons  = new TH1F("h_numNeutral","h_numNeutral",20,0,20);
	h_numNeuHadr  = new TH1F("h_numNeuHadr","h_numNeuHadr",20,0,20);
	h_numCharSignal = new TH1F("h_numCharSignal","h_numCharSignal",20,0,20);
	h_puCorr = new TH1F("h_puCorr","h_puCorr",200,0,20);
	h_muonloose = new TH1F("h_muonloose","h_muonloose",2,0,2);
	h_muontight  = new TH1F("h_muontight","h_muontight",2,0,2);
	h_evloose  = new TH1F("h_evloose","h_evloose",2,0,2);
	h_eloose  = new TH1F("h_eloose","h_eloose",2,0,2);
	h_emedium = new TH1F("h_emedium","h_emedium",2,0,2);
	h_etight = new TH1F("h_etight","h_etight",2,0,2);
	h_evtight = new TH1F("h_evtight","h_evtight",2,0,2);


	htau_Loose3Hits = new TH1F("htau_Loose3Hits","htau_Loose3Hits",2,0,2);

	htau_Medium3Hits = new TH1F("htau_Medium3Hits","htau_Medium3Hits",2,0,2);

	htau_Tight3Hits  = new TH1F("htau_Tight3Hits","htau_Tight3Hits",2,0,2);

	htau_RawIso  = new TH1F("htau_RawIso","htau_RawIso",100,0,40);

	htau_oldDMwoLTraw  = new TH1F("htau_oldDMwoLTraw","htau_oldDMwoLTraw",100,0,40);
	htau_newDMwoLTraw  = new TH1F("htau_newDMwoLTraw","htau_newDMwoLTraw",100,0,40);
	htau_oldDMwLT  = new TH1F("htau_oldDMwLT","htau_oldDMwLT",100,0,40);
	htau_newDMwLT = new TH1F("htau_newDMwLT","htau_newDMwLT",100,0,40);


	h_nbjets = new TH1F("h_nbjets","h_nbjets",10,0,10);
	h_nextraElectron = new TH1F("h_nextraElectron","h_nextraElectron",10,0,10);
	h_nextraMuon  = new TH1F("h_nextraMuon","h_nextraMuon",10,0,10);

	NTracks  = new TH1F("NTracks","NTracks",25,0,25);
	NTracksSignal  = new TH1F("NTracksSignal","NTracksSignal",25,0,25);
	NTracksIsolation  = new TH1F("NTracksIsolation","NTracksIsolation",25,0,25);


}
