// -*- C++ -*-
//
// Package:    ME0TimingAnalysis
// Class:      ME0TimingAnalysis
// 
/**\class ME0TimingAnalysis ME0TimingAnalysis.cc Analyzer/ME0TimingAnalysis/plugins/ME0TimingAnalysis.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  kaur amandeepkalsi
//         Created:  Tue, 22 Sep 2015 14:16:56 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoMuon/MuonIdentification/interface/ME0MuonSelector.h"
#include <DataFormats/MuonReco/interface/ME0Muon.h>
#include <DataFormats/MuonReco/interface/ME0MuonCollection.h>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
// class declaration
//
using namespace reco;
using namespace std;
class ME0TimingAnalysis : public edm::EDAnalyzer {
	public:
		explicit ME0TimingAnalysis(const edm::ParameterSet&);
		~ME0TimingAnalysis();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		// ----------member data ---------------------------
		int tmpindex;
		bool MatchedMuon(vector<int> me0muons, int recomuon) ;
		edm::Service<TFileService> fs;

		TH1F *hFillSignalMuontime,*hFillPUMuontime;
		TH1F *hFillSignalMuontimeErr,*hFillPUMuontimeErr;
		TH1F *hFillZMass , *hFillZGenMass, *hFillRecoEta;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ME0TimingAnalysis::ME0TimingAnalysis(const edm::ParameterSet& iConfig)

{
	//now do what ever initialization is needed

}


ME0TimingAnalysis::~ME0TimingAnalysis()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
ME0TimingAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	/*
	   vector<reco::ME0Muon>                 "me0SegmentMatching"        ""                "RECO"       
	   */
	edm::Handle <std::vector<ME0Muon> > OurMuons;
	iEvent.getByLabel <std::vector<ME0Muon> > ("me0SegmentMatching", OurMuons);

	edm::Handle <reco::GenParticleCollection> genparticles;
	iEvent.getByLabel("genParticles",genparticles);

	vector<unsigned int> indexmu;
	indexmu.clear();
	//============ gen particle collection  =============//

	for(unsigned int i = 0; i < genparticles->size();i++) {

		if(fabs(genparticles->at(i).eta()) > 2. && fabs(genparticles->at(i).eta()) < 3.) {
			if(abs(genparticles->at(i).pdgId()) == 13 && genparticles->at(i).status() == 1 && genparticles->at(i).numberOfMothers() > 0) { 
				if(fabs(genparticles->at(i).mother()->pdgId()) == 23) {
					indexmu.push_back(i); }
				else if(abs(genparticles->at(i).pdgId()) == abs(genparticles->at(i).mother()->pdgId())) {
					if(genparticles->at(i).mother()->numberOfMothers() > 0) {

						if(abs(genparticles->at(i).mother()->mother()->pdgId()) == 23) {
							indexmu.push_back(i); }  

					}
				}
			}
		}
	}
	if(indexmu.size() != 2) return;

	if(reco::deltaR(genparticles->at(indexmu.at(0)).eta(), genparticles->at(indexmu.at(0)).phi(),genparticles->at(indexmu.at(1)).eta(), genparticles->at(indexmu.at(1)).phi()) < 0.25) return;
	/*
	   for( unsigned int j = 0;  j < indexmu.size(); j++) {
	   for( unsigned int k = j+1;  k < indexmu.size(); k++) {
	   if()
	   }  
	   } 
	   */

	//double DRtmp = 0.25;
	std::vector<bool> IsMatched;
	std::vector<int> me0muons;
	me0muons.clear(); 
	tmpindex = -1; 
	for( unsigned int j = 0;  j < indexmu.size(); j++) {

		double DRtmp = 0.25;  

		for(unsigned int t = 0; t < OurMuons->size() ; t++) {

			if(int(t) == tmpindex) continue;
			double dr = reco::deltaR(genparticles->at(indexmu.at(j)).eta(), genparticles->at(indexmu.at(j)).phi(),OurMuons->at(t).eta(),OurMuons->at(t).phi());
			if(dr < DRtmp) {
				DRtmp = dr;
				tmpindex = t;
			}
		}   /////////////////////////////////////////
		me0muons.push_back(tmpindex); 
	}

	/// again filling
	if(me0muons.size() > 0) {hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time()); hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(0)).me0segment().timeErr()); } 
	if(me0muons.size() > 1) {hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(1)).me0segment().time()); hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(1)).me0segment().timeErr()); }
	if(me0muons.size() > 2) {hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(2)).me0segment().time()); hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(2)).me0segment().timeErr());}

	TLorentzVector muon1,muon2;
	for( unsigned int j = 0;  j < me0muons.size(); j++) {
		hFillRecoEta->Fill(OurMuons->at(me0muons.at(j)).eta());
		for( unsigned int ji = j+1; ji < me0muons.size(); ji++) { 

			//              hFillZMass->Fill((OurMuons->at(me0muons.at(j)).p()+ OurMuons->at(me0muons.at(ji)).p()).M());            
			//          if(reco::deltaR(OurMuons->at(me0muons.at(j)).eta(),OurMuons->at(me0muons.at(j)).phi(),OurMuons->at(me0muons.at(ji)).eta(),OurMuons->at(me0muons.at(ji)).phi()) < 0.3) continue;
			muon1.SetPtEtaPhiM(OurMuons->at(me0muons.at(j)).pt(),OurMuons->at(me0muons.at(j)).eta(),OurMuons->at(me0muons.at(j)).phi(),0);//OurMuons->at(me0muons.at(j)).energy());
			muon2.SetPtEtaPhiM(OurMuons->at(me0muons.at(ji)).pt(),OurMuons->at(me0muons.at(ji)).eta(),OurMuons->at(me0muons.at(ji)).phi(),0);//OurMuons->at(me0muons.at(ji)).energy());
			hFillZMass->Fill((muon1+muon2).M());
		}
	}
	for (unsigned int t = 0; t < OurMuons->size() ; t++){
		if(MatchedMuon(me0muons, int(t))) {
			if(indexmu.size() > 1) {
				TLorentzVector genmuon1, genmuon2;
				genmuon1.SetPtEtaPhiM(genparticles->at(indexmu.at(0)).pt(), genparticles->at(indexmu.at(0)).eta(),genparticles->at(indexmu.at(0)).phi(),genparticles->at(indexmu.at(0)).mass());

				genmuon2.SetPtEtaPhiM(genparticles->at(indexmu.at(1)).pt(), genparticles->at(indexmu.at(1)).eta(),genparticles->at(indexmu.at(1)).phi(),genparticles->at(indexmu.at(1)).mass());

				hFillZGenMass->Fill((genmuon1+genmuon2).M());


			}
		} else {
			hFillPUMuontime->Fill(OurMuons->at(t).me0segment().time()); 
			hFillPUMuontimeErr->Fill(OurMuons->at(t).me0segment().timeErr());
		}
	} 




#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
ME0TimingAnalysis::beginJob()
{

	hFillSignalMuontime = fs->make<TH1F>("hFillSignalMuontime","hFillSignalMuontime",1000,-150,-150);  
	hFillPUMuontime = fs->make<TH1F>("hFillPUMuontime","hFillPUMuontime",1000,-150,150);
	hFillSignalMuontimeErr = fs->make<TH1F>("hFillSignalMuontimeErr","hFillSignalMuontimeErr",1000,0,10);  
	hFillPUMuontimeErr = fs->make<TH1F>("hFillPUMuontimeErr","hFillPUMuontimeErr",1000,0,10);
	hFillZMass = fs->make<TH1F>("hFillZMass","hFillZMass",500,0,250); 
	hFillRecoEta = fs->make<TH1F>("hFillRecoEta","hFillZMass",500,-5,5);
	hFillZGenMass = fs->make<TH1F>("hFillZGenMass","hFillZGenMass",500,0,250);    
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
ME0TimingAnalysis::endJob() 
{
}

void
ME0TimingAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

bool ME0TimingAnalysis::MatchedMuon(vector<int> me0muons, int recomuon) {
	bool ismatched = false;
	for(unsigned int sd = 0; sd < me0muons.size(); sd++) {
		if(int(recomuon) == int(me0muons.at(sd))) {
			ismatched = true;
			break;
		}
	}
	return ismatched;
}
//define this as a plug-in
DEFINE_FWK_MODULE(ME0TimingAnalysis);
