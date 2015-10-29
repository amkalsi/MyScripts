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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/ParametersDefinerForTP.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"
//Associator for chi2: Including header files
//
// class declaration
//
using namespace reco;
using namespace std;
using namespace edm;
class ME0TimingAnalysis : public edm::EDAnalyzer {
	public:
		explicit ME0TimingAnalysis(const edm::ParameterSet&);
		~ME0TimingAnalysis();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		//		virtual void beginJob() override;
		void beginRun(edm::Run const&, edm::EventSetup const&) override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		void endRun(edm::Run const&, edm::EventSetup const&) override;
		//		virtual void endJob() override;
		// ----------member data ---------------------------
		int tmpindex;
		bool MatchedMuon(vector<int> me0muons, int recomuon) ;
		edm::Service<TFileService> fs;
		bool GenParticleSelection(const GenParticle* tp);
		TH1F *hFillSignalMuontime,*hFillPUMuontime;
		TH1F *hFillSignalMuontimeErr,*hFillPUMuontimeErr;
		TH1F *hFillZMass , *hFillZGenMass, *hFillRecoEta;
		TH1F *SignalMuonTime, *BGMuonTime;
		TH1F *hFillZMassInWindow,*hFillZMassOutWindow;
		bool UseAssociators;
		bool RejectEndcapMuons;
		TH1F *hFillGenMuonPtDen,*hFillGenMuonEtaDen,*hFillGenMuonPtNum,*hFillGenMuonEtaNum;
		TH1F *hFillDenRecPt, *hFillDenRecEta, *hFillNumRecPt, *hFillNumRecEta;


		std::vector<std::string> associators;
		std::vector<edm::InputTag> label;

		edm::ESHandle<MagneticField> theMF;
		std::vector<const TrackAssociatorBase*> associator;
		const TrackAssociatorByChi2* associatorByChi2;

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
	UseAssociators = iConfig.getParameter< bool >("UseAssociators");
	associators = iConfig.getParameter< std::vector<std::string> >("associators");
	label = iConfig.getParameter< std::vector<edm::InputTag> >("label");


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
	using namespace reco;

	Handle <std::vector<RecoChargedCandidate> > OurCandidates;
	iEvent.getByLabel <std::vector<RecoChargedCandidate> > ("me0MuonConverter", OurCandidates);

	edm::Handle <std::vector<ME0Muon> > OurMuons;
	iEvent.getByLabel <std::vector<ME0Muon> > ("me0SegmentMatching", OurMuons);

	edm::Handle <reco::GenParticleCollection> genparticles;
	iEvent.getByLabel("genParticles",genparticles);
	const GenParticleCollection genParticlesForChi2 = *(genparticles.product());

	edm::Handle<reco::TrackCollection> generalTracks;
	iEvent.getByLabel("generalTracks",generalTracks);

	Handle<ME0SegmentCollection> OurSegments;
	iEvent.getByLabel("me0Segments","",OurSegments);

	edm::ESHandle<ME0Geometry> me0Geom;
	iSetup.get<MuonGeometryRecord>().get(me0Geom);

	ESHandle<MagneticField> bField;
	iSetup.get<IdealMagneticFieldRecord>().get(bField);
	ESHandle<Propagator> shProp;
	iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", shProp);


	//Track Association by Chi2 or hits:

	//int w=0;

	for (unsigned int ww=0;ww<associators.size();ww++){
		associatorByChi2 = dynamic_cast<const TrackAssociatorByChi2*>(associator[ww]);

		if (associatorByChi2==0) continue;

		for (unsigned int www=0;www<label.size();www++){
			//
			//get collections from the event
			//
			edm::Handle<View<Track> >  trackCollection;

			reco::RecoToGenCollection recSimColl;
			reco::GenToRecoCollection simRecColl;
			unsigned int trackCollectionSize = 0;

			if(!iEvent.getByLabel(label[www], trackCollection)) {

				recSimColl.post_insert();
				simRecColl.post_insert();

			}

			else {

				trackCollectionSize = trackCollection->size();
				//associate tracks
				// problem here track collection has no entries
				recSimColl=associatorByChi2->associateRecoToGen(trackCollection,
						genparticles,
						&iEvent,
						&iSetup);

				simRecColl=associatorByChi2->associateGenToReco(trackCollection,
						genparticles,
						&iEvent,
						&iSetup);

			}

			/// start here
			for(GenParticleCollection::size_type i = 0; i < genParticlesForChi2.size(); i++) {

				double quality = 0.;

				GenParticleRef tpr(genparticles,i);
				GenParticle* tp= const_cast<GenParticle*>(tpr.get());
				TrackingParticle::Vector momentumTP;
				TrackingParticle::Point vertexTP;
				vertexTP  =  tp->vertex();
				momentumTP =  tp->momentum();

				///////TO DO -- apply quality cuts to muons
				if(!(GenParticleSelection(tp))) continue;
				if((fabs(tp->eta()) > 2.0) && (fabs(tp->eta()) < 2.8)) {
					hFillGenMuonPtDen->Fill(tp->pt());
					if(tp->pt() > 5 ) hFillGenMuonEtaDen->Fill(tp->eta());
                                         cout<<"total pt:"<<tp->pt()<<"\t"<<"total eta:"<<tp->eta()<<endl;
					std::vector<std::pair<RefToBase<Track>, double> > rt;
					if(simRecColl.find(tpr) != simRecColl.end()){
						rt = (std::vector<std::pair<RefToBase<Track>, double> >) simRecColl[tpr];
						if (rt.size()!=0) {
							RefToBase<Track> assoc_recoTrack = rt.begin()->first;

							quality = rt.begin()->second;
							cout<<quality;
                                                        cout<<"pt:"<<tp->pt()<<"\t"<<"eta:"<<tp->eta()<<endl;
							hFillGenMuonPtNum->Fill(tp->pt());
							if(tp->pt() > 5 ) hFillGenMuonEtaNum->Fill(tp->eta());
						}
					}
				}
			} // end of gen loop

			/////// tracking collection
			for(View<Track>::size_type i=0; i<trackCollectionSize; ++i){
				//bool Track_is_matched = false; 
				RefToBase<Track> track(trackCollection, i);
				if(fabs(track->eta()) <= 2) continue;
				if(fabs(track->eta()) >= 2.8) continue;
				hFillDenRecPt->Fill(track->pt());                         
				if(track->pt() > 5.) hFillDenRecEta->Fill(track->eta()); 
				//std::vector<std::pair<TrackingParticleRef, double> > tp;
				std::vector<std::pair<GenParticleRef, double> > tp;
				std::vector<std::pair<GenParticleRef, double> > tpforfake;
				//TrackingParticleRef tpr;
				GenParticleRef tpr;
				GenParticleRef tprforfake;

				//Check if the track is associated to any gen particle
				if(recSimColl.find(track) != recSimColl.end()){

					tp = recSimColl[track];
					if (tp.size()!=0) {
						//Track_is_matched = true;

						tpr = tp.begin()->first;
						GenParticle* tp1= const_cast<GenParticle*>(tpr.get());

						if(!(GenParticleSelection(tp1))) continue;
						if(fabs(tp1->eta()) <= 2) continue;
						if(fabs(tp1->eta()) >= 2.8) continue;                                                
						double assocChi2 = -(tp.begin()->second);
						cout<<assocChi2;

						//So this track is matched to a gen particle, lets get that gen particle now
						if (  (simRecColl.find(tpr) != simRecColl.end())    ){
							std::vector<std::pair<RefToBase<Track>, double> > rt;
							if  (simRecColl[tpr].size() > 0){
								rt=simRecColl[tpr];
								RefToBase<Track> bestrecotrackforeff = rt.begin()->first;
								//Only fill the efficiency histo if the track found matches up to a gen particle's best choice
								if (bestrecotrackforeff == track) {

								}
							}
						}
					}
				}
				//End checking of Efficient muons
				//For Fakes --------------------------------------------


				//Check if finding a track associated to a gen particle fails, or if there is no track in the collection at all

				if( (recSimColl.find(track) == recSimColl.end() ) || ( recSimColl[track].size() == 0  ) ){
					//So we've now determined that this track is not associated to any gen, and fill our histo of fakes:
					//					if ((track->pt() > FakeRatePtCut) && (TMath::Abs(track->eta()) < 2.8) ) {};// Chi2UnmatchedME0Muon_Eta->Fill(fabs(track->eta()));
					hFillNumRecPt->Fill(track->pt());
					if(track->pt() > 5.) hFillNumRecEta->Fill(track->eta());


				}

				//Its possible that the track is associated to a gen particle, but isn't the best match and would still fail
				//In that case, we go to the gen particle...
				else if (recSimColl[track].size() > 0){
					tpforfake = recSimColl[track];
					tprforfake=tpforfake.begin()->first;
					//We now have the gen particle, to check
					GenParticle* tp2= const_cast<GenParticle*>(tprforfake.get());

					if(!(GenParticleSelection(tp2))) continue;
					if(fabs(tp2->eta()) <= 2) continue;
					if(fabs(tp2->eta()) >= 2.8) continue;  
					//If for some crazy reason we can't find the gen particle, that means its a fake
					if (  (simRecColl.find(tprforfake) == simRecColl.end())  ||  (simRecColl[tprforfake].size() == 0)  ) {
						//						if ((track->pt() > FakeRatePtCut) && (TMath::Abs(track->eta()) < 2.8)) //  Chi2UnmatchedME0Muon_Eta->Fill(fabs(track->eta()));
						//	{};
						hFillNumRecPt->Fill(track->pt());
						if(track->pt() > 5.) hFillNumRecEta->Fill(track->eta());

					}
					//We can probably find the gen particle
					else if(simRecColl[tprforfake].size() > 0)  {
						//We can now access the best possible track for the gen particle that this track was matched to
						std::vector<std::pair<RefToBase<Track>, double> > rtforfake;
						rtforfake=simRecColl[tprforfake];

						RefToBase<Track> bestrecotrack = rtforfake.begin()->first;
						//if the best reco track is NOT the track that we're looking at, we know we have a fake, that was within the cut, but not the closest
						if (bestrecotrack != track) {
							//		if ( (track->pt() > FakeRatePtCut) && (TMath::Abs(track->eta()) < 2.8) ) {}; // Chi2UnmatchedME0Muon_Eta->Fill(fabs(track->eta()));
//							hFillNumRecPt->Fill(track->pt());
//							hFillNumRecEta->Fill(track->eta());

						}

					}
				}

				//End For Fakes --------------------------------------------


				//				if (TMath::Abs(track->eta()) < 3.0 && TMath::Abs(track->eta()) > 2.0) {  hFillDenRecPt->Fill(track->pt());                         
				//    hFillDenRecEta->Fill(track->eta());   
				//};// CheckME0Muon_Eta->Fill(fabs(track->eta()));

				//				NormChi2_h->Fill(track->normalizedChi2());
				//				NormChi2Prob_h->Fill(TMath::Prob(track->chi2(),(int)track->ndof()));
				//				NormChi2VsHits_h->Fill(track->numberOfValidHits(),track->normalizedChi2());
				//				chi2_vs_eta_h->Fill((track->eta()),track->normalizedChi2());

				//nhits_vs_eta_h->Fill((track->eta()),track->numberOfValidHits());


			}//END for(View<Track>::size_type i=0; i<trackCollectionSize; ++i){
		}//END for (unsigned int www=0;www<label.size();www++)
		}// END for (unsigned int www=0;www<label.size();www++)




		/////////////////////===============================================///////////////
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
		if(me0muons.size() > 0) {hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time()); hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(0)).me0segment().timeErr()); SignalMuonTime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time(), OurMuons->at(me0muons.at(0)).me0segment().timeErr()); }

		if(me0muons.size() > 1) {hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(1)).me0segment().time()); hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(1)).me0segment().timeErr()); SignalMuonTime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time(), OurMuons->at(me0muons.at(0)).me0segment().timeErr()); }
		if(me0muons.size() > 2) {hFillSignalMuontime->Fill(OurMuons->at(me0muons.at(2)).me0segment().time()); hFillSignalMuontimeErr->Fill(OurMuons->at(me0muons.at(2)).me0segment().timeErr());SignalMuonTime->Fill(OurMuons->at(me0muons.at(0)).me0segment().time(), OurMuons->at(me0muons.at(0)).me0segment().timeErr());}

		TLorentzVector muon1,muon2;
		for( unsigned int j = 0;  j < me0muons.size(); j++) {
			hFillRecoEta->Fill(OurMuons->at(me0muons.at(j)).eta());
			for( unsigned int ji = j+1; ji < me0muons.size(); ji++) { 

				//              hFillZMass->Fill((OurMuons->at(me0muons.at(j)).p()+ OurMuons->at(me0muons.at(ji)).p()).M());            
				//          if(reco::deltaR(OurMuons->at(me0muons.at(j)).eta(),OurMuons->at(me0muons.at(j)).phi(),OurMuons->at(me0muons.at(ji)).eta(),OurMuons->at(me0muons.at(ji)).phi()) < 0.3) continue;
				muon1.SetPtEtaPhiM(OurMuons->at(me0muons.at(j)).pt(),OurMuons->at(me0muons.at(j)).eta(),OurMuons->at(me0muons.at(j)).phi(),0);//OurMuons->at(me0muons.at(j)).energy());
				muon2.SetPtEtaPhiM(OurMuons->at(me0muons.at(ji)).pt(),OurMuons->at(me0muons.at(ji)).eta(),OurMuons->at(me0muons.at(ji)).phi(),0);//OurMuons->at(me0muons.at(ji)).energy());
				hFillZMass->Fill((muon1+muon2).M());
				if((OurMuons->at(me0muons.at(ji)).me0segment().time() >= 5.5  && OurMuons->at(me0muons.at(ji)).me0segment().time() <= 30.5) && (OurMuons->at(me0muons.at(j)).me0segment().time() >= 5.5  && OurMuons->at(me0muons.at(j)).me0segment().time() <= 30.5)) {

					hFillZMassInWindow->Fill((muon1+muon2).M());
				} else {
					hFillZMassOutWindow->Fill((muon1+muon2).M());

				}                     

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
				BGMuonTime->Fill(OurMuons->at(t).me0segment().time(),OurMuons->at(t).me0segment().timeErr()); 
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
	/*	void 
		ME0TimingAnalysis::beginJob()
		{

		hFillSignalMuontime = fs->make<TH1F>("hFillSignalMuontime","hFillSignalMuontime",1000,-300,300);  
		hFillPUMuontime = fs->make<TH1F>("hFillPUMuontime","hFillPUMuontime",1000,-300,300);
		hFillSignalMuontimeErr = fs->make<TH1F>("hFillSignalMuontimeErr","hFillSignalMuontimeErr",1000,0,10);  
		hFillPUMuontimeErr = fs->make<TH1F>("hFillPUMuontimeErr","hFillPUMuontimeErr",1000,0,10);

		SignalMuonTime  = fs->make<TH1F>("SignalMuonTime","SignalMuonTime",1000,-300,300);  
		BGMuonTime  = fs->make<TH1F>("BGMuonTime","BGMuonTime",1000,-300,300);  
		hFillZMass = fs->make<TH1F>("hFillZMass","hFillZMass",500,0,250); 
		hFillRecoEta = fs->make<TH1F>("hFillRecoEta","hFillRecoEta",500,-5,5);
		hFillZGenMass = fs->make<TH1F>("hFillZGenMass","hFillZGenMass",500,0,250);    

		hFillZMassInWindow = fs->make<TH1F>("hFillZMassInWindow","hFillZMassInWindow",500,0,250);
		hFillZMassOutWindow = fs->make<TH1F>("hFillZMassOutWindow","hFillZMassOutWindow",500,0,250);


		}
		*/
	// ------------ method called once each job just after ending the event loop  ------------
	/*	void 
		ME0TimingAnalysis::endJob() 
		{
		}
		*/

	void ME0TimingAnalysis::beginRun(Run const&, EventSetup const& setup) {

		hFillSignalMuontime = fs->make<TH1F>("hFillSignalMuontime","hFillSignalMuontime",1000,-300,300);  
		hFillPUMuontime = fs->make<TH1F>("hFillPUMuontime","hFillPUMuontime",1000,-300,300);
		hFillSignalMuontimeErr = fs->make<TH1F>("hFillSignalMuontimeErr","hFillSignalMuontimeErr",1000,0,10);  
		hFillPUMuontimeErr = fs->make<TH1F>("hFillPUMuontimeErr","hFillPUMuontimeErr",1000,0,10);

		SignalMuonTime  = fs->make<TH1F>("SignalMuonTime","SignalMuonTime",1000,-300,300);  
		BGMuonTime  = fs->make<TH1F>("BGMuonTime","BGMuonTime",1000,-300,300);  
		hFillZMass = fs->make<TH1F>("hFillZMass","hFillZMass",500,0,250); 
		hFillRecoEta = fs->make<TH1F>("hFillRecoEta","hFillRecoEta",500,-5,5);
		hFillZGenMass = fs->make<TH1F>("hFillZGenMass","hFillZGenMass",500,0,250);    

		hFillZMassInWindow = fs->make<TH1F>("hFillZMassInWindow","hFillZMassInWindow",500,0,250);
		hFillZMassOutWindow = fs->make<TH1F>("hFillZMassOutWindow","hFillZMassOutWindow",500,0,250);

		hFillGenMuonPtDen= fs->make<TH1F>("hFillGenMuonPtDen","hFillGenMuonPtDen",800,0,200);
		hFillGenMuonEtaDen = fs->make<TH1F>("hFillGenMuonEtaDen","hFillGenMuonEtaDen",100,-5,5);
		hFillGenMuonPtNum= fs->make<TH1F>("hFillGenMuonPtNum","hFillGenMuonPtNum",800,0,200);   
		hFillGenMuonEtaNum = fs->make<TH1F>("hFillGenMuonEtaNum","hFillGenMuonEtaNum",100,-5,5);
		hFillDenRecPt = fs->make<TH1F>("hFillDenRecPt","hFillDenRecPt",800,0,200);             
		hFillDenRecEta = fs->make<TH1F>("hFillDenRecEta","hFillDenRecEta",100,-5,5);
		hFillNumRecPt = fs->make<TH1F>("hFillNumRecPt","hFillNumRecPt",800,0,200);             
		hFillNumRecEta = fs->make<TH1F>("hFillNumRecEta","hFillNumRecEta",100,-5,5);

		if (UseAssociators) {                                                                
			edm::ESHandle<TrackAssociatorBase> theAssociator;
			for (unsigned int w=0;w<associators.size();w++) {
				setup.get<TrackAssociatorRecord>().get(associators[w],theAssociator);
				associator.push_back( theAssociator.product() );
			}
		}

	}


	void ME0TimingAnalysis::endRun(Run const&, EventSetup const& setup) {

	}

	void ME0TimingAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

	bool ME0TimingAnalysis::GenParticleSelection(const GenParticle* tp) {
		bool goodgen = false;
		if(!(fabs(tp->pdgId()) == 13)) return false;
		if(!(tp->status() == 1)) return false;
		if(!(tp->numberOfMothers() > 0)) return false; 
		if(fabs(tp->mother()->pdgId()) == 23) {
			goodgen = true;}
		else if(abs(tp->pdgId()) == abs(tp->mother()->pdgId())) {
			if(tp->mother()->numberOfMothers() > 0) {
				if(abs(tp->mother()->mother()->pdgId()) == 23) {
					goodgen = true; }
			}
		}  

		return goodgen;
	}
	//define this as a plug-in
	DEFINE_FWK_MODULE(ME0TimingAnalysis);
