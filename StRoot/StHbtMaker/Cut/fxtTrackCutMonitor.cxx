#include "StHbtMaker/Infrastructure/StHbtTrack.hh"
#include "StHbtMaker/Cut/fxtTrackCutMonitor.h"
#include <cstdio>
#include <string>
#include "StLorentzVector.hh"

#ifdef __ROOT__ 
ClassImp(fxtTrackCutMonitor)
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
fxtTrackCutMonitor::fxtTrackCutMonitor(const char* name,const double ParticleMass){ // default constructor
  monMass = ParticleMass;
  string s("fxtTrackCutMonitor");
  string n(name);
  mDCA= new StHbt1DHisto( (s+n+"mDCA").c_str(),"DCA",1000,0,5.); 
  mDCAGlobal= new StHbt1DHisto( (s+n+"mDCAGlobal").c_str(),"DCA Global",1000,0.,5.); 
  mNhits = new StHbt1DHisto( (s+n+"mNhits").c_str(),"nHits",50,0.,50.);
  mPt = new StHbt1DHisto( (s+n+"mPt").c_str(),"Pt",1300,0.,1.3);
  mChiSqr = new StHbt1DHisto( (s+n+"mChiSqr").c_str(),"ChiSqr",100,0.,10.);
  mNsigmaPion = new StHbt1DHisto( (s+n+"mNsigmaPion").c_str(),"NsigmaPion",70,-35.,35.);
  mNsigmaKaon = new StHbt1DHisto( (s+n+"mNsigmaKaon").c_str(),"NsigmaKaon",70,-35.,35.);
  mNsigmaProton = new StHbt1DHisto( (s+n+"mNsigmaProton").c_str(),"NsigmaProton",70,-35.,35.);
  mNsigmaElectron = new StHbt1DHisto( (s+n+"mNsigmaElectron").c_str(),"NsigmaElectron",70,-35.,35.);
  mPvsDedx= new StHbt2DHisto( (s+n+"mPvsDedx").c_str(),"Momentum (GeV/c) vs Energy loss (a.u.)",100,-3.,3.,50,0.,5e-6); 
  mPtvsDedx= new StHbt2DHisto( (s+n+"mPtvsDedx").c_str(),"Transverse Momentum (GeV/c) vs Energy loss (a.u.)",100,-3.,3.,50,0.,5e-6);
  mEtaVsPhi= new StHbt2DHisto( (s+n+"mEta").c_str(),"Pseudorapidity vs Phi",70,-1.2,1.2,180,0.,360.);
  mRapidity = new StHbt1DHisto( (s+n+"mRapidity").c_str(),"Rapidity",70,-1.2,1.2);
  mYvsPt = new StHbt2DHisto( (s+n+"mYvsPt").c_str(),"Rapidity vs Pt",130,0.,1.3,120,-1.2,1.2);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
fxtTrackCutMonitor::fxtTrackCutMonitor( const  fxtTrackCutMonitor& cutMoni)  {
  mDCA =new StHbt1DHisto(*(cutMoni.mDCA));
  mDCAGlobal =new StHbt1DHisto(*(cutMoni.mDCAGlobal));
  mNhits = new StHbt1DHisto(*(cutMoni.mNhits));
  mPt = new StHbt1DHisto(*(cutMoni.mPt));
  mChiSqr = new StHbt1DHisto(*(cutMoni.mChiSqr));
  mNsigmaPion = new StHbt1DHisto(*(cutMoni.mNsigmaPion));
  mNsigmaKaon = new StHbt1DHisto(*(cutMoni.mNsigmaKaon));
  mNsigmaProton = new StHbt1DHisto(*(cutMoni.mNsigmaProton));
  mNsigmaElectron = new StHbt1DHisto(*(cutMoni.mNsigmaElectron));
  mPvsDedx =new StHbt2DHisto(*(cutMoni.mPvsDedx));
  mPtvsDedx= new StHbt2DHisto(*(cutMoni.mPtvsDedx));
  mEtaVsPhi= new StHbt2DHisto(*(cutMoni.mEtaVsPhi));
  mRapidity= new StHbt1DHisto(*(cutMoni.mRapidity));
  mYvsPt= new StHbt2DHisto(*(cutMoni.mYvsPt));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
fxtTrackCutMonitor::~fxtTrackCutMonitor(){
  delete mDCA;
  delete mDCAGlobal;
  delete mNhits;
  delete mPt;
  delete mChiSqr;
  delete mNsigmaPion;
  delete mNsigmaKaon;
  delete mNsigmaProton;
  delete mNsigmaElectron;
  delete mPvsDedx ;
  delete mPtvsDedx;
  delete mEtaVsPhi;
  delete mRapidity;
  delete mYvsPt;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  fxtTrackCutMonitor::Fill(const StHbtTrack* track){

  float TEnergy = ::sqrt(track->P().mag2()+monMass*monMass);
  float TRapidity = 0.5*::log((TEnergy+track->P().z())/
			    (TEnergy-track->P().z()));
  float TEta = 0.5*::log((track->P().mag()+track->P().z())/(track->P().mag()-track->P().z()));
  float TPhi = atan2(track->P().y(),track->P().x());
  if(TPhi<0.0) { TPhi += 2.*TMath::Pi(); }
  TPhi *= 180./TMath::Pi();
  
  mDCA->Fill( track->DCAxy(), 1.); //using for DCA xyz
  mDCAGlobal->Fill( track->DCAxyGlobal(), 1.); //using for DCAxyz global
  mNhits->Fill(track->NHits(), 1.);
  mPt->Fill(track->Pt(), 1.);
  mChiSqr->Fill(track->ChiSquaredXY(), 1.); //use for ChiSqr variable (xyz or xy whichever it is)
  mNsigmaPion->Fill(track->NSigmaPion(), 1.);
  mNsigmaKaon->Fill(track->NSigmaKaon(), 1.);
  mNsigmaProton->Fill(track->NSigmaProton(), 1.);
  mNsigmaElectron->Fill(track->NSigmaElectron(), 1.);
  mPvsDedx->Fill( fabs(track->P().mag())*track->Charge(), track->dEdx(), 1.);
  mPtvsDedx->Fill( fabs(track->Pt())*track->Charge(), track->dEdx(), 1.);
  mEtaVsPhi->Fill(TEta,TPhi, 1.);
  mRapidity->Fill(TRapidity, 1.);
  mYvsPt->Fill(track->Pt(), TRapidity);
}


