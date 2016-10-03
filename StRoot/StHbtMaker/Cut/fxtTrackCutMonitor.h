#ifndef fxtTrackCutMonitor_hh
#define fxtTrackCutMonitor_hh

class StHbtTrack;

#include "StHbtMaker/Infrastructure/StHbtHisto.hh"
#include "StHbtMaker/Base/StHbtCutMonitor.hh"

class fxtTrackCutMonitor : public StHbtCutMonitor {

private:
  StHbt1DHisto* mDCA;
  StHbt1DHisto* mDCAGlobal;
  StHbt1DHisto* mNhits;
  StHbt1DHisto* mPt;
  StHbt1DHisto* mChiSqr;
  StHbt1DHisto* mNsigmaPion;
  StHbt1DHisto* mNsigmaKaon;
  StHbt1DHisto* mNsigmaProton;
  StHbt1DHisto* mNsigmaElectron;
  StHbt2DHisto* mPvsDedx;
  StHbt2DHisto* mPtvsDedx;
  StHbt2DHisto* mEtaVsPhi;
  StHbt1DHisto* mRapidity;
  StHbt2DHisto* mYvsPt;
  double monMass;
  
public:
  fxtTrackCutMonitor(const char*,const double);   // default constructor
  fxtTrackCutMonitor( const  fxtTrackCutMonitor& cutMoni);
  virtual ~fxtTrackCutMonitor();

  virtual void Fill(const StHbtTrack* track);

  // These dummy Fill() functions were introduced to remove a compiler
  //   warning related to overloaded base-class Fill() functions being 
  //   hidden by a single version of Fill() in this derived class
  void Fill(const StHbtParticleCollection* d) {;}
  void Fill(const StHbtEvent *d1, const StHbtParticleCollection* d2) {;}
  void Fill(const StHbtPair* d) {;}
  void Fill(const StHbtKink* d) {;}
  void Fill(const StHbtV0* d) {;}
  void Fill(const StHbtEvent* d) {;}
  
  StHbt1DHisto* DCA() { return mDCA; }
  StHbt1DHisto* DCAGlobal() { return mDCAGlobal; }
  StHbt1DHisto* Nhits() { return mNhits; }
  StHbt1DHisto* Pt() { return mPt; }
  StHbt1DHisto* ChiSqr() { return mChiSqr; }
  StHbt1DHisto* NsigmaPion() { return mNsigmaPion; }
  StHbt1DHisto* NsigmaKaon() { return mNsigmaKaon; }
  StHbt1DHisto* NsigmaProton() { return mNsigmaProton; }
  StHbt1DHisto* NsigmaElectron() { return mNsigmaElectron; }
  StHbt2DHisto* PvsDedx() { return mPvsDedx; }
  StHbt2DHisto* PtvsDedx() { return mPtvsDedx; }
  StHbt2DHisto* EtaVsPhi() { return mEtaVsPhi; }
  StHbt1DHisto* Rapidity() { return mRapidity; }
  StHbt2DHisto* YvsPt() { return mYvsPt; }
  
#ifdef __ROOT__  
 ClassDef(fxtTrackCutMonitor, 1)
#endif
};

#endif
