#ifndef fxtPairCutMonitor_hh
#define fxtPairCutMonitor_hh

class StHbtPair;

#include "StHbtMaker/Infrastructure/StHbtHisto.hh"
#include "StHbtMaker/Base/StHbtCutMonitor.hh"

class fxtPairCutMonitor : public StHbtCutMonitor {

private:
  StHbt1DHisto* mKt;
  StHbt1DHisto* mKtLowQ;
  StHbt1DHisto* mMt;
  StHbt1DHisto* mMtLowQ;
  StHbt1DHisto* mFractionOfMergedRow;
  StHbt1DHisto* mSplittingLevel;
  StHbt2DHisto* mPairRapidityVsEmissionAngle;

public:
  fxtPairCutMonitor(const char*);   // default constructor
  fxtPairCutMonitor( const  fxtPairCutMonitor& cutMoni);
  virtual ~fxtPairCutMonitor();

  virtual void Fill(const StHbtPair* pair);

  // These dummy Fill() functions were introduced to remove a compiler
  //   warning related to overloaded base-class Fill() functions being 
  //   hidden by a single version of Fill() in this derived class
  void Fill(const StHbtParticleCollection* d) {;}
  void Fill(const StHbtEvent *d1, const StHbtParticleCollection* d2) {;}
  void Fill(const StHbtTrack* d) {;}
  //void Fill(const StHbtPair* d) {;}
  void Fill(const StHbtKink* d) {;}
  void Fill(const StHbtV0* d) {;}
  void Fill(const StHbtEvent* d) {;}
  
  StHbt1DHisto* Kt() { return mKt; }
  StHbt1DHisto* KtLowQ() { return mKtLowQ; }
  StHbt1DHisto* Mt() { return mMt; }
  StHbt1DHisto* MtLowQ() { return mMtLowQ; }
  StHbt1DHisto* FractionOfMergedRow() { return mFractionOfMergedRow; }
  StHbt1DHisto* SplittingLevel() { return mSplittingLevel; }
  StHbt2DHisto* PairRapidityVsEmissionAngle() { return mPairRapidityVsEmissionAngle; }

#ifdef __ROOT__  
 ClassDef(fxtPairCutMonitor, 1)
#endif
};

#endif
