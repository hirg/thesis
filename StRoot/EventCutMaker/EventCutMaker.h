#ifndef EventCutMaker_hh
#define EventCutMaker_hh
#include "StMaker.h"
#include "TString.h"
#include "TArrayI.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "StRefMultCorr/StRefMultCorr.h"

//Leave this in
#ifndef St_NO_NAMESPACES
using std::string;
#endif

class TH2F;
class TFile;
class StMuDst;
class StMuDstMaker;
class StMuEvent;
class StMuTrack;

class EventCutMaker : public StMaker {

public:
  
  EventCutMaker(StMuDstMaker* muDst, StRefMultCorr* refmultCorr);
  ~EventCutMaker() {;}
  
  void Clear(Option_t *option="");
  Int_t Init();
  Int_t Make();
  Int_t Finish();
  
  void setMuDstMaker(StMuDstMaker*);
  void setRefmultCorrUtil(StRefMultCorr*);
  Bool_t eventPass();
  void SetFileName(TString fileName) {mFileName = fileName;}
  void SetZdcCut(const Int_t cut){ mZdcCut = cut;}
  void SetTofMatchCut(const Int_t cut){ mNTofMatchCut= cut;}
  void SetVrCut(const Float_t cut){ mVrCut = cut;}
  void SetVrCenter(const Float_t x, const Float_t y){ mVrCenter[0] = x; mVrCenter[1] = y;}
  void SetVzCut(const Float_t cut){ mVzCut = cut;}
  void SetCent9Cut(const UInt_t lo, const UInt_t high){ mCent9Cut[0] = lo; mCent9Cut[1] = high;}
  void SetCent16Cut(const UInt_t lo, const UInt_t high){ mCent16Cut[0] = lo; mCent16Cut[1] = high;}
  
private:
  
  bool acceptEvent(StMuEvent*);
    Bool_t checkTrigger(StMuEvent*,Bool_t);
    Bool_t checkVertex(StMuEvent*,Bool_t);
    Bool_t checkZdc(StMuEvent*,Bool_t);
    Bool_t checkRefmult(StMuEvent*,Bool_t);
    Bool_t checkNTofMatched(StMuDst*,Bool_t);
    Int_t nTofMatchedTracks(StMuDst*);
  
  StMuDstMaker* mMuDstMaker; //!
  TFile* mFile; //!
  TString mFileName; //!
  TString mLastMuFile; //!
  StMuEvent *mMuEvent;         //! pointer to Mu-DST Event array

  TH1I* mCutIndependent;
  TH1I* mCutCumulative;
  TH1I* mCentBin16All;
  TH1I* mCentBin16Pass;
  TH1I* mCentBin9All;
  TH1I* mCentBin9Pass;
  TH1I* mRefmult;
  TH1I* mNTofMatchedTracksAll;
  TH1I* mNTofMatchedTracksPass;
  TH1F* mCutVzAll;
  TH1F* mCutVzPass;
  TH2F* mCutVrAll;
  TH2F* mCutVrPass;
  
  int nEventsPassed;
  int nEventsFailed;
  int nBytes;
  int bin;

  StRefMultCorr* mRefmultCorrUtil;
  Double_t mZdcCoincidenceRate;

  Int_t mZdcCut;
  Int_t mNTofMatchCut;
  Float_t mVzCut;
  Float_t mVrCut;
  Float_t mVrCenter[2];
  UShort_t mCent9Cut[2];
  UShort_t mCent16Cut[2];
  Bool_t mEventPass;

  Float_t CalcDcaSigned(const StThreeVectorF vertex, const StPhysicalHelixD helix);

  ClassDef(EventCutMaker,1)
};

inline void EventCutMaker::setMuDstMaker(StMuDstMaker* f) {mMuDstMaker=f;}

inline void EventCutMaker::setRefmultCorrUtil(StRefMultCorr* f) {mRefmultCorrUtil=f;}

inline Bool_t EventCutMaker::eventPass() {return mEventPass;}

#endif

