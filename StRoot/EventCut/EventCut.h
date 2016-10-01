#ifndef EventCut_hh
#define EventCut_hh

#include <iostream>
#include <cmath>
#include "TTree.h"
#include "TEventList.h"
#include "TStopwatch.h"
#include "TLeaf.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#ifndef St_NO_NAMESPACES
using std::cout;
using std::endl;
#endif

class EventCut {

public:
  
    EventCut(TTree* tree);
    ~EventCut();

    void SetVxCut(const Float_t xLow, const Float_t xHigh ) { mVxCut[0] = xLow; mVxCut[1] = xHigh; }
    void SetVyCut(const Float_t yLow, const Float_t yHigh ) { mVyCut[0] = yLow; mVyCut[1] = yHigh; }
    void SetVzCut(const Float_t zLow, const Float_t zHigh ) { mVzCut[0] = zLow; mVzCut[1] = zHigh; }
    void SetVrCut(const Float_t rCut) { mVrCut = rCut; }
    void SetVrCenter(const Float_t x, const Float_t y) { mVrCenter[0] = x; mVrCenter[1] = y; }

    void SetZdcCut(const Float_t zdcCut) { mZdcCut = zdcCut; }
    void SetCent9Cut(const Float_t centLow, const Float_t centHigh) { mCent9Cut[0] = centLow; mCent9Cut[1] = centHigh; }
    void SetCent16Cut(const Float_t centLow, const Float_t centHigh) { mCent16Cut[0] = centLow; mCent16Cut[1] = centHigh; }
    void SetTriggers(const Float_t* triggers, Int_t nTriggers);
    void SetTree(TTree* tree) { mTree = tree; }

    TEventList* makeEventList();

private:

    Int_t mZdc[2];
    UShort_t mCent9;
    UShort_t mCent16;
    UShort_t mRefmult;
    Float_t mVx;
    Float_t mVy;
    Float_t mVz;
    Float_t mVr;
    Float_t mTrigger;

    Int_t mZdcCut;
    UShort_t mCent9Cut[2];
    UShort_t mCent16Cut[2];
    Double_t mRefmultCut[2];
    Double_t mVxCut[2];
    Double_t mVyCut[2];
    Double_t mVzCut[2];
    Double_t mVrCut;
    Float_t mTriggerCut[8];
    Int_t mNTriggers;
    Double_t mVrCenter[2];
    
    Bool_t mZdcPass;
    Bool_t mTriggerPass;
    Bool_t mVzPass;
    Bool_t mVrPass;
    Bool_t mHasVertexPass;
    Bool_t mRefmultPass;
    Bool_t mCent9Pass;
    Bool_t mCent16Pass;
    Bool_t mEventPass;

    TH1I* mAllCutsIndependentHist;
    TH1I* mAllCutsCumulativeHist;
    TH1F* mTriggerPassHist;
    TH1I* mRefmultAllHist;
    TH1I* mRefmultPassHist;
    TH1F* mVzAllHist;
    TH1F* mVzPassHist;
    TH2F* mVrAllHist;
    TH2F* mVrPassHist;
    TH2F* mZdcAllHist;
    TH2F* mZdcPassHist;

    Bool_t checkVz() { return (mVz >= mVzCut[0]) && (mVz <= mVzCut[1]); }
    Bool_t checkRefMult() { return (mRefmult >= mRefmultCut[0]) && (mRefmult <= mRefmultCut[1]); }
    Bool_t checkVr() { return (mVr <= mVrCut); }
    Bool_t checkHasVertex() { return ((mVx != 0) || (mVy != 0) || (mVz != 0)); }
    Bool_t checkZdc() { return (mZdc[0] <= mZdcCut) && (mZdc[1] <= mZdcCut); }
    Bool_t checkTrigger();

    void bookHistograms();
    void fillHistograms();
    void writeHistograms();

    TTree* mTree;
    TEventList* mEventList;

    ClassDef(EventCut,1)

};

#endif
