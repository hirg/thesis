#include "EventCutMaker.h"
#include "StEventTypes.h"

#include "TFile.h"
#include "TChain.h"
#include "StMessMgr.h"

#include "StMuDSTMaker/COMMON/StMuException.hh"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuDebug.h"
#include "StMuDSTMaker/COMMON/StMuCut.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"

#include "TClonesArray.h"
#include "TBranch.h"
#include "TTree.h"

#include "StRefMultCorr/StRefMultCorr.h"
#include "TH2.h"

ClassImp(EventCutMaker)

EventCutMaker::EventCutMaker(StMuDstMaker* muDst, StRefMultCorr* refmultCorr) : StMaker("EventCutMaker") {
	mFile = 0;
	mLastMuFile = "";
	nEventsPassed = nEventsFailed = 0;

	mMuDstMaker = muDst;
    mRefmultCorrUtil = refmultCorr;

}

Int_t EventCutMaker::Init() {
	
	mFile = new TFile("eventCuts.root", "RECREATE");
	mFile->SaveSelf();
	cout << "The output filename is eventCuts.root" << endl;
	
    mRefmultCorrUtil = new StRefMultCorr("refmult");
    mCutIndependent = new TH1I("mCutIndependent ","mCutIndependent ",9,-0.5,8.5);
    mCutCumulative = new TH1I("mCutCumulative ","mCutCumulative ",9,-0.5,8.5);
    mCentBin16All = new TH1I("mCentBin16All ","mCentBin16All ",16,-0.5,15.5);
    mCentBin16Pass = new TH1I("mCentBin16Pass ","mCentBin16Pass ",16,-0.5,15.5);
    mCentBin9All = new TH1I("mCentBin9All ","mCentBin9All ",9,-0.5,8.5);
    mCentBin9Pass = new TH1I("mCentBin9Pass ","mCentBin9Pass ",9,-0.5,8.5);
    mNTofMatchedTracksAll = new TH1I("mNTofMatchedTracksAll","mNTofMatchedTracksAll",1000,-0.5,999.5);
    mNTofMatchedTracksPass = new TH1I("mNTofMatchedTracksPass","mNTofMatchedTracksPass",1000,-0.5,999.5);
    mRefmult = new TH1I("mRefmult ","mRefmult ",1000,-0.5,999.5);
    mCutVzAll = new TH1F("mCutVzAll ","mCutVzAll ",200,-100,100);
    mCutVzPass = new TH1F("mCutVzPass ","mCutVzPass ",200,-100,100);
    mCutVrAll = new TH2F("mCutVrAll ","mCutVrAll ",200,-4,4,200,-4,4);
    mCutVrPass = new TH2F("mCutVrPass ","mCutVrPass ",200,-4,4,200,-4,4);

    const Char_t* binTitles[9] = {"All Events", "Trigger", "# ToF Matches", "Has Vertex", "V_R", "V_Z", "RefMult 9", "RefMult 16", "ZDC cut"};

    for(Int_t i = 0; i <= 8; i++)
    {
        mCutIndependent->GetXaxis()->SetBinLabel(i+1,binTitles[i]); 
        mCutCumulative->GetXaxis()->SetBinLabel(i+1,binTitles[i]); 
    }
	
	cout << "end of ECM::Init()\n";
	return StMaker::Init();
}

Int_t EventCutMaker::Make() {
	
	StMuDst* muDst = mMuDstMaker->muDst();
	mMuEvent = muDst->event();
	
	TString currentMuFile = mMuDstMaker->chain()->GetFile()->GetName();
	if(currentMuFile != mLastMuFile) {
        mRefmultCorrUtil->init(mMuEvent->runId());
		cout << "New File: " << currentMuFile.Data() << endl;
		mLastMuFile = currentMuFile;
	}

    mEventPass = ( !mRefmultCorrUtil->isBadRun(mMuEvent->runId()) );
    if( !mEventPass ) {return kStOK;} // Don't event bother with the rest if the run is bad
    
    mCutIndependent->Fill(0);
    mCutCumulative->Fill(0);

    if( !checkTrigger(mMuEvent,mEventPass)  ) {mEventPass = kFALSE;}
    if( !checkNTofMatched(muDst,mEventPass) ) {mEventPass = kFALSE;}
    if( !checkVertex(mMuEvent,mEventPass)   ) {mEventPass = kFALSE;}
    if( !checkRefmult(mMuEvent,mEventPass)  ) {mEventPass = kFALSE;}
    if( !checkZdc(mMuEvent,mEventPass)      ) {mEventPass = kFALSE;}

	if(!mEventPass) {nEventsFailed++;}
	else {nEventsPassed++;}

	mMuEvent->Clear();
    return kStOK;
}

Int_t EventCutMaker::Finish() {

	mFile->Write();
	mFile->Close();
	
	cout << "EventCutMaker::Finish()\n\n";
	cout << "\t  nEventsPassed: " << nEventsPassed << " events.\n";
	cout << "\t  nEventsFailed: " << nEventsFailed << " events.\n";
	cout << "Finish() ended : EventCutMaker.cxx\n";
	
	return kStOK;
}
	
//__________________________________________________________________//
void EventCutMaker::Clear(Option_t *opt) {
	StMaker::Clear();
}

Bool_t EventCutMaker::checkTrigger(StMuEvent* event, Bool_t prevPass) 
{

                        // 193 GeV U+U; 2012
	Bool_t triggerPass = (event->triggerIdCollection().nominal().isTrigger(400104) ||
				event->triggerIdCollection().nominal().isTrigger(400108) ||
				event->triggerIdCollection().nominal().isTrigger(400114) ||
				event->triggerIdCollection().nominal().isTrigger(400118) ||
				event->triggerIdCollection().nominal().isTrigger(400124) ||
				event->triggerIdCollection().nominal().isTrigger(400134) || 

                // 200 GeV Au+Au; 2011
	            event->triggerIdCollection().nominal().isTrigger(350003) ||
				event->triggerIdCollection().nominal().isTrigger(350013) ||
				event->triggerIdCollection().nominal().isTrigger(350023) ||
				event->triggerIdCollection().nominal().isTrigger(350033) ||
				event->triggerIdCollection().nominal().isTrigger(350043) || 

                // 39 GeV
	            event->triggerIdCollection().nominal().isTrigger(280001) ||
				event->triggerIdCollection().nominal().isTrigger(280002) ||

                // 27 GeV
	            event->triggerIdCollection().nominal().isTrigger(360001) ||
				event->triggerIdCollection().nominal().isTrigger(360002) ||

                // 19 GeV
	            event->triggerIdCollection().nominal().isTrigger(340001) ||
				event->triggerIdCollection().nominal().isTrigger(340011) ||
				event->triggerIdCollection().nominal().isTrigger(340021) ||

                // 11  GeV
	            event->triggerIdCollection().nominal().isTrigger(310004) ||
	            event->triggerIdCollection().nominal().isTrigger(310014) ||

                // 7  GeV
	            event->triggerIdCollection().nominal().isTrigger(310004) ||
	            event->triggerIdCollection().nominal().isTrigger(310014) ||

                // 14.5 GeV Au+Au; 2014
	            event->triggerIdCollection().nominal().isTrigger(440005) ||
				event->triggerIdCollection().nominal().isTrigger(440015) );
	
    if(triggerPass) { mCutIndependent->Fill(1); }
    if(triggerPass && prevPass) { mCutCumulative->Fill(1); }
	return triggerPass;
}

Bool_t EventCutMaker::checkNTofMatched(StMuDst* muDst, Bool_t prevPass)
{
    Bool_t tofPass = kTRUE;
    Int_t nTof = nTofMatchedTracks(muDst);

    mNTofMatchedTracksAll->Fill(nTof);

    if(nTof < mNTofMatchCut) { tofPass = kFALSE; }

    if(tofPass) {
        mNTofMatchedTracksPass->Fill(nTof);
        mCutIndependent->Fill(2);
    }


    if(tofPass && prevPass) { mCutCumulative->Fill(2); }

	return tofPass;
}



Bool_t EventCutMaker::checkVertex(StMuEvent* event, Bool_t prevPass) {

    Bool_t vertexZeroPass = kTRUE;
    Bool_t vertexRPass = kTRUE;
    Bool_t vertexZPass = kTRUE;
	Float_t vertexX = event->primaryVertexPosition().x();
	Float_t vertexY = event->primaryVertexPosition().y();
	Float_t vertexZ = event->primaryVertexPosition().z();
    Float_t Vr = sqrt( pow( (vertexX-mVrCenter[0]),2 ) + pow( (vertexY-mVrCenter[1]),2 ) );

    mCutVrAll->Fill(vertexX,vertexY);
    mCutVzAll->Fill(vertexZ);

	if (fabs(vertexX) < 1e-5 && fabs(vertexY) < 1e-5 && fabs(vertexZ) < 1e-5) { vertexZeroPass = kFALSE; }
	if (Vr >  mVrCut) { vertexRPass = kFALSE; }
	if (fabs(vertexZ) > mVzCut) { vertexZPass = kFALSE; }

    if(vertexZeroPass) { mCutIndependent->Fill(3);}

    if(vertexRPass) { 
        mCutIndependent->Fill(4);
        mCutVrPass->Fill(vertexX,vertexY);
    }
    
    if(vertexZPass) {
        mCutIndependent->Fill(5);
        mCutVzPass->Fill(vertexZ);
    }

    if(vertexZeroPass && prevPass) { mCutCumulative->Fill(3);}
    if(vertexZeroPass && vertexRPass && prevPass) { mCutCumulative->Fill(4);}
    if(vertexZeroPass && vertexRPass && vertexZPass && prevPass) { mCutCumulative->Fill(5);}

    return (vertexZeroPass && vertexRPass && vertexZPass);
}

Bool_t EventCutMaker::checkRefmult(StMuEvent* event, Bool_t prevPass) {

    Bool_t cent9Pass = kTRUE;
    Bool_t cent16Pass = kTRUE;
    mRefmultCorrUtil->initEvent(mMuEvent->refMult(), mMuEvent->primaryVertexPosition().z(), mMuEvent->runInfo().zdcCoincidenceRate());
	UShort_t cent9 = mRefmultCorrUtil->getCentralityBin9();
	UShort_t cent16 = mRefmultCorrUtil->getCentralityBin16();
    mRefmult->Fill(mRefmultCorrUtil->getRefMultCorr());

    const Double_t weight = mRefmultCorrUtil->getWeight();

    
    if(prevPass){
        mCentBin9All->Fill(cent9,weight);
        mCentBin16All->Fill(cent16,weight);
    }

    if( (cent9 < mCent9Cut[0]) || (cent9 > mCent9Cut[1]) )
    {
       cent9Pass = kFALSE; 
    }        

    if( (cent16 < mCent16Cut[0]) || (cent16 > mCent16Cut[1]) )
    {
       cent16Pass = kFALSE; 
    }        
        
    if(cent9Pass) {
        mCutIndependent->Fill(6); 
        if(prevPass) {mCentBin9Pass->Fill(cent9,weight);}
    }

    if(cent16Pass) {
        mCutIndependent->Fill(7);
        if(prevPass) {mCentBin16Pass->Fill(cent16,weight);}
    }

    if(cent9Pass && prevPass) { mCutCumulative->Fill(6); }
    if(cent9Pass && cent16Pass && prevPass) { mCutCumulative->Fill(7); }

    return (cent9Pass && cent16Pass);
}

Bool_t EventCutMaker::checkZdc(StMuEvent* event, Bool_t prevPass) {

    Bool_t zdcPass = kTRUE;
    Int_t zdcE = mMuEvent->zdcTriggerDetector().adc(4);
    Int_t zdcW = mMuEvent->zdcTriggerDetector().adc(0);

    if( (zdcE > mZdcCut) || (zdcW > mZdcCut)) 
    {
        zdcPass = kFALSE;
    }

    if(zdcPass ) { mCutIndependent->Fill(8); }
    if(zdcPass && prevPass) { mCutCumulative->Fill(8); }
    return zdcPass;
}

Int_t EventCutMaker::nTofMatchedTracks(StMuDst* muDst){

	Int_t nPrimary 	= muDst->primaryTracks()->GetEntries();
	Int_t nTofMatched = 0;
	for (int iNode = 0; iNode < nPrimary; iNode++ ){
		StMuTrack*	tPrimary 	= (StMuTrack*)muDst->primaryTracks(iNode);

		if ( !tPrimary ) continue;
		if ( tPrimary->vertexIndex() != 0 ) continue;

		StMuBTofPidTraits bTofPidTraits	= tPrimary->btofPidTraits();

		if ( bTofPidTraits.matchFlag() > 0 ) 
			nTofMatched++;
	}

	return nTofMatched;

}
