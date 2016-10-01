#include "EventCut.h"

ClassImp(EventCut)

EventCut::EventCut(TTree* tree)
{

    mTree = tree;

    mVxCut[0] = -999;
    mVxCut[1] = 999;
    mVyCut[0] = -999;
    mVyCut[1] = 999;
    mVzCut[0] = -999;
    mVzCut[1] = 999;
    mVrCut = 999;

    mRefmultCut[0] = -1;
    mRefmultCut[1] = 9999;
    mCent9Cut[0] = -1;
    mCent9Cut[1] = 99;
    mCent16Cut[0] = -1;
    mCent16Cut[1] = 99;

    mVrCenter[0] = 0.;
    mVrCenter[1] = 0.;

    mNTriggers = 0.;

    mEventList = new TEventList();

}

TEventList* EventCut::makeEventList()
{

    if(!mTree)
    {
        cout << "\n***** EventCut::makeEventList() *****\n";
        cout << "***** Warning! *****\n";
        cout << "EventCut has not been given a TTree/TChain\n";
        cout << "Finishing without doing anything.\n\n";

        return 0;
    }

    bookHistograms();
	//----------------- Read data using trees ---------------//

	TStopwatch* timer = new TStopwatch();
    Long64_t nTotal = mTree->GetEntries(); 

    mTree->SetBranchStatus("*",0);
    mTree->SetBranchStatus("MuEvent.mZdcTriggerDetector.mAdc*",1);
    mTree->SetBranchStatus("MuEvent.mTriggerIdCollection.mL1TriggerId.mId*",1);
    mTree->SetBranchStatus("PrimaryVertices.mPosition.mX1*",1);
    mTree->SetBranchStatus("PrimaryVertices.mPosition.mX2*",1);
    mTree->SetBranchStatus("PrimaryVertices.mPosition.mX3*",1);
    mTree->SetBranchStatus("MuEvent.mRefMultNeg*",1);
    mTree->SetBranchStatus("MuEvent.mRefMultPos*",1);
    mTree->SetBranchStatus("MuEvent.mRunInfo.mZdcCoincidenceRate*",1);
	Double_t overheadTime = timer->RealTime();

	timer->Start();

//    nTotal = 1000;
    for(Int_t i = 0; i <= (nTotal-1); i++)
    {
        mTree->GetEntry(i);
        mZdc[0] = mTree->GetLeaf("MuEvent.mZdcTriggerDetector.mAdc")->GetValue(4);
        mZdc[1] = mTree->GetLeaf("MuEvent.mZdcTriggerDetector.mAdc")->GetValue(0);
        mVx = mTree->GetLeaf("PrimaryVertices.mPosition.mX1")->GetValue();
        mVy = mTree->GetLeaf("PrimaryVertices.mPosition.mX2")->GetValue();
        mVz = mTree->GetLeaf("PrimaryVertices.mPosition.mX3")->GetValue();
        mRefmult = mTree->GetLeaf("MuEvent.mRefMultNeg")->GetValue() + mTree->GetLeaf("MuEvent.mRefMultPos")->GetValue();

        mVr = sqrt(mVx*mVx + mVy*mVy);

        mVzPass = checkVz();
        mVrPass = checkVr();
        mHasVertexPass = checkHasVertex();
        mZdcPass = checkZdc();
        mTriggerPass = checkTrigger();

        mEventPass = mVzPass && mVrPass && mHasVertexPass && mZdcPass && mTriggerPass;
        fillHistograms();

        if( mEventPass ) { mEventList->Enter(i); }

    }

	Double_t eventLoopTime = timer->RealTime();
    Long64_t nAdded = mEventList->GetN();
    Double_t percentAdded = (100. * nAdded) / (Double_t)nTotal;
    cout << "\n\n***** EventCut::makeEventList() Finished *****\n";
    cout << "Overhead time: " << overheadTime << "  Event loop time: " << eventLoopTime << endl;
    cout << "Processed " << nTotal << " entries. Ran at " << nTotal / eventLoopTime << " entries/sec.\n";
    cout << "Added " << nAdded << " entries (" << percentAdded << "%) to the event list.\n\n";

    delete timer;
    mTree->SetBranchStatus("*",1);

    writeHistograms();

    return mEventList;
}

void EventCut::SetTriggers(const Float_t* triggers, Int_t nTriggers)
{

    mNTriggers = nTriggers;

    if( mNTriggers > 8 )
    {
        cout << "\n***** EventCut::SetTriggers() *****\n";
        cout << "***** Warning! *****\n";
        cout << "EventCut will only use first eight (8) triggers you've specified.\n\n";

        mNTriggers = 8;
    }

    for(Int_t i = 0; i <= (mNTriggers - 1); i++)
    {
        mTriggerCut[i] = triggers[i];
    }
}

Bool_t EventCut::checkTrigger()
{

    Bool_t pass = kFALSE;
    Int_t triggerIndex = 0;
    Float_t eventTrigger = 1;

    do // Each trigger collection is 64 spaces long, but most of those are zero. This only checks the non-zero ones.
    {
        eventTrigger = mTree->GetLeaf("MuEvent.mTriggerIdCollection.mL1TriggerId.mId")->GetValue(triggerIndex);

        // Compare the trigger at index 'triggerIndex' to all the triggers in 'mTriggerCut[]`
        for(Int_t i = 0; i <= (mNTriggers - 1); i++)
        {
            if( mTriggerCut[i] == eventTrigger )
            {
                mTrigger = i;
                pass = kTRUE;        
                break;
            }
        }

        triggerIndex++;

    } while ( (eventTrigger != 0) && (pass == kFALSE) && (triggerIndex <= 63)) ;

    return pass;
}

void EventCut::bookHistograms()
{
    mAllCutsIndependentHist = new TH1I("mAllCutsIndependent", "Cuts - Independent", 7, 0.5, 7.5);
    mAllCutsCumulativeHist = new TH1I("mAllCutsCumulative", "Cuts - Cumulative", 7, 0.5, 7.5);

    mTriggerPassHist = new TH1F("mTriggerPass", "Triggers", mNTriggers, 0.5, (mNTriggers - 0.5));
    mRefmultAllHist = new TH1I("mRefmultAll", "Refmult - All events", 1000, -0.5, 999.5);
    mRefmultPassHist = new TH1I("mRefmultPass", "Refmult - Pass events", 1000, -0.5, 999.5);
    mZdcAllHist = new TH2F("mZdcAll", "Zdc: West vs. East - All events", 4050, 0.5, 4050.5, 4050, 0.5, 4050.5);
    mZdcPassHist = new TH2F("mZdcPass", "Zdc: West vs. East - Pass Cut", 4050, 0.5, 4050.5, 4050, 0.5, 4050.5);

    mVrAllHist = new TH2F("mVrAll", "V_{R} - All events", 300, -3, 3, 300, -3, 3);
    mVrPassHist = new TH2F("mVrPass", "V_{R} - Pass events", 300, -3, 3, 300, -3, 3);
    mVzAllHist = new TH1F("mVzAll", "V_{Z} - All events", 1000, -100, 100);
    mVzPassHist = new TH1F("mVzPass", "V_{Z} - Pass events", 1000, -100, 100);

    const Char_t* binTitles[6] = {"All Events", "Trigger", "Has Vertex", "V_{Z}", "V_{R}", "Zdc Cut"};

    for(Int_t i = 0; i <= 5; i++)
    {
        mAllCutsIndependentHist->GetXaxis()->SetBinLabel(i+1,binTitles[i]); 
        mAllCutsCumulativeHist->GetXaxis()->SetBinLabel(i+1,binTitles[i]); 
    }

    for(Int_t i = 0; i <= (mNTriggers - 1); i++)
    {
        TString label = TString::Itoa((Int_t)mTriggerCut[i],10);
        mTriggerPassHist->GetXaxis()->SetBinLabel(i+1,label.Data()); 
    }
}


void EventCut::fillHistograms()
{
    Bool_t pass = kTRUE;

    mAllCutsIndependentHist->Fill(1);
    mAllCutsCumulativeHist->Fill(1);

    pass = mTriggerPass && pass;
    if(mTriggerPass) { mAllCutsIndependentHist->Fill(2); }
    if(pass) { mAllCutsCumulativeHist->Fill(2); }

    pass = mHasVertexPass && pass;
    if(mHasVertexPass) { mAllCutsIndependentHist->Fill(3); }
    if(pass) { mAllCutsCumulativeHist->Fill(3); }

    pass = mVzPass && pass;
    if(mVzPass) { mAllCutsIndependentHist->Fill(4); }
    if(pass) { mAllCutsCumulativeHist->Fill(4); }

    pass = mVrPass && pass;
    if(mVrPass) { mAllCutsIndependentHist->Fill(5); }
    if(pass) { mAllCutsCumulativeHist->Fill(5); }

    pass = mZdcPass && pass;
    if(mZdcPass) { mAllCutsIndependentHist->Fill(6); }
    if(pass) { mAllCutsCumulativeHist->Fill(6); }

    mVzAllHist->Fill(mVz);
    mVrAllHist->Fill(mVx, mVy);
    mRefmultAllHist->Fill(mRefmult);
    mZdcAllHist->Fill(mZdc[0],mZdc[1]);

    if(mEventPass)
    {
        mVzPassHist->Fill(mVz); 
        mVrPassHist->Fill(mVx, mVy);
        mRefmultPassHist->Fill(mRefmult);
        mZdcPassHist->Fill(mZdc[0],mZdc[1]);
        mTriggerPassHist->Fill(mTrigger); 
    }

}

void EventCut::writeHistograms()
{
    TFile outputFile("eventCuts.root", "RECREATE");

    mAllCutsIndependentHist->Write();
    mAllCutsCumulativeHist->Write();
    mTriggerPassHist->Write();
    mRefmultAllHist->Write();
    mRefmultPassHist->Write();
    mVzAllHist->Write();
    mVzPassHist->Write();
    mVrAllHist->Write();
    mVrPassHist->Write();
    mZdcAllHist->Write();
    mZdcPassHist->Write();

    outputFile.Close();

}
