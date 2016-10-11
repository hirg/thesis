class StMuEvent;
class StMuDst;

Bool_t checkTrigger(StMuEvent* event, vector<int>* triggers); 

void QA(const TString fileList,
        const string configFile = "QA/19GeVQA.config",
        const TString outFile = "testOut.root",
		Int_t nEvents = 9999)
{

	TStopwatch* timer = new TStopwatch();

	//-----------------  Load Libraries ---------------//
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StRefMultCorr");
	gSystem->Load("ConfigReader");

    // --- Read config file --- // 
    ConfigReader cr(configFile);
    vector<int> triggers = cr.getVI("triggers", "193GeV");
    Float_t VzCut = cr.getF("vz");
    Float_t VrCut = cr.getF("vr");
    Int_t zdcCut = cr.getI("zdc");

	//----------------- Instantiate chain and MuDst reader ---------------//
	StChain* chain = new StChain("StChain");
	chain->SetDebug(0);
	StMuDebug::setLevel(0);

	StMuDstMaker* muMaker = new StMuDstMaker(0,0,"",fileList.Data(),
                                            "st:MuDst.root", 100, "MuDst");
    StRefMultCorr* refmultCorrUtil = new StRefMultCorr("refmult");

    muMaker->SetStatus("*",0);
    muMaker->SetStatus("MuEventAll",1);
    muMaker->SetStatus("BTofHeader",1);
    
    // muMaker->chain()->SetBranchStatus("*",0);
    // muMaker->chain()->SetBranchStatus("MuEvent*",1);
    // muMaker->chain()->SetBranchStatus("MuEvent*",1);
    
	//----------------------- Make Histograms, etc. -----------------------//
	TFile* fOut = new TFile(outFile.Data(), "RECREATE");

	TH2I* hZdcWestVsEast = new TH2I("zdcWestVsEast", "ZdcEast vs. ZdcWest",
            4050, -0.5, 4049.5, 4050, -0.5, 4049.5);
	TH2F* hVzDiffvsVzTPC = new TH2F("hVzDiffvsVzTPC",
            "V_{z,TPC} - V_{z,VPD} vs. V_{TPC}", 200, -50, 50, 200, -10, 10);
	TH1F* hTpcVz = new TH1F("hTpcVz", "V_{z,TPC}",200, -50, 50);
	TH1F* hVpdVz = new TH1F("hVpdVz", "V_{z,VPD}", 200, -50, 50);
	TH2F* hVyvsVx = new TH2F("hVyvsVx", "V_{y} vs. V_{x}", 
            100, -2.5, 2.5, 100, -2.5, 2.5);
	TH1I* hRefmult = new TH1I("hRefmult"," refmult", 800, -0.5, 800.5);
	TH1F* hRefmultCorr = new TH1F("hRefmultCorr","refmultCOrr", 800, -0.5, 799.5);
    TH2I* hTofmultVsRefmult = new TH2I("hTofmultVsRefmult",
            "ToF Multiplicity vs. Refmult", 800, -0.5, 799.5, 800, -0.5, 1599.5);

 
    TH1I* hAllCutsCumulative = new TH1I("hAllCutsCumulative",
            "Cuts - Cumulative", 7, 0.5, 7.5);
    const Char_t* binTitles[7] = {"All Events", "Good Run", "Trigger",
            "Has Vertex", "V_{Z,TPC}", "V_{R}", "1% Zdc"};

    // Set bin labels on cuts histogram
    for(Int_t i = 0; i <= 6; i++)
    { hAllCutsCumulative->GetXaxis()->SetBinLabel(i+1,binTitles[i]); }
	//--------------------- The STAR chain Event loop ---------------------//
	chain->Init();

	// Start  bookeeping stuff. Timer, event counter, etc.
	Int_t iReturn = 0;
    Int_t lastRunId = -1;
	Int_t percentCounter = 1;
	Int_t nEventsProcessed = 0;
    Float_t nEventsPerSecond = 0;
    Int_t totalEntries = muMaker->chain()->GetEntries();
    if( totalEntries < nEvents ) { nEvents = totalEntries; }

	// ------------ Event Loop --------------------- //
	Double_t overheadTime = timer->RealTime();
	timer->Continue();
	cout << "\n\n***** Starting Event Loop. " << overheadTime
         << " seconds elapsed so far. *****" << endl;

	for (Int_t iev = 0; iev < nEvents; iev++) {
        Bool_t goodRun = kFALSE, goodTrigger = kFALSE, goodZdc = kFALSE,
                goodVz = kFALSE, goodDeltaVz = kFALSE, goodVr = kFALSE;

		chain->Clear();
		iReturn = chain->Make(iev); 

		if (iReturn) {
            cout << "Return code: " << iReturn << endl;
			cout << "Exiting before we expected to!" << endl;
			break;
		}

		++nEventsProcessed;
        
		// Assign event quantities
		StMuEvent* event = muMaker->muDst()->event();
        Int_t runId = event->runId();
		Float_t Vx = event->primaryVertexPosition().x();
		Float_t Vy = event->primaryVertexPosition().y();
        Float_t Vr = sqrt( Vx*Vx + Vy*Vy );
		Float_t Vz = event->primaryVertexPosition().z();
		Float_t vpdVz = muMaker->muDst()->btofHeader()->vpdVz();
		Int_t zdcW = event->zdcTriggerDetector().adc(0);
		Int_t zdcE = event->zdcTriggerDetector().adc(4);
        Int_t zdcRate = event->runInfo().zdcCoincidenceRate();
        Int_t refmult = event->refMult();
        Int_t tofMult = event->btofTrayMultiplicity();

        hAllCutsCumulative->Fill(1); // All events

        //Check for good events
        if( lastRunId != runId ) {
            lastRunId = runId;
            refmultCorrUtil->init(runId);
        }

        // Is this a good run?
        goodRun = !(refmultCorrUtil->isBadRun(runId));
        if( goodRun ) { hAllCutsCumulative->Fill(2); }
        else { continue; }

        // Does the event have the right trigger?
        if( checkTrigger(event, &triggers) ) { hAllCutsCumulative->Fill(3); }
        else { continue; }

        // Does the event have a vertex?
        if( Vz != 0.0 ) { hAllCutsCumulative->Fill(4); }
        else { continue; }

        // Does the event pass the Vz cut?
        if( abs(Vz) <= VzCut ) { hAllCutsCumulative->Fill(5); }
        else { continue; }
        
        // Does the event pass the Vr cut?
        if( Vr <= VrCut ) { hAllCutsCumulative->Fill(6); }
        else { continue; }

        // Does the event pass the Zdc cut?
        goodZdc = (zdcE <= zdcCut) && (zdcW <= zdcCut);
        if(goodZdc) { hAllCutsCumulative->Fill(7); }
        else { continue; }

        refmultCorrUtil->initEvent(refmult, Vz, zdcRate);
        Double_t refmultCorr = refmultCorrUtil->getRefMultCorr();

		// Fill histograms
        hZdcWestVsEast->Fill(zdcE, zdcW);
        hVzDiffvsVzTPC->Fill( Vz, (Vz - vpdVz));
        hTpcVz->Fill(Vz);
        hVpdVz->Fill(vpdVz);
        hVyvsVx->Fill(Vx, Vy);
        hRefmult->Fill(refmult);
        hRefmultCorr->Fill(refmultCorr);
        hTofmultVsRefmult->Fill(refmult, tofMult);

        checkProgress(iev, nEvents, &percentCounter, timer);


	}  // Event loop

	Double_t eventLoopTime = timer->RealTime();
	cout << endl << "***** Finished Event Loop. " << eventLoopTime 
         << " seconds to process " << nEventsProcessed << " events. "
         << (Double_t)nEventsProcessed / eventLoopTime << " ev/s. *****\n\n";
	
	chain->Finish();

    fOut->Write();
    fOut->Close();

	delete chain;

}

Bool_t checkTrigger(StMuEvent* event, vector<int>* triggers)
{
	Bool_t hasTrigger = kFALSE;
    
    for(unsigned i = 0; i <= (triggers->size() - 1); i++)
    {
        if( event->triggerIdCollection().nominal().isTrigger(triggers->at(i)) )
        {
            hasTrigger = kTRUE;
        }
    }
	
	return hasTrigger;
}

void checkProgress(Int_t iev, Int_t nEvents, Int_t* percentCounter, TStopwatch* timer)
{
    Float_t progress = (Float_t)(iev * 100) / (Float_t)nEvents;	

    if(progress >= *percentCounter){
        Double_t time = timer->RealTime();
        timer->Continue();
        Float_t nEventsPerSecond = iev / time;

        cout << "\r" << *percentCounter << "% done. "
             << "Processing " << nEventsPerSecond << " ev/s.      "
             << flush;
        *percentCounter += 1;
    }
}
