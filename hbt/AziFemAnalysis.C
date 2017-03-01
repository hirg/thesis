// This macro performs an azimuthatlly sensitive femtoscopy analysis over
// pionDst files. In general, it defines/instantiates things from "smaller
// to larger", i.e., defined constants first, then event/track/pair cut
// objects, then Hbt analysis objects, then chains and makers.

#include <vector>

void AziFemAnalysis(const TString fileList = "hbt/AuAupions.list",
					   const TString outFile = "aziFemTestOut.root",
                       const Float_t vzLow = 999,
                       const Float_t vzHigh = -999,
					   const Int_t q2Bin = -1,
					   const Int_t multBin = 0,
					   const Int_t zdcBin = 0,
                       const string configFile = "hbt/femto.config",
                       const string speciesString = "AuAu",
					   const Int_t nEvents = 999
                       )
{
	TStopwatch* timer = new TStopwatch();
	//------------------- Load Shared Libraries ------------------//

	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StFlowMaker");
	gSystem->Load("StHbtMaker");
	gSystem->Load("PionDstMaker");
	gSystem->Load("ConfigReader");

    ConfigReader config(configFile);

	//------------------- Define Cut Values ------------------//

    // Event cuts
    const vector<Float_t> mult = config.getVF("mult", speciesString);
    const vector<Float_t> zdc = config.getVF("zdc", speciesString);
    const vector<Float_t> q2 = config.getVF("q2");
    const Float_t vx = config.getF("vx");
    const Float_t vy = config.getF("vy");
    const Float_t vz = config.getF("vz");

    // Track Cuts
    const Float_t piMass =  config.getF("piMass");
    const vector<Float_t> rapidity = config.getVF("eta");
    const vector<Float_t> nSigma = config.getVF("nSigma");
    const vector<Float_t> pt = config.getVF("pt");
    const vector<Float_t> nHitsTpc = config.getVF("nHitsTpc");
    const vector<Float_t> dcaGlobal = config.getVF("dcaGlobal");

    // Pair Cuts
    const vector<Float_t> kt = config.getVF("kt");
    const vector<Float_t> quality = config.getVF("quality");
    const Float_t maxFracMergedRows = config.getF("maxFracMergedRows");

    // Hbt Analysis parameters
    const Int_t ptSwitch = config.getI("ptSwitch");
    const Int_t nEventsToMix = config.getI("nEventsToMix");
    const Int_t nQbins = config.getI("nQbins");
    const vector<Float_t> qRange = config.getVF("qRange");
    const Int_t nPhiBins = config.getI("nPhiBins");
    const Int_t nEPBins = config.getI("nEPBins");
    const vector<Float_t> rpRange = config.getVF("rpRange");
    const Int_t nMultBins = config.getI("nMultBins");
    const vector<Float_t> multRange = config.getVF("multRange");
    const Int_t nVzBins = config.getI("nVzBins");
    // const vector<Float_t> vzRange = config.getVF("vzRange");
    vector<Float_t> vzRange;
    vzRange.push_back(vzLow);
    vzRange.push_back(vzHigh);

	//------------------- Instantiate Cut Objects ------------------//

    // We instantiate two sets of event cuts so that each can have it's own cut monitor,
    // which avoids double counting. However, we still only need to have one set of monitors,
    // since they're the same event cuts.
    fxtEventCut* eventCut[2];
    fxtEventCutMonitor* eventPass = new fxtEventCutMonitor("eventPass",""); 
    fxtEventCutMonitor* eventFail = new fxtEventCutMonitor("eventFail","");
    
    for(Int_t i = 0; i <= 1; i++)
    {
        eventCut[i] = new fxtEventCut();
        
        //Set mult bin
        if(multBin == -1) { 
            eventCut[i]->SetMult(mult.at(0), mult.at( mult.size() - 1 ) );
            Int_t sel = 1;
        } else {
            eventCut[i]->SetMult(mult.at(multBin), mult.at(multBin+1));
            Int_t sel = 0;
        }

        // Set q2 bin
        if(q2Bin == -1) {
            eventCut[i]->Setq2(q2.at(0), q2.at( q2.size() - 1 ));
        } else {
            eventCut[i]->Setq2(q2.at(q2Bin), q2.at(q2Bin+1));
        }

        eventCut[i]->SetVx(-1*vx, vx);
        eventCut[i]->SetVy(-1*vy, vy);
        eventCut[i]->SetVz(-1*vz, vz);
        eventCut[i]->SetZdc(zdc.at(zdcBin), zdc.at(zdcBin+1));
    }

    eventCut[0]->AddCutMonitor(eventPass,eventFail); 

    // Two sets of track cuts, one for pi+ and on for pi-
    franksTrackCut* trackCut[2];
    fxtTrackCutMonitor* trackPiPlusPass = new fxtTrackCutMonitor("_PiPlusPass_",piMass );
    fxtTrackCutMonitor* trackPiPlusFail = new fxtTrackCutMonitor("_PiPlusFail",piMass );
    fxtTrackCutMonitor* trackPiMinusPass = new fxtTrackCutMonitor("_PiMinusPass_",piMass );
    fxtTrackCutMonitor* trackPiMinusFail = new fxtTrackCutMonitor("_PiMinusFail",piMass );

    // i = 0 --> pi-minus; i = 1 --> pi-plus
    for(Int_t i = 0; i <=1; i++)
    {
        trackCut[i] = new franksTrackCut;
        trackCut[i]->SetMass(piMass);
        trackCut[i]->SetRapidity(rapidity[0],rapidity[1]);
        trackCut[i]->SetNSigmaElectron(nSigma[0],nSigma[1]);
        trackCut[i]->SetNSigmaPion(nSigma[0],nSigma[1]);
        trackCut[i]->SetNSigmaKaon(nSigma[0],nSigma[1]);
        trackCut[i]->SetNSigmaProton(nSigma[0],nSigma[1]);
        trackCut[i]->SetPt(pt[0],pt[1]);
        trackCut[i]->SetNHits(nHitsTpc[0],nHitsTpc[1]);
        trackCut[i]->SetDCAGlobal(dcaGlobal[0],dcaGlobal[1]);
        Int_t charge = i ? 1 : -1;
        trackCut[i]->SetCharge(charge);
    }

    trackCut[0]->AddCutMonitor(trackPiMinusPass,trackPiMinusFail);
    trackCut[1]->AddCutMonitor(trackPiPlusPass,trackPiPlusFail);

    // Pair Cuts
    kTPairCut* ktCut = new kTPairCut();
    ktCut->SetkTRange(kt[0],kt[1]);

	qualityPairCut* qualityCut = new qualityPairCut;
	qualityCut->SetQualityCut(quality[0],quality[1]);

	HitMergingPairCut* mergeCut = new HitMergingPairCut;
	mergeCut->setDefaultFullFieldMergingPar();
	mergeCut->setMaxFracOfMergedRow(maxFracMergedRows);

    ManyPairCuts* pairCut[2];
    pairCut[0] = new ManyPairCuts();
    pairCut[0]->AddPairCut(ktCut);
    pairCut[0]->AddPairCut(qualityCut);
    pairCut[0]->AddPairCut(mergeCut);
    pairCut[1] = new ManyPairCuts();
    pairCut[1]->AddPairCut(ktCut);
    pairCut[1]->AddPairCut(qualityCut);
    pairCut[1]->AddPairCut(mergeCut);

    fxtPairCutMonitor* pairPiPlusPass = new fxtPairCutMonitor("_PiPlusPass_" );
    fxtPairCutMonitor* pairPiPlusFail = new fxtPairCutMonitor("_PiPlusFail" );
    fxtPairCutMonitor* pairPiMinusPass = new fxtPairCutMonitor("_PiMinusPass_" );
    fxtPairCutMonitor* pairPiMinusFail = new fxtPairCutMonitor("_PiMinusFail" );

    pairCut[0]->AddCutMonitor(pairPiMinusPass,pairPiMinusFail);
    pairCut[1]->AddCutMonitor(pairPiPlusPass,pairPiPlusFail);

	//------------------- Instantiate Hbt Analyses and Correlation Functions ------------------//

    StHbtReactionPlaneAnalysis* azifemAnalysis[2];
    QoslCMSCorrFctnkT* qOslCF[2];
	QoslCMSCorrFctnRPkT* qOslRPCF[2];
    for(Int_t i = 0; i <=1; i++)
    {
        azifemAnalysis[i] = new StHbtReactionPlaneAnalysis(ptSwitch,nEPBins,rpRange[0],rpRange[1],nMultBins,multRange[0],multRange[1],nVzBins,vzRange[0],vzRange[1]);
        azifemAnalysis[i]->SetEventCut(eventCut[i]);
        azifemAnalysis[i]->SetFirstParticleCut(trackCut[i]);
        azifemAnalysis[i]->SetSecondParticleCut(trackCut[i]);
        azifemAnalysis[i]->SetPairCut(pairCut[i]);
        azifemAnalysis[i]->SetNumEventsToMix(nEventsToMix);

        TString title = "";
        if(i==0) {title += "PiMinus";}
        if(i==1) {title += "PiPlus";}
        qOslCF[i] = new QoslCMSCorrFctnkT(title.Data(),nQbins,qRange[0],qRange[1],nQbins,qRange[0],qRange[1],nQbins,qRange[0],qRange[1]);
        qOslRPCF[i] = new QoslCMSCorrFctnRPkT(title.Data(),nQbins,qRange[0],qRange[1],nQbins,qRange[0],qRange[1],nQbins,qRange[0],qRange[1],nPhiBins);
        azifemAnalysis[i]->AddCorrFctn(qOslRPCF[i]);
        azifemAnalysis[i]->AddCorrFctn(qOslCF[i]);
    }


	//------------------- Instantiate Chain and Makers ------------------//

	StChain* chain = new StChain("StChain");
	chain->SetDebug(0);

	StHbtPionDstReader* reader = new StHbtPionDstReader("",fileList.Data(),"root",99);
    reader->setSel(sel);
	StHbtMaker* hbtMaker = new StHbtMaker;
	hbtMaker->HbtManager()->SetEventReader(reader);

    for(Int_t i = 0; i <=1; i++)
    {
        hbtMaker->HbtManager()->AddAnalysis(azifemAnalysis[i]);
    }
        
	//--------------------- The STAR chain Event loop ---------------------//

	chain->Init();

	Double_t overheadTime = timer->RealTime();
	timer->Continue();
	cout << "\n\n***** Starting Event Loop. " << overheadTime
         << " seconds elapsed so far. *****" << endl;
  
    Int_t nEventsProcessed = 0;
	for (Int_t i = 0; i < nEvents; i++)
    {

		chain->Clear();
		Int_t makeReturn = chain->Make(i);
        ++nEventsProcessed;

		if (makeReturn)
        {
			cout << "Out of events!" << endl;
			break;
		}
	} 

    chain->Finish();

	Double_t eventLoopTime = timer->RealTime();
	cout << endl << "***** Finished Event Loop. " << eventLoopTime 
         << " seconds to process " << nEventsProcessed << " events. "
         << (Double_t)nEventsProcessed / eventLoopTime << " ev/s. *****\n\n";

// 	//--------------------- Write Output ---------------------//
	TFile histoOutput(outFile.Data(),"recreate");

    histoOutput.mkdir("EventCuts");
    histoOutput.cd("EventCuts");
    eventPass->VertexYvsVertexX()->Write();
    eventPass->VertexZ()->Write();
    eventPass->RefMult()->Write();
    eventFail->VertexYvsVertexX()->Write();
    eventFail->VertexZ()->Write();
    eventFail->RefMult()->Write();

    histoOutput.cd("../");
    histoOutput.mkdir("TrackCuts");
    histoOutput.cd("TrackCuts");
    trackPiPlusPass->DCA()->Write();
    trackPiPlusPass->DCAGlobal()->Write();
    trackPiPlusPass->Nhits()->Write();
    trackPiPlusPass->Pt()->Write();
    trackPiPlusPass->NsigmaPion()->Write();
    trackPiPlusPass->ChiSqr()->Write();
    trackPiPlusFail->DCA()->Write();
    trackPiPlusFail->DCAGlobal()->Write();
    trackPiPlusFail->Nhits()->Write();
    trackPiPlusFail->Pt()->Write();
    trackPiPlusFail->NsigmaPion()->Write();
    trackPiPlusFail->ChiSqr()->Write();
    trackPiMinusPass->DCA()->Write();
    trackPiMinusPass->DCAGlobal()->Write();
    trackPiMinusPass->Nhits()->Write();
    trackPiMinusPass->Pt()->Write();
    trackPiMinusPass->NsigmaPion()->Write();
    trackPiMinusPass->ChiSqr()->Write();
    trackPiMinusFail->DCA()->Write();
    trackPiMinusFail->DCAGlobal()->Write();
    trackPiMinusFail->Nhits()->Write();
    trackPiMinusFail->Pt()->Write();
    trackPiMinusFail->NsigmaPion()->Write();
    trackPiMinusFail->ChiSqr()->Write();

    histoOutput.cd("../");
    histoOutput.mkdir("PairCuts");
    histoOutput.cd("PairCuts");
    pairPiPlusPass->Kt()->Write();
    pairPiPlusPass->FractionOfMergedRow()->Write();
    pairPiPlusPass->SplittingLevel()->Write();
    pairPiPlusFail->Kt()->Write();
    pairPiPlusFail->FractionOfMergedRow()->Write();
    pairPiPlusFail->SplittingLevel()->Write();
    pairPiMinusPass->Kt()->Write();
    pairPiMinusPass->FractionOfMergedRow()->Write();
    pairPiMinusPass->SplittingLevel()->Write();
    pairPiMinusFail->Kt()->Write();
    pairPiMinusFail->FractionOfMergedRow()->Write();
    pairPiMinusFail->SplittingLevel()->Write();

    
    histoOutput.cd("../");
    for(Int_t i = 0; i <=1; i++)
    {

        for(int j=0; j<=4; j++) { //bins 0,1,2,(3) are typical kt bins for this version. Bin 4 is the ktIntegrated histo.

            for(int k=0; k<nPhiBins; k++) { 

				qOslRPCF[i]->Numerator3D(k,j)->Write();
				qOslRPCF[i]->Denominator3D(k,j)->Write();
				qOslRPCF[i]->CoulHisto3D(k,j)->Write();
				// qOslRPCF[i]->QinvHisto3D(k,j)->Write();
            } // End loop over phi bins

            qOslCF[i]->Numerator3D(j)->Write();
            qOslCF[i]->Denominator3D(j)->Write();
            qOslCF[i]->CoulHisto3D(j)->Write();
        } // End loop over kt bins
    } // End loop over pi+ / pi-

    histoOutput.Close();

    delete chain;
  
}
