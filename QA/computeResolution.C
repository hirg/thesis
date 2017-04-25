#include "src/azifemNamespace.cxx"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

Int_t getZdcBinLog(const Double_t zdc);
map<Int_t, Int_t> makeAuAuRunMap();
Int_t runYear(Int_t runId);
Int_t runDay(Int_t runId);

void computeResolution(	const TString inFile = "testPion.list", 
						const char* outFile = "testOut.root", 
                        const string configFile = "QA/computeResolution.config",
                        const Bool_t uuNotAuAu = kFALSE,
                        Int_t nEvents = 999999
                        )
{
	//------------------- Load Shared Libraries ------------------//

	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StFlowMaker");
	gSystem->Load("StHbtMaker");
	gSystem->Load("PionDstMaker");
	gSystem->Load("ConfigReader");

    ConfigReader config(configFile);
    const Int_t tofmultcut = config.getI("tofmultcut");
    const Int_t refmultcut = config.getI("refmultcut");
    const map<Int_t, Int_t> runMap = azifem::makeRunMap(uuNotAuAu);
    const Int_t nRuns = runMap.size();
    cout << runMap.size() << endl;

	//------------------- Set some variables ------------------//

	Float_t psiAMult = 0, psiBMult = 0, psiFullMult = 0;
	Float_t zdcE = 0, zdcW = 0, zdcHigher = 0;
	Float_t psiAq2 = 0, psiBq2 = 0, psiFullq2 = 0;
    Float_t resMult  = 0, resq2 = 0;
    Float_t res = 0, q2 = 0;
    Float_t Vx = 0, Vy = 0, Vz = 0;
    Int_t refmultBin = 0, refmult = 0, tofmult = 0; 
    Int_t runId = 0;
    Event* event = 0;

	TBranch* eventBranch = 0;

    TChain* chain = new TChain("PionDst","PionDst");
    cout << "Added ";
    cout << fillChain(chain,inFile.Data(),1000);
    cout << " files to chain." << endl;
	chain->SetBranchAddress("Event",&event); 

	//------------------- Create Histrograms, etc. ------------------//

	TFile* fOut = new TFile( outFile, "RECREATE" );
    Double_t pi = TMath::Pi();

    // Resolution TProfile's
	TProfile2D* subResSquaredRefmultExclusiveZdcCut = new TProfile2D(
        "subResSquaredRefmultExclusiveZdcCut","subResSquaredRefmultExclusiveZdcCut",
        5, -0.5, 4.5, 2, -0.5, 1.5);
	TProfile2D* subResSquaredq2ExclusiveZdcCut = new TProfile2D(
        "subResSquaredq2ExclusiveZdcCut","subResSquaredq2ExclusiveZdcCut",
        5, -0.5, 4.5, 2, -0.5, 1.5);
	TProfile2D* subResSquaredRefmultInclusiveZdcCut = new TProfile2D(
        "subResSquaredRefmultInclusiveZdcCut","subResSquaredRefmultInclusiveZdcCut",
        5, -0.5, 4.5, 2, -0.5, 1.5);
	TProfile2D* subResSquaredq2InclusiveZdcCut = new TProfile2D(
        "subResSquaredq2InclusiveZdcCut","subResSquaredq2InclusiveZdcCut",
        5, -0.5, 4.5, 2, -0.5, 1.5);

    // Vs. Runnumber
    TProfile* avgSubresVsRunnumber = new TProfile(
        "avgSubresVsRunnumber", "avgSubresVsRunnumber", 
        nRuns, 0, nRuns);
    TProfile* avgRefmultVsRunnumber = new TProfile(
        "avgRefmultVsRunnumber", "avgRefmultVsRunnumber", 
        nRuns, 0, nRuns);
    TH1I* nEventsPerRunnumber = new TH1I(
        "nEventsPerRunnumber", "nEventsPerRunnumber", 
        nRuns, 0, nRuns);

    // Refmult/q2 Dists
	TH1F* refmultDist = new TH1F("refmultDist","Refmult Distribution",1000,-0.5,999.5);
	TH1F* refmultDistInclusiveZdcCut[4];
	TH1F* refmultDistExclusiveZdcCut[4];
    TH1F* refmultDeltaPsiDist = new TH1F(
        "refmultDeltaPsiDist","#Delta #Psi_{2} Distribution - refmult",
        200, -1*pi, pi);

    TH1F* q2Dist = new TH1F("q2Dist","q2 Distribution",1000,0,5);
    TH1F* q2DeltaPsiDist = new TH1F(
        "q2DeltaPsiDist","#Delta #Psi_{2} Distribution - q2",
        200, -1*pi, pi);
    TH1F* q2DistExclusiveZdcCut[4];
    TH1F* q2DistInclusiveZdcCut[4];

    TH1F* zVertex = new TH1F("zVertex","V_{z}",500,-60,60);
	TH2F* rVertex = new TH2F("rVertex","V_{y} vs. V_{x}",200,-2,2,200,-2,2);
	TH2F* refmultq2 = new TH2F("refmultq2","Refmult vs. q2",1000,0,5,1000,-0.5,999.5);
	TH2F* zdcEvsZdcW = new TH2F("zdcEvsZdcW","East ZDC vs. West ZDC",700,0,2100,700,0,2100);
    TH2F* zdcBinVsVz = new TH2F("zdcBinVsVz", "Zdc Bin vs. V_{z}", 60, -30, 30, 4, -0.5, 3.5);

    TH1F* q2ByZdcLog[6];
    TH1F* refmultByZdcLog[6];

    TH2F* tofmultVsRefmult = new TH2F("tofmultVsRefmult", "ToF multiplicity vs. refmult",1000,-0.5,999.5,1000,-0.5,3999.5 );

    for(Int_t j = 0; j <= 3; j++)
    {

        TString histName = TString::Format("refmultDistExclusiveZdcCut_%d", j);
        TString histTitle = TString::Format("Refmult Distribution - Zdc Bin %d (Exclusive ZDC cut)", j);
        refmultDistExclusiveZdcCut[j] = new TH1F(histName.Data(), histTitle.Data(),1000,-0.5,999.5);

        histName = TString::Format("refmultDistInclusiveZdcCut_%d", j);
        histTitle = TString::Format("Refmult Distribution - Zdc Bin %d (Inclusive ZDC cut)", j);
        refmultDistInclusiveZdcCut[j] = new TH1F(histName.Data(), histTitle.Data(),1000,-0.5,999.5);

        histName = TString::Format("q2DistExclusiveZdcCut_%d", j);
        histTitle = TString::Format("q2 Distribution - Zdc Bin %d (Exclusive ZDC cut)", j);
        q2DistExclusiveZdcCut[j] = new TH1F(histName.Data(), histTitle.Data(),1000,0,5);

        histName = TString::Format("q2DistInclusiveZdcCut_%d", j);
        histTitle = TString::Format("q2 Distribution - Zdc Bin %d (Inclusive ZDC cut)", j);
        q2DistInclusiveZdcCut[j] = new TH1F(histName.Data(), histTitle.Data(),1000,0,5);

    }

    for(Int_t j = 0; j <= 5; j++)
    {

        TString q2HistName = TString::Format("q2ByZdcLog_%d", j);
        TString multHistName = TString::Format("refmultByZdcLog_%d", j);
        q2ByZdcLog[j] = new TH1F(q2HistName.Data(), q2HistName.Data(), 1000, 0, 5);
        refmultByZdcLog[j] = new TH1F(multHistName.Data(), multHistName.Data(), 1000, -0.5, 999.5);
    }

	//------------------- Event Loop ------------------//
    if(nEvents == 0) { nEvents = chain->GetEntries();}

	for (Long64_t i=0; i < nEvents; i++)
    {
        // event->Clear();
        Int_t makeReturn = chain->GetEntry(i);
		if (!makeReturn)
        {
			cout << "Out of events!" << endl;
			break;
		}


        zdcE = event->GetZDCe();
        zdcW = event->GetZDCw();
        zdcHigher = (zdcE > zdcW) ? zdcE : zdcW;
        Int_t zdcBinLinear = azifem::getZdcBin(zdcHigher, uuNotAuAu);
        Int_t zdcBinLog = getZdcBinLog(zdcHigher);

        Vx = event->GetVertexPos().x();
        Vy = event->GetVertexPos().y();
        Vz = event->GetVertexPos().z();

        psiFullMult = event->GetPsi(0,0);
        psiAMult = event->GetPsi(0,1);
        psiBMult = event->GetPsi(0,2);
        psiFullq2 = event->GetPsi(1,0);
        psiAq2 = event->GetPsi(1,1);
        psiBq2 = event->GetPsi(1,2);

		resMult = cos(2*(psiAMult - psiBMult));
		resq2 = cos(2*(psiAq2 - psiBq2));
        refmultDeltaPsiDist->Fill(psiAMult - psiBMult);
        q2DeltaPsiDist->Fill(psiAq2 - psiBq2);

        refmult = event->GetRefMultCorr();
        tofmult = event->GetTofMult();
        q2 = event->Getq2(1);
        runId = event->GetRunID();

        zVertex->Fill(Vz);
        rVertex->Fill(Vx,Vy);


        refmultDist->Fill(refmult);
        tofmultVsRefmult->Fill(refmult, tofmult);

        if(refmult >= 200)
        {
            if(zdcBinLinear >= 0) {
                refmultDistInclusiveZdcCut[zdcBinLinear]->Fill(refmult);
                q2DistInclusiveZdcCut[zdcBinLinear]->Fill(q2);
                refmultDistExclusiveZdcCut[zdcBinLinear]->Fill(refmult);
                q2DistExclusiveZdcCut[zdcBinLinear]->Fill(q2);

                subResSquaredRefmultInclusiveZdcCut->Fill(azifem::getRefmultBin(refmult, uuNotAuAu), zdcBinLinear, resMult);
                subResSquaredq2InclusiveZdcCut->Fill(azifem::getq2Bin(q2), zdcBinLinear, resq2);
                subResSquaredRefmultExclusiveZdcCut->Fill(azifem::getRefmultBin(refmult, uuNotAuAu), zdcBinLinear, resMult);
                subResSquaredq2ExclusiveZdcCut->Fill(azifem::getq2Bin(q2), zdcBinLinear, resq2);
            }

            if(zdcBinLinear == 0) {
                refmultDistInclusiveZdcCut[1]->Fill(refmult);
                q2DistInclusiveZdcCut[1]->Fill(q2);

                subResSquaredRefmultInclusiveZdcCut->Fill(azifem::getRefmultBin(refmult, uuNotAuAu), 1, resMult);
                subResSquaredq2InclusiveZdcCut->Fill(azifem::getq2Bin(q2), 1, resq2);
            }


            q2Dist->Fill(q2);
            q2ByZdcLog[zdcBinLog]->Fill(q2);
            refmultByZdcLog[zdcBinLog]->Fill(refmult);
        } // If refmult >= 200

        refmultq2->Fill(q2,refmult);
        zdcEvsZdcW->Fill(zdcW,zdcE);
        zdcBinVsVz->Fill(Vz, zdcBinLinear);


        avgRefmultVsRunnumber->Fill(runMap[runId], refmult);
        avgSubresVsRunnumber->Fill(runMap[runId], resMult);
        nEventsPerRunnumber->Fill(runMap[runId]);

	} // loop over events

	fOut->Write();
	fOut->Close();
}

Int_t fillChain(TChain* chain, const char* fileList, const Int_t maxFiles) {

    ifstream inputStream(fileList);

    char* temp;
    int count=0;
      
    while(inputStream.good())
    {
        temp = new char[200];
        inputStream.getline(temp,200);

        TString fileName(temp);
        if(fileName.Contains("root")) {
            chain->Add(temp);
            ++count;
        }
        delete temp;
        if (count>maxFiles) break;
    }   

  return count;
}

// Given a zdc value, return a logarithmically spaced bin number.
// Bins: 0-0.01%, 0-0.03%, 0-0.0625%, 0-0.125%, 0-0.25%, 0-0.5%, 0-1%
// Indices: 0, 1, 2, 3, 4, 5, 6
// Warning!!: These values are probably outdated!!
Int_t getZdcBinLog(const Double_t zdc)
{

    const Double_t zdcBinEdges[7] = {-1,274.5, 317.5, 370.5,437.5, 525.5, 10000}; // UU 193
    // const Double_t zdcBinEdges[7] = {-1,730.5,871.5,1027.5,1198.5,1390.5, 10000}; // AuAu200
    Int_t zdcBin = -1;

    for(Int_t i = 0; i <= 5; i++)
    {
        if( (zdc > zdcBinEdges[i]) && (zdc <= zdcBinEdges[i+1]) ) {zdcBin = i;}

    }

    return zdcBin;
}

// Given a runId of the form 12XXXXXX return runYear=12
Int_t runYear(Int_t runId)
{ return y = Int_t(runId / 1000000); }

// Given a runId of the form XX123XXX return runDay=123
Int_t runDay(Int_t runId)
{
    Int_t dayAndRun = runId % 1000000;
    Int_t day = Int_t (dayAndRun / 1000);

    return day;
}
