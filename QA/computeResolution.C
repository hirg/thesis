#include "../src/azifemNamespace_core.cxx"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

Int_t getCentBin(const Int_t rm);
Int_t getRefmultBin(const Double_t rm);
Int_t getq2Bin(const Double_t q2);
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

    TProfile* avgSubresVsRunnumber = new TProfile(
        "avgSubresVsRunnumber", "avgSubresVsRunnumber", 
        nRuns, 0, nRuns);
    TProfile* avgRefmultVsRunnumber = new TProfile(
        "avgRefmultVsRunnumber", "avgRefmultVsRunnumber", 
        nRuns, 0, nRuns);
    TH1I* nEventsPerRunnumber = new TH1I(
        "nEventsPerRunnumber", "nEventsPerRunnumber", 
        nRuns, 0, nRuns);
	TH1F* refmultDist = new TH1F("refmultDist","Refmult Distribution",1000,-0.5,999.5);
    TH1F* q2Dist = new TH1F("q2Dist","q2 Distribution",1000,0,5);
    TH1F* refmultDeltaPsiDist = new TH1F(
        "refmultDeltaPsiDist","#Delta #Psi_{2} Distribution - refmult",
        200, -1*pi, pi);
    TH1F* q2DeltaPsiDist = new TH1F(
        "q2DeltaPsiDist","#Delta #Psi_{2} Distribution - q2",
        200, -1*pi, pi);
	TH1F* refmultDistExclusiveZdcCut[2];
    TH1F* q2DistExclusiveZdcCut[2];
	TH1F* refmultDistInclusiveZdcCut[2];
    TH1F* q2DistInclusiveZdcCut[2];

	TH2F* rVertex = new TH2F("rVertex","V_{y} vs. V_{x}",200,-2,2,200,-2,2);
    TH1F* zVertex = new TH1F("zVertex","V_{z}",500,-60,60);
	TH2F* refmultq2 = new TH2F("refmultq2","Refmult vs. q2",1000,0,5,1000,-0.5,999.5);
	TH2F* zdcEvsZdcW = new TH2F("zdcEvsZdcW","East ZDC vs. West ZDC",700,0,2100,700,0,2100);

    TH1F* q2ByZdcLog[6];
    TH1F* refmultByZdcLog[6];

    TH2F* tofmultVsRefmult = new TH2F("tofmultVsRefmult", "ToF multiplicity vs. refmult",1000,-0.5,999.5,1000,-0.5,3999.5 );

    for(Int_t j = 0; j <= 2; j++)
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
        Int_t zdcBinLinear = getZdcBin(zdcHigher);
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

                subResSquaredRefmultInclusiveZdcCut->Fill(getRefmultBin(refmult), zdcBinLinear, resMult);
                subResSquaredq2InclusiveZdcCut->Fill(getq2Bin(q2), zdcBinLinear, resq2);
                subResSquaredRefmultExclusiveZdcCut->Fill(getRefmultBin(refmult), zdcBinLinear, resMult);
                subResSquaredq2ExclusiveZdcCut->Fill(getq2Bin(q2), zdcBinLinear, resq2);
            }

            if(zdcBinLinear == 0) {
                refmultDistInclusiveZdcCut[1]->Fill(refmult);
                q2DistInclusiveZdcCut[1]->Fill(q2);

                subResSquaredRefmultInclusiveZdcCut->Fill(getRefmultBin(refmult), 1, resMult);
                subResSquaredq2InclusiveZdcCut->Fill(getq2Bin(q2), 1, resq2);
            }


            q2Dist->Fill(q2);
            q2ByZdcLog[zdcBinLog]->Fill(q2);
            refmultByZdcLog[zdcBinLog]->Fill(refmult);
        } // If refmult >= 200

        refmultq2->Fill(q2,refmult);
        zdcEvsZdcW->Fill(zdcW,zdcE);


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

Int_t getCentBin(const Int_t rm)
{

    if((rm==0) || (rm==1)) { return -1; }
    else if ((rm==2) || (rm==3)) { return 8; }
    else if ((rm==4) || (rm==5)) { return 7; }
    else if ((rm==6) || (rm==7)) { return 6; }
    else if ((rm==8) || (rm==9)) { return 5; }
    else if ((rm==10) || (rm==11)) { return 4; }
    else if ((rm==12) || (rm==13)) { return 3; }
    else if (rm==14) { return 2; }
    else if (rm==15) { return 1; }
    else {return 0;}

}

Int_t getRefmultBin(const Double_t rm)
{

    const Float_t mult[6] = {200, 570, 589, 605, 623, 10000}; // UU 193
    // const Float_t mult[6] = {200, 501, 522, 538, 556, 10000}; // AuAu 200
    Int_t rmBin = -1;

    for(Int_t i = 0; i <= 5; i++)
    {
        if( (rm > mult[i]) && (rm <= mult[i+1]) ) {rmBin = i;}

    }

    return rmBin;

}

Int_t getq2Bin(const Double_t q2)
{

    const Double_t q2BinEdges[6] = {0,0.4925,0.7525,1.0075,1.3275,10}; 
    Int_t q2Bin = -1;

    for(Int_t i = 0; i <= 5; i++)
    {
        if( (q2 > q2BinEdges[i]) && (q2 <= q2BinEdges[i+1]) ) {q2Bin = i;}

    }

    return q2Bin;

}

Int_t getZdcBin(const Double_t zdc)
{

    const Double_t zdcBinEdges[3] = {-1,436.5,523.5}; // UU 193 - 0, 0.25, 0.5
    // const Double_t zdcBinEdges[3] = {-1, 1194.5, 1389.5}; // AuAu 200 - 0, 0.25, 0.5
    Int_t zdcBin = -1;

    for(Int_t i = 0; i <= 2; i++)
    {
        if( (zdc > zdcBinEdges[i]) && (zdc <= zdcBinEdges[i+1]) ) {zdcBin = i;}

    }

    return zdcBin;
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

// Return a map that can index runnumbers. Call with runIndex = runMap[runId]
map<Int_t, Int_t> makeAuAuRunMap()
{
    const Int_t runIds[] = {12126101, 12126104, 12126105, 12126106, 12126107, 12126108, 12127002, 12127003, 12127004, 12127005, 12127006, 12127008, 12127013, 12127017, 12127018, 12127019, 12127020, 12127021, 12127022, 12127023, 12127024, 12127030, 12127031, 12127032, 12127033, 12127034, 12127035, 12127038, 12127039, 12127040, 12127041, 12127042, 12127048, 12127049, 12127050, 12128001, 12128007, 12128008, 12128009, 12128010, 12128011, 12128012, 12128014, 12128015, 12128017, 12128018, 12128019, 12128020, 12128021, 12128024, 12128025, 12128029, 12128030, 12128032, 12128038, 12132008, 12132009, 12132010, 12132011, 12132013, 12132014, 12132017, 12132018, 12132019, 12132020, 12132021, 12132022, 12132023, 12132024, 12132025, 12132026, 12132032, 12132033, 12132034, 12132043, 12132044, 12132045, 12132046, 12132047, 12132048, 12132049, 12132050, 12132051, 12132052, 12132053, 12132054, 12132055, 12132056, 12132057, 12132058, 12132061, 12132062, 12132063, 12132064, 12132065, 12132066, 12132067, 12132068, 12132069, 12132070, 12132071, 12132072, 12133005, 12133006, 12133007, 12133008, 12133010, 12133011, 12133012, 12133013, 12133014, 12133016, 12133018, 12133019, 12133020, 12133021, 12133022, 12133023, 12133024, 12133025, 12133026, 12133027, 12133028, 12133038, 12133039, 12133040, 12133041, 12133052, 12133053, 12133054, 12133056, 12133057, 12133058, 12133059, 12133060, 12133061, 12134001, 12134002, 12134003, 12134005, 12134006, 12134007, 12134008, 12134009, 12134010, 12134011, 12134012, 12134013, 12134014, 12134015, 12134016, 12134017, 12134018, 12134023, 12134026, 12134028, 12134030, 12134031, 12134032, 12134033, 12134034, 12134035, 12134036, 12134037, 12134038, 12134040, 12134041, 12134044, 12134045, 12134046, 12134047, 12134048, 12134049, 12134050, 12134051, 12134052, 12134053, 12134055, 12134056, 12134057, 12134058, 12134059, 12134061, 12134062, 12134063, 12134064, 12134065, 12134066, 12134067, 12134068, 12135002, 12135003, 12135004, 12135005, 12135006, 12135007, 12135008, 12135009, 12135010, 12135011, 12135012, 12135013, 12135014, 12135019, 12135020, 12135021, 12135022, 12135023, 12135024, 12135025, 12135026, 12135027, 12135028, 12135029, 12135030, 12135031, 12135033, 12135034, 12135036, 12135037, 12135038, 12135039, 12135041, 12135042, 12135043, 12135044, 12135045, 12135046, 12135048, 12135049, 12135050, 12135051, 12135052, 12135053, 12135054, 12135055, 12135056, 12135057, 12135058, 12135059, 12135060, 12135061, 12136001, 12136002, 12136005, 12136006, 12136007, 12136014, 12136017, 12136018, 12136022, 12136023, 12136024, 12136025, 12136026, 12136027, 12136028, 12136029, 12136030, 12136031, 12136032, 12136034, 12136039, 12136040, 12136041, 12136043, 12136044, 12136045, 12136046, 12136047, 12136048, 12136050, 12136052, 12136053, 12136054, 12136055, 12136064, 12136065, 12136069, 12136070, 12136071, 12136074, 12136075, 12136076, 12136078, 12136079, 12136080, 12136081, 12136082, 12136083, 12136084, 12136085, 12136086, 12137003, 12137004, 12137005, 12137006, 12137007, 12137008, 12137009, 12137010, 12137011, 12137012, 12137013, 12137014, 12137015, 12137020, 12137021, 12137022, 12137023, 12137024, 12137025, 12137026, 12137027, 12137028, 12137029, 12137030, 12137033, 12137034, 12137035, 12137036, 12137037, 12137038, 12137039, 12137040, 12137041, 12137042, 12137043, 12138005, 12138007, 12138008, 12138009, 12138010, 12138011, 12138012, 12138013, 12138017, 12138021, 12138022, 12138023, 12138024, 12149028, 12149034, 12149035, 12149036, 12149037, 12149038, 12149039, 12149040, 12149044, 12149045, 12149046, 12149047, 12149048, 12149049, 12149050, 12149051, 12149058, 12149059, 12149060, 12149061, 12149062, 12149063, 12150001, 12150002, 12150003, 12150008, 12150009, 12150011, 12150012, 12150013, 12150014, 12150015, 12150016, 12150019, 12150020, 12150021, 12150026, 12150027, 12150034, 12150035, 12150036, 12150037, 12150038, 12150039, 12150040, 12150041, 12150047, 12150048, 12150050, 12150057, 12151001, 12151002, 12151003, 12151004, 12151005, 12151006, 12151009, 12151010, 12151011, 12151012, 12151013, 12151016, 12151017, 12151023, 12151032, 12151033, 12151034, 12151035, 12151036, 12151038, 12151039, 12151040, 12151041, 12151044, 12151045, 12151046, 12151056, 12151063, 12151064, 12151065, 12151066, 12151067, 12151068, 12151069, 12152005, 12152006, 12152007, 12152008, 12152009, 12152010, 12152011, 12152012, 12152013, 12152014, 12152015, 12152016, 12153002, 12153004, 12153007, 12153008, 12153009, 12153012, 12153013, 12153014, 12153015, 12153016, 12153017, 12153018, 12153021, 12153032, 12153033, 12153034, 12154001, 12154002, 12154003, 12154004, 12154005, 12154006, 12154007, 12154008, 12154009, 12154011, 12154012, 12154013, 12154014, 12154015, 12154016, 12154017, 12154018, 12154019, 12154020, 12154021, 12154038, 12154039, 12154043, 12154044, 12154045, 12154046, 12154047, 12154048, 12154065, 12154066, 12154067, 12154068, 12155001, 12155002, 12155004, 12155005, 12155006, 12155007, 12155008, 12155009, 12155011, 12155012, 12155013, 12155014, 12155015, 12155016, 12155017, 12155018, 12155019, 12155020, 12155021, 12155022, 12155023, 12155029, 12155030, 12155031, 12155036, 12155037, 12155038, 12155039, 12155040, 12155043, 12155044, 12155045, 12155046, 12155047, 12155048, 12155049, 12155050, 12155051, 12155054, 12155055, 12155056, 12155057, 12155058, 12155059, 12155060, 12155061, 12155062, 12155063, 12155064, 12155065, 12155066, 12156002, 12156003, 12156004, 12156005, 12156006, 12156008, 12156009, 12156010, 12156011, 12156014, 12156015, 12156016, 12156018, 12156019, 12156020, 12156021, 12156023, 12156029, 12156030, 12156031, 12156032, 12156034, 12156035, 12156036, 12156037, 12156038, 12156041, 12156042, 12156044, 12156045, 12156046, 12156047, 12156048, 12156049, 12156050, 12156051, 12156056, 12156057, 12156058, 12156059, 12156060, 12156061, 12156062, 12156063, 12156065, 12156066, 12157003, 12157009, 12157010, 12157011, 12157012, 12157013, 12157014, 12157015, 12157016, 12157017, 12157022, 12157023, 12157024, 12157025, 12157028, 12157030, 12157031, 12157038, 12157039, 12157040, 12157041, 12157042, 12157043, 12157044, 12157045, 12157046, 12157047, 12157048, 12157051, 12157052, 12158001, 12158002, 12158003, 12158004, 12158005, 12158006, 12158009, 12158010, 12158011, 12158015, 12158016, 12158017, 12158019, 12158020, 12158021, 12158025, 12158026, 12158027, 12158029, 12158040, 12158041, 12158051, 12158052, 12158053, 12158054, 12158056, 12158057, 12158058, 12158059, 12158060, 12158061, 12158066, 12158067, 12158068, 12158069, 12158070, 12158072, 12158073, 12158074, 12159003, 12159004, 12159005, 12159006, 12159012, 12159013, 12159018, 12159019, 12159020, 12159021, 12159022, 12159023, 12159024, 12159025, 12159026, 12159036, 12159039, 12159040, 12160001, 12160002, 12160003, 12160004, 12160005, 12160006, 12160007, 12160012, 12160013, 12160016, 12160017, 12160018, 12160019, 12160020, 12160021, 12160022, 12160023, 12160024, 12160025, 12161005, 12161006, 12161007, 12161009, 12161010, 12161011, 12161012, 12161013, 12161014, 12161015, 12161016, 12161017, 12161020, 12161021, 12161022, 12161024, 12161025, 12161026, 12161027, 12161028, 12161029, 12161030, 12161031, 12161032, 12161033, 12161050, 12161052, 12161053, 12161055, 12161056, 12161057, 12161058, 12161060, 12161061, 12162002, 12162003, 12162004, 12162005, 12162006, 12162007, 12162008, 12162009, 12162010, 12162011, 12162012, 12162014, 12162015, 12162016, 12162017, 12162018, 12162019, 12162020, 12162021, 12162027, 12162028, 12162029, 12162030, 12162031, 12162032, 12162033, 12162034, 12162035, 12162042, 12162043, 12162044, 12162045, 12162046, 12162055, 12162056, 12162057, 12162058, 12162059, 12163002, 12163003, 12163005, 12163006, 12163007, 12163008, 12163009, 12163010, 12163014, 12163015, 12163016, 12163017, 12163018, 12163019, 12163020, 12163021, 12163022, 12163024, 12163049, 12163050, 12163051, 12163052, 12163053, 12163055, 12163056, 12163057, 12163058, 12164002, 12164003, 12164004, 12164005, 12164008, 12164009, 12164010, 12164011, 12164019, 12164037, 12164040, 12164041, 12164042, 12164043, 12164044, 12164045, 12164046, 12164048, 12164050, 12164051, 12164053, 12164056, 12164057, 12164058, 12164063, 12164064, 12164065, 12164066, 12164067, 12164078, 12164079, 12164084, 12164085, 12164086, 12164088, 12164089, 12165001, 12165002, 12165003, 12165004, 12165005, 12165007, 12165008, 12165009, 12165010, 12165012, 12165013, 12165014, 12165015, 12165016, 12165020, 12165021, 12165022, 12165023, 12165024, 12165025, 12165026, 12165028, 12165029, 12165031, 12165037, 12165038, 12165042, 12166002, 12166003, 12166004, 12166005, 12166006, 12166008, 12166009, 12166010, 12166011, 12166051, 12166052, 12166054, 12166057, 12166058, 12166059, 12166060, 12166061, 12166062, 12167002, 12167003, 12167004, 12167005, 12167007, 12167008, 12167009, 12167010, 12167011, 12167012, 12167014, 12167015, 12167021, 12167024, 12167026, 12167029, 12167040, 12167041, 12167045, 12167046, 12167047, 12167048, 12167049, 12167052, 12168001, 12168002, 12168003, 12168004, 12168005, 12168006, 12168009, 12168010, 12168022, 12168054, 12168056, 12168057, 12168058, 12168059, 12168060, 12168061, 12168062, 12168063, 12168064, 12168077, 12168079, 12169013, 12169014, 12169015, 12169016, 12169017, 12169018, 12169019, 12169022, 12169023, 12169024, 12169025, 12169026, 12169027, 12169028, 12169029, 12169030, 12169031, 12169032, 12169033, 12169034, 12169041, 12169042, 12169043, 12169044, 12169045, 12169046, 12169047, 12169048, 12169049, 12169050, 12169052, 12169054, 12169055, 12169056, 12169057, 12169058, 12169059, 12169060, 12169061, 12169062, 12169063, 12169065, 12169066, 12169068, 12170001, 12170002, 12170003, 12170004, 12170005, 12170006, 12170007, 12170009, 12170010, 12170011, 12170012, 12170013, 12170014, 12170015, 12170016, 12170017, 12170018, 12170020, 12170025, 12170026, 12170027, 12170028, 12170029, 12170030, 12170031, 12170032, 12170033, 12170034, 12170035, 12170043, 12170044, 12170045, 12170046, 12170047, 12170048, 12170049, 12170050, 12170051, 12170054, 12170056, 12170057, 12170058, 12171001, 12171002, 12171003, 12171004, 12171007, 12171008, 12171009, 12171010, 12171011, 12171012, 12171013, 12171014, 12171015, 12171016 };
    Int_t nElements = sizeof(runIds) / sizeof(runIds[0]);

    map<Int_t, Int_t> runMap;

    for(unsigned i = 0; i <= (nElements - 1); i++)
    {
        runMap[runIds[i]] = i;
    }

    return runMap;
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
