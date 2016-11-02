#include <vector>
TH3D* histCopy(TH3* hist, TString nameTitle);
TH3D* correctedHistogram(TH3D* qInv, TH3D* denomHist, TString corrLabel, Float_t chargeRadius = 5.0);
void doPhiBin(TFile* inputFile, Int_t bin = 0);

void CoulombCorrect(const TString inputFileName =
                        "~/scratch/200GeV/testCopy.root",
                       )
{
	TStopwatch* timer = new TStopwatch();
	//------------------- Load Shared Libraries ------------------//

	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StFlowMaker");
	gSystem->Load("StHbtMaker");

	//------------------- Define Cut Values ------------------//
    //
    TFile* inputFile = new TFile( inputFileName.Data(), "UPDATE");
    if (inputFile->IsZombie()) { 
        cout << "Unable to open " << inputFileName.Data() <<". It is a zombie.\n";
        exit(-1);
    }

    Int_t nPhiBins = 8;
    for(Int_t i = 0; i <= (nPhiBins - 1); i++)
    { doPhiBin(inputFile, i); }
  
}

void doPhiBin(TFile* inputFile, Int_t bin)
{

    TString phiLabels[8] = {"0", "22", "45", "67", "90", "112", "135", "157"};
    TString qInvLabel = TString::Format("QinvPiPlus_phi%s_kt4",
                                        phiLabels[bin].Data());
    TString denLabel = TString::Format("DenPiPlus_phi%s_kt4",
                                        phiLabels[bin].Data());
    TString corrLabel = TString::Format("xCoulPiPlus_phi%s_kt4",
                                        phiLabels[bin].Data());
    TH3D* qInv = (TH3D*)inputFile->Get(qInvLabel.Data());
    TH3D* den = (TH3D*)inputFile->Get(denLabel.Data());
    TH3D* coul = correctedHistogram(qInv, den, corrLabel); 

    inputFile->cd();
    coul->Write();

}

TH3D* correctedHistogram(TH3D* qInv, TH3D* denomHist, TString corrLabel, Float_t chargeRadius)
{
    StHbtCoulomb* corr = new StHbtCoulomb();
    corr->SetRadius(chargeRadius);
    corr->SetChargeProduct(1.0);

    TH3D* returnHist = histCopy(denomHist, corrLabel.Data());

    Float_t mPion = 0.1396;
    Float_t fine_structure = 0.007297352;
    Float_t etaConstant = mPion * fine_structure;
    for(Int_t i = 0; i <= qInv->GetNbinsX(); i++){
        for(Int_t j = 0; j <= qInv->GetNbinsY(); j++){
            for(Int_t k = 0; k <= qInv->GetNbinsZ(); k++){
                Float_t eta = etaConstant / qInv->GetBinContent(i,j,k);
                Double_t denom = denomHist->GetBinContent(i,j,k);

                Float_t weight = corr->CoulombCorrect(eta, chargeRadius);
                returnHist->SetBinContent(i,j,k,denom*weight);

            } // z-bins
        } // y-bins
    } // x-bins

    return returnHist;
}

TH3D* histCopy(TH3* hist, TString nameTitle)
{

    Int_t nBins = hist->GetNbinsX();   
    Double_t lo = hist->GetXaxis()->GetXmin();
    Double_t hi = hist->GetXaxis()->GetXmax();

    TH3D* rHist = new TH3D(nameTitle.Data(),nameTitle.Data(),nBins,lo,hi,nBins,lo,hi,nBins,lo,hi);

    return rHist;

}
