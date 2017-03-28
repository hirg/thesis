#include "StHbtMaker/CorrFctn/QoslCMSCorrFctnkT.h"
#include <cstdio>
#include "StHbtMaker/Infrastructure/StHbtReactionPlaneAnalysis.h"
#include "PhysicalConstants.h"
#include "TString.h"
#include "TH3F.h"

#ifdef __ROOT__ 
ClassImp(QoslCMSCorrFctnkT)
#endif
//____________________________
QoslCMSCorrFctnkT::QoslCMSCorrFctnkT(char* title, 
               const int& nbinso, const float& QoLo, const float& QoHi,
			   const int& nbinss, const float& QsLo, const float& QsHi,
			   const int& nbinsl, const float& QlLo, const float& QlHi){

  mCorrection = new StHbtCoulomb;
  mCorrection->SetRadius(5.0);
  mCorrection->SetChargeProduct(1.0);

  qMax = QoHi;

  nKtBins = 5;
  nOverflowReal = 0;
  nOverflowMixed = 0;
  nEntriesReal = 0;
  nEntriesMixed = 0;
  

    for(int j=0; j<nKtBins; j++) {
        TString TitKt=Form("_kt%i",j);

        // set up numerator
        TString TitNum = "Num";
        TitNum += title;
        TitNum += TitKt.Data();
        mNumerator[j] = new TH3F(TitNum.Data(),TitNum.Data(),nbinso,QoLo,QoHi,
                    nbinss,QsLo,QsHi,
                    nbinsl,QlLo,QlHi);

        // set up denominator
        TString TitDen = "Den";
        TitDen += title;
        TitDen += TitKt.Data();
        mDenominator[j] = new TH3F(TitDen.Data(),TitDen.Data(),nbinso,QoLo,QoHi,
                      nbinss,QsLo,QsHi,
                      nbinsl,QlLo,QlHi);

        // set up coulomb correction
        TString TitCoul = "Coul";
        TitCoul += title;
        TitCoul += TitKt.Data();
        mCoulHisto[j] = new StHbt3DHisto(TitCoul.Data(),TitCoul.Data(),nbinso,QoLo,QoHi,nbinss,QsLo,QsHi,nbinsl,QlLo,QlHi);

        // set up Qinv
        TString TitQinv = "Qinv";
        TitQinv += title;
        TitQinv += TitKt.Data();
        mQinvHisto[j] = new StHbt3DHisto(TitQinv.Data(),TitQinv.Data(),nbinso,QoLo,QoHi,nbinss,QsLo,QsHi,nbinsl,QlLo,QlHi);

    }
  
}

//____________________________
QoslCMSCorrFctnkT::~QoslCMSCorrFctnkT(){
   for(int j=0; j<nKtBins; j++) {
    delete mNumerator[j];
    delete mDenominator[j];
    delete mCoulHisto[j];
    delete mQinvHisto[j];
   }
}
//_________________________
void QoslCMSCorrFctnkT::Finish(){

}

//____________________________
StHbtString QoslCMSCorrFctnkT::Report(){
  string stemp = "QoslCMS Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerators/denominators:\t%E/%E\n",nEntriesReal,nEntriesMixed);
  stemp += ctemp;
  sprintf(ctemp,"Number of entries overflowing in numerators/denominators:\t%E/%E\n",nOverflowReal,nOverflowMixed);
  stemp += ctemp;

  StHbtString returnThis = stemp;
  return returnThis;
}
//____________________________
void QoslCMSCorrFctnkT::AddRealPair(const StHbtPair* pair){

  int rpBin, ktBin;
  ktBin = GetKtBin(pair);
  if(ktBin<0) return;
  
  double Qo = pair->qOutCMS();
  double Qs = pair->qSideCMS();
  double Ql = pair->qLongCMS();

	if(fabs(Qo)>qMax || fabs(Qs)>qMax || fabs(Ql)>qMax) 
    {
        nOverflowReal++;
        return; 
    }

    nEntriesReal++;
  mNumerator[ktBin]->Fill(Qo,Qs,Ql);
  mNumerator[4]->Fill(Qo,Qs,Ql);
}
//____________________________
void QoslCMSCorrFctnkT::AddMixedPair(const StHbtPair* pair){

  int rpBin, ktBin;
  ktBin = GetKtBin(pair);
  if(ktBin<0) return;

  double weight = 1.0;
  if (mCorrection) weight = mCorrection->CoulombCorrect(pair);

  double Qo = pair->qOutCMS();
  double Qs = pair->qSideCMS();
  double Ql = pair->qLongCMS();
  double Qinv = fabs(pair->qInv());   

	if(fabs(Qo)>qMax || fabs(Qs)>qMax || fabs(Ql)>qMax) 
    {
        nOverflowMixed++;
        return; 
    }
    nEntriesMixed++;
  mDenominator[ktBin]->Fill(Qo,Qs,Ql);
  mDenominator[4]->Fill(Qo,Qs,Ql);
  mCoulHisto[ktBin]->Fill(Qo,Qs,Ql,weight);
  mCoulHisto[4]->Fill(Qo,Qs,Ql,weight);
  mQinvHisto[ktBin]->Fill(Qo,Qs,Ql,Qinv);
  mQinvHisto[4]->Fill(Qo,Qs,Ql,Qinv);

}
//____________________________
int QoslCMSCorrFctnkT::GetKtBin(const StHbtPair* pair) {

  double kT = fabs(pair->kT());
  int ktBin;

  if(kT<0.15 || kT>0.6) return -1;
  
  if(kT<0.25) 
    ktBin = 0;
  else if( kT >= 0.25 && kT < 0.35 ) 
    ktBin = 1;
  else if( kT >= 0.35 && kT < 0.45 ) 
    ktBin = 2;
  else if( kT >= 0.45 && kT <= 0.6 ) 
    ktBin = 3;

  return ktBin;

}
//____________________________
void QoslCMSCorrFctnkT::SetCorrection(StHbtCoulomb* coulomb) {
  mCorrection = coulomb;
}
