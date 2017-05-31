/***************************************************************************
 *
 * Author: John Campbell, OSU, campbell.1119@osu.edu
 *
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: StHbtMaker package
 *   2D correlation function: Qinv vs. Quality
 *
 **************************************************************************/

#include "TMath.h"
#include "StHbtMaker/CorrFctn/QualityvsQinv.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(QualityvsQinv)
#endif

//____________________________
  QualityvsQinv::QualityvsQinv(char* title, const int& nbinsX, const float& XLo, const float& XHi,
					   const int& nbinsY, const float& YLo, const float& YHi){
  // set up numerator
  char Tit[100];
  sprintf(Tit,"Quality_Num2D");
  strcat(Tit,title);
  mNumerator2D = new StHbt2DHisto(Tit,title,nbinsX,XLo,XHi,nbinsY,YLo,YHi);

  // set up denominator
  sprintf(Tit,"Quality_Den2D");
  strcat(Tit,title);
  mDenominator2D = new StHbt2DHisto(Tit,title,nbinsX,XLo,XHi,nbinsY,YLo,YHi);

  // set up ratio
  sprintf(Tit,"Quality_Rat2D");
  strcat(Tit,title);
  mRatio2D = new StHbt2DHisto(Tit,title,nbinsX,XLo,XHi,nbinsY,YLo,YHi);
  
  mNumerator2D->Sumw2();
  mDenominator2D->Sumw2();
  mRatio2D->Sumw2();
}
//____________________________
QualityvsQinv::~QualityvsQinv(){
  delete mNumerator2D;
  delete mDenominator2D;
  delete mRatio2D; 
}
//_________________________
void QualityvsQinv::Finish(){
  mRatio2D->Divide(mNumerator2D,mDenominator2D,1.0,1.0);
}


//____________________________
StHbtString QualityvsQinv::Report(){
  string stemp = "Qinv Correlation Function Report:\n";
  char ctemp[100];
  sprintf(ctemp,"Number of entries in numerator:\t%E\n",mNumerator2D->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in denominator:\t%E\n",mDenominator2D->GetEntries());
  stemp += ctemp;
  sprintf(ctemp,"Number of entries in ratio:\t%E\n",mRatio2D->GetEntries());
  stemp += ctemp;
  StHbtString returnThis = stemp;
  return returnThis;
}
//____________________________
void QualityvsQinv::AddRealPair(const StHbtPair* pair){
  double aveSep = pair->NominalTpcAverageSeparation();
  double fmr = pair->getFracOfMergedRow();
  double quality = pair->quality();

  Bool_t pass = ( (fmr >= fmrLo) &&
                  (fmr <= fmrHi) &&
                  (aveSep >= aveSepLo) &&
                  (aveSep <= aveSepHi) );

    if( pass ) { 
      double Qinv = fabs(pair->qInv());   // note - qInv() will be negative for identical pairs...
      mNumerator2D->Fill(Qinv,quality,1.0);
    }
}

//____________________________
void QualityvsQinv::AddMixedPair(const StHbtPair* pair){
  double aveSep = pair->NominalTpcAverageSeparation();
  double fmr = pair->getFracOfMergedRow();
  double quality = pair->quality();

  Bool_t pass = ( (fmr >= fmrLo) &&
                  (fmr <= fmrHi) &&
                  (aveSep >= aveSepLo) &&
                  (aveSep <= aveSepHi) );

    if( pass ) { 
      double Qinv = fabs(pair->qInv());   // note - qInv() will be negative for identical pairs...
      mDenominator2D->Fill(Qinv,quality,1.0);
    }
}
