#include "StHbtMaker/Infrastructure/StHbtPair.hh"
#include "StHbtMaker/Cut/fxtPairCutMonitor.h"
#include <cstdio>
#include <string>
#include "StLorentzVector.hh"

#ifdef __ROOT__ 
ClassImp(fxtPairCutMonitor)
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
fxtPairCutMonitor::fxtPairCutMonitor(const char* name){ // default constructor
  string s("fxtPairCutMonitor");
  string n(name);
  mKt= new StHbt1DHisto( (s+n+"mKt").c_str(),"Kt",1500,0.0,1.5); 
  mMt= new StHbt1DHisto( (s+n+"mMt").c_str(),"Mt",1500,0.0,1.5); 
  mKtLowQ= new StHbt1DHisto( (s+n+"mKtLowQ").c_str(),"KtLowQ",1500,0.0,1.5); 
  mMtLowQ= new StHbt1DHisto( (s+n+"mMtLowQ").c_str(),"MtLowQ",1500,0.0,1.5); 
  mFractionOfMergedRow= new StHbt1DHisto( (s+n+"mFractionOfMergedRow").c_str(),"Fraction of Merged Hits",100,0.,1.); 
  mSplittingLevel= new StHbt1DHisto( (s+n+"mSplittingLevel").c_str(),"Splitting Level",17,-0.6,1.1); 
  mPairRapidityVsEmissionAngle= new StHbt2DHisto( (s+n+"mPairRapidityVsEmissionAngle").c_str(),"Pair Rapidity vs EmissionAngle",180,0.,360, 200,-1.,1.);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
fxtPairCutMonitor::fxtPairCutMonitor( const  fxtPairCutMonitor& cutMoni)  {
  mKt =new StHbt1DHisto(*(cutMoni.mKt));
  mMt =new StHbt1DHisto(*(cutMoni.mMt));
  mKtLowQ =new StHbt1DHisto(*(cutMoni.mKtLowQ));
  mKtLowQ =new StHbt1DHisto(*(cutMoni.mMtLowQ));
  mFractionOfMergedRow =new StHbt1DHisto(*(cutMoni.mFractionOfMergedRow));
  mSplittingLevel =new StHbt1DHisto(*(cutMoni.mSplittingLevel));
  mPairRapidityVsEmissionAngle = new StHbt2DHisto(*(cutMoni.mPairRapidityVsEmissionAngle));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
fxtPairCutMonitor::~fxtPairCutMonitor(){
  delete mKt;
  delete mMt;
  delete mKtLowQ;
  delete mMtLowQ;
  delete mFractionOfMergedRow;
  delete mSplittingLevel ;
  delete mPairRapidityVsEmissionAngle;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  fxtPairCutMonitor::Fill(const StHbtPair* pair){
  Float_t mT = sqrt(pair->kT()*pair->kT() + (0.139*0.139)); // Hard coded pion mass!!
  mKt->Fill( pair->kT(), 1.);
  mMt->Fill( mT, 1.);
  if( pair->qInv() < 0.1 )
    {
      mKtLowQ->Fill( pair->kT(), 1.);
      mMtLowQ->Fill( mT, 1.);
    }
  mFractionOfMergedRow->Fill( pair->getFracOfMergedRow(), 1.);
  mSplittingLevel->Fill(pair->quality(), 1.);
  mPairRapidityVsEmissionAngle->Fill(pair->emissionAngle(), pair->rap(), 1.);
}


