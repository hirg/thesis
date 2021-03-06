#include "StHbtMaker/Cut/fxtEventCut.h"
#include <cstdio>

#ifdef __ROOT__
ClassImp(fxtEventCut)
#endif

fxtEventCut::fxtEventCut(){
  mNEventsPassed =  mNEventsFailed = 0;
  mVzPassed =  mVzFailed = 0;
  mVrPassed =  mVrFailed = 0;
  mq2Passed =  mq2Failed = 0;
  mMultPassed =  mMultFailed = 0;
  mZdcPassed =  mZdcFailed = 0;

  mq2[0] = -1;
  mq2[1] = 100;
  mMult[0] = 0;
  mMult[1] = 10000;
  mMinTofMatches = -2;
  mVx[0] = mVy[0] = mVz[0] = -500;
  mVx[1] = mVy[1] = mVz[1] = 500;
  mVr = 500;
  mZdc[0] = -1;
  mZdc[1] = 10000;
} 

Bool_t fxtEventCut::Pass(const StHbtEvent* event){
  Float_t mult = event->Refmult();
  Float_t q2 = event->q2();
  Int_t nTofMatches = event->NumberOfTofMatches();
  Float_t Vx = event->PrimVertPos().x();
  Float_t Vy = event->PrimVertPos().y();
  Float_t Vz = event->PrimVertPos().z();
  Float_t Vr = sqrt(Vx*Vx + Vy*Vy);
  Int_t zdcE = event->ZdcAdcEast();
  Int_t zdcW = event->ZdcAdcWest();
  Int_t zdcHigher = (zdcE > zdcW) ? zdcE : zdcW;

  Bool_t goodq2 = (q2 > mq2[0]) && (q2 < mq2[1]);
  Bool_t goodMult = (mult > mMult[0]) && (mult <= mMult[1]);
  Bool_t goodZdc = (zdcHigher > mZdc[0]) && (zdcHigher <= mZdc[1]);
  Bool_t goodVz = (Vz < mVz[1]) && (Vr < mVr);
  Bool_t goodVr = (Vr < mVr);

  Bool_t goodEvent =
    ((q2 > mq2[0]) && 
     (q2 < mq2[1]) && 
     (zdcHigher > mZdc[0]) && 
     (zdcHigher <= mZdc[1]) && 
     (mult > mMult[0]) && 
     (mult <= mMult[1]) && 
     (nTofMatches >= mMinTofMatches) && 
     (Vx > mVx[0]) &&
     (Vx < mVx[1]) &&
     (Vy > mVy[0]) &&
     (Vy < mVy[1]) &&
     (Vz > mVz[0]) &&
     (Vz < mVz[1]) &&
     (Vr < mVr));

  goodEvent ? mNEventsPassed++ : mNEventsFailed++ ;
  goodVz ? mVzPassed++ : mVzFailed++ ;
  goodVr ? mVrPassed++ : mVrFailed++ ;
  goodZdc ? mZdcPassed++ : mZdcFailed++ ;
  goodMult ? mMultPassed++ : mMultFailed++ ;
  goodq2 ? mq2Passed++ : mq2Failed++ ;
  return (goodEvent);
}
//------------------------------
StHbtString fxtEventCut::Report(){
  string Stemp;
  char Ctemp[300];
  sprintf(Ctemp,"\nMultiplicity:\t %d-%d",mMult[0],mMult[1]);
  Stemp = Ctemp;
  sprintf(Ctemp,"\nZdc:\t %d-%d",mZdc[0],mZdc[1]);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nq2:\t %f-%f",mq2[0],mq2[1]);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nMin Tof Matches:\t %d",mMinTofMatches);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nVertex X-position:\t %E-%E",mVx[0],mVx[1]);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nVertex Y-position:\t %E-%E",mVy[0],mVy[1]);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nVertex Z-position:\t %E-%E",mVz[0],mVz[1]);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nVertex R-position:\t %E",mVr);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nNumber of events which passed:\t%ld  Number which failed:\t%ld",mNEventsPassed,mNEventsFailed);
  Stemp += Ctemp;
  
  // Individual cuts
  sprintf(Ctemp,"\nNumber of events which passed Vz:\t%ld  Number which failed Vz:\t%ld",mVzPassed,mVzFailed);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nNumber of events which passed Vr:\t%ld  Number which failed Vr:\t%ld",mVrPassed,mVrFailed);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nNumber of events which passed q2:\t%ld  Number which failed q2:\t%ld",mq2Passed,mq2Failed);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nNumber of events which passed mult:\t%ld  Number which failed mult:\t%ld",mMultPassed,mMultFailed);
  Stemp += Ctemp;
  sprintf(Ctemp,"\nNumber of events which passed zdc:\t%ld  Number which failed zdc:\t%ld",mZdcPassed,mZdcFailed);
  Stemp += Ctemp;
  StHbtString returnThis = Stemp;
  return returnThis;
}
