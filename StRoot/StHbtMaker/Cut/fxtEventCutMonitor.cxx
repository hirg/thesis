#include "StHbtMaker/Infrastructure/StHbtEvent.hh"
#include "StHbtMaker/Infrastructure/StHbtTypes.hh"
#include "StHbtMaker/Cut/fxtEventCutMonitor.h"
#include <cstdio>
#include <cmath>

#ifdef __ROOT__ 
ClassImp(fxtEventCutMonitor)
#endif
fxtEventCutMonitor::fxtEventCutMonitor(){
  mVertexYvsVertexX = new StHbt2DHisto("VertexYvsVertexX", "VertexYvsVertexX", 600, -6.,6., 600, -6.,6.);
  mVertexZ = new StHbt1DHisto("VertexZ", "VertexZ", 1600, -80,80);
  mRefMult = new StHbt1DHisto("RefMult", "RefMult", 500, 0.,500);
  mNumberOfTofMatches = new StHbt1DHisto("NumberOfTofMatches", "NumberOfTofMatches", 500, 0.,500);
}
//------------------------------
fxtEventCutMonitor::fxtEventCutMonitor(const char* title1, const char* title2){
  char tit1[100];

  sprintf(tit1,"%s%s_VertexYvsVertexX",title1,title2);
  mVertexYvsVertexX = new StHbt2DHisto(tit1, "VertexYvsVertexX", 600, -6.,6., 600, -6.,6.);

  sprintf(tit1,"%s%s_VertexZ",title1,title2);
  mVertexZ = new StHbt1DHisto(tit1, "VertexZ", 1600, -80,80);

  sprintf(tit1,"%s%s_RefMult",title1,title2);
  mRefMult = new StHbt1DHisto(tit1, "RefMult", 500, 0., 500.);

  sprintf(tit1,"%s%s_NumberOfTofMatches",title1,title2);
  mNumberOfTofMatches= new StHbt1DHisto(tit1, "NumberOfTofMatches", 500, 0., 500.);

}
//------------------------------
fxtEventCutMonitor::~fxtEventCutMonitor(){
//  delete mScaler;
  delete mVertexYvsVertexX;
  delete mVertexZ;
  delete mRefMult;
  delete mNumberOfTofMatches;
}

//------------------------------
void fxtEventCutMonitor::Fill(const StHbtEvent* event){

  mVertexYvsVertexX->Fill( event->PrimVertPos().x(), event->PrimVertPos().y(), 1.);
  mVertexZ->Fill( event->PrimVertPos().z(), 1.);
  mRefMult->Fill( event->UncorrectedNumberOfPrimaries(), 1.);
  mNumberOfTofMatches->Fill( event->NumberOfTofMatches(), 1.);

}

//------------------------------
void fxtEventCutMonitor::Finish(){
  cout << " entries in histogram mVertexYvsVertexX : " << mVertexYvsVertexX->Integral() << endl;
  cout << " entries in histogram mVertexZ : " << mVertexZ->Integral() << endl;
  cout << " entries in histogram mRefMult : " << mRefMult->Integral() << endl;
  cout << " entries in histogram mNumberOfTofMatches : " << mNumberOfTofMatches->Integral() << endl;
}

//------------------------------
StHbtString fxtEventCutMonitor::Report(){
  string Stemp;
  char Ctemp[100];
  sprintf(Ctemp," fxtEventCutMonitor");
  Stemp=Ctemp;
  StHbtString returnThis = Stemp;
  return returnThis;
}


