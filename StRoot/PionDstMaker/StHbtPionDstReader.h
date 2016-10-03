#ifndef StHbtPionDstReader_hh
#define StHbtPionDstReader_hh

#include "StHbtMaker/Base/StHbtEventReader.hh"
#include "StHbtMaker/Infrastructure/StHbtEnumeration.hh"

#include "StMaker.h"
#include "StChain.h"
#include "St_DataSetIter.h"

#include "StEvent/StEnumerations.h"
#include "StStrangeMuDstMaker/StStrangeMuDstMaker.h"


#include <string>
#include <stdlib.h>
//#include <iterator.h>
//#include <algo.h>

class StHbtEvent;
class Event;
class StIOMaker;
class TFile; 
class TTree;
class TChain;

class StHbtPionDstReader : public StHbtEventReader {

private:
  string mCurrentFileName;
  string mDir;
  string mFile;
  string mFilter;

  int mMaxFiles;

  int mDebug;
  TChain*            mTChain; 
  TFile*             mCurrentFile;
  TTree*             mTTree;

  TH1D *fRPdistribution;

  string             mInputDir; 
  Event*             mEvent; 
  StHbtEvent*        mHbtEvent; 

  unsigned int       mEventIndex;
  int                 mSel;

  StHbtEvent* read();
  int initRead(string dir, string file, string filter, int mMaxFiles);
  int uninitRead();

  int fillChain(TChain* chain, const char* dir, const char* filter, const int maxFiles);
  int fillChain(TChain* chain, const char* list, const int maxFiles);

 protected:
  
  public:
  StHbtPionDstReader(const char* dirName, const char* fileName, 
		   const char* filter=".", int maxFiles=999);
 
 ~StHbtPionDstReader();
  
  StHbtEvent* ReturnHbtEvent();

  StHbtString Report();

  TH1D* ReturnRPdistribution() { return fRPdistribution;};
  int sel() {return mSel;}
  void setSel(Int_t sel) {mSel = sel;}
  
  void SetDebug(int);

  TArrayI shuffle_order(int nr);
  

  ClassDef(StHbtPionDstReader, 1)
};

inline void StHbtPionDstReader::SetDebug(int debug) {mDebug=debug;}
  
#endif
