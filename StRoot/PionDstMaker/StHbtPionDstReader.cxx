#include "StHbtPionDstReader.h"
#include <stdlib.h>
//#include <iterator.h>
//#include <algo.h>
#include "StChain.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include "StPhysicalHelixD.hh"

#include "SystemOfUnits.h"

#include "StIOMaker/StIOMaker.h"

#include "TVector3.h"
#include "TString.h"

#include <math.h>
#include <string>
#include <typeinfo>

#include "StHbtMaker/Infrastructure/StExceptions.hh"
#include "StHbtMaker/Infrastructure/StHbtTrackCollection.hh"
#include "StHbtMaker/Infrastructure/StHbtEvent.hh"

#include "StHbtMaker/Infrastructure/StHbtVector.hh"

#include "Event.h"
#include "Track.h"

#include "StarClassLibrary/StMemoryInfo.hh"

ClassImp(StHbtPionDstReader)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using std::random_shuffle;


//__________________
StHbtPionDstReader::StHbtPionDstReader( const char* dirName, const char* fileName, 
				   const char* filter, int maxFiles)
  :  mMaxFiles(maxFiles), mDebug(0), mCurrentFile(0), 
  mTTree(0) {
  if (mDebug) cout << "StHbtPionDstReader::StHbtPionDstReader(...)"<< endl;

  fRPdistribution = new TH1D("fRPdistribution","Reaction plane distribution",720,-2.*TMath::Pi(),2.*TMath::Pi());

  mDir = string(dirName);
  mFile = string(fileName);
  mFilter = string(filter);
  mReaderStatus = 0;  // "good"
  mEvent = new Event();
  if (mDebug) cout << "StHbtPionDstReader::StHbtPionDstReader(...) - leaving"<< endl;
}

//__________________
StHbtPionDstReader::~StHbtPionDstReader(){
  if (mCurrentFile) { mCurrentFile->Close(); delete mCurrentFile; mCurrentFile = 0;}
}

//__________________
StHbtString StHbtPionDstReader::Report(){
  StHbtString temp = "\n This is the StHbtPionDstReader\n";
  return temp;
}

//__________________
StHbtEvent* StHbtPionDstReader::ReturnHbtEvent(){
  if (mDebug) cout << "StHbtPionDstReader::ReturnHbtEvent()" << endl;

  StHbtEvent* hbtEvent = 0;

  try {
    hbtEvent = read();
  }
  catch(StExceptionEOF e) {
    e.print();
    mReaderStatus = 2;
    return 0;
  }
  catch(StException e) {
    e.print();
    mReaderStatus = 1;
    return 0;
  }
    
  if (!hbtEvent) cout << "StHbtPionDstReader::ReturnHbtEvent() - no hbtEvent" << endl;

  return hbtEvent;
}

//__________________
StHbtEvent* StHbtPionDstReader::read(){
    if (!mTChain) {
        try {
            cout << initRead(mDir,mFile,mFilter,mMaxFiles) << " files to analyse " << endl;
        }
        catch(StException e) {
            e.print();
            return 0;
        }
    }
 
    unsigned int nEvents = (unsigned int)mTChain->GetEntries();
    if (!nEvents) throw StException("StHbtPionDstReader::read() - no events to read ");

    mEvent->Clear();
    int iBytes= mTChain->GetEntry(mEventIndex++);

    if (nEvents<mEventIndex) throw StExceptionEOF("StHbtPionDstReader::read()");
    if (!iBytes) throw StException("StHbtPionDstReader::read() - no event ");

	StHbtEvent *hbtEvent = 0;

	if(mEvent) 
	{
	    hbtEvent = new StHbtEvent;
	    
	    hbtEvent->SetEventNumber(mEvent->GetEventID());
	    hbtEvent->SetRunNumber(mEvent->GetRunID()); //Set the runnumber for the event
	    hbtEvent->SetCtbMult(mEvent->GetCtbMult());
	    hbtEvent->SetZdcAdcEast(mEvent->GetZDCe());
	    hbtEvent->SetZdcAdcWest(mEvent->GetZDCw());
	    hbtEvent->SetNumberOfTpcHits(0);
	    hbtEvent->SetNumberOfTracks(mEvent->GetNumberOfTracks());
	    hbtEvent->SetNumberOfGoodTracks(mEvent->GetNumberOfGoodTracks());
	    hbtEvent->SetUncorrectedNumberOfPositivePrimaries(mEvent->GetRefMultPos());
	    hbtEvent->SetUncorrectedNumberOfNegativePrimaries(mEvent->GetRefMultNeg()); 
	    hbtEvent->SetUncorrectedNumberOfPrimaries(mEvent->GetRefMult());
	    hbtEvent->SetRefmult(mEvent->GetRefMultCorr());
	    hbtEvent->Setq2(mEvent->Getq2(1));
	    hbtEvent->SetReactionPlane(mEvent->GetPsi(mSel,0),0); 
	    hbtEvent->SetReactionPlane(0,1); // This is for pt-weighted RP, which we don't use.
	    hbtEvent->SetReactionPlaneError(0, 0);
	    hbtEvent->SetReactionPlaneSubEventDifference(0, 0);
	    hbtEvent->SetTriggerWord(0);
	    hbtEvent->SetTriggerActionWord(0);
	    hbtEvent->SetL3TriggerAlgorithm(0, 0);
	    hbtEvent->SetPrimVertPos(mEvent->GetVertexPos());  //Set the vertex for the event
	    
	    StHbtTrackCollection* mTrackCollection = hbtEvent->TrackCollection();

	    TClonesArray* tracks = mEvent->GetPionTracks(); 

	    if (tracks) 
        {
            int nTracks = tracks->GetEntries();
            
            TArrayI tab;
            tab=shuffle_order(nTracks);
            
            for ( int i=0; i<nTracks; i++) 
            {
                Track* track = (Track*) tracks->UncheckedAt(tab[i]);
                StHbtTrack* trackCopy = new StHbtTrack;
                trackCopy->SetHiddenInfo(0);
                
                trackCopy->SetCharge(track->GetCharge());
                trackCopy->SetNHits(track->GetNhits());
                
                //---- Set dummy values ----//
                
                trackCopy->SetNHitsDedx(0);
                if (track->GetCharge() < 0) {
                  trackCopy->SetNSigmaElectron(track->GetPidElectron());
                  trackCopy->SetNSigmaPion(track->GetPidPiMinus());
                  trackCopy->SetNSigmaKaon(track->GetPidKaonMinus());
                  trackCopy->SetNSigmaProton(track->GetPidAntiProton());
                } else if (track->GetCharge() > 0) {
                  trackCopy->SetNSigmaElectron(track->GetPidPositron());
                  trackCopy->SetNSigmaPion(track->GetPidPiPlus());
                  trackCopy->SetNSigmaKaon(track->GetPidKaonPlus());
                  trackCopy->SetNSigmaProton(track->GetPidProton());
                }
                trackCopy->SetPidProbElectron(0.);
                trackCopy->SetPidProbPion(0.);
                trackCopy->SetPidProbKaon(0.);
                trackCopy->SetPidProbProton(0.);
                trackCopy->SetdEdx(track->GetDedx());
                trackCopy->SetDCAxy(track->GetDca3().mag());
                trackCopy->SetDCAz(0.);
                trackCopy->SetDCAxyGlobal(track->GetDcaGlobal3().mag());
                trackCopy->SetDCAzGlobal(0.);
                trackCopy->SetChiSquaredXY(track->GetChi2());
                trackCopy->SetChiSquaredZ(0.);
                
                //---- Set momentum ----//
            
                Float_t px = track->P().x();
                Float_t py = track->P().y();
                Float_t pz = track->P().z();
                
                StHbtThreeVector v(px,py,pz);
                
                trackCopy->SetP(v);
                trackCopy->SetPt(sqrt(px*px+py*py));
                
            // NEED TO SET HELIX!!
                
                const StThreeVectorD p((double)px,(double)py,(double)pz);
                const StThreeVectorD origin((double)mEvent->GetVertexPos().x(),(double)mEvent->GetVertexPos().y(),(double)mEvent->GetVertexPos().z());
                
                StPhysicalHelixD helix(p,origin,(double)(mEvent->GetMagField())*kilogauss,(double)(track->GetCharge())); 
                
                trackCopy->SetHelix(helix);
            
                trackCopy->SetTopologyMap(0,track->GetTopologyMap().data(0));
                trackCopy->SetTopologyMap(1,track->GetTopologyMap().data(1));
                
                mTrackCollection->push_back(trackCopy);
            }
        }
	}

	return hbtEvent; 
}

int StHbtPionDstReader::initRead(string dir, string file, string filter, int mMaxFiles){
  mEventIndex =0;
  mTChain = new TChain("PionDst","PionDst");
  cout<<"StHbtPionDstReader::initReader()"<<endl;
  
  int nFiles =0;
  if (file!="") { // if a filename was given
    if( strstr(file.c_str(),".lis") || strstr(file.c_str(),".list") ) { // if a file list is specified
	cout<<"StHbtPionDstReader::initReader() list = "<<file.c_str()<<endl;
      try {
	nFiles = fillChain(mTChain, (dir+file).c_str(), mMaxFiles);
	cout<<"StHbtPionDstReader::InitReader() = "<< (dir+file).c_str()<<endl;
	
      }
      catch(StException e) {
	throw e;
      }
    }
    else { // a single file was specified
      mTChain->Add((dir+file).c_str());
	cout<<"StHbtPionDstReader::InitReader() 2 = "<< (dir+file).c_str()<<endl;
      nFiles++;
    }
  }
  else {
    try {
      nFiles = fillChain(mTChain,dir.c_str(), filter.c_str(), mMaxFiles);
    }
    catch(StException e) {
      throw e;
    }
  }  
	
	mTChain->SetBranchAddress("Event",&mEvent); 
  return nFiles;
}


int StHbtPionDstReader::uninitRead(){
  if (mEvent) delete mEvent;
  if (mTChain) delete mTChain;
  mEvent = 0;
  mTChain = 0;
  return 0;
}

int StHbtPionDstReader::fillChain(TChain* chain, const char* fileList, const int maxFiles) {
  ifstream* inputStream = new ifstream;
  cout<<"StHbtPionDstReader::fillChain() filelist = "<<fileList<<endl;
  
  inputStream->open(fileList);
  if (!(inputStream)) throw StException("StHbtPionDstReader::fillChain(string dir) - can not open directory");
  char* temp;
  int count=0;
  if (mDebug>1) cout << " StHbtPionDstReader::fillChain(...)- inputStream->good() : " << inputStream->good() << endl;
  for (;inputStream->good();) {
    temp = new char[200];
    inputStream->getline(temp,200);

		TString fileName(temp);
		if(fileName.Contains("root")) {
       chain->Add(temp);
       ++count;
		}
 	  delete temp;
    if (count>maxFiles) break;
  }   
  delete inputStream;
  if (mDebug) cout << "StHbtPionDstReader::(string dir)(string dir) - Added " << count << " files to the chain" << endl;
  return count;
}

int StHbtPionDstReader::fillChain(TChain* chain, const char* dir, const char* filter, const int maxFiles) {
  // read directory
  void *pDir = gSystem->OpenDirectory(dir);
  if(!pDir) throw StException("StHbtPionDstReader::fillChain(string dir) - can not open directory");
  // now find the files that end in the specified searchString
  const char* fileName(0);
  int count(0);
  while((fileName = gSystem->GetDirEntry(pDir))){
    if(strcmp(fileName,".")==0 || strcmp(fileName,"..")==0) continue;
    if(strstr(fileName,filter) ) { // found a match
      char* fullFile = gSystem->ConcatFileName(dir,fileName);
      // add it to the chain
      cout << "StHbtPionDstReader::fillChain(string dir) - Adding " << fullFile << " to the chain" << endl;
      chain->Add(fullFile);
      delete fullFile;
      ++count;
      if (count>maxFiles) break;
    }   
  }
  cout << "StHbtPionDstReader::(string dir)(string dir) - Added " << count << " files to the chain" << endl;
  return count;
}

TArrayI StHbtPionDstReader::shuffle_order(int nr)
{
    TArrayI tablica;
    tablica.Set(nr);
    int i=0;
    int test=0;
    int r;
    while (i<nr)
    {
        r=int((double)nr*rand()/(double)RAND_MAX);
        test=0;
        for (int j=0; j<i; j++)
            if (r==tablica[j]) test=1;
        if (test==1)
            continue;
        else
        {
            tablica[i] = r;
            i++;
        }
    }
    return tablica;
}

