#ifndef TRACK_H
#define TRACK_H
#include <string.h>
#include <math.h>
#include "Rtypes.h"
#include "StObject.h"
//#include "../StFlowMaker/StFlowConstants.h"
#include "StTrackTopologyMap.h"
#include "StThreeVectorD.hh"
//#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"

class Track : public StObject {

public:

  Track(double Px, double Py, double Pz);
  Track() { maxInt = 32.;}
  virtual ~Track() {};

  
  Float_t       GetPidPiPlus()       const;
  Float_t       GetPidPiMinus()      const;
  Float_t       GetPidProton()       const;
  Float_t       GetPidKaonMinus()    const;
  Float_t       GetPidKaonPlus()     const;
  Float_t       GetPidAntiProton()   const;
  Float_t       GetPidDeuteron()     const;
  Float_t       GetPidAntiDeuteron() const;
  Float_t       GetPidElectron()     const;
  Float_t       GetPidPositron()     const;

  const Char_t* GetPid()         const;
  Float_t       GetPhi()         const;
  Float_t       GetPhiGlobal()   const;
  Float_t       GetEta()         const;
  Float_t       GetEtaGlobal()   const;
  Float_t       GetZFirstPoint() const;
  Float_t       GetZLastPoint()  const;
  Float_t       GetDedx()        const;
  Float_t       GetPt()          const;
  Float_t       GetPtGlobal()    const;
  Float_t       GetP()           const;
  Float_t       GetPGlobal()     const;
  Float_t       GetY()           const;
  Short_t       GetCharge()      const;
  Float_t       GetDca()         const;
  StThreeVectorD GetDca3()  const;
  Float_t       GetDcaSigned()   const;
  Float_t       GetDcaGlobal()   const;
  Float_t       GetChi2()        const;
  Int_t         GetFitPts()      const;  // contains fit points in TPC xor FTPC only (SVT and/or SSD hits subtracted)
  Int_t         GetMaxPts()      const;  // contains possible hits in TPC xor FTPC only (SVT and/or SSD hits subtracted)
  Int_t         GetNhits()       const;  // contains ALL hits on the track (TPC + SVT + SSD + FTPC east + FTPC west)
  Int_t         GetNdedxPts()    const;
  Float_t       GetTrackLength() const;
//  Int_t Select(Int_t harmonic, Int_t selection, Int_t subevent= -1) const;
  Int_t         GetMostLikelihoodPID()    const; 
  Float_t       GetMostLikelihoodProb()   const; 
  Int_t         GetExtrapTag()            const;
  Float_t       GetElectronPositronProb() const;
  Float_t       GetPionPlusMinusProb()    const; 
  Float_t       GetKaonPlusMinusProb()    const; 
  Float_t       GetProtonPbarProb()       const;
  StThreeVectorD GetDcaGlobal3()          const; 
  const StTrackTopologyMap& GetTopologyMap() const;
  Int_t  GetFlag() const;
  Int_t  GetIndex2Global() const;
  Float_t GetZFirstPointX() const;
  Float_t GetZFirstPointY() const;
  StPhysicalHelixD GetHelix();

  Float_t GetMsquared() const;
  void SetMsquared(Float_t);
  
  void SetPidPiPlus(Float_t);
  void SetPidPiMinus(Float_t);
  void SetPidProton(Float_t);
  void SetPidKaonMinus(Float_t);
  void SetPidKaonPlus(Float_t);
  void SetPidAntiProton(Float_t);
  void SetPidDeuteron(Float_t);
  void SetPidAntiDeuteron(Float_t);
  void SetPidElectron(Float_t);
  void SetPidPositron(Float_t);
  void SetPid(const Char_t*);
  void SetPhi(Float_t);
  void SetPhiGlobal(Float_t);
  void SetEta(Float_t);
  void SetEtaGlobal(Float_t);
  void SetZFirstPoint(Float_t);
  void SetZLastPoint(Float_t);
  void SetDedx(Float_t);
  void SetPt(Float_t);
  void SetPtGlobal(Float_t);
  void SetCharge(Short_t);
  void SetDca(Float_t);
  void SetDca3(StThreeVectorD);
  void SetDcaSigned(Float_t);
  void SetDcaGlobal(Float_t);
  void SetChi2(Float_t);
  void SetFitPts(Int_t);
  void SetMaxPts(Int_t);
  void SetNhits(Int_t);
  void SetNdedxPts(Int_t);
  void SetTrackLength(Float_t);
//  void SetSelect(Int_t harmonic, Int_t selection);
//  void SetSubevent(Int_t harmonic, Int_t selection, Int_t subevent);
  void SetMostLikelihoodPID(Int_t); 
  void SetMostLikelihoodProb(Float_t); 
  void SetExtrapTag(Int_t); 
  void SetElectronPositronProb(Float_t);
  void SetPionPlusMinusProb(Float_t);
  void SetKaonPlusMinusProb(Float_t);
  void SetProtonPbarProb(Float_t);
  void SetDcaGlobal3(StThreeVectorD gdca3);
  void SetTopologyMap(StTrackTopologyMap map);

  void SetFlag(Int_t);
  void SetIndex2Global(Int_t);
  void SetZFirstPointX(Float_t);
  void SetZFirstPointY(Float_t);
  void SetHelix(StPhysicalHelixD);

  void SetP(double Px, double Py, double Pz) { fP[0] = Px; fP[1] = Py; fP[2] = Pz; }
  void SetP(StThreeVectorD v) { fP = v; }
  StThreeVectorD& P() { return fP; }

//  Float_t maxInt;
  
private:

  StThreeVectorD fP;
  Float_t mMsquared;
  
  Int_t   mFlag;
  Int_t   mIndex2Global;
  Float_t mZFirstPointX;
  Float_t mZFirstPointY;
  StPhysicalHelixD mHelix;
  
  Int_t   mPidPiPlus;
  Int_t   mPidPiMinus;
  Int_t   mPidProton;
  Int_t   mPidKaonPlus;
  Int_t   mPidKaonMinus;
  Int_t   mPidAntiProton;
  Int_t   mPidDeuteron;
  Int_t   mPidAntiDeuteron;
  Int_t   mPidElectron;
  Int_t   mPidPositron;
  Char_t  mPid[10];
  Float_t mPhi;
  Float_t mPhiGlobal;
  Float_t mEta;
  Float_t mEtaGlobal;
  Float_t mZFirstPoint;
  Float_t mZLastPoint;
  Float_t mDedx;
  Float_t mPt;
  Float_t mPtGlobal;
  Short_t mCharge;
  Float_t mDca;
  Float_t mDcaSigned;
  Float_t mDcaGlobal;
  Float_t mChi2;
  Int_t   mFitPts; // contains fit points in TPC xor FTPC only (SVT and/or SSD hits subtracted)
  Int_t   mMaxPts; // contains possible hits in TPC xor FTPC only (SVT and/or SSD hits subtracted)
  Int_t   mNhits;  // contains ALL hits on the track (TPC + SVT + SSD + FTPC east + FTPC west)
  Int_t   mNdedxPts;
  Float_t mTrackLength;
//  Int_t   mSelection;
//  Short_t mSubevent[Flow::nHars][Flow::nSels];
    Float_t maxInt;  //declared globally in Track.cxx
  Int_t   mMostLikelihoodPID;
  Float_t mMostLikelihoodProb;
  Int_t   mExtrapTag; //merging area tag.
  Float_t mElectronPositronProb;
  Float_t mPionPlusMinusProb;  
  Float_t mKaonPlusMinusProb;  
  Float_t mProtonPbarProb;
  StThreeVectorD mDcaGlobal3;  
  StThreeVectorD mDca3;
  StTrackTopologyMap mTopology; //Not set in FillFromMuDst but needed for Hbt
  				//SO SET IT TOO LIKE IN FillFromEvent()

  ClassDef(Track,1)
};
inline Float_t Track::GetMsquared() const { return mMsquared; }
inline void Track::SetMsquared(Float_t ms) { mMsquared = ms; }

//Define the Get functions for the Track
inline Float_t  Track::GetPidPiPlus()    const { return mPidPiPlus/1000.; }
inline Float_t  Track::GetPidPiMinus()   const { return mPidPiMinus/1000.; }
inline Float_t  Track::GetPidProton()    const { return mPidProton/1000.; }
inline Float_t  Track::GetPidKaonMinus() const { return mPidKaonMinus/1000.; }
inline Float_t  Track::GetPidKaonPlus()  const { return mPidKaonPlus/1000.; }
inline Float_t  Track::GetPidAntiProton() const { return mPidAntiProton/1000.; }
inline Float_t  Track::GetPidDeuteron()  const { return mPidDeuteron/1000.; }
inline Float_t  Track::GetPidAntiDeuteron() const { return mPidAntiDeuteron/1000.; }
inline Float_t  Track::GetPidElectron()  const { return mPidElectron/1000.; }
inline Float_t  Track::GetPidPositron()  const { return mPidPositron/1000.; }
inline const Char_t* Track::GetPid()     const { return mPid; }
inline Float_t  Track::GetPhi()          const { return mPhi; }   
inline Float_t  Track::GetPhiGlobal()    const { return mPhiGlobal; }   
inline Float_t  Track::GetEta()          const { return mEta; }     
inline Float_t  Track::GetEtaGlobal()    const { return mEtaGlobal; }     
inline Float_t  Track::GetZFirstPoint()  const { return mZFirstPoint; }     
inline Float_t  Track::GetZLastPoint()   const { return mZLastPoint; }     
inline Float_t  Track::GetDedx()         const { return mDedx; }     
inline Float_t  Track::GetPt()           const { return mPt; }
inline Float_t  Track::GetPtGlobal()     const { return mPtGlobal; }
inline Short_t  Track::GetCharge()       const { return mCharge; }   
inline Float_t  Track::GetDca()          const { return mDca; }
inline Float_t  Track::GetDcaSigned()    const { return mDcaSigned; }
inline Float_t  Track::GetDcaGlobal()    const { return mDcaGlobal; }
inline Float_t  Track::GetChi2()         const { return mChi2; } 
inline Int_t    Track::GetFitPts()       const { return mFitPts; }  
inline Int_t    Track::GetMaxPts()       const { return mMaxPts; }  
inline Int_t    Track::GetNhits()        const { return mNhits; }
inline Int_t    Track::GetNdedxPts()     const { return mNdedxPts; }  
inline Float_t  Track::GetTrackLength()  const { return mTrackLength; }  
inline Int_t    Track::GetMostLikelihoodPID() const
{ return mMostLikelihoodPID; } 
inline Float_t  Track::GetMostLikelihoodProb() const 
{ return mMostLikelihoodProb; } 
inline Int_t    Track::GetExtrapTag()    const { return mExtrapTag; } 
inline Float_t  Track::GetElectronPositronProb() const { return mElectronPositronProb; }
inline Float_t  Track::GetPionPlusMinusProb() const { return mPionPlusMinusProb; }
inline Float_t  Track::GetKaonPlusMinusProb() const { return mKaonPlusMinusProb; }
inline Float_t  Track::GetProtonPbarProb() const { return mProtonPbarProb; }
inline StThreeVectorD Track::GetDcaGlobal3() const { return mDcaGlobal3; }
inline StThreeVectorD Track::GetDca3() const { return mDca3; }
inline const StTrackTopologyMap& Track::GetTopologyMap() const { return mTopology; }


inline Float_t Track::GetP()             const { 
  float momentum = mPt/::sqrt(1-(tanh(mEta)*tanh(mEta)));
  return momentum; }

inline Float_t Track::GetPGlobal()       const { 
  float momentum = mPtGlobal/::sqrt(1-(tanh(mEtaGlobal)*tanh(mEtaGlobal)));
  return momentum; }

inline Float_t Track::GetY()             const { 
  float M = 0.139; 
  if (strcmp(mPid, "none") == 0)          { M = 0.139; }
  else if (strcmp(mPid, "pi+") == 0)      { M = 0.139; }
  else if (strcmp(mPid, "pi-") == 0)      { M = 0.139; }
  else if (strcmp(mPid, "pr+") == 0)      { M = 0.938; }
  else if (strcmp(mPid, "pr-") == 0)      { M = 0.938; }
  else if (strcmp(mPid, "k+")  == 0)      { M = 0.494; }
  else if (strcmp(mPid, "k-")  == 0)      { M = 0.494; }
  else if (strcmp(mPid, "d+")  == 0)      { M = 1.876; }
  else if (strcmp(mPid, "d-")  == 0)      { M = 1.876; }
  else if (strcmp(mPid, "e-")  == 0)      { M = 0.0005; }
  else if (strcmp(mPid, "e+")  == 0)      { M = 0.0005; }
  double Pz = ::sqrt(this->GetP()*this->GetP() - mPt*mPt); 
  if (mEta < 0) { Pz = -Pz; }
  double E = ::sqrt(this->GetP()*this->GetP() + M*M);
  float rapidity = 0.5*::log((E + Pz)/(E - Pz));
  return rapidity;
}

inline Int_t Track::GetFlag() const { return mFlag; }

inline Int_t Track::GetIndex2Global() const { return mIndex2Global; }  

inline Float_t Track::GetZFirstPointX() const { return mZFirstPointX; } 

inline Float_t Track::GetZFirstPointY() const { return mZFirstPointY; } 

inline StPhysicalHelixD Track::GetHelix() { return mHelix; }


inline Track::Track(double Px, double Py, double Pz) {
  fP[0] = Px;
  fP[1] = Py;
  fP[2] = Pz;
  maxInt = 32.;
}


//Define the Set functions for the Track
inline void Track::SetMostLikelihoodPID(Int_t val) {
         mMostLikelihoodPID=val; } 

inline void Track::SetMostLikelihoodProb(Float_t val) {
         mMostLikelihoodProb=val; } 

inline void Track::SetExtrapTag(Int_t val) {
         mExtrapTag=val; } 

inline void Track::SetElectronPositronProb(Float_t val) {
  mElectronPositronProb = val; }

inline void Track::SetPionPlusMinusProb(Float_t val) {
  mPionPlusMinusProb = val; }

inline void Track::SetKaonPlusMinusProb(Float_t val) {
  mKaonPlusMinusProb = val; }

inline void Track::SetProtonPbarProb(Float_t val) {
  mProtonPbarProb = val; }

inline void Track::SetPidPiPlus(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidPiPlus = (Int_t)(pid*1000.); }

inline void Track::SetPidPiMinus(Float_t pid) {
  if (fabs(pid) > maxInt) pid = maxInt; mPidPiMinus = (Int_t)(pid*1000.); }

inline void Track::SetPidProton(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidProton = (Int_t)(pid*1000.); }

inline void Track::SetPidKaonMinus(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidKaonMinus = (Int_t)(pid*1000.); }

inline void Track::SetPidKaonPlus(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidKaonPlus = (Int_t)(pid*1000.); }

inline void Track::SetPidAntiProton(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidAntiProton = (Int_t)(pid*1000.); }

inline void Track::SetPidDeuteron(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidDeuteron = (Int_t)(pid*1000.); }

inline void Track::SetPidAntiDeuteron(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidAntiDeuteron = (Int_t)(pid*1000.); }

inline void Track::SetPidElectron(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidElectron = (Int_t)(pid*1000.); }

inline void Track::SetPidPositron(Float_t pid)  {
  if (fabs(pid) > maxInt) pid = maxInt; mPidPositron = (Int_t)(pid*1000.); }

inline void Track::SetPid(const Char_t* pid)  { strncpy(mPid, pid, 9);
                                                         mPid[9] = '\0'; }
inline void Track::SetPhi(Float_t phi)        { mPhi = phi; }      

inline void Track::SetPhiGlobal(Float_t gphi) { mPhiGlobal = gphi; }   

inline void Track::SetEta(Float_t eta)        { mEta = eta; }       

inline void Track::SetEtaGlobal(Float_t geta) { mEtaGlobal = geta; } 

inline void Track::SetZFirstPoint(Float_t zFirst) { mZFirstPoint = zFirst; } 

inline void Track::SetZLastPoint(Float_t zLast) { mZLastPoint = zLast; } 

inline void Track::SetDedx(Float_t dedx)      { mDedx = dedx; }       

inline void Track::SetPt(Float_t pt)          { mPt = pt; }              

inline void Track::SetPtGlobal(Float_t gpt)   { mPtGlobal = gpt; }

inline void Track::SetCharge(Short_t charge)  { mCharge = charge; }     

inline void Track::SetDca(Float_t dca)        { mDca = dca; }

inline void Track::SetDcaSigned(Float_t sdca) { mDcaSigned = sdca; }

inline void Track::SetDcaGlobal(Float_t gdca) { mDcaGlobal = gdca; }

inline void Track::SetChi2(Float_t chi2)      { mChi2 = chi2; }

inline void Track::SetFitPts(Int_t fitPts)    { mFitPts = fitPts; }

inline void Track::SetMaxPts(Int_t maxPts)    { mMaxPts = maxPts; }

inline void Track::SetNhits(Int_t nhits)      { mNhits = nhits; }

inline void Track::SetNdedxPts(Int_t ndedxPts) { mNdedxPts = ndedxPts; }

inline void Track::SetTrackLength(Float_t tl) { mTrackLength = tl; }

inline void Track::SetDcaGlobal3(StThreeVectorD gdca3) { mDcaGlobal3 = gdca3; }

inline void Track::SetDca3(StThreeVectorD pdca3) { mDca3 = pdca3; }

inline void Track::SetFlag(Int_t flag) { mFlag = flag; }

inline void Track::SetIndex2Global(Int_t i2g) { mIndex2Global = i2g; }  

inline void Track::SetZFirstPointX(Float_t zFirstX) { mZFirstPointX = zFirstX; } 

inline void Track::SetZFirstPointY(Float_t zFirstY) { mZFirstPointY = zFirstY; } 

inline void Track::SetHelix(StPhysicalHelixD h) { mHelix = h; }

//NOTE:  This is different than in the Track class used for making PionDst's.
inline void Track::SetTopologyMap(StTrackTopologyMap map) { mTopology = map; }

#endif
