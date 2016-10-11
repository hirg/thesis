#ifndef fxtEventCut_hh
#define fxtEventCut_hh

#include "StHbtMaker/Base/StHbtEventCut.h"

class fxtEventCut : public StHbtEventCut {

public:

  fxtEventCut();
  fxtEventCut(fxtEventCut&);
  //~fxtEventCut();

  void SetMult(const Int_t lo,const Int_t hi);
  void SetZdc(const Int_t lo,const Int_t hi);
  void Setq2(const Float_t lo,const Float_t hi);
  void SetMinTofMatches(const Int_t min);
  void SetVx(const Float_t lo, const Float_t hi);
  void SetVy(const Float_t lo, const Float_t hi);
  void SetVz(const Float_t lo, const Float_t hi);
  void SetVr(const Float_t Vr);
  Int_t NEventsPassed();
  Int_t NEventsFailed();

  virtual StHbtString Report();
  virtual Bool_t Pass(const StHbtEvent*);

  fxtEventCut* Clone();

private:

  Int_t mMult[2];
  Int_t mZdc[2];
  Float_t mq2[2];
  Int_t mMinTofMatches;
  Float_t mVx[2];
  Float_t mVy[2];
  Float_t mVz[2];
  Float_t mVr;

  long mNEventsPassed;
  long mNEventsFailed;
  long mVzPassed;
  long mVrPassed;
  long mMultPassed;
  long mq2Passed;
  long mZdcPassed;
  long mVzFailed;
  long mVrFailed;
  long mMultFailed;
  long mq2Failed;
  long mZdcFailed;

#ifdef __ROOT__
  ClassDef(fxtEventCut, 1)
#endif

};

inline void fxtEventCut::SetMult(const Int_t lo, const Int_t hi){mMult[0]=lo; mMult[1]=hi;}
inline void fxtEventCut::SetZdc(const Int_t lo, const Int_t hi){mZdc[0]=lo; mZdc[1]=hi;}
inline void fxtEventCut::Setq2(const Float_t lo, const Float_t hi){mq2[0]=lo; mq2[1]=hi;}
inline void fxtEventCut::SetMinTofMatches(const Int_t min){mMinTofMatches=min;}
inline void fxtEventCut::SetVx(const Float_t lo, const Float_t hi){mVx[0]=lo; mVx[1]=hi;}
inline void fxtEventCut::SetVy(const Float_t lo, const Float_t hi){mVy[0]=lo; mVy[1]=hi;}
inline void fxtEventCut::SetVz(const Float_t lo, const Float_t hi){mVz[0]=lo; mVz[1]=hi;}
inline void fxtEventCut::SetVr(const Float_t Vr){mVr=Vr;}
inline Int_t  fxtEventCut::NEventsPassed() {return mNEventsPassed;}
inline Int_t  fxtEventCut::NEventsFailed() {return mNEventsFailed;}
inline fxtEventCut* fxtEventCut::Clone() { fxtEventCut* c = new fxtEventCut(*this); return c;}
inline fxtEventCut::fxtEventCut(fxtEventCut& c) : StHbtEventCut(c) {
  mq2[0] = c.mq2[0];
  mq2[1] = c.mq2[1];
  mMult[0] = c.mMult[0];
  mMult[1] = c.mMult[1];
  mMinTofMatches = c.mMinTofMatches;
  mVx[0] = c.mVx[0];
  mVx[1] = c.mVx[1];
  mVy[0] = c.mVy[0];
  mVy[1] = c.mVy[1];
  mVz[0] = c.mVz[0];
  mVz[1] = c.mVz[1];
  mVr = c.mVr;
  mNEventsPassed = 0;
  mNEventsFailed = 0;
}


#endif
