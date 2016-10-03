#ifndef QoslCMSCorrFctnkT_hh
#define QoslCMSCorrFctnkT_hh

#include "StHbtMaker/Base/StHbtCorrFctn.hh"
#include "StHbtMaker/Base/StHbtPairCut.h"
#include "StHbtMaker/Infrastructure/StHbtCoulomb.h"

class TH3S;

class QoslCMSCorrFctnkT : public StHbtCorrFctn {
public:
  QoslCMSCorrFctnkT(char* title,
           const Int_t& nbinso, const Float_t& QoLo, const Float_t& QoHi,
	       const Int_t& nbinss, const Float_t& QsLo, const Float_t& QsHi,
	       const Int_t& nbinsl, const Float_t& QlLo, const Float_t& QlHi);

  virtual ~QoslCMSCorrFctnkT();

  virtual StHbtString Report();
  virtual void AddRealPair(const StHbtPair*);
  virtual void AddMixedPair(const StHbtPair*);
  void SetCorrection(StHbtCoulomb*);

  virtual void Finish();

  Float_t qMax;
  Int_t nKtBins;
  Double_t nOverflowReal;
  Double_t nOverflowMixed;
  Double_t nEntriesReal;
  Double_t nEntriesMixed;
  
  TH3S* Numerator3D(const Int_t& ktBin);
  TH3S* Denominator3D(const Int_t& ktBin);
  StHbt3DHisto* QinvHisto3D(const Int_t& ktBin);
  StHbt3DHisto* CoulHisto3D(const Int_t& ktBin);

private:
  TH3S* mNumerator[5];
  TH3S* mDenominator[5];
  StHbt3DHisto* mQinvHisto[5];
  StHbt3DHisto* mCoulHisto[5];

  StHbtCoulomb* mCorrection; //!

  Int_t GetKtBin(const StHbtPair*);
  
#ifdef __ROOT__ 
  ClassDef(QoslCMSCorrFctnkT, 1)
#endif
};

inline  TH3S* QoslCMSCorrFctnkT::Numerator3D(const Int_t& ktBin){return mNumerator[ktBin];}
inline  TH3S* QoslCMSCorrFctnkT::Denominator3D(const Int_t& ktBin){return mDenominator[ktBin];}
inline  StHbt3DHisto* QoslCMSCorrFctnkT::QinvHisto3D(const Int_t& ktBin){return mQinvHisto[ktBin];}
inline  StHbt3DHisto* QoslCMSCorrFctnkT::CoulHisto3D(const Int_t& ktBin){return mCoulHisto[ktBin];}

#endif
