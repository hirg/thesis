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

#ifndef QualityvsQinv_hh
#define QualityvsQinv_hh

#include "StHbtMaker/Base/StHbtCorrFctn.hh"

class QualityvsQinv : public StHbtCorrFctn {
public:
  QualityvsQinv(char* title, const int& nbinsX, const float& XLo, const float& XHi,
		      const int& nbinsY, const float& YLo, const float& YHi);
  virtual ~QualityvsQinv();

  virtual StHbtString Report();
  virtual void AddRealPair(const StHbtPair*);
  virtual void AddMixedPair(const StHbtPair*);

  virtual void Finish();

  StHbt2DHisto* Numerator2D();
  StHbt2DHisto* Denominator2D();
  StHbt2DHisto* Ratio2D();

private:
  StHbt2DHisto* mNumerator2D;
  StHbt2DHisto* mDenominator2D;
  StHbt2DHisto* mRatio2D;

#ifdef __ROOT__ 
  ClassDef(QualityvsQinv, 1)
#endif
};

inline  StHbt2DHisto* QualityvsQinv::Numerator2D(){return mNumerator2D;}
inline  StHbt2DHisto* QualityvsQinv::Denominator2D(){return mDenominator2D;}
inline  StHbt2DHisto* QualityvsQinv::Ratio2D(){return mRatio2D;}

#endif
