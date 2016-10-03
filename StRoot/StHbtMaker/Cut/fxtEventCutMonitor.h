#ifndef fxtEventCutMonitor_hh
#define fxtEventCutMonitor_hh

#include "StHbtMaker/Infrastructure/StHbtTypes.hh"
#include "StHbtMaker/Infrastructure/StHbtEvent.hh"
#include "StHbtMaker/Base/StHbtCutMonitor.hh"
#include "TProfile.h"

class fxtEventCutMonitor : public StHbtCutMonitor{

private:
  StHbt2DHisto*  mVertexYvsVertexX;  
  StHbt1DHisto*  mVertexZ; 
  StHbt1DHisto*  mRefMult;
  StHbt1DHisto*  mNumberOfTofMatches;

public:
  fxtEventCutMonitor();
  fxtEventCutMonitor(const char* TitCutMoni, const char* title);
  virtual ~fxtEventCutMonitor();


  virtual StHbtString Report(); 
  virtual void Fill(const StHbtEvent* event);
  virtual void Finish();

  // These dummy Fill() functions were introduced to remove a compiler
  //   warning related to overloaded base-class Fill() functions being 
  //   hidden by a single version of Fill() in this derived class
  void Fill(const StHbtParticleCollection* d) {;}
  void Fill(const StHbtEvent *d1, const StHbtParticleCollection* d2) {;}
  void Fill(const StHbtPair* d) {;}
  void Fill(const StHbtKink* d) {;}
  void Fill(const StHbtV0* d) {;}
  void Fill(const StHbtTrack* d) {;}

  StHbt2DHisto* VertexYvsVertexX() {return mVertexYvsVertexX;}
  StHbt1DHisto* VertexZ() {return mVertexZ;}
  StHbt1DHisto* RefMult() {return mRefMult;}
  StHbt1DHisto* NumberOfTofMatches() {return mNumberOfTofMatches;}
  
#ifdef __ROOT__ 
 ClassDef(fxtEventCutMonitor, 1)
#endif
};

#endif

