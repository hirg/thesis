#include "Event.h"

#include "TClonesArray.h"

ClassImp(Event)


Event::Event() {
// Create an Event object.

  fPionTracks = new TClonesArray("Track", 2000);
  mNPionTracks = 0;

}


//_____________________________________________________________________________
Track* Event::AddPionTrack() {
//// This function adds a Track to the TClonesArray for Event
//
//  //_____________________________________________________________________________                                                                           
  TClonesArray &tracks = *fPionTracks;
  Track* track = new(tracks[mNPionTracks++]) Track;
  return track;
}

//_____________________________________________________________________________
void Event::Clear() {

   fPionTracks->Clear("C"); 
   mNPionTracks = 0;

}


//_____________________________________________________________________________
Event::~Event() {}
