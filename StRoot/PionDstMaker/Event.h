#include "Track.h"
#include "StThreeVectorF.hh"
#include "StZdcTriggerDetector.h"

class Event : public TObject {

public:

	Event();
	~Event();

	// Get Methods
	TClonesArray*   GetPionTracks() const;
	Int_t           GetEventID() const;
	Int_t           GetRunID() const;
	UInt_t          GetRefMult() const;
	Double_t        GetRefMultCorr() const;
	Double_t        GetRefMultCorrWeight() const;
	Int_t           GetCent9() const;
	Int_t           GetCent16() const;
	UInt_t          GetRefMultPos() const;
	UInt_t          GetRefMultNeg() const;
    Int_t           GetNPionTracks() const;  
    Int_t           GetTofMult() const;  
	StThreeVectorF  GetVertexPos() const;
	Float_t		    GetZDCe() const;
	Float_t		    GetZDCw() const;
	Double_t		GetCtbMult() const;
    Int_t           GetNumberOfTracks() const;  
    Int_t           GetNumberOfGoodTracks() const;  
	Double_t		GetMagField() const;
    Float_t         GetPsi(Int_t sel, Int_t sub) const;   
    Float_t         Getq2(Int_t sel) const;  

	// Set Methods
	void SetEventID(const Int_t&);
	void SetRunID(const Int_t&);
	void SetRefMult(const Int_t&);
	void SetRefMultCorr(const Double_t&);
	void SetRefMultCorrWeight(const Double_t&);
	void SetCent9(const Int_t&);
	void SetCent16(const Int_t&);
	void SetRefMultPos(const Int_t&);
	void SetRefMultNeg(const Int_t&);
	void SetTofMult(const UInt_t&);
	void SetVertexPos(const StThreeVectorF&);
	void SetZDCe(const Float_t zdce);
	void SetZDCw(const Float_t zdcw);
    void SetCtbMultiplicity(const Double_t);
	void SetNumberOfTracks(const Int_t);
	void SetNumberOfGoodTracks(const Int_t);
    void SetMagField(const Double_t);
	void SetPsi(Float_t psi, Int_t sel, Int_t sub); 
    void SetEPMult(const UInt_t& epMult, Int_t sel, Int_t sub);
	void Setq2(Float_t q2, Int_t sel);   


    Track* AddPionTrack();
	void Clear();

private:

	Int_t		            mEventID;			    // ID of the event
	Int_t		            mRunID;				    // ID of the run
	UInt_t		            mRefMult;               // Reference multiplicity
    Double_t                mRefMultCorr;           // Corrected reference multiplicity
    Double_t                mRefMultCorrWeight;     // Refmult weight for histograms
    Int_t                   mCent9;                 // Centrality bin, 9 indices
    Int_t                   mCent16;                // Centrality bin, 16 indices
	UInt_t		            mRefMultPos;            // Reference multiplicity of positive tracks
	UInt_t		            mRefMultNeg;            // Reference multiplicity of negative tracks
	UInt_t		            mNPionTracks;           // Number of pion tracks
	UInt_t                  mTofMult;               // ToF multiplicity
	StThreeVectorF          mVertexPos;			    // primary vertex position
	Float_t		            mZDCe;				    // ZDC east
	Float_t		            mZDCw;				    // ZDC west
    Double_t                mCtbMultiplicity;       // Central Trigger Barrel Multiplicity
    Int_t                   mNumberOfTracks;        // Total number of TPC tracks
    Int_t                   mNumberOfGoodTracks;    // Total number of "good" tracks
    Double_t                mMagField;              // Magnetic field magnitude
	Float_t                 mPsi[2][3];             // Event plane, indexed by [selection][subEvent]
	Float_t                 mq2[2];                 // q2, indexed by [selection]
	TClonesArray*           fPionTracks;            // Array of primary tracks

	ClassDef(Event,1)
};

//Define the Get methods
inline TClonesArray* Event::GetPionTracks() const  {return fPionTracks;}

inline Int_t Event::GetEventID() const { return mEventID; }

inline Int_t Event::GetRunID() const { return mRunID; }

inline UInt_t Event::GetRefMult() const { return mRefMult; }

inline Double_t Event::GetRefMultCorr() const { return mRefMultCorr; }

inline Double_t Event::GetRefMultCorrWeight() const { return mRefMultCorrWeight; }

inline Int_t Event::GetCent9() const { return mCent9; }

inline Int_t Event::GetCent16() const { return mCent16; }

inline UInt_t Event::GetRefMultPos() const { return mRefMultPos; }

inline UInt_t Event::GetRefMultNeg() const { return mRefMultNeg; }

inline Int_t  Event::GetNPionTracks() const { return mNPionTracks; }  

inline Int_t  Event::GetTofMult() const { return mTofMult; }  

inline StThreeVectorF Event::GetVertexPos() const { return mVertexPos; }

inline Float_t  Event::GetZDCe() const { return mZDCe; }

inline Float_t  Event::GetZDCw() const { return mZDCw; }  

inline Double_t  Event::GetCtbMult() const { return mCtbMultiplicity; }  

inline Int_t  Event::GetNumberOfTracks() const { return mNumberOfTracks; }  

inline Int_t  Event::GetNumberOfGoodTracks() const { return mNumberOfGoodTracks; }  

inline Double_t  Event::GetMagField() const { return mMagField; }  

inline Float_t  Event::GetPsi(Int_t sel, Int_t sub) const { return mPsi[sel][sub]; }  

inline Float_t  Event::Getq2(Int_t sel) const { return mq2[sel]; }  

//Define the Set methods
inline void Event::SetEventID(const Int_t& id) { mEventID = id; }

inline void Event::SetRunID(const Int_t& id) { mRunID = id; }

inline void Event::SetRefMult(const Int_t& mult) { mRefMult = mult; }

inline void Event::SetRefMultCorr(const Double_t& mult) { mRefMultCorr = mult; }

inline void Event::SetRefMultCorrWeight(const Double_t& weight) { mRefMultCorrWeight = weight; }

inline void Event::SetCent9(const Int_t& cent) { mCent9 = cent; }

inline void Event::SetCent16(const Int_t& cent) { mCent16 = cent; }

inline void Event::SetRefMultPos(const Int_t& mult) { mRefMultPos = mult; }

inline void Event::SetRefMultNeg(const Int_t& mult) { mRefMultNeg = mult; }

inline void Event::SetTofMult(const UInt_t& tm) { mTofMult = tm; }

inline void Event::SetVertexPos(const StThreeVectorF& vertexPos) {
	mVertexPos = vertexPos; }

inline void Event::SetZDCe(const Float_t zdce) { mZDCe = zdce; }

inline void Event::SetZDCw(const Float_t zdcw) { mZDCw = zdcw; }

inline void Event::SetCtbMultiplicity(const Double_t ctb) { mCtbMultiplicity = ctb; }

inline void Event::SetNumberOfTracks(const Int_t numTracks) { mNumberOfTracks = numTracks; }

inline void Event::SetNumberOfGoodTracks(const Int_t numTracks) { mNumberOfGoodTracks = numTracks; }

inline void Event::SetMagField(const Double_t mag) { mMagField = mag; }

inline void Event::SetPsi(Float_t psi, Int_t sel, Int_t sub) { mPsi[sel][sub] = psi; }  

inline void Event::Setq2(Float_t q2, Int_t sel) { mq2[sel] = q2; }  
