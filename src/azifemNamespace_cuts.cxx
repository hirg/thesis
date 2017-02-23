#ifndef AZIFEMNAMESPACE_CUTS_CXX
#define AZIFEMNAMESPACE_CUTS_CXX

#include "azifemNamespace_cuts.h"

Int_t azifem::getZdcBin(const Float_t zdc, const Bool_t uuNotAuAu)
{

    Int_t zdcBin = -1;
    Float_t zdcLow = 999999;
    Float_t zdcHigh = -1;

    for(Int_t i = 0; i <= 4; i++)
    {
        if(uuNotAuAu) {
            zdcLow = UUzdcBinEdges[i];
            zdcHigh = UUzdcBinEdges[i+1];
        } else {
            zdcLow = AuAuzdcBinEdges[i];
            zdcHigh = AuAuzdcBinEdges[i+1];
        }

        if( (zdc > zdcLow) && (zdc <= zdcHigh) ) { zdcBin = i; }

    }

    return zdcBin;
}

Int_t azifem::getRefmultBin(const Float_t rm, const Bool_t uuNotAuAu)
{

    Int_t rmBin = -1;
    Float_t rmLow = 999999;
    Float_t rmHigh = -1;

    for(Int_t i = 0; i <= 5; i++)
    {
        if(uuNotAuAu) {
            rmLow = UUmult[i];
            rmHigh = UUmult[i+1];
        } else {
            rmLow = AuAumult[i];
            rmHigh = AuAumult[i+1];
        }

        if( (rm > rmLow) && (rm <= rmHigh) ) { rmBin = i; }

    }

    return rmBin;
}

Int_t azifem::getq2Bin(const Double_t q2)
{

    Int_t q2Bin = -1;

    for(Int_t i = 0; i <= 5; i++)
    {
        if( (q2 > q2BinEdges[i]) && (q2 <= q2BinEdges[i+1]) ) {q2Bin = i;}

    }

    return q2Bin;

}

#endif
