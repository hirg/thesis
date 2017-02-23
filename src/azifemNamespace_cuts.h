#ifndef AZIFEMNAMESPACE_CUTS_H
#define AZIFEMNAMESPACE_CUTS_H

namespace azifem {
    // The five bin edges define four ZDC bins:
    //      0-.25%, 0.25-0.50%, 0.50%-0.75%, and 0.75-1.0%
    const Float_t UUzdcBinEdges[5] = {-1, 436.5, 523.5, 592.5, 999999};
    const Float_t AuAuzdcBinEdges[5] = {-1, 1194.5, 1389.5, 1521.5, 999999};

    const Float_t q2BinEdges[6] = {0, 0.4925, 0.7525, 1.0075, 1.3275, 10}; 
    const Float_t UUmult[6] = {200, 570, 589, 605, 623, 10000}; // UU 193
    const Float_t AuAumult[6] = {200, 501, 522, 538, 556, 10000}; // AuAu 200

    Int_t getZdcBin(const Float_t zdc, const Bool_t uuNotAuAu);
    Int_t getq2Bin(const Double_t q2);
}

#endif
