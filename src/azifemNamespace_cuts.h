#ifndef AZIFEMNAMESPACE_CUTS_H
#define AZIFEMNAMESPACE_CUTS_H

namespace azifem {
    const Float_t UUzdcBinEdges[3] = {-1, 436.5, 523.5}; // UU193 - 0%, 0.25%, 0.5%
    const Float_t AuAuzdcBinEdges[3] = {-1, 1194.5, 1389.5}; // AuAu 200 - 0, 0.25, 0.5
    const Float_t q2BinEdges[6] = {0,0.4925,0.7525,1.0075,1.3275,10}; 
    const Float_t UUmult[6] = {200, 570, 589, 605, 623, 10000}; // UU 193
    const Float_t AuAumult[6] = {200, 501, 522, 538, 556, 10000}; // AuAu 200

    Int_t getZdcBin(const Float_t zdc, const Bool_t uuNotAuAu);
}

#endif
