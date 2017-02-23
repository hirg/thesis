#ifndef AZIFEMNAMESPACE_RUNID_H
#define AZIFEMNAMESPACE_RUNID_H
#include <map>

namespace azifem {

    // This is the list of Au+Au run ID's. You can recreate this list with:
    //      {QUERY} | awk '{print $1}' | sort | uniq
    //  where {QUERY} is: 
    //      get_file_list.pl -onefile -keys 'runnumber' -delim " " -cond
    //      "production=P11id,trgsetupname=AuAu200_production_2011,
    //      runnumber[]{RANGE},filetype=daq_reco_MuDst,filename~st_physics,storage!=HPSS"
    //      -limit "0"
    //
    //  where {RANGE} is either of:
    //      12126101-12138024
    //      12149028-12171016
    const Int_t AuAuRunIds[];

    // This is the list of U+U 1% run ID's. You can recreate this list with:
    //      {QUERY} | awk '{print $1}' | sort | uniq
    //  where {QUERY} is: 
    //      get_file_list.pl -onefile -keys 'runnumber' -delim " " -cond 
    //      "filetype=daq_reco_mudst,trgsetupname=UU_production_2012,
    //      production=P12id,filename~centralpro,tpx=1,sanity=1,storage!=HPSS"
    //      -limit "0"
    const Int_t UURunIds[];

    // These are the lengths of the arrays defined above. Yes, it's bad that 
    // these are hard-coded, but these numbers should never change, and using
    // sizeof(array) inside a namespace just returns the size of one element
    const Int_t nAuAuRuns = 1013;
    const Int_t nUURuns = 722;

    map<Int_t, Int_t> makeRunMap(Bool_t uuNotAuAu = kTRUE);
}

#endif
