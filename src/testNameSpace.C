#include "azifemNamespace_core.cxx"
#include <map>
// #include "../src/azifem_namespace.cxx"
void testNameSpace()
{
    Int_t a = 3;
    azifem::makeRunMap(kTRUE);
    azifem::makeRunMap(kFALSE);
    // const map<Int_t, Int_t> runMap = azifem::makeRunMap();
    // const map<Int_t, Int_t> auMap = azifem::makeRunMap(0);
    // cout << auMap[12127003] << endl;
    // cout << runMap[13136015] << endl;
    // cout << azifem::UURunIds[3] << endl;

}
