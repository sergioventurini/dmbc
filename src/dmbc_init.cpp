// dmbc_init.cpp

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "dmbc.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(DMBC_MCMC, 21),
    CALLDEF(RELABEL, 15),
    CALLDEF(PACK_PAR, 7),
    {NULL, NULL, 0}
};

extern "C" void attribute_visible R_init_dmbc(DllInfo* dll){
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
