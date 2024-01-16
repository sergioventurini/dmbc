#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "dmbc.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
  CALLDEF(dmbc_mcmc, 22),
  CALLDEF(dmbc_pack_par, 7),
  CALLDEF(dmbc_relabel, 16),
  {NULL, NULL, 0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_dmbc(DllInfo *dll) {
  R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
