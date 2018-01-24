#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP CEmc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CPSicalc(SEXP, SEXP, SEXP, SEXP);
extern SEXP CWTcMC(SEXP, SEXP, SEXP);
extern SEXP CWTdMC(SEXP, SEXP, SEXP);
extern SEXP MCprocedure(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"CEmc",        (DL_FUNC) &CEmc,        6},
    {"CPSicalc",    (DL_FUNC) &CPSicalc,    4},
    {"CWTcMC",      (DL_FUNC) &CWTcMC,      3},
    {"CWTdMC",      (DL_FUNC) &CWTdMC,      3},
    {"MCprocedure", (DL_FUNC) &MCprocedure, 4},
    {NULL, NULL, 0}
};

void R_init_RInSp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
