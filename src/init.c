#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _MixPLNClust_aitkens_accel(SEXP, SEXP, SEXP);
extern SEXP _MixPLNClust_compute_F_matrices(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MixPLNClust_compute_iOsigO(SEXP, SEXP, SEXP);
extern SEXP _MixPLNClust_compute_pi_g(SEXP, SEXP);
extern SEXP _MixPLNClust_create_lib_mat_full(SEXP, SEXP);
extern SEXP _MixPLNClust_invert_matrices(SEXP, SEXP, SEXP);
extern SEXP _MixPLNClust_PD_check(SEXP, SEXP);
extern SEXP _MixPLNClust_update_g_params(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MixPLNClust_update_mu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MixPLNClust_update_sig(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MixPLNClust_update_sig2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_MixPLNClust_aitkens_accel",       (DL_FUNC) &_MixPLNClust_aitkens_accel,        3},
    {"_MixPLNClust_compute_F_matrices",  (DL_FUNC) &_MixPLNClust_compute_F_matrices,   8},
    {"_MixPLNClust_compute_iOsigO",      (DL_FUNC) &_MixPLNClust_compute_iOsigO,       3},
    {"_MixPLNClust_compute_pi_g",        (DL_FUNC) &_MixPLNClust_compute_pi_g,         2},
    {"_MixPLNClust_create_lib_mat_full", (DL_FUNC) &_MixPLNClust_create_lib_mat_full,  2},
    {"_MixPLNClust_invert_matrices",     (DL_FUNC) &_MixPLNClust_invert_matrices,      3},
    {"_MixPLNClust_PD_check",            (DL_FUNC) &_MixPLNClust_PD_check,             2},
    {"_MixPLNClust_update_g_params",     (DL_FUNC) &_MixPLNClust_update_g_params,     14},
    {"_MixPLNClust_update_mu",           (DL_FUNC) &_MixPLNClust_update_mu,            8},
    {"_MixPLNClust_update_sig",          (DL_FUNC) &_MixPLNClust_update_sig,          12},
    {"_MixPLNClust_update_sig2",         (DL_FUNC) &_MixPLNClust_update_sig2,         12},
    {NULL, NULL, 0}
};

void R_init_MixPLNClust(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}