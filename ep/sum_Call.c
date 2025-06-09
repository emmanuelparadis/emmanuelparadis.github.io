#include <R.h>
#include <Rinternals.h>

SEXP sum_Call(SEXP x)
{
    double *xp, s = 0.;
    SEXP ans;
    PROTECT(x = coerceVector(x, REALSXP));
    xp = REAL(x);
    PROTECT(ans = allocVector(REALSXP, 1));
    for (int i = 0; i < LENGTH(x); i++) s += xp[i];
    REAL(ans)[0] = s;
    UNPROTECT(2);
    return ans;
}
