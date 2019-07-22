#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(gintp)(double *grad, double *gvar, double *bt2, double *beta2, double *rho, double *omega, int *npar, int *link, int *maxy, double *x2, int *y2, double *theta2, double *work2, int *n, double *li, double *ls, double *epsabs, double *epsrel, int *key, int *limit);
extern void F77_NAME(gintp0)(double *grad, double *gvar, double *bt2, double *beta2, double *omega, int *npar, int *link, int *maxy, double *x2, int *y2, double *theta2, double *work2, int *n, double *li, double *ls, double *epsabs, double *epsrel, int *key, int *limit);
extern void F77_NAME(intp)(double *logL, double *bt2, double *beta2, double *rho, double *omega, int *npar, int *link, int *maxy, double *x2, int *y2, double *theta2, double *work2, int *n, double *li, double *ls, double *epsabs, double *epsrel, int *key, int *limit);
extern void F77_NAME(intp0)(double *logL, double *bt2, double *beta2,  double *omega, int *npar, int *link, int *maxy, double *x2, int *y2, double *theta2, double *work2, int *n, double *li, double *ls, double *epsabs, double *epsrel, int *key, int *limit);
extern void F77_NAME(pssgrd)(double *grad, double *beta, double *rho, int *npar, double *x, int *y, double *theta, double *work, int *n, double *f, int *link);
extern void F77_NAME(pssgrd0)(double *grad, double *beta, int *npar, double *x, int *y, double *theta, double *work, int *n, int *link);
extern void F77_NAME(psslik)(double *logL, double *beta, double *rho, int *np, double *x, int *y, double *theta, double *work, int *n, double *fact, int *link);
extern void F77_NAME(psslik0)(double *logL, double *beta, int *np, double *x, int *y, double *theta, double *work, int *n, double *fact, int *link);

static const R_FortranMethodDef FortranEntries[] = {
    {"gintp",   (DL_FUNC) &F77_NAME(gintp),   20},
    {"gintp0",  (DL_FUNC) &F77_NAME(gintp0),  19},
    {"intp",    (DL_FUNC) &F77_NAME(intp),    19},
    {"intp0",   (DL_FUNC) &F77_NAME(intp0),   18},
    {"pssgrd",  (DL_FUNC) &F77_NAME(pssgrd),  11},
    {"pssgrd0", (DL_FUNC) &F77_NAME(pssgrd0),  9},
    {"psslik",  (DL_FUNC) &F77_NAME(psslik),  11},
    {"psslik0", (DL_FUNC) &F77_NAME(psslik0), 10},
    {NULL, NULL, 0}
};

void R_init_cold(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
