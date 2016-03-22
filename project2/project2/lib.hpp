    /*
     * The definition module lib.h
     * for the library function common for all C programs.
     */

     // Standard ANSI-C++ include files 


#include <iostream>
#include <new>
//#include <cstdio>
//#include <cstdlib>
#include <cmath>
#include <cstring>

static float sqrarg;
#define SQR(a)    ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)          // Required by pythag
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))                        // Required by tqli

void tqli(double *, double *, int, double **);
double pythag(double, double);                                          // Required by tqli



