    /*
    ** The library program  module
    **             lib.cpp
    **


void tqli(double d[], double e[], int n, double **z)
    ** determine eigenvalues and eigenvectors of a real symmetric
    ** tri-diagonal matrix, or a real, symmetric matrix previously
    ** reduced by function tred2[] to tri-diagonal form. On input,
    ** d[] contains the diagonal element and e[] the sub-diagonal
    ** of the tri-diagonal matrix. On output d[] contains the
    ** eigenvalues and  e[] is destroyed. If eigenvectors are
    ** desired z[][] on input contains the identity matrix. If
    ** eigenvectors of a matrix reduced by tred2() are required,
    ** then z[][] on input is the matrix output from tred2().
    ** On output, the k'th column returns the normalized eigenvector
    ** corresponding to d[k]. 
    ** The function is modified from the version in Numerical recipe.

double pythag(double a, double b)
    ** The function is modified from the version in Numerical recipe.

void jacobi(double** a, double* d, double** v, int n, int& nrot)
    ** Computes the eigenvalues and eigenvectors of the square symmetric matrix
    ** a  by use of the Jacobi method.
    ** It puts the eigenvalues in d and eigenvectors in v.
    ** n is an integer denoting the size of a
    **nrot keeps track of the number of rotations
    ** The function is as in the Numerical recipe

void jacobi_rot(double** a, double s, double tau, int i, int j, int k, int l)
	 ** A helping function for jacobi making the actual rotations
	 ** a is the matrix to be rotated, s is sine of the rotation
	 ** angle, tau is s/(1 + c) where c is cosine of the angle.
	 ** The integers i-l denotes the matrix element to be
	 ** rotated.
     */

#include "lib.hpp"

    /*
    ** The function
    **                 tqli()
    ** determine eigenvalues and eigenvectors of a real symmetric
    ** tri-diagonal matrix, or a real, symmetric matrix previously
    ** reduced by function tred2[] to tri-diagonal form. On input,
    ** d[] contains the diagonal element and e[] the sub-diagonal
    ** of the tri-diagonal matrix. On output d[] contains the
    ** eigenvalues and  e[] is destroyed. If eigenvectors are
    ** desired z[][] on input contains the identity matrix. If
    ** eigenvectors of a matrix reduced by tred2() are required,
    ** then z[][] on input is the matrix output from tred2().
    ** On output, the k'th column returns the normalized eigenvector
    ** corresponding to d[k]. 
    ** The function is modified from the version in Numerical recipe.
    */

void tqli(double *d, double *e, int n, double **z)
{
    int    m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    
    for(i = 1; i < n; i++)
        e[i-1] = e[i];
        e[n]   = 0.0;
        for(l = 0; l < n; l++) {
            iter = 0;
            do {
                for(m = l; m < n-1; m++) {
                    dd = fabs(d[m]) + fabs(d[m+1]);
                    if((double)(fabs(e[m])+dd) == dd) break;
                }
                if(m != l) {
                    if(iter++ == 30) {
                        printf("\n\nToo many iterations in tqli.\n");
                        exit(1);
                    }
                    g = (d[l+1] - d[l])/(2.0 * e[l]);
                    r = pythag(g,1.0);
                    g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
                    s = c = 1.0;
                    p = 0.0;
                    for(i = m-1; i >= l; i--) {
                        f      = s * e[i];
                        b      = c*e[i];
                        e[i+1] = (r=pythag(f,g));
                        if(r == 0.0) {
                            d[i+1] -= p;
                            e[m]    = 0.0;
                        break;
                    }
                    s      = f/r;
                    c      = g/r;
                    g      = d[i+1] - p;
                    r      = (d[i] - g) * s + 2.0 * c * b;
                    d[i+1] = g + (p = s * r);
                    g      = c * r - b;
                    for(k = 0; k < n; k++) {
                        f         = z[k][i+1];
                        z[k][i+1] = s * z[k][i] + c * f;
                        z[k][i]   = c * z[k][i] - s * f;
                    } /* end k-loop */
                } /* end i-loop */
                if(r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l]  = g;
                e[m]  = 0.0;
            } /* end if-loop for m != 1 */
        } while(m != l);
    } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */


double pythag(double a, double b)
{
  double absa = fabs(a);
  double absb = fabs(b);
    if (absa > absb) {
        return absa * sqrt(1.0 + SQR(absb/absa));
    } else {
        if (absb == 0.0) {
            return 0.0;
        } else {
            return absb * sqrt(1.0 + SQR(absa/absb));
        }
    }
}
// End: function pythag(), (C) Copr. 1986-92 Numerical Recipes Software )%.


void orthogonalityCheck ( ) {
    
}