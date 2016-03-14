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
    register int   m,l,iter,i,k;
    double         s,r,p,g,f,dd,c,b;
    
    for(i = 1; i < n; i++) e[i-1] = e[i];
    e[n] = 0.0;
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


/*
	 ** The function 
         **           jacobi_rot()
	 ** A helping function for jacobi making the actual rotations
	 ** a is the matrix to be rotated, s is sine of the rotation
	 ** angle, tau is s/(1 + c) where c is cosine of the angle.
	 ** The integers i-l denotes the matrix element to be
	 ** rotated.
	 **
         */ 

inline void jacobi_rot(double** a, double s, double tau, int i, int j, int k, int l){
  double g,h;

  g = a[i][j];
  h = a[k][l];
  a[i][j] = g - s * (h + g*tau);
  a[k][l] = h + s * (g - h*tau);

}//End function jacobi_rot

       /*
       ** The function 
       **              jacobi()
       ** Computes the eigenvalues and eigenvectors of the square symmetric matrix
       ** A  by use of the Jacobi method.
       ** It puts the eigenvalues in d and eigenvectors in v.
       ** n is an integer denoting the size of A
       ** The function is as in the Numerical recipe
       */
void jacobi(double** a, double* d, double** v, int n, int &nrot){
  int i,j, ip, iq;
  double tresh, theta, tau, t, sm, s, h, g, c;
  
  double* b = new double[n];
  double* z = new double[n];
  for(ip = 0; ip < n; ip++){
    for(iq = 0; iq < n; iq++){
      v[ip][iq] = 0.0;          //Initializing v to the identity matrix
      v[ip][ip] = 1.0;
    }
  }

  for(ip = 0; ip <n; ip++){
    b[ip] = d[ip] = a[ip][ip];  //Initializing d and b to the diagonal of a
    z[ip] = 0.0;                //z will accumulate terms of the form
                                //t*a[ip][iq]
  }


  nrot = 0;
  for(i = 1; i <= 50; i++){
    sm = 0.0;
    for(ip = 0; ip < n - 1; ip++){
      for(iq = ip + 1; iq < n; iq++){
	sm += fabs(a[ip][iq]);   //Sum magnitude of off-diagonal elements
      }
    }
    if(sm == 0.0){
      return;                  //The normal return at convergence
    }
    if(i < 4){
      tresh = 0.2 * sm/(n*n);  //On the first four sweeps
    }else{
      tresh = 0.0;             //... thereafter
    }
    for(ip = 0; ip < n-1; ip++){
      for(iq = ip + 1; iq < n; iq++){
	g = 100.0*fabs(a[ip][iq]);
	//After four sweeps we skip the rotation if the off-diagonal element is small
	if(i >4 && (fabs(d[ip]) + g) == fabs(d[ip]) 
	   && (fabs(d[iq]) + g) == fabs(d[iq])){
	  a[ip][iq] = 0.0;
	}else if(fabs(a[ip][iq]) > tresh){
	  h = d[iq] - d[ip];
	  if((fabs(h) + g) == fabs(h)){
	    t = (a[ip][iq])/h;
	  }else{
	    theta = 0.5*h/(a[ip][iq]);
	    t = 1.0/(fabs(theta) + sqrt(1.0 + theta*theta));
	    if(theta < 0.0){
	      t = -t;
	    }
	  }
	  c = 1.0/sqrt(1 + t*t);
	  s = t*c;
	  tau = s/(1.0 + c);
	  h = t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq] = 0.0;
	  for(j = 0; j < ip; j++){
	    jacobi_rot(a, s, tau, j, ip, j, iq); // Rotations for 0 <= j < ip
	  }
	  for(j = ip + 1; j < iq; j++){
	    jacobi_rot(a, s, tau, ip, j, j, iq); // Rotations for ip < j < iq
	  }
	  for(j = iq + 1; j < n; j++){
	    jacobi_rot(a, s, tau, ip, j, iq, j); // Rotations for q < j < n
	  }
	  for(j = 0; j < n; j++){
	    jacobi_rot(v, s, tau, j, ip, j, iq); //Updating v
	  }
	  nrot++;
	}
      }
    }
    for(ip = 0; ip < n; ip ++){
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  printf("\n\nToo many iterations in routine jacobi.\n");
  exit(1); 
}//End function jacobi()