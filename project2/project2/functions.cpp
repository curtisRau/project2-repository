//
//  functions.cpp
//  project2
//
//  Created by Curtis Rau on 3/4/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#include "functions.hpp"
#include <iostream>         // For std::cout
#include "math.h"           // For sqrt function.


namespace function {
    
    // A function to print a matrix for debugging purposes
    void printMatrix (double** A, unsigned int m, unsigned int n) {
        for (unsigned int i = 0; i<m; i++) {
            for (unsigned int j = 0; j<n; j++) {
                std::cout << A[i][j] << "\t\t\t";
            }
            std::cout << "\r";
        }
    }
    
    // A function to print a vector d
    void printVector (double* d, unsigned int n) {
        for (unsigned int i = 0; i<n; i++) {
            std::cout << d[i] << " ";
        }
    }
    
    // Generate a Constant Vector of length "length", initialized
    // with all values of value "val".  This function is vectorizable.
    double* generateConstantVector (unsigned int length, double val) {
        double* vec = new double[length];
        for (unsigned int i = 0; i<length; i++) {
            vec[i] = val;
        }
        return vec;
    }
    
    // Generate Tridiagonal Matrix with Constant Arguments Exact:
    // Returns a n by n tridiagonal matrix.  The sub-diagonal elements are
    // value a, diagonal elements are value b, super-diagonal elements are
    // value c.  All other elements are 0.  This function is non-vectorizable,
    // due to the if else statements in the for loop, so it will be slow, but
    // it generates a matrix exactly how it advertizes.
    double** genTridiagMatConstArgsExact (unsigned int n, double a, double b, double c) {
        double** A = new double*[n];
        for(unsigned int i = 0; i < n; i++) {
            
            A[i] = new double[n];
            for (unsigned int j = 0; j<n; j++) {
                
                if (j == (i-1)) {
                    A[i][j] = a;
                } else if (j == i) {
                    A[i][j] = b;
                } else if (j == (i+1)) {
                    A[i][j] = c;
                } else {
                    A[i][j] = 0.0;
                }
            }
        }
        return A;
    }
    
    // Generate Tridiagonal Matrix with Constant Arguments Fast:
    // Returns a n by n tridiagonal matrix.  The sub-diagonal elements are
    // value a, diagonal elements are value b, super-diagonal elements are
    // value c.  All other elements are roughly 0.  This function assumes
    // values of an array are initialized to 0, which they are sometimes not.
    // This function is fast though, being vectorizable, and could even be
    // parallelized (with some tweeking?).
    double** genTridiagMatConstArgsFast (unsigned int n, double a, double b, double c) {
        double** A = new double*[n];
        A[0]       = new double[n];
        A[0][0]    = b;
        A[0][1]    = c;
        
        for(unsigned int i = 1; i < n-1; i++) {
            A[i]      = new double[n];
            A[i][i-1] = a;
            A[i][i]   = b;
            A[i][i+1] = c;
        }
        
        A[n-1]      = new double[n];
        A[n-1][n-2] = a;
        A[n-1][n-1] = b;
        
        return A;
    }
    
    // Generate Tridiagonal Matrix with Vector Arguments Exact:
    // Returns a n by n tridiagonal matrix.  The sub-diagonal elements are
    // value a, diagonal elements are value b, super-diagonal elements are
    // value c.  All other elements are 0.  This function is non-vectorizable,
    // due to the if else statements in the for loop, so it will be slow, but
    // it generates a matrix exactly how it advertizes.
    double** genTridiagMatVectArgsExact (unsigned int n, double* a, double* b, double* c) {
        double** A = new double*[n];
        for(unsigned int i = 0; i < n; i++) {
            
            A[i] = new double[n];
            for (unsigned int j = 0; j<n; j++) {
                
                if (j == (i-1)) {
                    A[i][j] = a[i];
                } else if (j == i) {
                    A[i][j] = b[i];
                } else if (j == (i+1)) {
                    A[i][j] = c[i];
                } else {
                    A[i][j] = 0.0;
                }
            }
        }
        return A;
    }
    
    // Generate Tridiagonal Matrix with Vector Arguments Fast:
    // Returns a n by n tridiagonal matrix.  The sub-diagonal elements are
    // value a, diagonal elements are value b, super-diagonal elements are
    // value c.  All other elements are 0.  This function is vectorizable.
    double** genTridiagMatVectArgsFast (unsigned int n, double* a, double* b, double* c) {
        double** A = new double*[n];
        A[0]       = new double[n];
        A[0][0]    = b[0];
        A[0][1]    = c[0];
        
        for(unsigned int i = 1; i < n-1; i++) {
            A[i]      = new double[n];
            A[i][i-1] = a[i];
            A[i][i]   = b[i];
            A[i][i+1] = c[i];
        }
        
        A[n-1]      = new double[n];
        A[n-1][n-2] = a[n-2];
        A[n-1][n-1] = b[n-1];
        
        return A;
    }
    
    // Functions specific to Jacobi's Method
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double frobeniusNorm (double** A, unsigned int m, unsigned int n) {
        double sum = 0.0;
        for (unsigned int i = 0; i<m; i++) {         // The order of this sum will matter.  row-column indexed?
            for (unsigned int j = 0; j<n; j++) {
                sum += A[i][j] * A[i][j];
            }
        }
        return sqrt(sum);
    }
    
    //
    void jacobiRotation (double** A, unsigned int matrixSize, unsigned int i, unsigned int j, float theta) {
        double s = sin(theta);
        double c = cos(theta);
        
        A[i][i] = (c*c) * A[i][i] - (2*s*c) * A[i][j] + (s*s) * A[j][j];
        A[j][j] = (s*s) * A[i][i] + (2*s*c) * A[i][j] + (c*c) * A[j][j];
        A[i][j] = (c*c - s*s) * A[i][j] + (s*c) * (A[i][i] - A[j][j]);
        A[j][i] = A[i][j];
        
        for (unsigned int k = 0; (k < matrixSize) && (k != i) && (k != j); k++) {
            A[i][k] = c * A[i][k] - s * A[j][k];
            A[k][i] = A[i][k];
            A[j][k] = s * A[i][k] + c * A[j][k];
            A[k][j] = A[j][k];
        }
    }
    
    // This function is not vectorizable as stands, but probably could be.
    double off (double** A, unsigned int matrixSize) {
        // The square root of the sum of the squares of the off diagonal elements.
        double sum = 0.0;
        for (unsigned int i=0; i < matrixSize; i++) {
            for (unsigned int j=0; (j<matrixSize && i!=j); j++) {
                sum += (A[i][j]) * (A[i][j]);
            }
        }
        return sqrt(sum);
    }
    
    // This function determines the indicies of the largest absolute valued, off diagonal, element of a matrix.
    // A = the matrix
    // m is the number of rows of A, n is the number of columns of A.
    // p is the row, q is the column of the largest absolute valued element.
    //unsigned int x;
    //    unsigned int y;
    //    unsigned int *p = &x;
    //    unsigned int *q = &y;
    void maxOffDiagnalElement (double** A, unsigned int matrixSize, double* value, unsigned int* p, unsigned int* q) {
        *value = 0.0;
        for (unsigned int i = 0; i<matrixSize; i++) {
            for (unsigned int j = 0; (j<matrixSize) && (j != i); j++) {
                if ( fabs(A[i][j]) > *value ) {
                    *value = fabs(A[i][j]);
                    *p = i;
                    *q = j;
                }
            }
        }
    }
    
    
    // A function that returns the smallest on-diagonal matrix element.
    // This is used to extract the smallest eigenvalue.
    double minDiagonalElement (double** A, unsigned int matrixSize) {
        double min = A[0][0];
        for (unsigned int i = 1; i < matrixSize; i++) {
            if (A[i][i] < min) {
                min = A[i][i];
            }
        }
        return min;
    }
}
