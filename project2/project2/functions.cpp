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
#include <fstream>              // for working with files.

namespace function {
    
    // --description--
    void plotVector (double* vec, unsigned int xMin, unsigned int xMax, unsigned int xScale, double yScale) {
        for (unsigned int i = xMin; i < xMax; i += xScale) {
            for (unsigned int j = 0; j < yScale * vec[i] * vec[i]; j++) {           // Plots the square of the vec[i]
                std::cout << "-";
            }
            std::cout << std::endl;
        }
    }
    
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
        std::cout << "\n";
    }

    // A function to print a diagonals of matrix A
    void printDiagonals (double** A, unsigned int n) {
        for (unsigned int i = 0; i<n; i++) {
            std::cout << A[i][i] << " ";
        }
        std::cout << "\n";
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
    
    // Computational Overhead:
    // -- Additions/Subtractions: 7 + (N-2)*2
    // -- Multiplications       : 
    // -- Memory Read           :
    // -- Memory Write          :
    void jacobiRotation (double** A, unsigned int matrixSize, unsigned int i, unsigned int j, float theta) {
        if (i != j) {
            double s = sin(theta);                                                      // Read = 1; Write = 1;
            double c = cos(theta);                                                      // Read = 1; Write = 1;
            
            A[i][i] = (c*c) * A[i][i] - (2.0*s*c) * A[i][j] + (s*s) * A[j][j];            // Mult = 7; Add = 2; Read =
            A[j][j] = (s*s) * A[i][i] + (2.0*s*c) * A[i][j] + (c*c) * A[j][j];            // Mult = 7; Add = 2
            A[i][j] = (c*c - s*s) * A[i][j] + (s*c) * (A[i][i] - A[j][j]);              // Mult = 5; Add = 3
            A[j][i] = A[i][j];
            
            for (unsigned int k = 0; k < matrixSize; k++) {   // Number of executions = matrixSize - 2
                if ((k != i) && (k != j)) {
                    A[i][k] = c * A[i][k] - s * A[j][k];                                    // Mult = 2; Add = 1;
                    A[k][i] = A[i][k];
                    A[j][k] = s * A[i][k] + c * A[j][k];                                    // Mult = 2; Add = 1;
                    A[k][j] = A[j][k];
                }
            }
        }
    }
    
    // Computational Overhead:
    // -- Additions/Subtractions: 7 + (N-2)*2
    // -- Multiplications       :
    // -- Memory Read           :
    // -- Memory Write          :
    void jacobiRotationSC (double** A, unsigned int matrixSize, unsigned int i, unsigned int j, double s, double c) {
        if (i != j) {
            
            double AII = A[i][i];
            double AJJ = A[j][j];
            
            A[i][i] = (c*c) * AII - (2.0*s*c) * A[i][j] + (s*s) * AJJ;            // Mult = 7; Add = 2; Read =
            A[j][j] = (s*s) * AII + (2.0*s*c) * A[i][j] + (c*c) * AJJ;            // Mult = 7; Add = 2
            A[i][j] = (c*c - s*s) * A[i][j] + (s*c) * (AII - AJJ);              // Mult = 5; Add = 3
            A[j][i] = A[i][j];
            
            double AIK;
            for (unsigned int k = 0; k < matrixSize; k++) {   // Number of executions = matrixSize - 2
                if ((k != i) && (k != j)) {
                    AIK = A[i][k];
                    A[i][k] = c * AIK - s * A[j][k];                                    // Mult = 2; Add = 1;
                    A[k][i] = A[i][k];
                    A[j][k] = s * AIK + c * A[j][k];                                    // Mult = 2; Add = 1;
                    A[k][j] = A[j][k];
                }
            }
        }
    }
    
    // This function is not vectorizable as stands, but probably could be.
    double off (double** A, unsigned int matrixSize) {
        // The square root of the sum of the squares of the off diagonal elements.
        double sum = 0.0;
        for (unsigned int i=0; i < matrixSize; i++) {
            for (unsigned int j=0; j < matrixSize; j++) {
                if (i!=j) {
                    sum += (A[i][j]) * (A[i][j]);
                }
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
    
    // THE MATRIX IS SYMMETRIC SO WE CAN LOOP OVER JUST THE UPPER OR LOWER HALF.
    // THE MATRIX IS SYMMETRIC SO WE CAN LOOP OVER JUST THE UPPER OR LOWER HALF.
    // THE MATRIX IS SYMMETRIC SO WE CAN LOOP OVER JUST THE UPPER OR LOWER HALF.
    // THE MATRIX IS SYMMETRIC SO WE CAN LOOP OVER JUST THE UPPER OR LOWER HALF.
    // THE MATRIX IS SYMMETRIC SO WE CAN LOOP OVER JUST THE UPPER OR LOWER HALF.
    void maxOffDiagnalElement (double** A, unsigned int matrixSize, double* value, unsigned int* p, unsigned int* q) {
        *value = 0.0;
        for (unsigned int i = 0; i < matrixSize; i++) {
            for (unsigned int j = 0; j < matrixSize; j++) {// THE MATRIX IS SYMMETRIC SO WE CAN LOOP OVER JUST THE UPPER OR LOWER HALF.
                if (j != i) {// THE MATRIX IS SYMMETRIC SO WE CAN LOOP OVER JUST THE UPPER OR LOWER HALF.
                    if ( fabs(A[i][j]) > *value ) {// THE MATRIX IS SYMMETRIC SO WE CAN LOOP OVER JUST THE UPPER OR LOWER HALF.
                        *value = fabs(A[i][j]);
                        *p = i;
                        *q = j;
                    }
                }
            }
        }
    }
    
    
    // A function that returns the smallest on-diagonal matrix element.
    // This is used to extract the smallest eigenvalue.
    double minDiagonalElement (double** A, unsigned int matrixSize) {
        double min = A[0][0];
        for (unsigned int i = 1; i < matrixSize; i++) {
            if (fabs(A[i][i]) < min) {
                min = fabs(A[i][i]);
            }
        }
        return min;
    }
    
    
    // A function to returen the "numOfElem2Return" smallest values from a one dimentional
    // array "vec" of size "vecSize".
    double* minVectorElements (double* vec, unsigned int vecSize, unsigned int numOfElem2Return) {
        double* minElemsVec = new double [numOfElem2Return];
        double* XX = new double[numOfElem2Return];
        
        minElemsVec[0] = vec[0];
        XX[0]          = vec[0];
        for (unsigned int i = 0; i < vecSize; i++) {
            for (unsigned int j = 0; j < numOfElem2Return; j++) {
                
                if (vec[i] < minElemsVec[j]) {
                    minElemsVec[j] = vec[i];
                    for (NULL; j < (numOfElem2Return - 1); j++) {
                        minElemsVec[j+1] = XX[j];                           // Shift forward
                        XX[j] = minElemsVec[j];                             // Copy minElemsVec to XX
                    }
                    XX[numOfElem2Return] = minElemsVec[numOfElem2Return];   // Copy minElemsVec to XX
                }
            }
        }
        delete [] XX;
        return minElemsVec;
    }
    
    void saveMatrix4Mathematica (const char* filename, double** matrix, unsigned int matrixSizeM, unsigned int matrixSizeN) {
        std::ofstream outputFile;
        outputFile.open(filename, std::ios::out | std::ios::trunc);         // Open a file for output and overwrite current content if it exists.
        
        if (outputFile.is_open()) {                                         // If the file is open...
            for (unsigned int i = 0; i < (matrixSizeM - 1); i++) {
                outputFile << matrix[i][0];
                for (unsigned int j = 1; j < matrixSizeN; j++) {
                    outputFile << "\t" << matrix[i][j];
                }
                outputFile << "\r";
            }
            outputFile << matrix[matrixSizeM - 1][0];
            for (unsigned int j = 1; j < matrixSizeN; j++) {
                outputFile << "," << matrix[matrixSizeM - 1][j];
            }
        } else {
            std::cout << "File '" << filename << "' did not open /r";
        }
        outputFile.close();
    }
    
    // Can be vectorized.
    void saveArray4Mathematica (const char* filename, double* array, unsigned int arraySize) {
        std::ofstream outputFile;
        outputFile.open(filename, std::ios::out | std::ios::trunc);         // Open a file for output and overwrite current content if it exists.
        
        if (outputFile.is_open()) {                                         // If the file is open...
            outputFile << array[0];
            for (unsigned int i = 0; i < arraySize; i++) {
                outputFile << "\t" << array[i];
            }
        } else {
            std::cout << "File '" << filename << "' did not open /r";
        }
        outputFile.close();
    }
}
