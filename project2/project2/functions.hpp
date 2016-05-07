//
//  functions.hpp
//  project2
//
//  Created by Curtis Rau on 3/4/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>

#endif

namespace function {
    void     plotVector     (double* vec, unsigned int xMin, unsigned int xMax, unsigned int xScale, double yScale);
    void     printMatrix    (double** A, unsigned int m, unsigned int n);
    void     printVector    (double* d, unsigned int n);
    void     printDiagonals (double** A, unsigned int n);
    
    double*  generateConstantVector      (unsigned int length, double val);
    double** genTridiagMatConstArgsExact (unsigned int n, double a, double b, double c);
    double** genTridiagMatConstArgsFast  (unsigned int n, double a, double b, double c);
    double** genTridiagMatVectArgsExact  (unsigned int n, double* a, double* b, double* c);
    double** genTridiagMatVectArgsFast   (unsigned int n, double* a, double* b, double* c);
    
    double   frobeniusNorm (double** A, unsigned int m, unsigned int n);
    void     jacobiRotation   (double** A, unsigned int matrixSize, unsigned int i, unsigned int j, float theta);
    void     jacobiRotationSC (double** A, unsigned int matrixSize, unsigned int i, unsigned int j, double s, double c);
    double   off (double** A, unsigned int matrixSize);
    void     maxOffDiagnalElement (double** A, unsigned int matrixSize, double* value, unsigned int* p, unsigned int* q);
    double   minDiagonalElement (double** A, unsigned int matrixSize);
    double*  threeMinDiagonalElements (double** A, unsigned int matrixSize);
    double*  minVectorElements (double* vec, unsigned int vecSize, unsigned int numOfElem2Return);
    
    void     saveMatrix4Mathematica (const char* filename, double** matrix, unsigned int matrixSizeM, unsigned int matrixSizeN);
    void     saveArray4Mathematica (const char* filename, double* array, unsigned int arraySize);
    
    double   vectorDotProduct (double* u, double* v, unsigned int vectorLength);
    void     transposeMatrix (double** A, unsigned int matrixSize);
    
    void     reorderSolution     (double** eigenvectors, double* eigenvalues, unsigned int numberOfEigenvalues);
    
    void     deradializeSolution (double* u, unsigned int vectorLength, double r0, double dr);
    void     squareVector        (double* v, unsigned int vectorLength);
    void     vectorNormalize     (double* v, unsigned int vectorLength, double dx = 1.0);
}
