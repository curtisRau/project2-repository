//
//  main.cpp
//  project2
//
//  Created by Curtis Rau on 3/4/16.
//  Copyright © 2016 Curtis Rau. All rights reserved.
//

#include <iostream>
#include "functions.hpp"
#include <math.h>           // for atan
#include "time.h"
#include "lib.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//const double beta = 1.0;             // beta may not be equal to one.  This is a constant in the equation we are trying to solve.

// This is the potential
double V (double rho) {
    return rho * rho;           // + beta / rho;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char * argv[]) {
    
    unsigned int N    = 10;
    //double       rho0 = 0.00000000000000001;        // The starting position, probably 0.0, but 1/0 encountered.
    double       h    = 1;                  // The step length
    double       h2   = h*h;                // Step Length Squared;
    
    // Generate the A matrix which the Jacobi Method will diagonalize
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    double* a = function::generateConstantVector(N-1, -1/h2);
    double* c = function::generateConstantVector(N-1, -1/h2);
    double* b = new double[N];
    b[0]=2/h2;
    for (int i = 1; i < N; i++) {
        b[i] = 2/h2 + V(i*h);
    }
    function::printVector(b,N);
    // Passing vector arguments instead of making calls to the functions
    // that generate the vector arguments allows the following function
    // to be vectorizable.
    double** A = function::genTridiagMatVectArgsExact(N, a, b, c);
    
    //delete [] a;
    //delete [] b;
    delete [] c;
    
    // Implement the Jacobi Method
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    unsigned int x;
    unsigned int y;
    unsigned int* p = &x;
    unsigned int* q = &y;
    double theta = 100;
    unsigned int maxRecursion = 100;        // Maximum number of times for loop will run.
    double minTheta = .000000000000000000000000000000000000000001;

    for (unsigned int i = 0; ( theta*theta > minTheta); i++) {
        function::indiciesOfMaxOffDiagnalElement(A, N, N, p, q);
        theta = atan(
                     (2*A[*p][*q]) / (A[*q][*q] - A[*p][*p])
                     ) / 2.0;
        function::jacobiRotation(A, N, *p, *q, theta);
        std::cout<<theta<<endl;
        //function::printMatrix(A, N, N);
    }

    function::printMatrix(A, N, N);
    //make identity matrix for tqli
    double* ones = function::generateConstantVector(N, 1);
    double* zeros = function::generateConstantVector(N-1, 0);
    double** z  = function::genTridiagMatVectArgsExact(N, zeros, ones, zeros);

    tqli(b,a,N,z); //householder method

    function::printDiagonals(A,N);
    function::printVector(b,N);



    delete [] a;
    delete [] b;
    return 0;
}
