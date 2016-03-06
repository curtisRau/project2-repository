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
#include "lib.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//const double beta = 1.0;             // beta may not be equal to one.  This is a constant in the equation we are trying to solve.

// This is the potential
double V (double rho) {
    return rho * rho;           // + beta / rho;
}

//// Convert the eigenvalue to energy
//double eigenvalue2Energy (double eigenvalue) {
//    double hbar = ;
//    double omega = ;
//    return hbar * omega / 2.0;
//}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char * argv[]) {
    
    // Why does the smallest eigenvalue get more negitive as N gets bigger?
    unsigned int N      = 1000;                       // The number of ??? should this be N-2 or something?
    double       rhoMin = 0.0;                      // The starting position.
    double       rhoMax = 5000;
    double       h      = (rhoMax - rhoMin) / N;    // The step length
    double       h2     = h*h;                      // Step Length Squared;
    
    // Generate the A matrix which the Jacobi Method will diagonalize
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    double* a = function::generateConstantVector(N-1, -1.0/h2);
    double* c = function::generateConstantVector(N-1, -1.0/h2);
    double* b = new double[N];

    b[0]=0;
    for (int i = 1; i < N; i++) {
        b[i] = 2.0/h2 + V(rhoMin + i*h);
    }


    
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
    std::cout << " ----------- Jacobi Method ----------- " << std::endl;
    
    unsigned int x;
    unsigned int y;
    unsigned int* p = &x;
    unsigned int* q = &y;
    double theta;
    unsigned int maxRecursion = 5000;        // Maximum number of times for loop will run.
    double z;
    double* maxValue = &z;
    double tolerance = 0.0001;
    unsigned int numberOfItterations = 0;
    
    const clock_t begin_time = clock();
    
    function::maxOffDiagnalElement(A, N, maxValue, p, q);
    for (unsigned int* i = &numberOfItterations; (*i < maxRecursion) && (*maxValue > tolerance); *i += 1) {
        function::maxOffDiagnalElement(A, N, maxValue, p, q);
        theta = atan(
                     (2*A[*p][*q]) / (A[*q][*q] - A[*p][*p])
                     ) / 2.0;
        function::jacobiRotation(A, N, *p, *q, theta);
    }
    
    std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
    std::cout << "Number of itterations performed = " << numberOfItterations << std::endl;
    std::cout << "Lagrest Off Diagonal Element = " << *maxValue << std::endl;
    std::cout << "Smallest eigenvalue = " << function::minDiagonalElement(A, N) << std::endl;

    
    // Implement Householder Algorithm
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << " ----------- Householder Algorithm ----------- " << std::endl;

    //make identity matrix for tqli
    double* ones = function::generateConstantVector(N, 1);
    double* zeros = function::generateConstantVector(N-1, 0);
    double** I  = function::genTridiagMatVectArgsExact(N, zeros, ones, zeros);

    tqli(b,a,N,I); //householder method

    std::cout << "Smallest eigenvalue = " << function::minDiagonalElement(A, N) << std::endl;

    // Still need to delete the A matrix and I matrix
    delete [] a;
    delete [] b;
    delete [] ones;
    delete [] zeros;
    return 0;
}
