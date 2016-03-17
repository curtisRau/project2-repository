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
#include <string>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//const double beta = 1.0;             // beta may not be equal to one.  This is a constant in the equation we are trying to solve.

// This is the potential
double V (double rho) {
    return rho * rho;           // + beta / rho;
}

double Vc (double rho, double omega) {
    return omega*omega*rho*rho + 1/rho;           // + beta / rho;
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
    
    unsigned int N      = 200;                           // The Matrix Size.  Nstep = N + 1.  Npoints = N + 2.
    double       rhoMin = 0.0;                          // The starting position.
    double       rhoMax = 5.0;
    double       h      = (rhoMax - rhoMin) / (N + 1);  // The step length
    double       h2     = h*h;                          // Step Length Squared;

    clock_t begin_time;             // Variable for keeping track of computation times.


    // Implement the Jacobi Method
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (true) {
        std::cout << " ----------- Jacobi Method ----------- " << std::endl;

        unsigned int maxRecursion = 50000;          // Maximum number of times "for" loop will run.
        double       tolerance    = 0.1;         // When all off diagonal matrix elements are < this, the matrix is
                                                    // considered diagonalized.

        begin_time = clock();                       // Start the clock.

        // Generate the A matrix which the Jacobi Method will diagonalize
        double* a = function::generateConstantVector(N-1, -1.0/h2);
        double* c = function::generateConstantVector(N-1, -1.0/h2);
        
        double* b = new double [N];
        for (int i = 0; i < N; i++) {
            b[i] = (2.0 / h2) + V(rhoMin + (i+1)*h);
        }

        // Passing vector arguments instead of making calls to the functions
        // that generate the vector arguments allows the following function
        // to be vectorizable.
        double** A = function::genTridiagMatVectArgsExact(N, a, b, c);

        delete [] a;
        delete [] b;
        delete [] c;

        // Variables for Jacobi's Metod "for" loop
        double theta;
        unsigned int x;
        unsigned int y;
        unsigned int* p = &x;
        unsigned int* q = &y;
        double z;
        double* maxValue = &z;
        unsigned int numberOfItterations = 0;
        
        
        double tau;
        double tan;
        double sin;
        double cos;
        
//        function::maxOffDiagnalElement(A, N, maxValue, p, q);
//        for (unsigned int* i = &numberOfItterations; (*i < maxRecursion) && (*maxValue > tolerance); *i += 1) {
//            function::maxOffDiagnalElement(A, N, maxValue, p, q);
//            theta = atan(
//                         (2.0 * A[*p][*q]) / (A[*q][*q] - A[*p][*p])
//                         ) / 2.0;
//            function::jacobiRotation(A, N, *p, *q, theta);
//        }
        
        function::maxOffDiagnalElement(A, N, maxValue, p, q);
        for (unsigned int* i = &numberOfItterations; (*i < maxRecursion) && (*maxValue > tolerance); *i += 1) {
            function::maxOffDiagnalElement(A, N, maxValue, p, q);
            
            tau = (A[*q][*q] - A[*p][*p]) / (2.0 * A[*p][*q]);
            if (tau >= 0.0) {
                tan =  1.0 / ( tau + sqrt(1.0 + tau * tau));
            } else {
                tan = -1.0 / (-tau + sqrt(1.0 + tau * tau));
            }
            
            cos = 1.0 / sqrt(1.0 + tan * tan);
            sin = tan * cos;
            
            function::jacobiRotationSC(A, N, *p, *q, sin, cos);
        }
        
        std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
        std::cout << "Number of itterations performed = " << numberOfItterations << std::endl;
        std::cout << "Lagrest Off Diagonal Element = " << *maxValue << std::endl;
        std::cout << "Smallest eigenvalue = " << function::minDiagonalElement(A, N) << std::endl;

        // Deallocate memory for "A" matrix.
        for (unsigned int i = 0; i<N; i++) {
            delete [] A[i];
        }
        delete [] A;

    }


    // Implement Householder Algorithm
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (true) {
        std::cout << " ----------- Householder Algorithm ----------- " << std::endl;

        begin_time = clock();                                   // Start the clock.

        //make identity matrix for tqli
        double* ones  = function::generateConstantVector(N, 1);
        double* zeros = function::generateConstantVector(N-1, 0);
        double** I    = function::genTridiagMatVectArgsExact(N, zeros, ones, zeros);

        // "ones" and "zeros" no longer needed.
        delete [] ones;
        delete [] zeros;

        double* a = function::generateConstantVector(N-1, -1.0/h2);

        double* b = new double [N];
        for (int i = 0; i < N; i++) {
            b[i] = (2.0 / h2) + V(rhoMin + (i+1)*h);
        }

        tqli(b,a,N,I); // householder method

        delete [] a;    // "a" no longer needed.
        std::cout.precision(17);
        std::cout << "Total computation time [s] = "                 << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
        std::cout << "Smallest eigenvalue (should be 3) = \t\t"      << function::minVectorElements(b, N, 3)[0]          << std::endl;
        std::cout << "Second smallest eigenvalue (should be 7) = \t" << function::minVectorElements(b, N, 3)[1]          << std::endl;
        std::cout << "Third smallest eigenvalue (should be 11) = \t" << function::minVectorElements(b, N, 3)[2]          << std::endl;

        //function::printVector(I[1], N);
        //function::plotVector(I[1], 0, 1000, 20, 10000);

        // Save output
        if (true) {
            function::saveMatrix4Mathematica("/Volumes/userFilesPartition/Users/curtisrau/Documents/School/Physics/PHY480ComputationalPhysics/Project2/project2-repository/dataOut/solutionMatrix.csv", I, N, N);
            function::saveArray4Mathematica("/Volumes/userFilesPartition/Users/curtisrau/Documents/School/Physics/PHY480ComputationalPhysics/Project2/project2-repository/dataOut/eigenvalueArray.csv", b, N);
        }

        // Deallocate Memory
        delete [] b;
        for (unsigned int i = 0; i<N; i++) {
            delete [] I[i];
        }
        delete [] I;

    }

//part c
    if (false) {
        double* omega = new double [4];
        omega[0]=.01;
        omega[1]=0.5;
        omega[2]=1;
        omega[3]=5;

    for (int r = 0; r < 4; r++) {
    std::cout<<"-----------------\nOmega = "<<omega[r]<<std::endl;

        // Implement the Jacobi Method
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if (true) {
            std::cout << " ----------- Jacobi Method ----------- " << std::endl;

            unsigned int maxRecursion = 50000;          // Maximum number of times "for" loop will run.
            double       tolerance    = 0.1;         // When all off diagonal matrix elements are < this, the matrix is
                                                        // considered diagonalized.

            begin_time = clock();                       // Start the clock.

            // Generate the A matrix which the Jacobi Method will diagonalize
            double* a = function::generateConstantVector(N-1, -1.0/h2);
            double* c = function::generateConstantVector(N-1, -1.0/h2);

            double* b = new double [N];
            for (int i = 0; i < N; i++) {
                b[i] = (2.0 / h2) + Vc(rhoMin + (i+1)*h,omega[r]);
            }

            // Passing vector arguments instead of making calls to the functions
            // that generate the vector arguments allows the following function
            // to be vectorizable.
            double** A = function::genTridiagMatVectArgsExact(N, a, b, c);

            delete [] a;
            delete [] b;
            delete [] c;

            // Variables for Jacobi's Metod "for" loop
            double theta;
            unsigned int x;
            unsigned int y;
            unsigned int* p = &x;
            unsigned int* q = &y;
            double z;
            double* maxValue = &z;
            unsigned int numberOfItterations = 0;


            double tau;
            double tan;
            double sin;
            double cos;

    //        function::maxOffDiagnalElement(A, N, maxValue, p, q);
    //        for (unsigned int* i = &numberOfItterations; (*i < maxRecursion) && (*maxValue > tolerance); *i += 1) {
    //            function::maxOffDiagnalElement(A, N, maxValue, p, q);
    //            theta = atan(
    //                         (2.0 * A[*p][*q]) / (A[*q][*q] - A[*p][*p])
    //                         ) / 2.0;
    //            function::jacobiRotation(A, N, *p, *q, theta);
    //        }

            function::maxOffDiagnalElement(A, N, maxValue, p, q);
            for (unsigned int* i = &numberOfItterations; (*i < maxRecursion) && (*maxValue > tolerance); *i += 1) {
                function::maxOffDiagnalElement(A, N, maxValue, p, q);

                tau = (A[*q][*q] - A[*p][*p]) / (2.0 * A[*p][*q]);
                if (tau >= 0.0) {
                    tan =  1.0 / ( tau + sqrt(1.0 + tau * tau));
                } else {
                    tan = -1.0 / (-tau + sqrt(1.0 + tau * tau));
                }

                cos = 1.0 / sqrt(1.0 + tan * tan);
                sin = tan * cos;

                function::jacobiRotationSC(A, N, *p, *q, sin, cos);
            }

            std::cout << "Total computation time [s] = " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
            std::cout << "Number of itterations performed = " << numberOfItterations << std::endl;
            std::cout << "Lagrest Off Diagonal Element = " << *maxValue << std::endl;
            std::cout << "Smallest eigenvalue = " << function::minDiagonalElement(A, N) << std::endl;

            // Deallocate memory for "A" matrix.
            for (unsigned int i = 0; i<N; i++) {
                delete [] A[i];
            }
            delete [] A;

        }
    }
    
    }

 //part d
    if (true){
        double* omega = new double [4];
        omega[0]=.01;
        omega[1]=0.5;
        omega[2]=1;
        omega[3]=5;
        int count = 1;
    for (int r = 0; r < 4; r++) {
        count++;
        std::cout<<"-----------------\nOmega = "<<omega[r]<<std::endl;

            std::cout << " ----------- Householder Algorithm ----------- " << std::endl;

            begin_time = clock();                                   // Start the clock.

            //make identity matrix for tqli
            double* ones  = function::generateConstantVector(N, 1);
            double* zeros = function::generateConstantVector(N-1, 0);
            double** I    = function::genTridiagMatVectArgsExact(N, zeros, ones, zeros);

            // "ones" and "zeros" no longer needed.
            delete [] ones;
            delete [] zeros;

            double* a = function::generateConstantVector(N-1, -1.0/h2);

            double* b = new double [N];
            for (int i = 0; i < N; i++) {
                b[i] = (2.0 / h2) + Vc(rhoMin + (i+1)*h,omega[r]);
            }

            tqli(b,a,N,I); // householder method

            delete [] a;    // "a" no longer needed.
            std::cout.precision(17);
            std::cout << "Total computation time [s] = "                 << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
            std::cout << "Smallest eigenvalue (should be 3) = \t\t"      << function::minVectorElements(b, N, 3)[0]          << std::endl;
            std::cout << "Second smallest eigenvalue (should be 7) = \t" << function::minVectorElements(b, N, 3)[1]          << std::endl;
            std::cout << "Third smallest eigenvalue (should be 11) = \t" << function::minVectorElements(b, N, 3)[2]          << std::endl;

            //function::printVector(I[1], N);
            //function::plotVector(I[1], 0, 1000, 20, 10000);

            // Save output
            if (true) {
                std::string r_str = std::to_string(count);
                std::string matrixname = "matrix" + r_str + ".txt";
                std::string arrayname = "array" + r_str + ".txt";
                function::saveMatrix4Mathematica(matrixname, I, N, N);
                function::saveArray4Mathematica(arrayname, b, N);
            }

            // Deallocate Memory
            delete [] b;
            for (unsigned int i = 0; i<N; i++) {
                delete [] I[i];
            }
            delete [] I;

        }


    }//end of part d
    return 0;
}
