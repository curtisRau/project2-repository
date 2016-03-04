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
#include <iostream>     // For std::cout

#endif

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
    
}
