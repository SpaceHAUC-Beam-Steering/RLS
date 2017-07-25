/*
    UMass Lowell SPACEHAUC Beam Steering

    RLS.cpp     

    GNU GPLv3 License

*/


#include <cstdio>
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Matrix>
#include <eigen3/Eigen/MatrixBase>

using namespace std;

int[] RLS(int lambda, int delta, int sysorder){

    /*i
    Input:
    @param lambda           forgetting factor, dim 1x1   
    @param delta            initial value, P(0)=delta^(-1) * (Identity Matrix), dim 1x1
    @param sysorder         filter order, dim 1x1
        

    Output:
    @error                  a priority estimation error, dim nx1
    @weights                final filter coefficients, dim sysorderx1
    */


    // Data Parameters

    int NUM_PTS     = 1000;                 // Number of points to generate
    long double FREQ     = 12e9;                 // Frequency of fundamental tone
    int FILT_ORDER  = 16;                   // Filter order
    int FILT        = rand(FILT_ORDER,1);   // Filter coefficients, FILT_ORDERx1 matrix
    int N_VAR       = 1;                    // White Noise variance
    int SNR         = -20;                  // Signal to noise ratio of tone

    int[] gen_values = gen_Signal(NUM_PTS, FREQ, FILT, N_VAR, SNR);

    int received = gen_values[0];
    int desired  = gen_values[1];
    int noise    = gen_values[2];

            
}


