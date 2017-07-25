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

    /*
    Input:
    @param lambda           forgetting factor, dim 1x1   
    @param delta            initial value, P(0)=delta^(-1) * (Identity Matrix), dim 1x1
    @param sysorder         filter order, dim 1x1
        

    Output:
    @error                  a priority estimation error, dim nx1
    @weights                final filter coefficients, dim sysorderx1
    */


    // Data Parameters
    
    int E_W[1] = 0;
    
    int NUM_PTS          = 1000;                 // Number of points to generate
    long double FREQ     = 12e9;                 // Frequency of fundamental tone
    int FILT_ORDER       = 16;                   // Filter order
    int FILT             = rand(FILT_ORDER,1);   // Filter coefficients, FILT_ORDERx1 matrix
    int N_VAR            = 1;                    // White Noise variance
    int SNR              = -20;                  // Signal to noise ratio of tone

    int[] gen_values = gen_Signal(NUM_PTS, FREQ, FILT, N_VAR, SNR);

    int received = gen_values[0];
    int desired  = gen_values[1];
    int noise    = gen_values[2];
                    
}

int[] gen_Signal(int numPts, long double freq, int filt, int nVar, int SNR){

    /*

        Produces a 1x3 array with indexes:
        [0] - received      input signal, dim 1x1
        [1] - desired       desired signal, dim 1x1
        [2] - noise         noise 
        
        Input:
        @param numPTs       number of points to generate 
        @param freq         frequency of fundamental tone
        @param filt         filter coefficients
        @param nVar         white noise variance
        @param SNR          signal to noise ratio of tone

    */
    

    // Generate time values
    time_values = linspace(0,1,numPts);

    // Generate tone


    // Wavevector from angle of arrival 
    float elevation = PI / 3;
    float azimuth   = PI / 3;

    // Generate noise


    // Adjust SNR of tone

    
    // Add noise to signal

    return
}

