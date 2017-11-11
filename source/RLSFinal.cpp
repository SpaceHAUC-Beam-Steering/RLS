/*
    UMass Lowell SPACEHAUC Beam Steering
    RLSFinal.cpp     
    GNU GPLv3 License
*/

#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>
#include <complex>
#include "sigpack-1.2.3/sigpack/sigpack.h"
#include <cstdlib>

#ifndef M_PI
#define M_PI 3.1416 //(3.14159265358979323846)
#endif

using namespace std;
using namespace arma;
using namespace sp;

typedef complex<double> cdouble;

class GenSignalReturn {
public:
  mat received;
  mat desired;
  mat noise;

  GenSignalReturn(mat received, mat desired, mat noise) {
    this->received = received;
    this->desired = desired;
    this->noise = noise;
  }
};

class RLSReturn {
public:
  mat error;
  mat weights;

  RLSReturn(mat error, mat weights) {
    this->error = error;
    this->weights = weights;
  }
};

GenSignalReturn genSignal(double num_points, double frequency, mat filt,
			  double n_var, double SNR) {
    /*
      Produces object containing three members:
      received      input signal, dim 1x1
      desired       desired signal, dim 1x1
      noise         noise 
      
      Input:
      @param num_points       number of points to generate 
      @param frequency         frequency of fundamental tone
      @param filt         filter coefficients
      @param n_var         white noise variance
      @param SNR          signal to noise ratio of tone
    */
  cout << "genSignal function call" << endl;

  cdouble i(0, 1);

  mat t = linspace(0, 1, num_points);

  cout << "Mean t: " << mean(t) << endl;

  mat desired = sin(2*M_PI*frequency*t);

  cout << "Frequency: " << frequency << endl;

  cout << "Intermediate: " << mean(2*M_PI*frequency*t);

  cout << "Original desired: " << mean(desired) << endl;
  cout << "Original desired var: " << var(vectorise(desired)) << endl;

  
  // Testing purposes
 
  // double elevation = M_PI/3;
  // double azimuth = M_PI/3;
 
  double elevation = genRandAngle();
  double azimuth   = genRandAngle();

  cout << "complex math" << endl;
  
  std::vector<complex<double> > g;
  for(double m = 0; m < 4; m++) {
    cdouble u = exp(i*M_PI*(m-1));
    cdouble a = u*sin(elevation)*cos(azimuth);
    g.push_back(a);
  }
  cx_rowvec g_arma = cx_rowvec(g);

  std::vector<complex<double> > h;
  for(double p = 0; p < 4; p++) {
    cdouble v = exp(i*M_PI*(p-1));
    cdouble b = v * sin(elevation)*sin(azimuth);
    h.push_back(b);
  }
  cx_rowvec h_arma = cx_rowvec(h);

  cout << "noise calculations" << endl;

  // Generate noise
  FIR_filt<double, double, double> G;
  G.set_coeffs(filt);
  mat noise = sqrt(n_var)*randn(num_points,1);
  cout << "Mean noise: " << vectorise(mean(noise)) << endl;
  mat addnoise = G.filter(noise);

  cout << "SNR calculations" << endl;
  desired = (desired) / sqrt(var(vectorise(desired))/(pow(10,(SNR/10))*var(vectorise(noise))));
  // Display calculated SNR
  cout << "Calculated SNR = " << 10*log10(var(vectorise(desired))/var(vectorise(noise))) << endl;

  // Add noise to signal - beware of the real function!

  cx_mat result = desired*g_arma*h_arma.t();

  cout << "Mean addnoise: " << mean(vectorise(addnoise)) << endl;
  cout << "Variance addnoise: " << var(vectorise(addnoise)) << endl;
  
  mat received = real(result) + addnoise;

  cout << "genSignal return" << endl;

  cout << "Received mean value: " << mean(received) << endl;
  cout << "Desired mean value: " << mean(desired) << endl;
  cout << "Original desired var: " << var(vectorise(desired)) << endl;
  //cout << noise << endl;
  
  return GenSignalReturn(received, desired, noise);
}


RLSReturn RLS(double lambda, double delta, double sysorder) {
  /*
  Input:
  @param lambda           forgetting factor, dim 1x1   
  @param delta            initial value, P(0)=delta^(-1) * (Identity Matrix), dim 1x1
  @param sysorder         filter order, dim 1x1
     
  Output:
  Object containing members .error and .weights
  @error                  a priority estimation error, dim nx1
  @weights                final filter coefficients, dim sysorderx1
  */
  cout << "RLS function call" << endl;

  arma_rng::set_seed_random();
  
  int num_points = 1000;
  double frequency = 9e9;
  int filter_order = 16;
  mat filter = randu<mat>(filter_order, 1);//should be rand(filtord, 1);
  int n_var = 1;
  int SNR = -20;
  
  GenSignalReturn gsr = genSignal(num_points, frequency, filter, n_var, SNR);
  mat received = gsr.received;
  mat desired = gsr.desired;
  mat noise = gsr.noise;

  cout << "matrices initialize" << endl;
  mat P = eye(sysorder, sysorder) * 1.0/delta;
  cout << "weights" << endl;
  mat weights = zeros<mat>(sysorder,1);
  cout << "received_vec" << endl;
  vec received_vec = vectorise(received);
  cout << "error" << endl;
  mat error = desired*0;
  cout << "lambda1 created from lambda=" << lambda << endl;
  double lambda1 = pow(lambda, -1);
  cout << "len" << endl;
  int len = received_vec.n_elem;
  cout << "sysorder: " << sysorder << endl;
  cout << "len: " << len << endl;
  for(int n = sysorder-1; n < len; n++) {
    
    cout << "process input" << endl;
    cout << "Mean/var noise" << mean(noise) << ", " << var(noise) << endl;
    // Must do indexing manually: input = noise(n:-1:n-sysorder+1);
    mat input = zeros<mat>(sysorder, 1);
    for(int i = n - sysorder+1; i < n+1; i++) {
      input(i) = noise(i);
    }
    cout << "Mean/var input: " << mean(input) << ", " << var(input) << endl;
    mat tmp = input;
    cout << "flipping" << endl;
    input = flipud(input);
    cout << "testing for reversal of input" << endl;
    if(tmp(0) != input(input.n_elem-1)) {
      cout << "Input may not be reversed!" << endl;
      cout << "Final element: " << input.n_elem-1 << endl;
      cout << "Values: " << tmp(0) << ", " << input(input.n_elem-1) << endl;
    }

    cout << "Calculate K" << endl;
    cout << "Value of denominator: " << as_scalar(1+lambda1*input.t()*P*input) << endl;
    cout << "lambda1*input.t()*P: " << lambda1 << endl;
    mat K = (lambda1*P*input)/as_scalar(1+lambda1*input.t()*P*input); //ahh!
    cout << "Mean/var K: " << mean(K) << ", " << var(K) << endl;
    cout << "calculate output" << endl;
    
    mat output = weights.t() * input; // there was another '* here!
    error(n) = received_vec(n) - as_scalar(output);
    weights = weights + K * error(n); //error(n) instead
    P=(lambda1*P)-(lambda1*K*input.t()*P); // another '* here
  }

  cout << "Weights mean: " << mean(weights) << endl;
  cout << "Weights var: " << var(vectorise(weights)) << endl;

  cout << "calculate error statistics" << endl;
  
  mat e = zeros<mat>(1000-16, 1); //error(16:1000);
  cout << e.size() << endl;
  cout << error.size() << endl;
  cout << error(0) << endl;
  for(int index = 0; index < 1000-16; index++) {
    double tmp = error(15+index);
    e(index) = tmp;
  }
  cout << "linspace" << endl;
  mat q = linspace(1, 985, 985);

  mat e2 = abs(pow(e,2));

  cout << "plot" << endl;

  // Calculate SNR improvement
  double SNRi = 10*log10(var(vectorise(received))/var(vectorise(error)));
  cout << SNRi << "dB SNR Improvement" << endl;

  cout << "RLS return" << endl;

  return RLSReturn(error, weights);
}

double genRandAngle(){
    double randNum = ((double) rand() / (RAND_MAX)); 
    return randNum * M_PI;
}

int main(int argc, char* argv[]) {
  // Driver program.. need parameters to test this function with.
  RLS(0.98, 100.0, 16.0);
  return 0;
}

