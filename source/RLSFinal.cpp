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
#include "sig/sigpack/sigpack.h"
#include <cstdlib>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace arma;
using namespace sp;

typedef complex<double> cdouble;

class GenSignalReturn {
	/*Matrices returned by the genSignal function.
	These matrices are named exactly as found in the
	Matlab program.*/
public:
  mat received;
  mat desired;
  mat noise;

  GenSignalReturn(mat received, mat desired, mat noise) {
	  /*Create return object from matrices.*/
    this->received = received;
    this->desired = desired;
    this->noise = noise;
  }
};

class RLSReturn {
	/*Matrices returned by the RLS() function. Matrices
	are named as they appear in the Matlab project.*/
public:
  mat error;
  mat weights;

  RLSReturn(mat error, mat weights) {
    this->error = error;
    this->weights = weights;
  }
};

double genRandAngle(){
	/*Generate a random angle between 0 and pi.*/
    double randNum = ((double) rand() / (RAND_MAX));
    return randNum * M_PI;
}

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

  cdouble i(0, 1);

  mat t = linspace(0, 1, num_points);

  mat desired = sin(2*M_PI*frequency*t);
 
  double elevation = genRandAngle();
  double azimuth   = genRandAngle();
  
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

  // Generate noise
  FIR_filt<double, double, double> G;
  G.set_coeffs(filt);
  mat noise = sqrt(n_var) * randn(num_points, 1);
  mat addnoise = G.filter(noise);
 
  desired = (desired) / sqrt(var(vectorise(desired))/(pow(10,(SNR/10))*var(vectorise(noise))));

  // Add noise to signal - beware of the real function!

  cx_mat result = desired*g_arma*h_arma.t();
  
  mat received = real(result) + addnoise;
  
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

  arma_rng::set_seed_random();
  
  // parameters!
  int num_points = 1000;
  double frequency = 12e9;
  int filter_order = 16;
  mat filter = randu<mat>(filter_order, 1);
  int n_var = 1;
  int SNR = -20;
  
  // Generate some sort of signal - wow!
  GenSignalReturn gsr = genSignal(num_points, frequency, filter, n_var, SNR);
  mat received = gsr.received;
  mat desired = gsr.desired;
  mat noise = gsr.noise;

  mat P = eye(sysorder, sysorder) * 1.0/delta;
  mat weights = zeros<mat>(sysorder,1);
  vec received_vec = vectorise(received);
  mat error = desired*0;
  double lambda1 = pow(lambda, -1);
  int len = received_vec.n_elem;

  // Big calculation here
  int number_of_loop_executions = 0;
  for(int n = sysorder-1; n < len; n++) {
    ++number_of_loop_executions;

    mat input = zeros<mat>(sysorder, 1);
    for(int i = n; i > n-sysorder; i--) {
      //cout << "i: " << i << endl;
      input(n - i) = noise(i);
    }

    mat K = (lambda1*P*input)/as_scalar(1+lambda1*input.t()*P*input); //ahh!
	  
    mat output = weights.t() * input; // there was another '* here!

    error(n) = received_vec(n) - as_scalar(output);

    weights = weights + K * error(n); //error(n) instead

    P=(lambda1*P)-(lambda1*K*input.t()*P); // another '* here
  }
  
  mat e = zeros<mat>(1000-16, 1);

  for(int index = 0; index < 1000-16; index++) {
    double tmp = error(15+index);
    e(index) = tmp;
  }

  mat q = linspace(1, 985, 985);

  mat e2 = abs(pow(e,2));

  // Calculate SNR improvement
  double SNRi = 10*log10(var(vectorise(received))/var(vectorise(error)));
  cout << SNRi << "dB SNR Improvement" << endl;

  return RLSReturn(error, weights);
}


int main(int argc, char* argv[]) {
  
  srand (time(NULL));
  // Driver program for RLS() function.
  RLS(0.98, 100.0, 16.0);
  return 0;
}

