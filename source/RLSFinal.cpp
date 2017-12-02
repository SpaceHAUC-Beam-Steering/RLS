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
//#include "sigpack-1.2.2/sigpack/sigpack.h"
#include "sig/sigpack/sigpack.h"
#include <cstdlib>
#include <time.h>

#ifndef M_PI
#define M_PI 3.1416 //(3.14159265358979323846)
#endif

using namespace std;
using namespace arma;
using namespace sp;

// Random generator written by Fezvez
// https://stackoverflow.com/questions/5829499/
// how-to-let-boostrandom-and-matlab-produce-the-same-random-numbers

// uint32_t xor128(void) {
//   static uint32_t x = 123456789;
//   static uint32_t y = 362436069;
//   static uint32_t z = 521288629;
//   static uint32_t w = 88675123;
//   uint32_t t;

//   t = x ^ (x << 11);
//   x = y; y = z; z = w;
//   return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
// }

// double urand() {
//   double max_uint = (double)numeric_limits<uint32_t>::max();
//   uint32_t random_integer = xor128();
//   return random_integer / max_uint;
// }

// double nrand() {
//   double random_number = urand();
//   //cout << "random_number: " << random_number << endl;
//   double log_rand = log(random_number);
//   //cout << "log_rand: " << log_rand << endl;
//   double square_root = pow(-2 * log_rand, 0.5);
//   //cout << "square_root: " << square_root << endl;
//   return square_root * sin(2 * M_PI * urand());
// }

// mat urand_matrix(int num_rows, int num_cols) {
//   /*Generate a rows by cols matrix filled with
//    random uniformly generated numbers on interval [0,1].*/
//   mat matrix(num_rows, num_cols, fill::zeros);
//   for(int row_index=0; row_index<num_rows; row_index++) {
//     for(int col_index=0; col_index<num_cols; col_index++) {
//       matrix(row_index, col_index) = urand();
//     }
//   }
//   return matrix;
// }

// mat nrand_matrix(int rows, int cols) {
//   mat matrix(rows, cols, fill::zeros);
//   for(int row_index=0; row_index<rows; row_index++) {
//     for(int col_index=0; col_index<cols; col_index++) {
//       matrix(row_index, col_index) = nrand();
//     }
//   }
//   return matrix;
// }

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

double genRandAngle(){
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
  //cout << "genSignal function call" << endl;

  cdouble i(0, 1);

  mat t = linspace(0, 1, num_points);

  //cout << "Mean t: " << mean(t) << endl;

  mat desired = sin(2*M_PI*frequency*t);

  //cout << "Frequency: " << frequency << endl;

  //cout << "Intermediate: " << mean(2*M_PI*frequency*t);

  //cout << "Original desired: " << mean(desired) << endl;
  //cout << "Original desired var: " << var(vectorise(desired)) << endl;

  
  // Testing purposes
 
  //double elevation = M_PI/3;
  //double azimuth = M_PI/3;
 
  double elevation = genRandAngle();
  double azimuth   = genRandAngle();

  //cout << "elevation: " << elevation << endl;
  //cout << "azimuth:   " << azimuth << endl;
  //cout << "complex math" << endl;
  
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

  //cout << "noise calculations" << endl;

  // Generate noise
  //mat temp = nrand_matrix(10, 1);
  //cout << temp << endl;
  //cout << "Temp mean: " << mean(temp) << endl;
  FIR_filt<double, double, double> G;
  G.set_coeffs(filt);
  //mat noise = sqrt(n_var) * nrand_matrix(num_points, 1); //randn(num_points,1);
  mat noise = sqrt(n_var) * randn(num_points, 1);
  //cout << "Mean noise: " << mean(noise) << endl;
  mat addnoise = G.filter(noise);
  //cout << "Mean addnoise: " << mean(addnoise) << endl;
  
  //cout << "SNR calculations" << endl;
  desired = (desired) / sqrt(var(vectorise(desired))/(pow(10,(SNR/10))*var(vectorise(noise))));
  // Display calculated SNR
  //cout << "Calculated SNR = " << 10*log10(var(vectorise(desired))/var(vectorise(noise))) << endl;

  // Add noise to signal - beware of the real function!

  cx_mat result = desired*g_arma*h_arma.t();

  //cout << "Variance addnoise: " << var(vectorise(addnoise)) << endl;
  
  mat received = real(result) + addnoise;

  //cout << "genSignal return" << endl;

  //cout << "Received mean value: " << mean(received) << endl;
  //cout << "Desired mean value: " << mean(desired) << endl;
  //cout << "Desired var: " << var(vectorise(desired)) << endl;
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
  //cout << "RLS function call" << endl;

  arma_rng::set_seed_random();
  
  int num_points = 1000;
  double frequency = 12e9;
  int filter_order = 16;
  //mat filter = urand_matrix(filter_order, 1); //
  mat filter = randu<mat>(filter_order, 1); //should be rand(filtord, 1);
  int n_var = 1;
  int SNR = -20;

  //cout << "Mean filter: " << mean(filter) << endl;
  //cout << "Var filter: " << var(vectorise(filter)) << endl;
  
  GenSignalReturn gsr = genSignal(num_points, frequency, filter, n_var, SNR);
  mat received = gsr.received;
  mat desired = gsr.desired;
  mat noise = gsr.noise;

  //cout << "matrices initialize" << endl;
  mat P = eye(sysorder, sysorder) * 1.0/delta;
  //cout << "weights" << endl;
  mat weights = zeros<mat>(sysorder,1);
  //cout << "received_vec" << endl;
  vec received_vec = vectorise(received);
  //cout << "error" << endl;
  mat error = desired*0;
  //cout << "lambda1 created from lambda=" << lambda << endl;
  double lambda1 = pow(lambda, -1);
  //cout << "len" << endl;
  int len = received_vec.n_elem;
  //cout << "sysorder: " << sysorder << endl;
  //cout << "len: " << len << endl;
  //cout << "First and last noise: " << noise(0) << "," << noise(size(noise)[0]-1) << endl;
  int number_of_loop_executions = 0;
  for(int n = sysorder-1; n < len; n++) {
    //cout << "Iteration: " << n << endl;
    ++number_of_loop_executions;
    //cout << "process input" << endl;
    // Must do indexing manually: input = noise(n:-1:n-sysorder+1);
    //cout << "Noise size: " << size(noise) << endl;
    //cout << "Indices: " << size(noise)[0]-1 << "," << 1 << endl;

 
    mat input = zeros<mat>(sysorder, 1);
    for(int i = n; i > n-sysorder; i--) {
      //cout << "i: " << i << endl;
      input(n - i) = noise(i);
    }
    //cout << "Mean/var noise" << mean(noise) << ", " << var(noise) << endl;
    //cout << "Mean/var input: " << mean(input) << ", " << var(input) << endl;
    //cout << "First and last input: " << input(0) << "," << input(size(input)[0]-1) << endl;

    //cout << "Calculate K" << endl;
    //cout << "Value of denominator: " << as_scalar(1+lambda1*input.t()*P*input) << endl;
    mat K = (lambda1*P*input)/as_scalar(1+lambda1*input.t()*P*input); //ahh!
    //cout << "Mean/var K: " << mean(K) << ", " << var(vectorise(K)) << endl;
    mat output = weights.t() * input; // there was another '* here!
    //cout << "Mean/var K: " << mean(K) << ", " << var(vectorise(K)) << endl;
    error(n) = received_vec(n) - as_scalar(output);
    //cout << "Mean/var error: " << mean(error) << ", " << var(vectorise(error)) << endl;
    weights = weights + K * error(n); //error(n) instead
    //cout << "Mean/var weights: " << mean(weights) << ", " << var(vectorise(weights)) << endl;
    P=(lambda1*P)-(lambda1*K*input.t()*P); // another '* here
    //cout << "Mean/var P: " << mean(vectorise(P)) << ", " << var(vectorise(P)) << endl;
    //if(number_of_loop_executions == 2) {
    //  break;
    //}
  }
  //cout << "Number of loop executions: " << number_of_loop_executions << endl;

  //cout << "Weights mean: " << mean(weights) << endl;
  //cout << "Weights var: " << var(vectorise(weights)) << endl;

  //cout << "calculate error statistics" << endl;
  
  mat e = zeros<mat>(1000-16, 1); //error(16:1000);
  //cout << e.size() << endl;
  //cout << error.size() << endl;
  //cout << error(0) << endl;
  for(int index = 0; index < 1000-16; index++) {
    double tmp = error(15+index);
    e(index) = tmp;
  }
  //cout << "linspace" << endl;
  mat q = linspace(1, 985, 985);

  mat e2 = abs(pow(e,2));

  //cout << "plot" << endl;

  // Calculate SNR improvement
  double SNRi = 10*log10(var(vectorise(received))/var(vectorise(error)));
  cout << SNRi << "dB SNR Improvement" << endl;
  /*
  cout << "uint32 max: " << numeric_limits<uint32_t>::max() << endl;
  cout << "RLS return" << endl;

  mat test_matrix = nrand_matrix(4, 3);
  cout << test_matrix << endl;

  for(int index = 0; index < 10; index++) {
    cout << nrand() << endl;
  }

  mat generated_numbers = nrand_matrix(1000, 1);
  cout << "Variance: " << var(vectorise(generated_numbers)) << endl;
  */
  
  return RLSReturn(error, weights);
}


int main(int argc, char* argv[]) {
  
  srand (time(NULL));
  // Driver program.. need parameters to test this function with.
  RLS(0.98, 100.0, 16.0);
  return 0;
}

