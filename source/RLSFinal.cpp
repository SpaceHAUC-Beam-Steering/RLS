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

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

using namespace std;
using namespace arma;

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

  cdouble i(0, 1);

  mat t = linspace(0, 1, num_points);

  mat desired = sin(2*M_PI*frequency*t);

  double elevation = M_PI/3;
  double azimuth = M_PI/3;
  
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
  mat noise = sqrt(n_var)*randn(num_points,1);
  //double addnoise = 1; //filter(filt, 1, noise);
  //freqz(filt, 1, num_points);

  desired = (desired) / sqrt(var(desired)/(pow(10,(SNR/10))*var(noise)));
  // Display calculated SNR
  cout << "Calculated SNR = " << 5.0 << endl; // num2str(10*log10(var(desired)/var(noise)))])

  // Add noise to signal - beware of the real function!
  
  mat received = real(desired*g_arma*h_arma.t()); // + addnoise;
  
  return GenSignalReturn(received, desired, noise);
}


RLSReturn RLS(int lambda, int delta, int sysorder) {
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
  
  int num_points = 1000;
  double frequency = 12e9;
  int filter_order = 16;
  mat filter = randu<mat>(filter_order, 1);//should be rand(filtord, 1);
  int n_var = 1;
  int SNR = -20;
  
  GenSignalReturn gsr = genSignal(num_points, frequency, filter, n_var, SNR);
  mat received = gsr.received;
  mat desired = gsr.desired;
  mat noise = gsr.noise;

  mat P = eye(sysorder, sysorder) * 1.0/delta;
  //double P = 0; //eye(sysorder)* 1/delta;
  mat weights = zeros<mat>(sysorder,1);
  vec received_vec = vectorise(received);
  mat error = desired*0;
  double lambda1 = 1/lambda;
  int len = received_vec.n_elem;

  for(int n = sysorder; n < len; n++) {

    // Must do indexing manually: input = noise(n:-1:n-sysorder+1);
    int size = sysorder + 1;
    mat input = zeros<mat>(size, 1);
    for(int i = 0; i < n - sysorder; i++) {
      input(i, 1) = noise(i, 1);
    }
    

    mat K = (lambda1*P*input)/(1+lambda1*input.t()*P*input); //ahh! there was a '* here somewhere...
    mat output = weights.t() * input; // there was another '* here!
    error(n) = received_vec(n) - output(0);
    weights = weights + K * error(n); //error(n) instead
    P=(lambda1*P)-(lambda1*K*input.t()*P); // another '* here
  }
  
  double e = 0; //error(16:1000);
  double q = -1; //linspace(1, 985, 985);

  double e2 = 0; //abs(e.^2);

  // Plot here!

  double t = 0; //linspace(0, length(received)/num_points, length(received));

  //figure;
  //plot(t, desired, t, error);
  //title('Comparison of Filtered Signal to Original Signal');
  //figure;
  //plot(q, e2);

  // Calculate SNR improvement
  double SNRi = 10; //10*log10(var(received)/var(error));
  cout << SNRi << "dB SNR Improvement" << endl;

  return RLSReturn(error, weights);
}

int main(int argc, char* argv[]) {
  // Driver program.. need parameters to test this functio with.
  RLS(1,.1,1);
  return 0;
}

