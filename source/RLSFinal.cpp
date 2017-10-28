/*
    UMass Lowell SPACEHAUC Beam Steering
    RLS.cpp     
    GNU GPLv3 License
*/

#include <iostream>
#include <armadillo>
#include <cmath>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

using namespace std;
using namespace arma;

class GenSignalReturn {
public:
  double received;
  double desired;
  double noise;

  GenSignalReturn(double received, double desired, double noise) {
    this->received = received;
    this->desired = desired;
    this->noise = noise;
  }
};

class RLSReturn {
public:
  double error;
  double weights;

  RLSReturn(double error, double weights) {
    this->error = error;
    this->weights = weights;
  }
};

GenSignalReturn genSignal(double num_points, double frequency, double filter,
			  double n_var, double SNR) {
    /*
      Produces object containing three members:
      received      input signal, dim 1x1
      desired       desired signal, dim 1x1
      noise         noise 
      
      Input:
      @param num_points       number of points to generate 
      @param frequency         frequency of fundamental tone
      @param filter         filter coefficients
      @param n_var         white noise variance
      @param SNR          signal to noise ratio of tone
    */

  double t = 0;//linspace(0, 1, num_points);

  double desired = 0; //sin(2*pi*frequency*t);

  double elevation = M_PI/3;
  double azimuth = M_PI/3;
  
  double g = 0; // [];
  for(int m = 1; m < 4; m++) {
    double u = 0;//exp(1i*pi*(m-1);
    double a = u*sin(elevation)*cos(azimuth);
    double g = 0; //[g a];
  }

  double h = 0; //[];
  for(int p = 0; p < 4; p++) {
    double v = 0;//exp(1i*pi*(p-1));
    double b = v * sin(elevation)*sin(azimuth);
    h = 0; //[h b];
  }

  // Generate noise
  double noise = sqrt(n_var)*1;//sqrt(nVar)*randn(num_points,1);
  double addnoise = 1; //filter(filter, 1, noise);
  //freqz(filt, 1, num_points);

  desired = (desired); // (desired) / sqrt(var(desired)/(10^(SNR/10)*var(noise)));
  // Display calculated SNR
  cout << "Calculated SNR = " << 5.0 << endl; // num2str(10*log10(var(desired)/var(noise)))])

  // Add noise to signal
  double received = 0; //(desired*g*h.') + addnoise;
  
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
  int filter = -1; //should be rand(filtord, 1);
  int n_var = 1;
  int SNR = -20;
  
  GenSignalReturn gsr = genSignal(num_points, frequency, filter, n_var, SNR);
  double received = gsr.received;
  double desired = gsr.desired;
  double oise = gsr.noise;

  double P = 0; //eye(sysorder)* 1/delta;
  double weights = 0; //zeros(sysorder,1);
  received = received; //received(:);
  double error = desired*0;
  double lambda1 = 1/lambda;
  int len = 1; //length(received);

  for(int n = sysorder; n < len; n++) {
    double input = 5; //noise(n:-1:n-sysorder+1);

    double K = (lambda1*P*input)/(1+lambda1*input*P*input); //ahh! there was a '* here somewhere...
    double output = weights * input; // there was another '* here!
    //error(n) = received(n)-output;
    weights = weights + K * error; //error(n) instead
    P=(lambda1*P)-(lambda1*K*input*P); // another '* here
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

