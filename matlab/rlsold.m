%function [error,weights]=rls3(lambda,sysorder) 
%% Recursive Least Squares
%  Call: 
%  Input arguments: 
%   lambda = forgetting factor, dim 1x1 
%   M = filter length, dim 1x1 
%   received = input signal, dim Nx1 
%   desired = desired signal, dim Nx1
%   delta = initial value, P(0)=delta^-1*I, dim 1x1 
% Output arguments:
%   error = a priori estimation error, dim Nx1
%   weights = final filter coefficients, dim Mx1
% inital values
%w = 0, P = delta*I
%% run with 
   %lambda = .95 ;
    sysorder = 16;
%% received and desired signals
%see{wavevector.mfile}

%(wrong way to represent it i think)
frequency = 1e9;
time = linspace(0,1e-7,100);
desired = (1:100)'*sin(2*pi*frequency);
%received = sin(2*pi*(frequency)*time*3);%(abitrary constant to make it different));
%M x N planar array, w/ spacing wavelength/2

theta = 30;
phi = 20;
%for x elements
x = [];
for n=1:4
    e = exp(1i*pi*(n-1));
    a = e * sin(theta)*cos(phi);
    x = [x a];
end
%for y elements
y = [];
for n=1:4
    e = exp(1i*pi*(n-1));
    a = e * sin(theta)* sin(phi);
    y = [y a];
end
received = (desired*x*y');%+noise;
received = real(received);


%% intitalize algorithm
delta = 100* var(received);
P=eye(sysorder)/delta;
weights=zeros(sysorder,1); 
MSE = [];
received = received(:);
desired = desired(:);
error = desired;
inputlength = length(received);
%% start algorithm
for n=sysorder:inputlength 
    uvec=(received(n:-1:n-sysorder+1));
    K =(lambda^(-1)*P*uvec)/(1+lambda^(-1)*uvec'*P*uvec); 
    error = desired(n)-weights'*uvec;
    weights=weights+K*conj(error); 
    P=lambda^(-1)*P-lambda^(-1)*K*uvec'*P;
    e = error^2;
    MSE = [MSE e];
end
plot(MSE)


