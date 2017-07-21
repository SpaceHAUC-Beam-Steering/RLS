%RLS beamsteering algorithm
%Jake Hogan
function[weights,error] = RLS(lambda,delta,sysorder)
%% Recursive Least Squares
% Input arguments: 
%   lambda = forgetting factor, dim 1x1 
%   sysorder = filter order, dim 1x1 
%   received = input signal, dim Nx1 
%   desired = desired signal, dim Nx1
%   delta = initial value, P(0)=delta^-1*(Identity Matrix), dim 1x1 
% Output arguments:
%   error = a priori estimation error, dim nx1
%   weights = final filter coefficients, dim sysorderx1

%% create signal    
    % Data Parameters
    numPts  = 1000;             % number of points to generate
    freq    = 12e9;              % frequency of fundamental tone
    filtord = 16;                % filter order
    filt    = rand(filtord, 1); % filter coefficients
    nVar    = 1;                % white noise variance
    SNR     = -20;              % signal to noise ratio of tone
    
    [received,desired,noise] = genSignal(numPts, freq, filt, nVar, SNR);
    
%% intitalize algorithm
%lambda = ;
%sysorder = ;
%delta = ;
P=eye(sysorder)* 1/delta; % identity matrix with delta as diagonals
weights=zeros(sysorder,1); % zeroing out initial weights
received = received(:);% makes the signal column vectors
error = desired*0;
lambda1 = (lambda)^-1;% inverse of forgetting factor
len = length(received);

%% begin
for n = sysorder:len
    input = noise(n:-1:n-sysorder+1);% corrects dimensions

    K =(lambda1*P*input)/(1+lambda1*input'*P*input); % calculates Gain matrix
    
    output = weights' * input ; % calculates array output w/ ccurrent weights
    
    error(n) = received(n)-output; % a priori error of the signal
       
    weights = weights + K*error(n); % updates weights
    
    P=(lambda1*P)-(lambda1*K*input'*P); % calculates new correlation matrix
    
    
end

e = error(16:1000);
q = linspace(1,985,985);

e2 = abs(e.^2);



%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------

% Plot filter results
t = linspace(0,length(received)/numPts,length(received));
%figure;
%plot(t,received,t,error);

%title('Comparison of Filtered Signal to Reference Input');


% Plot comparison of results to original signal 
    figure;
    plot(t,desired,t,error);
   
    title('Comparison of Filtered Signal to Original Signal');
%plot learning curve
figure;
plot(q, e2);





% Calculate SNR improvement
SNRi    = 10*log10(var(received)/var(error));

disp([num2str(SNRi) 'dB SNR Improvement'])



return

function[received,desired,noise] = genSignal(numPts, freq, filt, nVar, SNR)

    % Generate time values
    t = linspace(0,1,numPts)';
    
    
    % Generate tone
    desired = sin(2*pi*freq*t);
    
    %wavevector from angle of arrival
    elevation = pi/3;
    azimuth = pi/3;
    g = [];
    for m=1:4
        u = exp(1i*pi*(m-1));
        a = u*sin(elevation)*cos(azimuth);
        g = [g a];
    end
    %for y elements
    h = [];
    for p=1:4
        v = exp(1i*pi*(p-1));
        b = v * sin(elevation)* sin(azimuth);
        h = [h b];
    end
    
    % Generate noise
    noise = sqrt(nVar)*randn(numPts,1);
    addnoise = filter(filt, 1, noise);
    freqz(filt,1,numPts)
    
    % Adjust SNR of tone
    desired = (desired)/sqrt(var(desired)/(10^(SNR/10)*var(noise)));
    disp(['Calculated SNR = ' num2str(10*log10(var(desired)/var(noise)))])
    
    % Add noise to signal
    received = (desired*g*h.') + addnoise ;
    
return
