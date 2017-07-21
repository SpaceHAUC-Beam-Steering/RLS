%M x N planar array, w/ spacing wavelength/2
function[received] = wavevector(M,N,theta,phi,desired)
%for x elements
x = [];
for n=1:M
    e = exp(1i*pi*(n-1));
    a = e * sin(theta)*cos(phi);
    x = [x a];
end
%for y elements
y = [];
for n=1:N
    e = exp(1i*pi*(n-1));
    a = e * sin(theta)* sin(phi);
    y = [y a];
end
received = desired*x*y';
end


