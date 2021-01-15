function [Noise11, Noise22] = Noise(nSteps, i)
DIM = [nSteps+1,1];
BETA = 0;   % White noise
%BETA = -1; % Brown noise
%BETA = -2; % Pink noise
Noise1 = spatialPattern(DIM,BETA);
Noise2 = spatialPattern(DIM,BETA);
Noise3 = spatialPattern(DIM,BETA);
Noise4 = spatialPattern(DIM,BETA);
Noi1 = mean(Noise1);
Noi2 = mean(Noise2);
Noi3 = mean(Noise3);
Noi4 = mean(Noise4);

Noise11 = [Noise1(i,1)-Noi1;Noise2(i,1)-Noi2];
Noise22 = [Noise3(i,1)-Noi3;Noise4(i,1)-Noi4];
end

function x = spatialPattern(DIM,BETA)

u = [(0:floor(DIM(1)/2)) -(ceil(DIM(1)/2)-1:-1:1)]'/DIM(1);
% Reproduce these frequencies along ever row
u = repmat(u,1,DIM(2));
% v is the set of frequencies along the second dimension.  For a square
% region it will be the transpose of u
v = [(0:floor(DIM(2)/2)) -(ceil(DIM(2)/2)-1:-1:1)]/DIM(2);
% Reproduce these frequencies along ever column
v = repmat(v,DIM(1),1);

% Generate the power spectrum
S_f = (u.^2 + v.^2).^(BETA/2);

% Set any infinities to zero
S_f(S_f==inf) = 0;

% Generate a grid of random phase shifts
phi = rand(DIM);

% Inverse Fourier transform to obtain the the spatial pattern
x = ifft2(S_f.^0.5 .* (cos(2*pi*phi)+1i*sin(2*pi*phi)));

% Pick just the real component
x = real(x);
end