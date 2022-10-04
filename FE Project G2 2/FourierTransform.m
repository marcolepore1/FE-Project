function I = FourierTransform(f, moneyness, M, dz)
%
% Function that computes the fourier transform with only two input
% parameters p
%
% INPUTS 
% f:                  Function to compute the Fourier transform
% moneyness:          grid of values in which evaluate
% 
%
% OUTPUT
% I:                  Value of the integral with FFT method
%
%


%% Defintion of Parameters
N  = 2^M;
dx = 2*pi/(N * dz);
x1 = - 0.5*(N-1)*dx;
x  = [x1:dx:-x1];
z1 = - 0.5.*(N-1)*dz;
z  = [z1:dz:-z1];

%% Computing the Fourier Transform
fj        = f(x).*exp(-1i*dx*[0:N-1]*z1);
prefactor = dx*exp(-1i*x1*z);
I_FFT     = real(prefactor.*(fft(fj)));
I         = interp1(z,I_FFT,moneyness);

end
