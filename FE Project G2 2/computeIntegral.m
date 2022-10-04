function I = computeIntegral(f, moneyness, modelParams, numericalParams, method)
%
% INPUT
% f:                  Function to compute the Fourier transform
% moneyness:          grid of values in which evaluate
% modelParams:        model parameters, e.g. empty
% numericalParams:    parameters of the numerical method as a struct, e.g.
%                     x1, xN, dx, z1, zN, dz, M for fft
%                     x1, xN for quadrature
% method:             flag 1 for fft, 0 for quadrature
%
% OUTPUT
% I:            Value of the integral
%
%

z = [numericalParams.z1:numericalParams.dz:numericalParams.zn];
I = zeros(length(moneyness),1);

if method == 1 || nargin < 5

    N         = 2^numericalParams.M;
    x         = [numericalParams.x1:numericalParams.dx:numericalParams.xn];
    fj        = f(z).*exp(-1i*numericalParams.dz*[0:N-1]*numericalParams.x1);
    prefactor = numericalParams.dz*exp(-1i*numericalParams.z1*x);
    I_FFT     = real(prefactor.*(fft(fj)));
%     hold on
%     plot(z, I_FFT,'LineWidth', 2)
    I         = interp1(x,I_FFT,moneyness);

elseif method == 0
    g = @(z,x) f(x).*exp(- 1i*z.*x);
    integralfun = @(z) quadgk(@(x) g(z,x), numericalParams.z1, numericalParams.zn);
%     hold on
%     plot(z, real(arrayfun(integralfun,z)),'--' , 'LineWidth', 4)
    I = real(arrayfun(integralfun,moneyness));
    
end 
   
   

end