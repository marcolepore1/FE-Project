function [prices] = CallPricesNMVMFFT( forward, discount, moneyness, timeToMaturity, params,numericalMethodParameters, alpha)
%
% Compute the prices of Calls with FFT method
%
% INPUT 
%
% forward:                    forward value at settlement date
% discount:                   1 year discount
% moneyness:                  grid of the moneyness
% timeToMaturity:             time to maturity
% params:                     sigma, eta, k of LTS/ATS process
% numericalMethodParameters:  parameters of the numerical method as a struct x1, xN, dx, z1, zN, dz, M for fft
% alpha:                      stability parameter (0.5 -> NIG, 0 -> VG)
%
% OUTPUT
% prices:                     Call option prices
%
% CALLS
% computeIntegral
%

%% model params
sigma = params(1);
eta   = params(2);
k     = params(3);

% loglaplace exponent
LaplaceExpTS   = @(x) (timeToMaturity/k) * ((1-alpha)/alpha) * (1 - (1 + (x.*k.*sigma^2)./(1-alpha)).^alpha);

% char function
phi            = @(x) exp(-1i.*x.*LaplaceExpTS(eta)).*exp(LaplaceExpTS((x.^2 + 1i*(1+2.*eta).*x)./2));

% integrand function
fTS            = @(v) 1./(v.^2 + 0.25);
LevyFunctionTS = @(v) 1/(2*pi).*phi(-v-1i/2).*fTS(v);
Integrals      =  computeIntegral(LevyFunctionTS, moneyness, [], numericalMethodParameters, 1);
prices         = ((1 - exp(-moneyness./2).* Integrals)*discount*forward)';

end


