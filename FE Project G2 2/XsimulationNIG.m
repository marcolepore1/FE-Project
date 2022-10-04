function [X] = XsimulationNIG(reset_date, NIGparams, Nsim, prevDate)
%
% Function that simulates NIG process increment f(s,t) where s=prevDate
% t=reset_date
%
% INPUT
% reset_date:  second date of the increment
% prevDate:    first date
% NIGparams:   struct with NIG parameters
% Nsim:        number of simulations to be performed
%
% OUTPUT
% X:           simulations performed of f(s,t)            
%
%


%% Day conventions
IBDaycount = 3;

%% Parameters
alpha = 0.5;
k     = unique(NIGparams.k);
eta   = unique(NIGparams.eta);
sigma = unique(NIGparams.sigma);

TTM_1 = yearfrac(prevDate, reset_date, IBDaycount);
g1    = randn(Nsim, 1); % generation of a std normal

%% NTS simulation
G1      = random('InverseGaussian', 1, TTM_1/k, Nsim, 1); 
exp_lap = @(w,TTM) TTM/k *  (1 - alpha)/alpha * (1-(1+(k*w*sigma^2)/(1-alpha)).^alpha); 
X       = -exp_lap(eta,TTM_1) - (1/2+eta)*sigma^2*TTM_1*G1 + sqrt(G1).*sigma*sqrt(TTM_1).*g1; 


end