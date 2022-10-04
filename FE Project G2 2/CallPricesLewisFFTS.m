function [price, IC95, SD] = CallPricesLewisFFTS(forward, discount, moneyness, timeToMaturity, params, M, Nsim, flagspline)
% 
% Computes the prices of call options with a grid of moneyness exploiting
% Lewis FFT simulation of the underlying prices up to maturity date,
% exploiting spline or linear according to flagspline
%
% INPUT
% forward:          initial price of the forward F(t0, t) t=maturity
% discount:         discount @maturity
% moneyness:        grid of the moneyness of the options
% timeToMaturity:   yearfraction up to maturity of the options (same for all one)
% params:           ATS parameters @maturity
% M:                discretization parameter of FFT grid
% Nsim:             number of simulations of the underlying with LewisFFTS
% flagspline:       1 -> simulate with spline, 2 -> simulate with linear
%
% OUTPUT
% price:            grid of prices of the options
% IC95:             confidence interval with confidence level 95%
% SD:               standard deviation of MC
%
% CALLS
% simulateATS with FFT method (no SUM!)
% 

%% Setting data
rng(3)                              % to compare results with linear and spline
sigma_t = params(1);
eta_t   = params(2);
k_t     = params(3);

%% simulation of the underlying
[X]     = simulateATS(M, Nsim, flagspline, timeToMaturity, 0, sigma_t, params(1), eta_t, params(2), k_t, params(3), 0);
St      = forward.*exp(X); %forward=F_0   St==  S_T=FT_T
K       = forward.*exp(-moneyness); % grid of strikes


%% Pricing with a Montecarlo-like approach
price   = zeros(length(K), 1);
errors  = zeros(length(K), 1);
MCdev   = zeros(length(K), 1);

for i = 1:length(K)
    
    payoff    = max(St - K(i), 0);
    errors(i) = std(payoff)/sqrt(Nsim);
    MCdev(i)  = std(payoff);
    price(i)  = discount.*mean(payoff);

end

%% confidence intervals & SD
SD   = mean(MCdev); 
IC95 = [price - 1.96.*errors, price + 1.96.*errors];



end