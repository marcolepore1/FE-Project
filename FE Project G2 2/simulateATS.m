function [X] = simulateATS(M, Nsim, flagspline, t, s, sigma_t, sigma_s, eta_t, eta_s, k_t, k_s, flag)
%
% Function that calls different methods to invert the CDF according to
% flagspline
%
% INPUT 
% self-explainatory variables (as simulateATSFFT, simulateATSSUM inputs)
% flag:     0 -> FFT  1 -> the sum formula evaluating point by point
%
% OUTPUT
% X:        simulated increments of the ATS process (f(s,t))
%
% CALLS
% simulateATSFFT, simulateATSSUM
%
% NB: simulateATSSUM is really slow and implemented just for completeness

% decide which method to be used 
switch flag

    case 0 
        
        [X] = simulateATSFFT(M, Nsim, flagspline, t, s,sigma_t,sigma_s,eta_t,eta_s,k_t,k_s);

    case 1
        % not efficient
        [X] = simulateATSSUM(M, Nsim, flagspline, t, s,sigma_t,sigma_s,eta_t,eta_s,k_t,k_s);

end

end

