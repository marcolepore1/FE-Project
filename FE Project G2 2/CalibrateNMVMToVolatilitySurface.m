function [sigma, eta, k, MSE, MAPE, mktskew, calskew] = CalibrateNMVMToVolatilitySurface(setDate, discount, Forward, strikes, nCall, prices, maturity, alpha, x00, flag)
%
% Function that calls calibrations of NIG parameters at a
% certain maturity using corresponding options according to flag 
%  
% INPUT
% setDate:   settlement date considered
% discount:  discount factor @maturity
% Forward:   forward price with expiry @maturity (F(t0,t))
% strikes:   strikes of the options
% nCall:     number of calls between options provided
% prices:    prices of options passed as input
% maturity:  maturity of the options
% alpha:     parameter of LTS process (0.5 -> NIG, 0 -> VG)
% x00:       starting point of the calibration
% flag:      0 -> calibration only on calls, 1 -> calibration only 
%            on puts, 2 -> calibration on both (always exploited)
%
% OUTPUT
% sigma:     calibrated model parameter which represents the avarage volatility
% k:         calibrated model parameter which represents 'volatility of the volatility'
% eta:       calibrated model parameter which represents the skewness of the volatility curve
% MSE:       Mean Squared Error of the performed calibration
% MAPE:      Mean Absolute Percentage Error of the performed calibration
% mktSkew:   implied skew of the market (difference between last & first impl vols from market prices)
% calSkew:   reproduced implied skew  (difference between last & first impl vols from calibrated prices)
%
% CALLS 
% fmincon    for the constrained minimization of the distance (NIG constraint only)
% CallPricesNMVMFFT
% 


% decide which calibration to be performed
switch flag
    
    case 0
        [sigma, eta, k, MSE, MAPE] = CalibrateNMVMToVolatilitySurfaceCall(setDate, discount, Forward, strikes, prices, maturity, alpha, x00);
    
    case 1
        [sigma, eta, k, MSE, MAPE] = CalibrateNMVMToVolatilitySurfacePut(setDate, discount, Forward, strikes, prices, maturity, alpha, x00);

    case 2
        [sigma, eta, k, MSE, MAPE, mktskew, calskew] = CalibrateNMVMToVolatilitySurfaceCallPut(setDate, discount, Forward, strikes, nCall, prices, maturity, alpha, x00);

end


end