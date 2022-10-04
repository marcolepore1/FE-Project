function [sigma, eta, k, MSE, MAPE, mktskew, calskew] = CalibrateNMVMToVolatilitySurfaceConst(setDate, discount, Forward, strikes, nCall, prices, maturity, prevTTM, alpha, x00, flag)
%
% Function that calls one of the calibrations of the parameters of
% ATS NIG at a certain maturity with the constraints coming from previous
% parameters
%  
% INPUT
% setDate:   settlement date considered
% discount:  discount factor @maturity
% Forward:   forward price with expiry @maturity (F(t0,t))
% strikes:   strikes of the options
% nCall:     number of calls between options provided
% prices:    prices of options passed as input
% maturity:  maturity of the options
% prevTTM:   yearfraction between the setDate & previous maturity
% alpha:     parameter of ATS process (0.5 -> NIG, 0 -> VG)
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
% CalibrateNMVMToVolatilitySurfaceConstCall,
% CalibrateNMVMToVolatilitySurfaceConstPut, CalibrateNMVMToVolatilitySurfaceConstCallPut


% decide which calibration to be performed
switch flag
    
    case 0
        
        [sigma, eta, k, MSE, MAPE] = CalibrateNMVMToVolatilitySurfaceConstCall(setDate, discount, Forward, strikes, prices, maturity, prevTTM, alpha, x00);

    case 1

        [sigma, eta, k, MSE, MAPE] = CalibrateNMVMToVolatilitySurfaceConstPut(setDate, discount, Forward, strikes, prices, maturity, prevTTM, alpha, x00);

    case 2

        [sigma, eta, k, MSE, MAPE, mktskew, calskew] = CalibrateNMVMToVolatilitySurfaceConstCallPut(setDate, discount, Forward, strikes, nCall, prices, maturity, prevTTM, alpha, x00);

end

end