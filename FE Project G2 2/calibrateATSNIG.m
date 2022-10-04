function [eta, k, sigma, MSE, MAPE] = calibrateATSNIG(setDate, optionTable, Forwards, discountCurve, flag)
%
% Function that calls calibration of ATS NIG model, using
% different options according to flag
% 
% INPUT
% setDate:          settlement date considered
% optionTable:      struct with option data ([ASK, BID, STRIKES, MATURITY, OPENING, FLAG])
%                   ASK     -> askprice  BID      -> bidprice  OPENING -> opening price 
%                   STRIKES -> strike    MATURITY -> maturity  FLAG    -> 'C' if call, 'P' if put                
% Forwards:         initial forward prices from synthetic fwd bootstrap 
% discountCurve:    struct with [Dates, discounts]
% flag:             0 -> calibration only on calls, 1 -> calibration only 
%                   on puts, 2 -> calibration on both (always exploited)
%
% OUTPUT
% eta:              vector of skew @each maturity (ATS NIG time-dependent)
% k:                vector of vol-of-vol @each maturity (ATS NIG time-dependent)
% sigma:            vector of volatility @each maturity (ATS NIG time-dependent)
% MSE:              vector of Mean Squared Error @each maturity 
% MAPE:             vector of Mean Absolute Percentage Error @each maturity 
%
% CALLS:
% calibrateATSNIGcall, calibrateATSNIGput, calibrateATSNIGcallput
%

% decide on which options calibrate
switch flag
    case 0
        [eta, k, sigma, MSE, MAPE] = calibrateATSNIGcall(setDate, optionTable, Forwards, discountCurve);

    case 1
        [eta, k, sigma, MSE, MAPE] = calibrateATSNIGput(setDate, optionTable, Forwards, discountCurve);
    
    case 2
        [eta, k, sigma, MSE, MAPE] = calibrateATSNIGcallput(setDate, optionTable, Forwards, discountCurve);

end
