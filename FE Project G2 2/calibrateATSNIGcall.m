function [eta, k, sigma, MSE, MAPE] = calibrateATSNIGcall(setDate, optionTable, Forwards, discountCurve)
%
% Function that calibrates time dependent parameters of ATSNIG from data
% provided in optionTable
%
% INPUT
% setDate:          settlement date considered
% optionTable:      struct with option data ([ASK, BID, STRIKES, MATURITY, OPENING, FLAG])
%                   ASK     -> askprice  BID      -> bidprice  OPENING -> opening price 
%                   STRIKES -> strike    MATURITY -> maturity  FLAG    -> 'C' if call, 'P' if put                
% Forwards:         initial forward prices from synthetic fwd bootstrap 
% discountCurve:    struct with [Dates, discounts]
%
% OUTPUT
% eta:              vector of skew @each maturity (ATS NIG time-dependent)
% k:                vector of vol-of-vol @each maturity (ATS NIG time-dependent)
% sigma:            vector of volatility @each maturity (ATS NIG time-dependent)
% MSE:              vector of Mean Squared Error @each maturity 
% MAPE:             vector of Mean Absolute Percentage Error @each maturity 
%
% CALLS
% CalibrateNMVMToVolatilitySurface, CalibrateNMVMToVolatilitySurfaceConst
%

% table with maturities, strikes, 
options            = table;
options.STRIKES    = optionTable.STRIKES;
options.MATURITIES = optionTable.MATURITIES;
options.PRICES     = 0.5*(optionTable.ASK+optionTable.BID);

% data
IBDaycount = 3;
maturities = discountCurve.dates;
discounts  = discountCurve.discounts;
alpha      = 1/2;
sigma      = zeros(length(maturities),1);
eta        = zeros(length(maturities),1);
k          = zeros(length(maturities),1);
MSE        = zeros(length(maturities), 1);
MAPE       = zeros(length(maturities), 1);
TTM        = yearfrac(setDate, maturities, IBDaycount);
mktskew    = zeros(length(maturities),1);
Calskew    = zeros(length(maturities),1);


%% calibration
x00                      = [0.3; 1; TTM(1)];
idx                      = find(~((optionTable.MATURITIES - maturities(1))));
OTMcalls                 = options.PRICES(idx(1):2:idx(end-1)).*(Forwards(1) < options.STRIKES(idx(1):2:idx(end-1)));     
OTMcalls                 = OTMcalls(find(OTMcalls));
OTMstrikes               = options.STRIKES(idx(1):2:idx(end-1)).*(Forwards(1) < options.STRIKES(idx(1):2:idx(end-1)));
OTMstrikes               = OTMstrikes(find(OTMstrikes));
[sigma(1), eta(1), k(1), MSE(1), MAPE(1)] = CalibrateNMVMToVolatilitySurface(setDate, discounts(1), Forwards(1), OTMstrikes, length(OTMcalls),...
    OTMcalls, maturities(1), alpha, x00, 0);



for i=2:length(maturities)
    idx             = find(~((optionTable.MATURITIES - maturities(i))));
    OTMcalls        = options.PRICES(idx(1):2:idx(end-1)).*(Forwards(i) < options.STRIKES(idx(1):2:idx(end-1)));     
    OTMcalls        = OTMcalls(find(OTMcalls));
    OTMstrikes      = options.STRIKES(idx(1):2:idx(end-1)).*(Forwards(i) < options.STRIKES(idx(1):2:idx(end-1)));
    OTMstrikes      = OTMstrikes(find(OTMstrikes));
    x00             = [sigma(i-1), eta(i-1), k(i-1)]';
    [sigma(i), eta(i), k(i), MSE(i), MAPE(i)] = CalibrateNMVMToVolatilitySurfaceConst(setDate, discounts(i), Forwards(i), OTMstrikes, length(OTMcalls),...
        OTMcalls, maturities(i), TTM(i-1), alpha, x00, 0);
%     [sigma(i), eta(i), k(i)] = CalibrateNMVMToVolatilitySurface(setDate, discounts(i), Forwards(i), OTMstrikes,...
%         OTMcalls, maturities(i), alpha, x00);
end

%% constraints
% verifying the constraints are satisfied with the calibrated parameters
g1 = 0.5 + eta - sqrt((0.5 + eta).^2 + 2*(1-alpha)./(sigma.^2.*k));
g2 = -0.5 - eta - sqrt((0.5 + eta).^2 + 2*(1-alpha)./(sigma.^2.*k));
g3 = (TTM.^(1/alpha).*sigma.^2)./(k.^((1-alpha)/alpha)).*sqrt((0.5 + eta).^2 + 2*(1-alpha)./(sigma.^2.*k));

% plot of the constraints
figure()
subplot(3,1,1)
plot(TTM, g1, '-', 'Marker', '*', 'LineWidth',2)
grid on
title('g1')
xlabel('t')
subplot(3,1,2)
plot(TTM, g2, '-', 'Marker', '*', 'LineWidth',2)
grid on
title('g2')
xlabel('t')
subplot(3,1,3)
plot(TTM, g3, '-', 'Marker', '*', 'LineWidth',2)
grid on
title('g3')
xlabel('t')





end
