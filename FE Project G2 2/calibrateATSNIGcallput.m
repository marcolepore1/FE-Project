function [eta, k, sigma, MSE, MAPE] = calibrateATSNIGcallput(setDate, optionTable, Forwards, discountCurve)
%
% Function that calibrates time dependent parameters of ATSNIG from data
% provided in optionTable (only OTM options are used)
%
% INPUT
% setDate:         settlement date considered
% optionTable:     table with data of call & put options
% Forwards:        F(t0, T) @options maturities
% discountCurve:   struct with [maturities, discounts]
%
% OUTPUT
% eta:             skew
% k:               vol-of-vol
% sigma:           volatility
% MSE:             vector of Mean Squared error wrt prices @maturities
% MAPE:            vector of Mean absolute % error wrt prices @maturities                
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
% k0                       = [0.0270, 0.1075, 0.2023, 0.3231, 0.4192, 0.4899, 0.6142, 0.7090, 0.8826, 1.4, 1.8]';
% eta0                     = [29.6742, 21.9358, 19.2525, 16.7753, 15.2209, 14.5, 14.2, 13.5, 13.0, 12.5, 11.0]';
% vol0                     = [0.1270, 0.1209, 0.1176, 0.1162, 0.1168, 0.1145, 0.1164, 0.1146, 0.1098, 0.1088, 0.08]';
% x00                      = [0.3; 1; TTM(1)];
x00                      = [0.2; 1; TTM(1)];
idx                      = find(~((optionTable.MATURITIES - maturities(1))));
OTMcalls                 = options.PRICES(idx(1):2:idx(end-1)).*(Forwards(1) < options.STRIKES(idx(1):2:idx(end-1)));     
OTMcalls                 = OTMcalls(find(OTMcalls));
OTMstrikesCall           = options.STRIKES(idx(1):2:idx(end-1)).*(Forwards(1) < options.STRIKES(idx(1):2:idx(end-1)));
OTMstrikesCall           = OTMstrikesCall(find(OTMstrikesCall));
OTMputs                  = options.PRICES(idx(2):2:idx(end)).*(Forwards(1) > options.STRIKES(idx(2):2:idx(end)));     
OTMputs                  = OTMputs(find(OTMputs));
OTMstrikesPut            = options.STRIKES(idx(2):2:idx(end)).*(Forwards(1) > options.STRIKES(idx(2):2:idx(end)));
OTMstrikesPut            = OTMstrikesPut(find(OTMstrikesPut));
OTMoptions               = [OTMcalls; OTMputs];
OTMstrikes               = [OTMstrikesCall; OTMstrikesPut];
[sigma(1), eta(1), k(1), MSE(1), MAPE(1), mktskew(1), Calskew(1)] = CalibrateNMVMToVolatilitySurface(setDate, discounts(1), Forwards(1), OTMstrikes, length(OTMcalls),...
    OTMoptions, maturities(1), alpha, x00, 2);



for i=2:length(maturities)
    idx            = find(~((optionTable.MATURITIES - maturities(i))));
    OTMcalls       = options.PRICES(idx(1):2:idx(end-1)).*(Forwards(i) < options.STRIKES(idx(1):2:idx(end-1)));     
    OTMcalls       = OTMcalls(find(OTMcalls));
    OTMstrikesCall = options.STRIKES(idx(1):2:idx(end-1)).*(Forwards(i) < options.STRIKES(idx(1):2:idx(end-1)));
    OTMstrikesCall = OTMstrikesCall(find(OTMstrikesCall));
    OTMputs        = options.PRICES(idx(2):2:idx(end)).*(Forwards(i) > options.STRIKES(idx(2):2:idx(end)));     
    OTMputs        = OTMputs(find(OTMputs));
    OTMstrikesPut  = options.STRIKES(idx(2):2:idx(end)).*(Forwards(i) > options.STRIKES(idx(2):2:idx(end)));
    OTMstrikesPut  = OTMstrikesPut(find(OTMstrikesPut));
    OTMoptions     = [OTMcalls; OTMputs];
    OTMstrikes     = [OTMstrikesCall; OTMstrikesPut];
%     x002            = [vol0(i-1), eta0(i-1), k0(i-1)]';
    x00            = [sigma(i-1), eta(i-1), k(i-1)]';
    [sigma(i), eta(i), k(i), MSE(i), MAPE(i), mktskew(i), Calskew(i)] = CalibrateNMVMToVolatilitySurfaceConst(setDate, discounts(i), Forwards(i), OTMstrikes, length(OTMcalls),...
        OTMoptions, maturities(i), TTM(i-1), alpha, x00, 2);
%     [sigma(i), eta(i), k(i), MSE(i), MAPE(i), mktskew(i), Calskew(i)] = CalibrateNMVMToVolatilitySurfaceConstCallPutN(setDate, discounts(i), Forwards(i), OTMstrikes, length(OTMcalls),...
%         OTMoptions, maturities(i), TTM(i-1), alpha, x00, x002);
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

%% Skew reproduction
figure()
plot(TTM, mktskew, '--', 'Marker', 'diamond', 'LineWidth', 2)
hold on 
plot(TTM, Calskew, '--', 'Marker', 'diamond', 'LineWidth', 2)
grid on
legend('mktskew', 'Calskew')
title('volatility skews')


end