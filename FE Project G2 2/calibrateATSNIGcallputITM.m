function [eta, k, sigma, MSE, MAPE] = calibrateATSNIGcallputITM(setDate, optionTable, Forwards, discountCurve)
%
% Function that calibrates time dependent parameters of ATSNIG from data
% provided in optionTable (only ITM are used)
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
% x00                      = [0.3; 1; TTM(1)];
x00                      = [0.1; 1; TTM(1)];
idx                      = find(~((optionTable.MATURITIES - maturities(1))));
ITMcalls                 = options.PRICES(idx(1):2:idx(end-1)).*(Forwards(1) > options.STRIKES(idx(1):2:idx(end-1)));     
ITMcalls                 = ITMcalls(find(ITMcalls));
ITMstrikesCall           = options.STRIKES(idx(1):2:idx(end-1)).*(Forwards(1) > options.STRIKES(idx(1):2:idx(end-1)));
ITMstrikesCall           = ITMstrikesCall(find(ITMstrikesCall));
ITMputs                  = options.PRICES(idx(2):2:idx(end)).*(Forwards(1) < options.STRIKES(idx(2):2:idx(end)));     
ITMputs                  = ITMputs(find(ITMputs));
ITMstrikesPut            = options.STRIKES(idx(2):2:idx(end)).*(Forwards(1) < options.STRIKES(idx(2):2:idx(end)));
ITMstrikesPut            = ITMstrikesPut(find(ITMstrikesPut));
ITMoptions               = [ITMcalls; ITMputs];
ITMstrikes               = [ITMstrikesCall; ITMstrikesPut];
[sigma(1), eta(1), k(1), MSE(1), MAPE(1)] = CalibrateNMVMToVolatilitySurface(setDate, discounts(1), Forwards(1), ITMstrikes, length(ITMcalls),...
    ITMoptions, maturities(1), alpha, x00, 2);


% calibration with constraints
for i=2:length(maturities)
    idx            = find(~((optionTable.MATURITIES - maturities(i))));
    ITMcalls       = options.PRICES(idx(1):2:idx(end-1)).*(Forwards(i) > options.STRIKES(idx(1):2:idx(end-1)));     
    ITMcalls       = ITMcalls(find(ITMcalls));
    ITMstrikesCall = options.STRIKES(idx(1):2:idx(end-1)).*(Forwards(i) > options.STRIKES(idx(1):2:idx(end-1)));
    ITMstrikesCall = ITMstrikesCall(find(ITMstrikesCall));
    ITMputs        = options.PRICES(idx(2):2:idx(end)).*(Forwards(i) < options.STRIKES(idx(2):2:idx(end)));     
    ITMputs        = ITMputs(find(ITMputs));
    ITMstrikesPut  = options.STRIKES(idx(2):2:idx(end)).*(Forwards(i) < options.STRIKES(idx(2):2:idx(end)));
    ITMstrikesPut  = ITMstrikesPut(find(ITMstrikesPut));
    ITMoptions     = [ITMcalls; ITMputs];
    ITMstrikes     = [ITMstrikesCall; ITMstrikesPut];
%     x00            = [0.3, 5, TTM(i)]
    x00            = [sigma(i-1), eta(i-1), k(i-1)]';
    [sigma(i), eta(i), k(i), MSE(i), MAPE(i), mktskew(i), Calskew(i)] = CalibrateNMVMToVolatilitySurfaceConst(setDate, discounts(i), Forwards(i), ITMstrikes, length(ITMcalls),...
        ITMoptions, maturities(i), TTM(i-1), alpha, x00, 2);
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
% figure()
% plot(TTM, mktskew, '--', 'Marker', 'diamond', 'LineWidth', 2)
% hold on 
% plot(TTM, Calskew, '--', 'Marker', 'diamond', 'LineWidth', 2)
% grid on
% legend('mktskew', 'Calskew')
% title('volatility skews')


end