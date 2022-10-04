function [eta, k, sigma, MSE, MAPE] = calibrateNIGcallput(setDate, optionTable, Forwards, discountCurve)
%
% Function that calibrates NIG model on the whole volatility surface
% whose data are stored in optionTable
%
% INPUT
% setDate:           settlement date of the contracts
% optionTable:       table with data of call & put options
% Forwards:          F(t0, T) T @options maturities
% discountCurve:     struct with [maturities, discounts]
%
% OUTPUT
% eta:               skew
% k:                 vol-of-vol
% sigma:             volatility
% MSE:               vector with Mean Squared errors @maturities wrt prices
% MAPE:              vector with Mean Absolute % errors @maturities wrt prices               
%
% FUNCTIONS CALLED:
% FFTparameters
% CallPricesNMVMFFT
% funNIG

%% useful data
IBDaycount = 3;
maturities = discountCurve.dates;
discounts  = discountCurve.discounts;
alpha      = 1/2;
M          = 15;


%% contracts table
% table with maturities, strikes & prices
options            = table;
options.STRIKES    = optionTable.STRIKES;
options.MATURITIES = optionTable.MATURITIES;
options.PRICES     = 0.5*(optionTable.ASK+optionTable.BID);

%% preallocate data
TTM        = yearfrac(setDate, maturities, IBDaycount);
nCalls         = [0; zeros(length(maturities),1)];
nPuts          = [0; zeros(length(maturities),1)];
idxCalls       = 0;
idxPuts        = 0;
OTMcalloptions = nan(length(options.PRICES), 6);
OTMputoptions  = nan(length(options.PRICES), 6);
rates          = zeros(length(TTM), 1);

% parameters for FFT
Params = FFTparameters(M, 0.0025, 1);

%% calibration
for i=1:length(maturities)
    idx                           = find(~((optionTable.MATURITIES - maturities(i))));
    OTMcalls                      = options.PRICES(idx(1):2:idx(end-1)).*(Forwards(i) < options.STRIKES(idx(1):2:idx(end-1)));     
    OTMcalls                      = OTMcalls(find(OTMcalls));
    OTMstrikesCall                = options.STRIKES(idx(1):2:idx(end-1)).*(Forwards(i) < options.STRIKES(idx(1):2:idx(end-1)));
    OTMstrikesCall                = OTMstrikesCall(find(OTMstrikesCall));
    OTMputs                       = options.PRICES(idx(2):2:idx(end)).*(Forwards(i) > options.STRIKES(idx(2):2:idx(end)));     
    OTMputs                       = OTMputs(find(OTMputs));
    OTMstrikesPut                 = options.STRIKES(idx(2):2:idx(end)).*(Forwards(i) > options.STRIKES(idx(2):2:idx(end)));
    OTMstrikesPut                 = OTMstrikesPut(find(OTMstrikesPut));
    moneynessCall                 = log(Forwards(i)./OTMstrikesCall);
    moneynessPut                  = log(Forwards(i)./OTMstrikesPut);
    nCalls(i+1)                   = length(OTMcalls);
    nPuts(i+1)                    = length(OTMputs);
    OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),:) = [OTMcalls, moneynessCall, OTMstrikesCall, repmat([Forwards(i), discounts(i), TTM(i)], nCalls(i+1),1)];
    OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),:)   = [OTMputs, moneynessPut, OTMstrikesPut, repmat([Forwards(i), discounts(i), TTM(i)], nPuts(i+1),1)];
    idxCalls                      = idxCalls+nCalls(i+1);
    idxPuts                       = idxPuts+nPuts(i+1);
    rates(i)                      = -log(discounts(i))/TTM(i);
end

% reduce the matrices to the non-nan values
OTMcalloptions = OTMcalloptions(find(~isnan(OTMcalloptions(:,1))), :);
OTMputoptions  = OTMputoptions(find(~isnan(OTMputoptions(:,1))), :);

%% Constraints of NIG model
function [c,ceq] = constraints (cal)
  w   = 1/(2*cal(3).*cal(1).^2);  
  c   = [-w-cal(2),...
      -cal(1)];
  ceq = [];
end

%% calibration
x0      = [0.1, 1, 0.1]';
LB      = [0 , -5, 0]';
UB      = [1, 30, 5]';
options = optimoptions('fmincon', 'Display','off');
cal     = fmincon(@(cal) funNIG(Forwards, discounts, OTMcalloptions, OTMputoptions, TTM, cal, Params, alpha), x0, [],[], [], [], LB, UB, @(cal)constraints(cal), options);

sigma = cal(1);
eta   = cal(2);
k     = cal(3);

%% checking results (prices with calibrated parameters)
CalPricesCall = zeros(length(OTMcalloptions), 1);
CalPricesPut  = zeros(length(OTMputoptions), 1);
MSE           = zeros(length(maturities), 1);
MAPE          = zeros(length(maturities), 1);
idxCalls      = 0;
idxPuts       = 0;

for i=1:length(maturities)
    CalPricesCall(idxCalls+1:idxCalls+nCalls(i+1)) = real(CallPricesNMVMFFT(Forwards(i), discounts(i), OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),2), TTM(i), cal, Params, alpha));
    CalPricesPut(idxPuts+1:idxPuts+nPuts(i+1))     = real(CallPricesNMVMFFT(Forwards(i), discounts(i), OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),2), TTM(i), cal, Params, alpha)) - discounts(i).*(Forwards(i) - OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1), 3))';
    CallVols                                       = blkimpv(Forwards(i), OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),3), rates(i), TTM(i), CalPricesCall(idxCalls+1:idxCalls+nCalls(i+1)));
    PutVols                                        = blkimpv(Forwards(i), OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),3), rates(i), TTM(i), CalPricesPut(idxPuts+1:idxPuts+nPuts(i+1)), "Class", "put");
    mktvolsCall                                    = blkimpv(Forwards(i), OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),3), rates(i), TTM(i), OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),1));
    mktvolsPut                                     = blkimpv(Forwards(i), OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),3), rates(i), TTM(i), OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),1), "Class","put");
    strikes                                        = [OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),3); OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),3)];
    mktvols                                        = [mktvolsCall; mktvolsPut];
    CalVols                                        = [CallVols; PutVols];
    MAPE(i)                                        = mean(100*abs([CalPricesCall(idxCalls+1:idxCalls+nCalls(i+1)); CalPricesPut(idxPuts+1:idxPuts+nPuts(i+1))] ...
        - [OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),1); OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),1)])./[OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),1); OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),1)]);
    MSE(i)                                         = mean(abs([CalPricesCall(idxCalls+1:idxCalls+nCalls(i+1)); CalPricesPut(idxPuts+1:idxPuts+nPuts(i+1))] ...
        - [OTMcalloptions(idxCalls+1:idxCalls+nCalls(i+1),1); OTMputoptions(idxPuts+1:idxPuts+nPuts(i+1),1)]).^2)';

    % plot of implied vols of calibrated prices vs mkt implied vols
    moneyness = log(strikes./Forwards(i));
    figure()
    plot(moneyness, CalVols, '+', 'MarkerSize', 6)
    hold on
    plot(moneyness, mktvols, 'square', 'MarkerSize', 6)
    grid on
    t = text(min(moneyness), min(mktvols), ['\bf MSEPrices : ', num2str(MSE(i))], 'Color', 'k');  % MSE
    t.FontSize = 10;
    title('Smile NIG@', num2str(datestr(maturities(i))))
    legend('CalVol', 'mktVols')

    idxCalls                                       = idxCalls+nCalls(i+1);
    idxPuts                                        = idxPuts+nPuts(i+1);
end


end
