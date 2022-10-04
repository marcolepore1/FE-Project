function [sigma, eta, k, MSE, MAPE, mktSkew, calSkew] = CalibrateNMVMToVolatilitySurfaceCallPut(setDate, discount, Forward, strikes, nCall, prices, maturity, alpha, x00)
%
% Function that calibrates the parameters of an ATS NIG
% process at a certain maturity, making use of OTM calls & puts whose prices are
% stored in prices. 
% It performs also a plot of implied volatility smile obtained and computes the most common errors
% MSE & MAPE
%  
% INPUT
% setDate:   settlement date considered
% discount:  discount factor @maturity
% Forward:   forward price with expiry @maturity (F(t0,t))
% strikes:   strikes of the options
% nCall:     number of calls between options provided
% prices:    prices of options passed as input
% maturity:  maturity of the options
% alpha:     parameter of ATS process (0.5 -> NIG, 0 -> VG)
% x00:       starting point of the calibration
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


%% Displaying data
IBdaycount = 3;
TTM        = yearfrac(setDate, maturity, IBdaycount);
moneyness  = log(Forward./strikes);
rate       = -log(discount)/TTM;

%% Model Prices
% Calibrating  parameters for the FFT
I_res     = 2*pi*exp(-sign(moneyness)*0.5.*moneyness);
M         = 15;
options   = optimset('TolFun',1e-5);
x0        = 0.0025;
LB        = eps;
UB        = 0.01;
fTS       = @(v) 1./(v.^2 + 0.25);
ff        = @(x) 1./(x.^2 + 1/4);
% dz        = 0.0025; % find dz such that FFT replicates the residual integral
dz        = lsqnonlin(@(dz) abs( FourierTransform(fTS, moneyness,  M, dz)- I_res), x0, LB, UB, options); % find dz such that FFT replicates the residual integral
Params    = FFTparameters(M, dz, 1);
I         = computeIntegral(ff, moneyness, [], Params, 1);
MdlPricesCall = @(cal) CallPricesNMVMFFT(Forward, discount, moneyness(1:nCall), TTM, cal, Params, alpha);
MdlPricesPut  = @(cal) CallPricesNMVMFFT(Forward, discount, moneyness(nCall+1:end), TTM, cal, Params, alpha) - discount.*(Forward - strikes(nCall+1:end))'; 
MdlPrices     = @(cal) [MdlPricesCall(cal), MdlPricesPut(cal)];
errorFFTINT = abs(I - I_res);

% 
% figure()
% plot(moneyness, I_res, '*b', 'LineWidth', 2)
% hold on
% plot(moneyness, I, '+g', 'LineWidth', 2)
% grid on
% text(moneyness - 0.0002, min(I, I_res) , num2str(abs(I - I_res)),'Color', 'g', 'FontSize', 8); % error for each point 
% tINT = text(min(moneyness), I(2), ['\bf Error : ', num2str(sum(errorFFTINT))], 'Color', 'g');  % total error
% tINT.FontSize = 13;
% legend('Residuals', 'FFT')
% hold off

% NIG constraints (passed to fmincon)
% Matlab help : easy way to implement a non linear constraint 
function [c,ceq] = constraints (cal)
  w   = 1/(2*cal(3).*cal(1).^2);  
  c   = [-w-cal(2),...
      -cal(1)];
  ceq = [];
end

% starting points, lower bounds & upper bounds
% x00 = [0.3; 3; TTM];     % sigma, eta, k
% LB = [0.17; 8; 0.2];
LB = [0.001; -5; 0.001];   % sigma, eta, k
% UB = [1; 20; 5];
UB = [1; 50; 5];

aux = @(cal) sum(abs(prices-MdlPrices(cal)').^2); % function that we want to minimize the weigths are negligible since they are equal to 1

% implied volatility calibration
% aux             = @(cal) sqrt(sum(abs(arrayfun(@(i) blkimpv(Forward, strikes(i), rate, TTM, prices(i)) ...
%     -blkimpv(Forward, strikes(i), rate, TTM, MdlPrices(cal))').^2));


%% Calibration & Results
const = @constraints;

% trying with multistart (no particular improvements)
% opts = optimoptions(@fmincon,'Algorithm','sqp');
% problem = createOptimProblem('fmincon','objective',...
%     aux,'x0', x00,'lb', LB,'ub', UB, 'nonlcon', const, 'options',opts);
% ms    = MultiStart;
% [cal, f] = run(ms,problem,20);

% minimization
cal   = fmincon(@(cal) aux(cal),x00,[],[],[],[],LB,UB,const);

% calibrated parameters
sigma  = cal(1);
eta    = cal(2);
k      = cal(3);

%% Compute Model Prices with calibrated parameters
CalPricesCall = real(CallPricesNMVMFFT(Forward, discount, moneyness(1:nCall), TTM, cal, Params, alpha));
CalPricesPut  = real(CallPricesNMVMFFT(Forward, discount, moneyness(nCall+1:end), TTM, cal, Params, alpha)) - discount.*(Forward - strikes(nCall+1:end))';
CalPrices     = [CalPricesCall, CalPricesPut];

%% Black Volatility obtain from the Calibrated Model Prices 
mktvolsCall   = blkimpv(Forward, strikes(1:nCall), rate, TTM, prices(1:nCall));
CalVolCall    = blkimpv(Forward, strikes(1:nCall), rate, TTM, CalPricesCall');
mktvolsPut    = blkimpv(Forward, strikes(nCall+1:end), rate, TTM, prices(nCall+1:end), 'Class','put');
CalVolPut     = blkimpv(Forward, strikes(nCall+1:end), rate, TTM, CalPricesPut', 'Class', 'put');
mktvols       = [mktvolsCall; mktvolsPut];
CalVol        = [CalVolCall; CalVolPut];

% calibration errors (on prices & volatilities)
PriceSquaredErrors = abs(prices - CalPrices').^2;
PercError          = abs(prices - CalPrices')./prices;
VolErrors          = abs(mktvols - CalVol);
TotVolError        = sum(VolErrors);
MAPE               = 100*mean(PercError);
MSE                = mean(PriceSquaredErrors);
moneyness          = log(strikes./Forward);

%% Skew comparison
mktSkew = mktvolsCall(1) - mktvolsPut(end);
calSkew = CalVolCall(1) - CalVolPut(end);

% plot
figure()
plot(moneyness, CalVol,'+', 'MarkerSize', 5)
hold on 
plot(moneyness, mktvols, 'square', 'MarkerSize',5)
grid on
t = text(min(moneyness), min(mktvols), ['\bf MSEPrices : ', num2str(MSE)], 'Color', 'k');  % MSE
t.FontSize = 10;
legend('Calvols', 'mktvols')
xlabel('Moneyness')
ylabel('Volatilities')
title('Smile ATS @', num2str(datestr(maturity)))

% grid on
% text(moneyness - 0.0002, min(I, I_res) , num2str(abs(I - I_res)),'Color', 'g', 'FontSize', 8); % error for each point 
% tINT = text(min(moneyness), I(2), ['\bf Error : ', num2str(sum(errorFFTINT))], 'Color', 'g');  % total error
% tINT.FontSize = 13;
% legend('Residuals', 'FFT')

  

end