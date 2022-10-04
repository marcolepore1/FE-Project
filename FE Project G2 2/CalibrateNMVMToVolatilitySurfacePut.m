function [sigma, eta, k, MSE, MAPE] = CalibrateNMVMToVolatilitySurfacePut(setDate, discount, Forward, strikes, prices, maturity, alpha, x00)
%
% Function that calibrates the parameters of an ATS NIG
% process at a certain maturity, making use of OTM puts only, whose prices are
% stored in prices. 
% It performs also a plot of implied volatility smile obtained and computes the most common errors
% MSE & MAPE
%  
% INPUT
% setDate:   settlement date considered
% discount:  discount factor @maturity
% Forward:   forward price with expiry @maturity (F(t0,t))
% strikes:   strikes of the options
% prices:    prices of options passed as input
% maturity:  maturity of the options
% alpha:     parameter of LTS process (for example 0.5 -> NIG, 0 -> VG)
% x00:       starting point of the calibration
%
% OUTPUT
% sigma:     calibrated model parameter which represents the avarage volatility
% k:         calibrated model parameter which represents 'volatility of the volatility'
% eta:       calibrated model parameter which represents the skewness of the volatility curve
% MSE:       Mean Squared Error of the performed calibration
% MAPE:      Mean Absolute Percentage Error of the performed calibration
%
% CALLS 
% fmincon    for the constrained minimization of the distance (NIG constraint only)
% CallPricesNMVMFFT
% 
%% Set data
IBdaycount = 3;
TTM        = yearfrac(setDate, maturity, IBdaycount);
moneyness  = log(Forward./strikes);
rate       = -log(discount)/TTM;

%% Model Prices (tempered stable with ɑ = 1/3)
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
dz        = lsqnonlin(@(dz) abs(arrayfun(@(i) FourierTransform(fTS, moneyness(i),  M, dz)- I_res(i), [1:length(moneyness)]')), x0, LB, UB, options); % find dz such that FFT replicates the residual integral
Params    = FFTparameters(M, dz, 1);
I         = computeIntegral(ff, moneyness, [], Params, 1);
MdlPrices = @(cal) CallPricesNMVMFFT(Forward, discount, moneyness, TTM, cal, Params, alpha) - discount.*(Forward - strikes)';
errorFFTINT = abs(I - I_res); % precision of FFT calibration

% Uncomment for a plot of FFT calibration
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

% Matlab help : easy way to implement a non linear constraint 
function [c,ceq] = constraints (cal)
  w   = 1/(2*cal(3).*cal(1).^2);  
  c   = [-w-cal(2),...
      -cal(1)];
  ceq = [];
end

%% Calibration & Results
% constraints of the calibration
% x00 = [0.3; 3; TTM]; % sigma, eta, k
LB = [0.001; -5; 0.001];
UB = [1; 50; 5];

% x00             = [0.1, 3, TTM]';  % sigma, eta, k
% LB              = [0.01, 1, eps]';
% UB              = [0.5, 20, 1.5]';
aux   = @(cal) sum(abs(prices-MdlPrices(cal)').^2); % function that we want to minimize the weigths are negligible since they are equal to 1
% aux             = @(cal) sqrt(sum(abs(arrayfun(@(i) blkimpv(Forward, strikes(i), rate, TTM, prices(i)) ...
%     -blkimpv(Forward, strikes(i), rate, TTM, MdlPrices(cal))').^2));

const = @constraints;
cal   = fmincon(@(cal) aux(cal),x00,[],[],[],[],LB,UB,const);

% cal    =  ...
%      lsqnonlin(@(cal) abs(arrayfun(@(i) blkimpv(Forward, strikes(i), rate, TTM, prices(i)) ...
%      - blkimpv(Forward, strikes(i), rate, TTM, CallPricesNMVMFFT(Forward, discount, moneyness(i), TTM, cal, Params, alpha)), [1:length(moneyness)]')), x00,LB,UB);

sigma  = cal(1);
eta    = cal(2);
k      = cal(3);

%% Try with lsqnonlin ( save time but no constraints )

% x0                = [0.25,1,5]';
% cal               = lsqnonlin(@(cal) arrayfun(@(i) abs(MktPrices(i)- CallPricesNMVMFFT(F0, discount, moneyness(i), TTM, cal(1), cal(3), cal(2), Params,alpha)), [1:length(MktPrices)]'), x0);
% sigma             = cal(1);
% eta               = cal(2);
% k                 = cal(3);

%% Compute Model Prices with calibrated parameters
CalPrices = real(CallPricesNMVMFFT(Forward, discount, moneyness, TTM, cal, Params, alpha)) - discount.*(Forward - strikes)';

%% Black Volatility obtain from the Calibrated Model Prices 
% CalPrices = CallPricesNMVMFFT(Forward, discount, moneyness, TTM, cal, Params, alpha);
mktvols   = blkimpv(Forward, strikes, rate, TTM, prices, 'Class', 'put');
CalVol    = blkimpv(Forward, strikes, rate, TTM, CalPrices', 'Class', 'put');

% calibration errors (on prices & volatilities)
PriceSquaredErrors = abs(prices - CalPrices').^2;
PercError          = abs(prices - CalPrices')./prices;
VolErrors          = abs(mktvols - CalVol);
TotVolError        = sum(VolErrors);
MAPE               = 100*mean(PercError);
MSE                = mean(PriceSquaredErrors);
moneyness          = log(strikes./Forward);


%% plot
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
  

end



