function [sigma, eta, k, MSE, MAPE, mktSkew, calSkew] = CalibrateNMVMToVolatilitySurfaceConstCallPut(setDate, discount, Forward, strikes, nCall, prices, maturity, prevTTM, alpha, x00)
%
% Function that calibrates the parameters of an ATS NIG
% process at a certain maturity, making use of OTM calls & puts whose prices are
% stored in prices. 
% It performs also a plot of implied volatility smile obtained and computes the most common errors
% MSE & MAPE. 
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
% fmincon    for the constrained minimazion of the distance (NIG constraints and the one from previous parameters)
% CallPricesNMVMFFT
% 


%% Displaying datas 
IBdaycount  = 3;
TTM         = yearfrac(setDate, maturity, IBdaycount);
moneyness   = log(Forward./strikes);
rate        = -log(discount)/TTM;

% to prove with weights
% mktvegaCall = blsvega(S0, strikes(1:nCall), rate, TTM, mktvolsCall, dividends);
% mktvegaPut  = blsvega(S0, strikes(nCall+1:end), rate, TTM, mktvolsPut, dividends);
% mktvega     = [mktvegaCall; mktvegaPut];

 
%% Model Prices 
% Calibrating  parameters for the FFT
I_res         = 2*pi*exp(-sign(moneyness)*0.5.*moneyness);
M             = 15;
options       = optimset('TolFun',1e-5, 'Display', 'off');
x0            = 0.0025;
LB            = eps;
UB            = 0.01;
fTS           = @(v) 1./(v.^2 + 0.25);
ff            = @(x) 1./(x.^2 + 1/4);
% dz          = 0.0025; 
dz            = lsqnonlin(@(dz) abs( FourierTransform(fTS, moneyness,  M, dz)- I_res), x0, LB, UB, options); % find dz such that FFT replicates the residual integral
Params        = FFTparameters(M, dz, 1);
I             = computeIntegral(ff, moneyness, [], Params, 1);
MdlPricesCall = @(cal) CallPricesNMVMFFT(Forward, discount, moneyness(1:nCall), TTM, cal, Params, alpha);
MdlPricesPut  = @(cal) CallPricesNMVMFFT(Forward, discount, moneyness(nCall+1:end), TTM, cal, Params, alpha) - discount.*(Forward - strikes(nCall+1:end))'; 
MdlPrices     = @(cal, moneyness) [MdlPricesCall(cal), MdlPricesPut(cal)];
errorFFTINT   = abs(I - I_res);

%% Precision of the FFT calibration
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

%% constraints

% g1, g2 & g3 from previous parameters
g1 = 0.5 + x00(2) - sqrt((0.5 + x00(2))^2 + 2*(1-alpha)/(x00(1)^2*x00(3)));
g2 = -0.5 - x00(2) - sqrt((0.5 + x00(2))^2 + 2*(1-alpha)/(x00(1)^2*x00(3)));
g3 = (prevTTM^(1/alpha)*x00(1)^2)/(x00(3)^((1-alpha)/alpha))*sqrt((0.5 + x00(2))^2 + 2*(1-alpha)/(x00(1)^2*x00(3)));

% Matlab help : easy way to implement a non linear constraint
% constraints of NIG & of ATS non increasing g1, g2 & g3
function [c,ceq] = constraints (cal, g1, g2, g3, TTM)
  w   = 1/(2*cal(3).*cal(1).^2);  
  c   = [-w-cal(2),...
      -cal(1), ...
      g1 - (0.5 + cal(2) - sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3)))), ...
      g2 - (-0.5 - cal(2) - sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3)))), ...
      g3 - ((TTM^(1/0.5).*cal(1)^2)./(cal(3)^((1-0.5)/0.5))).*sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3))) ]';    % eta>-w & sigma>0

  ceq = [];
end

%% Calibration
LB      = [0.01; -5; 0.01];  % sigma, eta, k
% LB      = [eps; 10; 0.01];
% UB      = [1; 20; 5];
UB      = [0.5; 50; 5];
% UB      = [1; 50; 5];
options = optimoptions('fmincon', 'Display','off');

% penalty algorithm
% w           = @(cal) 1/(2*cal(3).*cal(1).^2);
% constraints = @(cal) (g1 > (0.5 + cal(2) - sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3)))))*1e16 + ...
%     (g2 > (-0.5 - cal(2) - sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3)))))*1e16 + ...
%     (g3 > ((TTM^(1/0.5).*cal(1)^2)./(cal(3)^((1-0.5)/0.5))).*sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3))))*1e16+...
%     (cal(2) < - w(cal))*1e16 + (cal(1)<0)*1e16;
% cal = lsqnonlin(aux, x00', LB, UB);

% Trying with weights
% aux = @(cal) sum((1./(mktvega).^2).*abs(prices-MdlPrices(cal)').^2);

% minimizaion of sum of squared distances
aux = @(cal) sum(abs(prices-MdlPrices(cal)').^2);
cal = fmincon(@(cal) aux(cal),x00,[],[],[],[],LB,UB, @(cal)constraints(cal, g1, g2, g3, TTM), options);

% calibrated parameters
sigma  = cal(1);
eta    = cal(2);
k      = cal(3);

%% Compute Model Prices with calibrated parameters
CalPricesCall = real(CallPricesNMVMFFT(Forward, discount, moneyness(1:nCall), TTM, cal, Params, alpha));
CalPricesPut  = real(CallPricesNMVMFFT(Forward, discount, moneyness(nCall+1:end), TTM, cal, Params, alpha)) - discount.*(Forward - strikes(nCall+1:end))';
CalPrices     = [CalPricesCall, CalPricesPut];

%% Black Volatility obtained from the Calibrated Model Prices 
mktvolsCall   = blkimpv(Forward, strikes(1:nCall), rate, TTM, prices(1:nCall));
CalVolCall    = blkimpv(Forward, strikes(1:nCall), rate, TTM, CalPricesCall');
mktvolsPut    = blkimpv(Forward, strikes(nCall+1:end), rate, TTM, prices(nCall+1:end), 'Class','put');
CalVolPut     = blkimpv(Forward, strikes(nCall+1:end), rate, TTM, CalPricesPut', 'Class', 'put');
mktvols       = [mktvolsCall; mktvolsPut];
CalVol        = [CalVolCall; CalVolPut];


%% Skew comparison
mktSkew = mktvolsCall(1) - mktvolsPut(end);
calSkew = CalVolCall(1) - CalVolPut(end);


%% errors (on prices & volatilities)
PriceSquaredErrors = abs(prices - CalPrices').^2;
VolErrors          = abs(mktvols - CalVol);
TotVolError        = sum(VolErrors);
PercErrors         = abs(prices - CalPrices')./prices;
MSE                = mean(PriceSquaredErrors);
MAPE               = 100*mean(PercErrors);


%% plot
moneyness  = log(strikes./Forward);
figure()
plot(moneyness, CalVol,'+', 'MarkerSize',5)
hold on 
plot(moneyness, mktvols, 'square', 'MarkerSize',5)
grid on
t = text(min(moneyness), min(mktvols), ['\bf MSEPrices : ', num2str(MSE)], 'Color', 'k');  % MSE
t.FontSize = 10;
xlabel('moneyness')
ylabel('volatilities')
legend('Calvols', 'mktvols')
title('Smile ATS @', num2str(datestr(maturity)))
% text(strikes - 0.002, min(CalVol, mktvols) , num2str(VolErrors), 'Color', 'k', 'FontSize', 6); 
% tVol           = text(strikes(end)-0.001, mktvols(1), ['\bf Total Error : ', num2str(TotVolError)], 'Color', 'r');  
% tVol.FontSize  = 5;
%

end

