function [sigma, eta, k, MSE, MAPE] = CalibrateNMVMToVolatilitySurfaceConstCall(setDate, discount, Forward, strikes, prices, maturity, prevTTM, alpha, x00)
%
% Function that calibrates the parameters of an ATS NIG
% process at a certain maturity, making use of OTM calls only whose prices are
% stored in prices. 
% It performs also a plot of implied volatility smile obtained and computes the most common errors
% MSE & MAPE. 
%  
% INPUT
% setDate:   settlement date considered
% discount:  discount factor @maturity
% Forward:   forward price with expiry @maturity (F(t0,t))
% strikes:   strikes of the options
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
MdlPrices = @(cal) CallPricesNMVMFFT(Forward, discount, moneyness, TTM, cal, Params, alpha);
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

% g1, g2 & g3 with previous params
g1 = 0.5 + x00(2) - sqrt((0.5 + x00(2))^2 + 2*(1-alpha)/(x00(1)^2*x00(3)));
g2 = -0.5 - x00(2) - sqrt((0.5 + x00(2))^2 + 2*(1-alpha)/(x00(1)^2*x00(3)));
g3 = (prevTTM^(1/alpha)*x00(1)^2)/(x00(3)^((1-alpha)/alpha))*sqrt((0.5 + x00(2))^2 + 2*(1-alpha)/(x00(1)^2*x00(3)));

% constraints of NIG e non-increasing g1, g2 & g3
% Matlab help : easy way to implement a non linear constraint 
function [c,ceq] = constraints (cal, g1, g2, g3, TTM)
  w   = 1/(2*cal(3).*cal(1).^2);  
  c   = [-w-cal(2),...
      -cal(1), ...
      g1 - (0.5 + cal(2) - sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3)))), ...
      g2 - (-0.5 - cal(2) - sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3)))), ...
      g3 - ((TTM^(1/0.5).*cal(1)^2)./(cal(3)^((1-0.5)/0.5))).*sqrt((0.5 + cal(2)).^2 + 2*(1-0.5)./(cal(1).^2.*cal(3))) ]';    % eta>-w & sigma>0

  ceq = [];
end

%% Calibration & Results
% constraints of the calibration
LB  = [0.001; -5; 0.001];
UB  = [1; 20; 5];
aux = @(cal) sum(abs(prices-MdlPrices(cal)').^2); % function that we want to minimize the weigths are negligible since they are equal to 1
% aux             = @(cal) sqrt(sum(abs(arrayfun(@(i) blkimpv(Forward, strikes(i), rate, TTM, prices(i)) ...
%     -blkimpv(Forward, strikes(i), rate, TTM, MdlPrices(cal))').^2));

% opts     = optimoptions(@fmincon,'Algorithm','sqp');
% problem  = createOptimProblem('fmincon','objective',...
%     aux,'x0', x00,'lb', LB,'ub', UB, 'nonlcon', @(cal) const(g1, g2, g3), 'options',opts);
% ms       = MultiStart;
% [cal, f] = run(ms,problem,20);

cal = fmincon(@(cal) aux(cal),x00,[],[],[],[],LB,UB, @(cal)constraints(cal, g1, g2, g3, TTM));

% calibrated parameters
sigma  = cal(1);
eta    = cal(2);
k      = cal(3);


%% Compute Model Prices with calibrated parameters
CalPrices = real(CallPricesNMVMFFT(Forward, discount, moneyness, TTM, cal, Params, alpha));

%% Black Volatility obtain from the Calibrated Model Prices 
mktvols   = blkimpv(Forward, strikes, rate, TTM, prices);
CalVol    = blkimpv(Forward, strikes, rate, TTM, CalPrices');

% calibration errors (on prices & volatilities)
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
% text(strikes - 0.002, min(CalVol, mktvols) , num2str(VolErrors), 'Color', 'k', 'FontSize', 6); 
% grid on
% tVol           = text(strikes(end)-0.001, mktvols(1), ['\bf Total Error : ', num2str(TotVolError)], 'Color', 'r');  
% tVol.FontSize  = 5;
% legend('Calvols', 'mktvols')
%  

end

