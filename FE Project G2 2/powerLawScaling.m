function [etabar, kbar, beta, delta, theta, sigmabar, gamma] = powerLawScaling(sigma, eta, k, setDate, dates)
%
% Function that fits a Power-law ATS in order to estimate the parameters of
% the ATS
%
% INPUTS:
% sigma:             volatility
% eta:               skew
% k:                 vol-of-vol
% setDate:           settlement date
% dates:             maturities
%
% OUTPUTS:
% etabar             parameter
% kbar               parameter
% beta               scaling parameters
% delta              scaling parameters
% theta              scaling time
%
% CALLS
% fitlm 


%% Data
IBDaycount = 3;
TTM        = yearfrac(setDate, dates, IBDaycount);
% theta    = TTM.*sigma.^2;   % for time change
khat       = k;
% khat     = k.*sigma.^2;     % for time change
etahat     = eta;
sigmahat   = sigma;

%% Fitting through linear regression

% With time change k
% model1          = fitlm(log(theta), log(khat));
% kbar            = exp(model1.Coefficients.Estimate(1));
% beta            = model1.Coefficients.Estimate(2);
% regressionline1 = model1.Coefficients.Estimate(1) + beta*log(theta);
% figure()
% plot(log(theta), log(khat), '*', 'LineWidth', 3)
% hold on
% plot(log(theta), regressionline1, '-r', 'LineWidth', 2)
% xlabel('ln(θ)')
% ylabel('ln(k_{θ})')
% legend('Calibrated', 'Fitted')

% with 'normal' time k model
model1          = fitlm(log(TTM), log(khat));
kbar            = exp(model1.Coefficients.Estimate(1));
beta            = model1.Coefficients.Estimate(2);
regressionline1 = model1.Coefficients.Estimate(1) + beta*log(TTM);

% plot of the regression line
figure()
plot(log(TTM), log(khat), '*', 'LineWidth', 3)
hold on
plot(log(TTM), regressionline1, '-r', 'LineWidth', 2)
xlabel('ln(t)')
ylabel('ln(k_{t})')
grid on
legend('Calibrated', 'Fitted')


% With time change eta
% model2 = fitlm(log(theta), log(etahat));
% etabar = exp(model2.Coefficients.Estimate(1));
% delta  = model2.Coefficients.Estimate(2);
% regressionline2 = model2.Coefficients.Estimate(1) + delta*log(theta);
% figure()
% plot(log(theta), log(etahat), '*', 'LineWidth', 3)
% hold on
% plot(log(theta), regressionline2, '-r', 'LineWidth', 2)
% xlabel('ln(θ)')
% ylabel('ln(η_{θ})')
% legend('Calibrated', 'Fitted')

% with 'normal' time eta model
model2          = fitlm(log(TTM), log(etahat));
etabar          = exp(model2.Coefficients.Estimate(1));
delta           = model2.Coefficients.Estimate(2);
regressionline2 = model2.Coefficients.Estimate(1) + delta*log(TTM);

% plot of the regression line for eta
figure()
plot(log(TTM), log(etahat), '*', 'LineWidth', 3)
hold on
plot(log(TTM), regressionline2, '-r', 'LineWidth', 2)
xlabel('ln(t)')
ylabel('ln(η_{t})')
grid on
legend('Calibrated', 'Fitted')

% model for sigma (just to find a sigmabar consistent for later purposes)
model3   = fitlm(log(TTM), log(sigmahat));
sigmabar = exp(model3.Coefficients.Estimate(1));
gamma  = model3.Coefficients.Estimate(2);
regressionline3 = model3.Coefficients.Estimate(1) + gamma*log(TTM);
theta = TTM;

% figure()
% plot(log(TTM), log(sigmahat), '*', 'LineWidth', 3)
% hold on
% plot(log(TTM), regressionline3, '-r', 'LineWidth', 2)
% xlabel('ln(t)')
% ylabel('ln(sigma_{t})')
% grid on
% legend('Calibrated', 'Fitted')




end