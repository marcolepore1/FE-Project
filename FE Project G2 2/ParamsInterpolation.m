function [sigmaResetDates, etaResetDates, kResetDates] = ParamsInterpolation(setDate, Dates, sigma, eta, k, ResetDates)
%
% Function which performs the interpolation of ATS parameters @ResetDates
% Confronts also the result of linear & spline interpolation with a plot
%
% INPUT
% setDate:         Settlement date 
% Dates:           Dates @ which parameters of ATS are calibrated
% sigma:           Volatility of ATS
% eta:             Skew of ATS
% k:               Vol-of-vol of ATS
% ResetDates:      Resettlements dates at which interpolation is performed
%
% OUTPUT
% sigmaResetDates: sigma interpolated @ reset dates
% etaResetDates:   eta interpolated @ reset dates
% kResetDates:     k interpolated @ reset dates
%

%% Yearfractions
IBDaycount      = 3;
TTMDates        = yearfrac(setDate, Dates, IBDaycount);
TTMResetDates   = yearfrac(setDate, ResetDates, IBDaycount);

%% interpolation comparison
sigmaResetDatesL  = interp1(TTMDates, sigma, TTMResetDates, 'linear');  %%NB: Caso parametri LTS non interpola , passo i valori cost in ingresso 
etaResetDatesL    = interp1(TTMDates, eta, TTMResetDates, 'linear');
kResetDatesL      = interp1(TTMDates, k, TTMResetDates, 'linear');

sigmaResetDatesS  = interp1(TTMDates, sigma, TTMResetDates, 'spline');
etaResetDatesS    = interp1(TTMDates, eta, TTMResetDates, 'spline');
kResetDatesS      = interp1(TTMDates, k, TTMResetDates, 'spline');

% pp                = polyfit(TTMDates, eta, 8);
% etaResetDatesP    = polyval(pp, TTMResetDates);

% comparison
sigmaDiff = sum(abs(sigmaResetDatesL - sigmaResetDatesS));
etaDiff   = sum(abs(etaResetDatesL - etaResetDatesS));
kDiff     = sum(abs(kResetDatesL - kResetDatesS));
% add the differences

%% plot of the behavior (only if initial params are non-constant)
if (eta(2) ~= eta(1))
    figure();
    plot(Dates, eta, 'LineStyle','-', 'Marker', 'square','Color','k', 'LineWidth',2)
    hold on
    plot(ResetDates, etaResetDatesL, 'LineStyle','-', 'Marker', 'square', 'Color','b', 'LineWidth',2)
    hold on
    plot(ResetDates, etaResetDatesS, 'LineStyle','-', 'Marker', 'diamond', 'Color','g', 'LineWidth',2)
    grid on 
    gx = gca;
    gx.XTick = ResetDates;
    gx.XTickLabelRotation = 30;
    datetick('x', 1)
    title('eta')
    legend('eta', 'LinearInterp', 'SplineInterp')
    
    figure()
    plot(Dates, sigma, 'LineStyle','-', 'Marker', 'square', 'Color','k', 'LineWidth',2)
    hold on
    plot(ResetDates, sigmaResetDatesL, 'LineStyle','-', 'Marker', 'square', 'Color','b', 'LineWidth',2)
    hold on
    plot(ResetDates, sigmaResetDatesS, 'LineStyle','-', 'Marker', 'diamond', 'Color','g', 'LineWidth',2)
    grid on
    gx = gca;
    gx.XTick = ResetDates;
    gx.XTickLabelRotation = 30;
    datetick('x', 1)
    title('sigma')
    legend('sigma', 'LinearInterp', 'SplineInterp')
    
    figure()
    plot(Dates, k, 'LineStyle','-', 'Marker', 'square', 'Color','k', 'LineWidth',2)
    hold on
    plot(ResetDates, kResetDatesL, 'LineStyle','-', 'Marker', 'square', 'Color','b', 'LineWidth',2)
    hold on 
    plot(ResetDates, kResetDatesS, 'LineStyle','-', 'Marker', 'diamond', 'Color','g', 'LineWidth',2)
    grid on
    gx = gca;
    gx.XTick = ResetDates;
    gx.XTickLabelRotation = 30;
    datetick('x', 1)
    legend('k', 'LinearInterp', 'SplineInterp')
    title('k')
end

%% output
sigmaResetDates = sigmaResetDatesL;
etaResetDates   = etaResetDatesS;
kResetDates     = kResetDatesL;


end