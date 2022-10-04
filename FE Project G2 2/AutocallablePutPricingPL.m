function [PutPrice, IC] = AutocallablePutPricingPL(OptionParams, UnderlyingParams, setDate, discountCurve, M, Nsim, flagspline)
%
% Function that prices a sprint autocallable put option simulating the
% underlying as an ATS process
%
% INPUT
%
%

%% Day conventions
IBDaycount   = 3;

%% Option & underlying data
ResetDates   = OptionParams.ResetDates;
PaymentDates = OptionParams.PaymentDates;
strike       = OptionParams.strike;
S0           = OptionParams.S0;

sigma        = UnderlyingParams.sigma;
eta          = UnderlyingParams.eta;
k            = UnderlyingParams.k;

Discounts    = discountCurve.discounts;
Dates        = discountCurve.dates;
dividends    = UnderlyingParams.dividends;

%% Interpolating parameters & rates
zRates            = 0.01*zeroRates(setDate, Dates, Discounts);
ResetRates        = interp1(Dates, zRates, ResetDates);
ResetDiscounts    = exp(-ResetRates.*yearfrac(setDate, ResetDates, IBDaycount));
FwdResetDiscounts = ResetDiscounts(2:end)./ResetDiscounts(1:end-1);
FwdResetRates     = [ResetRates(1); -log(FwdResetDiscounts)./yearfrac(ResetDates(1:end-1), ResetDates(2:end), IBDaycount)];

PaymentRates      = interp1(Dates, zRates, PaymentDates);
PaymentDiscounts  = exp(-PaymentRates.*yearfrac(setDate, PaymentDates, IBDaycount));

sigmaResetDates   = [0.5; sigma];
etaResetDates     = [0.5; eta];
kResetDates       = [0.5; k];


dResetDates       = interp1(Dates, dividends, ResetDates);
S                 = zeros(Nsim, length(ResetDates)+1);
S(:,1)            = S0.*ones(length(S),1); 
F                 = zeros(Nsim, length(ResetDates)+1);
X                 = zeros(Nsim, length(ResetDates));
DiscPayoff        = zeros(length(PaymentDates), 1);

%% simulation of the underlying and pricing 
deltaResetDates  = [0; yearfrac(setDate, ResetDates, IBDaycount)];
Payoff           = zeros(Nsim, length(ResetDates));




for i = 1:length(ResetDates)
    X(:,i)       = simulateATS(M, Nsim, flagspline, deltaResetDates(i+1), deltaResetDates(i), sigmaResetDates(i+1), sigmaResetDates(i)...
                        ,etaResetDates(i+1), etaResetDates(i), kResetDates(i+1), kResetDates(i), 0);            
    F(:,i)        = S(:, i).*exp((FwdResetRates(i)-dResetDates(i)).*(deltaResetDates(i+1))); %Ft_i-1_t_i
    S(:, i+1)     = F(:,i).*exp(X(:,i));
    Payoff(:, i)  = 1/S0*(max(strike - S(:,i+1), 0).*(strike < min(S(:,1:i), [], 2)).*1/S0.*max(S(:, 1:i), [], 2));
    DiscPayoff(i) = PaymentDiscounts(i)*mean(1/S0*(max(strike - S(:,i+1), 0).*(strike < min(S(:,1:i), [], 2)).*1/S0.*max(S(:, 1:i), [], 2)));
end


% output
TotalPayoff = sum(Payoff, 2);
Error       = std(TotalPayoff)/sqrt(Nsim);
PutPrice    = sum(DiscPayoff);
IC          = [PutPrice - 1.96.*Error, PutPrice + 1.96.*Error];
end

