function [PutPrice, IC,S] = AutocallablePutPricing(OptionParams, UnderlyingParams, setDate, discountCurve, M, Nsim, flagspline)
%
% Function that prices a sprint autocallable put option simulating the
% underlying as an ATS process or as an LTS (depending on the parameters passed as input)
%
% INPUT
% OptionParams:        struct with data of the option (ResetDates, PaymentDates, strike, Forwards, S0)
% UnderlyingParams:    struct with parameters of the underlying
% setDate:             settlement date
% discountCurve:       struct with [dates, discounts]
% M:                   grid of the moneyness parameter (N=2^M) used in the FFT Lewis simulation
% Nsim:                number of simulations performed with Lewis FFT
% flagspline:          1 for spline 2 for linear interpolation in Lewis FFT
% 
% OUTPUT
% PutPrice:            price of the contract obtained
% IC:                  confidence interval of the price
% 
% CALLS
% ParamsInterpolation, XsimulationNIG, simulateATS


%% Day conventions
IBDaycount   = 3;
SwapDaycount = 2;


%% Option & underlying data
ResetDates   = OptionParams.ResetDates;
PaymentDates = OptionParams.PaymentDates;
strike       = OptionParams.strike;
Forwards     = OptionParams.Forwards;
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

dResetDates       = interp1(Dates, dividends, ResetDates);  % linear interpolation of dividends
S                 = zeros(Nsim, length(ResetDates)+1);
S(:,1)            = S0.*ones(length(S),1); 
F                 = zeros(Nsim, length(ResetDates)+1);
X                 = zeros(Nsim, length(ResetDates));
DiscPayoff        = zeros(length(PaymentDates), 1);

% interpolation of the model parameters
[sigmaResetDates, etaResetDates, kResetDates] = ParamsInterpolation(setDate, Dates, sigma, eta, k, ResetDates);

sigmaResetDates   = [1; sigmaResetDates]; % add one just for later purposes (not important)
etaResetDates     = [1; etaResetDates];
kResetDates       = [1; kResetDates];

% flag to decide how to simulate the underlying (LTS of ATS)
simulationFlag    = (all(sigmaResetDates(3:end) ~= sigmaResetDates(2)) || all(etaResetDates(3:end) ~= etaResetDates(2)) || all(kResetDates(3:end) ~= kResetDates(2)));


%% simulation of the underlying and pricing 
deltaResetDates  = [0; yearfrac(setDate, ResetDates, IBDaycount)];
prevDates        = [setDate; ResetDates(1:end-1)];
Payoff           = zeros(Nsim, length(ResetDates));

switch simulationFlag
    
    case 0 % if the model is a std LTS
        LTSparams.sigma = sigmaResetDates(2:end);
        LTSparams.eta   = etaResetDates(2:end);
        LTSparams.k     = kResetDates(2:end);
        
        % NIG simulation
        for i = 1:length(ResetDates)
            X(:,i)        = XsimulationNIG(ResetDates(i), LTSparams, Nsim, prevDates(i));
            F(:,i)        = S(:, i).*exp((FwdResetRates(i)-dResetDates(i)).*deltaResetDates(i+1)); %Ft_i-1_t_i
            S(:, i+1)     = F(:,i).*exp(X(:,i));
            Payoff(:, i)  = 1/S0*(max(strike - S(:,i+1), 0).*(strike < min(S(:,1:i), [], 2)).*1/S0.*max(S(:, 1:i), [], 2));
            DiscPayoff(i) = PaymentDiscounts(i)*mean(1/S0*(max(strike - S(:,i+1), 0).*(strike < min(S(:,1:i), [], 2)).*1/S0.*max(S(:, 1:i), [], 2)));
        end
        
    case 1 % if the model is an ATS

        % ATS simulation
        for i = 1:length(ResetDates)
            X(:,i)        = simulateATS(M, Nsim, flagspline, deltaResetDates(i+1), deltaResetDates(i), sigmaResetDates(i+1), sigmaResetDates(i)...
                                ,etaResetDates(i+1), etaResetDates(i), kResetDates(i+1), kResetDates(i), 0);            
            F(:,i)        = S(:, i).*exp((FwdResetRates(i)-dResetDates(i)).*(deltaResetDates(i+1))); %Ft_i-1_t_i
            S(:, i+1)     = F(:,i).*exp(X(:,i));
            Payoff(:, i)  = 1/S0*(max(strike - S(:,i+1), 0).*(strike < min(S(:,1:i), [], 2)).*1/S0.*max(S(:, 1:i), [], 2));
            DiscPayoff(i) = PaymentDiscounts(i)*mean(1/S0*(max(strike - S(:,i+1), 0).*(strike < min(S(:,1:i), [], 2)).*1/S0.*max(S(:, 1:i), [], 2)));
        end
        
end


% figure
% plot([setDate;ResetDates],S)
% gx=gca;
% gx.XTick=[setDate;ResetDates];
% gx.XTickLabelRotation = 30;
% datetick('x',1,'keepticks')
% grid on
% title('Simulation of undelying')

%% output
TotalPayoff = sum(Payoff, 2);
Error       = std(TotalPayoff)/sqrt(Nsim);
PutPrice    = sum(DiscPayoff);
IC          = [PutPrice - 1.96.*Error, PutPrice + 1.96.*Error];




end



 







