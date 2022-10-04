function [Dates, Discounts, Forwards] = bootstrapMarketDiscountslm(optionTable)
%
% Function which computes: 
% i)  Implied market discount factors at the maturities of the options provided in optionTable 
%     with linear regression
% ii) Initial forward prices F(t0, t), where t stands for the maturity
%     of every option in optionTable
% Discounts & fwds are computed exploiting synthetic fwds technique making
% use of the put-call parity
%
% INPUT
% optionTable:  table with provided options data
% 
% OUTPUT
% Dates:        dates @ which discounts are computed
% Discounts:    implied discounts computed
% Forwards:     forward prices computed F(t0, t)
%
%

%% setting data
Dates       = unique(optionTable.MATURITIES);
nSynthFwd   = zeros(length(Dates),1);
Discounts   = zeros(length(Dates),1);
Forwards    = zeros(length(Dates),1);
CI          = zeros(length(Dates),2);
R2          = zeros(length(Dates),1);

%% fitting regression lines  exploiting synthetic forward technique 
for i=1:length(Dates)
    idx          = find(~(optionTable.MATURITIES - Dates(i)));          % indexes of options with ith maturity
    nSynthFwd(i) = length(idx/2);
    Gbid         = optionTable.BID(idx(1):2:idx(end-1)) - optionTable.ASK(idx(2):2:idx(end));
    Gask         = optionTable.ASK(idx(1):2:idx(end-1)) - optionTable.BID(idx(2):2:idx(end));
    G            = 0.5*(Gbid + Gask);
    K            = unique(optionTable.STRIKES(idx));

    % fitting linear model
    model        = fitlm(K,G);                                          
    CIresults    = coefCI(model);
    Discounts(i) = -model.Coefficients.Estimate(2,1);
    Forwards(i)  = model.Coefficients.Estimate(1,1)/Discounts(i);
    CI(i,:)      = CIresults(1,:)./Discounts(i);
    R2(i)        = model.Rsquared.Ordinary;
end

% plot of the discounts
figure()
plot(Dates, Discounts, '-*','LineWidth',2)
grid on
gx = gca;
gx.XTick = Dates;
gx.XTickLabelRotation = 30;
datetick('x', 1)
ylabel('Implied Discounts')
xlabel('Option Maturities')
title('Bootstrapped discounts')

% plot of the Forwards
figure()
plot(Dates, Forwards, '-*','LineWidth',2)
hold on
plot(Dates, CI, '-*','LineWidth',2)
grid on
gx = gca;
gx.XTick = Dates;
gx.XTickLabelRotation = 30;
datetick('x', 1)
ylabel('Forwards & CI')
xlabel('Option Maturities')
title('Forwards')

% display at output
fprintf('The minimum R2 of regressions is: %1.8f \n', min(R2, [], 1));

end