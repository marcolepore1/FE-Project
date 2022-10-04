function [Dates, Discounts, Forwards] = bootstrapMarketDiscounts(optionTable)
%
% Function which computed the implied market discount factors at the
% maturities of the options provided in optionTable
%
% INPUT
% optionTable:              table with given options
% 
% OUTPUT
% dates:                    dates @ which discounts are computed
% discounts:                implied discounts computed
%

% Gask        = optionTable.ASK(1:2:end-1) - optionTable.ASK(2:2:end);
% Gbid        = optionTable.BID(1:2:end-1) - optionTable.BID(2:2:end);
Dates       = unique(optionTable.MATURITIES);
% nMaturities = length(dates);
nSynthFwd   = zeros(length(Dates),1);
Discounts   = zeros(length(Dates),1);
Forwards    = zeros(length(Dates),1);


for i=1:length(Dates)
    idx          = find(~(optionTable.MATURITIES - Dates(i)));
    nSynthFwd(i) = length(idx/2);
    Gbid         = optionTable.BID(idx(1):2:idx(end-1)) - optionTable.ASK(idx(2):2:idx(end));
    Gask         = optionTable.ASK(idx(1):2:idx(end-1)) - optionTable.BID(idx(2):2:idx(end));
    G            = 0.5*(Gbid + Gask);
    Ghat         = mean(G)';
    K            = unique(optionTable.STRIKES(idx));
    Khat         = mean(K)';
    diffK        = K - Khat;
    diffG        = G - Ghat;
    Discounts(i) = -(diffK'*diffG)./(diffK'*diffK);
    FakeForwards = G./Discounts(i) + K;
    plot(K, FakeForwards)
    axis([1000 4000 2500 3000])
    hold on
    Forwards(i)  = mean(FakeForwards)';
    
end










% try with vectorial code
% arrayfun(@(i) find(~(optionTable.MATURITIES - Maturities(i))), [1:length(Maturities)]', 'UniformOutput', false)

















end
