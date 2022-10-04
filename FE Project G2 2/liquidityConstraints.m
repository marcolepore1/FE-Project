function optionTableCleaned = liquidityConstraints(optionTable)
%
% Function that applies the two liquidity constraints: 
% i)  removes options with price less than 10% the mini-
%     mum difference in the grid of market strikes;
% ii) removes options with bid-ask spread over bid bigger than 60%.
%
% INPUT
% optionTable:         Table which contains all options provided
%                     
% OUTPUT
% optionTableCleaned:  Table with call and put that satisfy liq constraints
% 
%

%% first constraint
strikesDiff = abs(optionTable.STRIKES(2:end)-optionTable.STRIKES(1:end-1));
idx         = find(strikesDiff);
strikesDiff = strikesDiff(idx);
minDiff     = min(strikesDiff);
idx         = find(0.5*(optionTable.ASK + optionTable.BID) >= 0.1*minDiff);
optionTable = optionTable(idx, :);

%% second constraint
bidAskSpread = abs(optionTable.ASK - optionTable.BID);
idx          = find(bidAskSpread./optionTable.BID <= 0.6);
optionTable  = optionTable(idx, :);

% removing if not both C & P are present
idx = zeros(length(optionTable.STRIKES),1);
for i=2:length(optionTable.STRIKES)-1
    if (optionTable.STRIKES(i) ~= optionTable.STRIKES(i+1)) && (optionTable.STRIKES(i) ~= optionTable.STRIKES(i-1))
        idx(i) = i;
    end
end
idx2        = find(~idx);
optionTable = optionTable(idx2, :);

% the first one & last one have to be removed
optionTableCleaned  = optionTable(2:end-1, :);

end
