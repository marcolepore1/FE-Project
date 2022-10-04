function zRates = zeroRates(setDate, dates, discounts)
% computes the zero rates from the discount factors for a given set of
% dates
% 
% INPUT
% dates:       dates where the discounts act (between the setDate and the date)
% discounts:   given discount factors for those dates
%
% OUTPUT
% zRates:      zero rates @dates in percentage
% 

%% data
IBDayCount = 3;

%% zero_rates
zRates = 100*(-log(discounts)./yearfrac(setDate, dates, IBDayCount));

end