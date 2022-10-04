function d = funNIG(Forwards, discounts, OTMcalloptions, OTMputoptions, TTM, cal, Params, alpha)
%
% Function that computes the distance (Euclidean one) between NIG model
% prices and market ones. It is needed in the calibration of NIG
%
% INPUT
% Forwards:           Forward prices
% discounts:          discount factors
% OTMcalloptions:     vector with OTM call options prices for diff maturities
% OTMputoptions:      vector with OTM put options prices for diff maturities
% TTM:                yearfractions to maturities of all options
% cal:                parameters to be calibrated
% Params:             FFT parameters
% alpha:              LTS parameter
%
% OUTPUT
% d:                  Euclidean distance
% 
% CALLS
% CallPricesNMVMFFT
%


for i = 1:length(TTM)
    index1                = (OTMcalloptions(:, 6) == TTM(i));
    index2                = (OTMputoptions(:, 6) == TTM(i));
    MdlPricesCall(index1) = CallPricesNMVMFFT(Forwards(i), discounts(i), OTMcalloptions(index1, 2), TTM(i), cal, Params, alpha);
    MdlPricesPut(index2)  = CallPricesNMVMFFT(Forwards(i), discounts(i), OTMputoptions(index2, 2), TTM(i), cal, Params, alpha)'- discounts(i).*(Forwards(i) - OTMputoptions(index2, 3));
end
allPrices = [MdlPricesCall, MdlPricesPut]';
mktPrices = [OTMcalloptions(:,1); OTMputoptions(:,1)];
d         = sum(abs(allPrices - mktPrices).^2);

end
