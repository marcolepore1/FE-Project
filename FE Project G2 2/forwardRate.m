function fwdR = forwardRate(dateT, dateS, discountT, discountS)
% 
% Function that computes the forward rate between two dates r(s, t)
%
% INPUT
% dateT:        second date (t)
% dateS:        first date (s)
% discountT:    discount @second date
% discountS:    discounts @first date
%
% OUTPUT
% fwdR:         forward rate r_s_t between the provided dates in input;
%

% Daycount
IBDaycount  = 3;

% Computation of rate & discount
FwdDiscount = discountT /discountS;
fwdR        =  -log(FwdDiscount)/yearfrac(dateS, dateT, IBDaycount);


end