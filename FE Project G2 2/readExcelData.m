function [OptionStruct] = readExcelData( filename, formatData)
%
% Reads data from excel
%  
%
% INPUT
% filename: excel file name where data are stored
% formatData: data format in Excel
% 
% OUTPUT
% OptionStruct: struct with settlement, expiries, ask price, bid price,
%               strike price, open interest, type of option
%

%% Dates from Excel

%Settlement date
[~, settlement] = xlsread(filename, 1, 'B3');
%Date conversion
OptionStruct.settlement = datenum(settlement, formatData);

ask_price         = xlsread(filename, 1, 'C3:C2053');
bid_price         = xlsread(filename, 1, 'D3:D2053');
strike_price      = xlsread(filename, 1, 'E3:E2053');
open_price        = xlsread(filename, 1, 'J3:J2053');
excercise_price   = xlsread(filename, 1, 'G3:G2053');
[~, maturities]   = xlsread(filename, 1, 'L3:L2053');
[~, type_option]  = xlsread(filename, 1, 'M3:M2053');

OptionStruct.maturities   = datenum(maturities, formatData);
OptionStruct.askPrices    = ask_price;
OptionStruct.bidPrices    = bid_price;
OptionStruct.strikes      = strike_price;
OptionStruct.opening      = open_price;
OptionStruct.flag         = type_option;


end % readExcelData