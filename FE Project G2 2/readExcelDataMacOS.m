function [optionData] = readExcelDataMacOS(filename)
% 
% Reads data from excel (MacOS version)
% It reads data from excel
%
% INPUT
% filename: excel file name where data are stored
% 
% OUTPUTS:
% OptionStruct: struct with settlement, expiries, ask price, bid price,
%               strike price, open interest, type of option
%

%% Dates from Excel

% Define cell-to-array function
cell2arr = @(c) reshape(c, size(c));

% Import all content
allcontent = readcell(filename);

% Settlement date
settlement            = allcontent{3,1};            % 07-06-2019 
optionData.settlement = datenum(settlement);

% Maturity dates of the options
maturities            = cell2arr(allcontent(3:end, 12));
optionData.maturities = datenum(maturities);

% Bid-Ask-Opening prices
bidPrices            = cell2arr(allcontent(3:end, 4));
askPrices            = cell2arr(allcontent(3:end, 3));
openingPrices        = cell2arr(allcontent(3:end, 10));
highPrices           = cell2arr(allcontent(3:end, 8));
optionData.bidPrices = cell2mat(bidPrices);
optionData.askPrices = cell2mat(askPrices);
optionData.opening   = cell2mat(openingPrices);
optionData.high      = cell2mat(highPrices);


% Call or put
flag                 = cell2arr(allcontent(3:end, 13));
optionData.flag      = flag;

% Strike
strikes              = cell2arr(allcontent(3:end, 5));
optionData.strikes   = cell2mat(strikes);


end % readExcelDataOSX