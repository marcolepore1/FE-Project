function t = dateMoveVec(startdate, datepart, num, businessdayconvention, market)
% Replaces DATEMOVE that works and only if num is a number
%
% Moves a given date forward or backward by a given number of time units
% according to a given business days convention on a given market
% 
% INPUT:
% 	- startdate               Start date [integer]
% 	- datepart                Time unit for the movement: 1='d', 2='w', 3='m', 4='y'
% 	- num                     Number of time units of the movement
% 	- businessdayconvention	  Business Day Convention: 'F', 'P', 'MF', 'MP', 'U'
% 	- market                  Market: target [vector of holidays]
% 
% Uses: busdate, isbusday, dateAdd
%
% NB:
%  - Either startdate or num can be arrays (not both)
%
% RondPoint 2012
% Last Modified: 27.01.2012 R. Baviera, A. Cassaro 

if (isempty(num))
    t=[];
    return
end

if ischar(startdate)
    startdate = datenum(startdate) ;
end    

% %the case @datepart = 'd' must be treated separately and no adjustment is made, i.e. only following business day convention
if strcmp(datepart,'d')
    if (length(startdate) == 1) % una data di partenza, piu' jump
        numberOfJumps = length(num);
        direc = sign(num);
        t = startdate*ones(size(direc)); 
        counter = abs(num);
        for i = 1:numberOfJumps
            while ( counter(i) > 0 )
                t(i) = busdate(t(i),direc(i),market);
                counter(i) = counter(i) - 1;
            end
        end
    else % un solo jump, piu' date di partenza
        numberOfStart = length(startdate);
        direc = sign(num);
        t = startdate; 
        counter = abs(num)*ones(numberOfStart);
        for i = 1:numberOfStart
            while ( counter(i) > 0 )
                t(i) = busdate(t(i),direc,market);
                counter(i) = counter(i) - 1;
            end
        end
        
    end
else
    t = dateAdd(startdate,num,datepart);

% if unadjusted, return
    if ( strcmp(businessdayconvention, 'U') )
        return
    end

%no adjustment if all t are business day
    if isbusday(t,market)
       return
    end

    hol = ~isbusday(t, market);
    thol = t(hol==1); % date holiday da aggiustare
    monthbeforeconv = month(thol);

    if ( strcmp(businessdayconvention,'F') || strcmp(businessdayconvention,'MF') )
        direc = 1;
    else
        direc = -1;
    end

%following or preceding business day convention
    thol = busdate(thol,direc,market);

    if ( strcmp(datepart,'m') || strcmp(datepart,'y') )
%modified following or preceding business day convention
        if ( strcmp(businessdayconvention,'MF') || strcmp(businessdayconvention,'MP') )
            for i = 1:length(thol)
                while ( month(thol(i)) ~= monthbeforeconv(i) )
                    thol(i) = busdate(thol(i),-direc,market);
                end
            end
        end
    end

    t(hol==1) = thol;

end

end % function dateMoveVec