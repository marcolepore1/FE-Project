function fdate = dateAdd(startdate, num, datepart)
% Aggiunge ad una data un num di d, w, m, y 
%
% RondPoint 2012
% Last Modified: 27.01.2012 R. Baviera, A. Cassaro 

s = size(num);
if (s(1)==1)
    num = num';
    n = s(2);
else
    n = s(1);
end

if ischar(startdate)
    startdate = datenum(startdate);
end

startdate = repmat(startdate,n,1);

switch datepart
    case 'd'
        fdate = startdate + num;
    case 'w'
        fdate = startdate + 7*num;
    case 'm'
        t = datevec(startdate);
        t(:,1) = t(:,1)+dividi12(t(:,2)+num);
        t(:,2) = resto12(t(:,2)+num);
        t(:,3) = min(t(:,3),eomday(t(:,1),t(:,2)));

        fdate = datenum(t);
    case 'y'
        t = datevec(startdate);
        t(:,1)=t(:,1) + num;

        fdate = datenum(t);
    otherwise
end

end % function dateadd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = dividi12(x)
d = floor((x-1)/12) ;
end % function dividi12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = resto12(x)
r = mod(x-1,12)+1 ;
end % function resto12