function [flag] = powerLawConditions(beta, delta, alpha)
%
% Function that checks the conditions on the scaling parameters of the power-law
% for the existence of the ATS and displays a message at output
%
% INPUTS:
% beta:        scaling parameter
% delta:       scaling parameter
% alpha:       stability parameter            
% 
% OUTPUTS:
% flag:        true or false if conditions are satisfied or  not
%


%% check of the conditions
flag = false;

if 0 <= beta & beta <= 1/(1-alpha*0.5) & -min(beta, (1-beta*(1-alpha))/alpha) < delta & delta <= 0
    
    flag = true;
    disp('Power-law scaling parameters satisfy the conditions')

else 
    
    disp('Power-law scaling parameters do not satisfy the conditions')

end

end