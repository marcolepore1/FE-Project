function [flag1, flag2] = ATSassumptions(k, eta, sigma, TTM, alpha)
% 
% Function that checks Assumptions 1 and 2 for the ATS:
% 1- pplus and pminus non increasing in t
% 2- |phi_s_t(u-ia)| < Bexp{-bu^w} for sufficiently large u
% 
% INPUTS:
% k:                 vol-of-vol
% eta:               skew
% sigma:             volatility
% TTM:               time to maturity
% alpha:             parameter of LTS/ATS (0.5 -> NIG, 0 -> VG)
%
% OUTPUTS:
% flag1:             true or false, first assumption
% flag2:             true or false, second assumption
%
%


%% Assumption 1
g1         = 0.5 + eta - sqrt((0.5 + eta).^2 + 2*(1-alpha)./(sigma.^2.*k));
g2         = -0.5 - eta - sqrt((0.5 + eta).^2 + 2*(1-alpha)./(sigma.^2.*k));
pplus      = -g2 - 1;  
pminus     = -g1;

figure()
plot(TTM, pplus, '-*', 'LineWidth', 3);
hold on 
plot(TTM, pminus, '-+', 'LineWidth', 3);
grid on
title('Non-increasing pplus & pminus')
legend('pplus', 'pminus')

% pplus and pminus non increasing in t
flag1 = (pplus(1:end-1) >= pplus(2:end)).*(pminus(1:end-1) >= pminus(2:end));

if all(flag1)
    flag1 = true;
    disp('First Assumption is satisfied');
else
    flag1 = false;
    disp('First Assumption is not satisfied');
end

%% Assumption 2 

a              = 0.5*(pplus +1);
uu             = linspace(0, 10000, 1000);
flag2          = true;

for i = 1:length(a)-1
    
    LaplaceExp_t   = @(x) (TTM(i+1)./k(i+1)) .* ((1-alpha)/alpha) .* (1 - (1 + (x.*k(i+1).*sigma(i+1).^2)./(1-alpha)).^alpha);
    LaplaceExp_s   = @(x) (TTM(i)./k(i)) .* ((1-alpha)/alpha) .* (1 - (1 + (x.*k(i).*sigma(i).^2)./(1-alpha)).^alpha);
    
    % phi_s_t = phi_t / phi_s
    phi_s_t        = @(x) exp(-1i.*x.*(LaplaceExp_t(eta(i+1)) - LaplaceExp_s(eta(i)))).*exp(LaplaceExp_t((x.^2 + 1i.*(1+2.*eta(i+1)).*x)./2)-LaplaceExp_s((x.^2 + 1i.*(1+2.*eta(i)).*x)./2));
    
    % parameters B, b, omega
    b              = (1-alpha)^(1-alpha)./(2^(alpha)*alpha).*(TTM(i+1)./k(i+1).^(1-alpha).*sigma(i+1).^(2*alpha) - TTM(i)./k(i).^(1-alpha).*sigma(i).^(2*alpha)) - 0.0001;
    omega          = 2*alpha - 0.0001;
    B              = exp((1-alpha)/alpha.*(TTM(i+1)./k(i+1) - TTM(i)./k(i)) + 1);
    
    % |phi_s_t(u-ia)| < Bexp{-bu^w} for sufficiently large u
    limit          = B.*exp(-b.*uu.^omega);
    yy             = abs(phi_s_t(uu - 1i.*a(i+1)));
    flag2          = flag2*(yy(end) <= limit(end));

%     figure()
%     plot(uu, yy, 'LineWidth',2)
%     hold on
%     plot(uu, limit, 'LineWidth',2)
%     legend('phi', 'limit')
end

% output
if flag2
    disp('Second Assumption is satisfied')
else
    disp('Second Assumption is not satisfied')
end





end
