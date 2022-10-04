function [X] = simulateATSSUM(M, Nsim, flagspline, t, s, sigma_t, sigma_s, eta_t, eta_s, k_t, k_s)
%
% Function that simulates an increment of the ATS process via CDF
% inversion, using the approximation of the integral with the sum 
% PAY ATTENTION! Really slow since the CDF is evaluated point by point of
% the grid
%
% INPUT 
% M:          discretiazion param for grid
% Nsim:       number of simulations 
% flagspline: 1 for spline 2 for linear
% t:          second date                   % s: first date
% sigma_t:    volatility @second date       % sigma_s: volatility @first date
% eta_t:      skew @second date             % eta_s: skew @first date
% k_t:        vol-of-vol @second date       % k_s: vol-of-vol @first date
%
% OUTPUT
% X:          Nsim increments simulated f(s,t)
%
% CALLS
% 


%% Assumption 2 parameters
alpha       = 0.5;
IBDaycount  = 3;
b           = (1-alpha)^(1-alpha)./(2^(alpha)*alpha).*(t./k_t.^(1-alpha).*sigma_t.^(2*alpha) - s./k_s.^(1-alpha).*sigma_s.^(2*alpha)) - 1e-4;
omega       = 2*alpha-1e-4;
B           = exp((1-alpha)/alpha.*(t./k_t - s./k_s) + 1);
g2          = -0.5 - eta_t - sqrt((0.5 + eta_t).^2 + 2*(1-alpha)./(sigma_t.^2.*k_t));
pplus       = -g2 - 1;  
a           = 0.5*(pplus +1);
limit       = @(u) B.*exp(-b.*u.^(omega));

% characteristic function of the increment of the process (between s & t)
LaplaceExp_t   = @(x) (t./k_t) .* ((1-alpha)/alpha) .* (1 - (1 + (x.*k_t.*sigma_t.^2)./(1-alpha)).^alpha);
LaplaceExp_s   = @(x) (s./k_s) .* ((1-alpha)/alpha) .* (1 - (1 + (x.*k_s.*sigma_s.^2)./(1-alpha)).^alpha);
phi_s_t        = @(x) exp(-1i.*x.*(LaplaceExp_t(eta_t) - LaplaceExp_s(eta_s))).*exp(LaplaceExp_t((x.^2 + 1i.*(1+2.*eta_t).*x)./2)-LaplaceExp_s((x.^2 + 1i.*(1+2.*eta_s).*x)./2));

% plot to check the assumption 2 is satisfied with chosen parameters
% uu             = linspace(0, 100, 1000);
% figure
% plot(uu, abs(phi_s_t(uu - 1i*a)))
% hold on 
% plot(uu, limit(uu))
% legend('phi' , 'Bexp(-bu^w)')

%% FFT parameters
N           = 2^M;
h           = (2.*pi.*a./(b*N.^(omega))).^(1/(omega+1));
delta       = t-s;
D           = 5;
% D           = 1.7;
gamma       = (2*pi)/(N*h);
X0          = -gamma*(N-1)*0.5;
Xk          = -X0;
z1          = -0.5*h*(N-1);
zn          = -z1;
xx          = X0:gamma:Xk;


% trying with the formula in the paper
Phatsum    = @(x) 1-(exp(-a.*x)/pi).*h.*sum(real(exp(-1i.*([0:N-1]+0.5).*h.*x).*phi_s_t(([0:N-1]+0.5).*h-1i.*a)./(1i.*([0:N-1]+0.5).*h+a)));
Phatsumres = arrayfun(@(i) Phatsum(xx(i)), [1:length(xx)]');


%% CDF on truncated domain
point       = D*sqrt(delta);
idx         = find(xx >= point);

% if last point of the grid is less than D*sqrt(delta) keep the last point
if isempty(idx)
    idx = length(xx);
end

Xk = xx(idx(1)).*(abs(xx(idx(1)) - point) < abs(xx(idx(1)-1) - point))+ xx(idx(1)-1).*(abs(xx(idx(1)) - point) > abs(xx(idx(1)-1) - point));
X0 = -Xk;
xx = X0:gamma:Xk;

% trying with the formula in the paper
Phatsum    = @(x) 1-(exp(-a.*x)/pi).*h.*sum(real(exp(-1i.*([0:N-1]'+0.5).*h.*x).*phi_s_t(([0:N-1]'+0.5).*h-1i.*a)./(1i.*([0:N-1]'+0.5).*h+a)));
Phatsumres = arrayfun(@(i) Phatsum(xx(i)), [1:length(xx)]');

% put to 0 values before the CDF is increasing
ind1 = find(Phatsumres>1);
Phatsumres(ind1) = 1;
diff = Phatsumres(2:end) - Phatsumres(1:end-1);
for i = 1:length(Phatsumres)-1
    Phatsumres(i) = Phatsumres(i).*(all(diff(i:end)>=0)).*(Phatsumres(i)>=0);
end


%% simulating
% X1 = min(Phat(find(Phat)));
% U  = X1 + (Phat(end) - X1).*rand(Nsim, 1);
U  = rand(Nsim, 1);

switch flagspline
    
    case 1 % spline
        
        [Phatunique, idxunique, ~] = unique(Phatsumres, 'stable'); 
        X   = interp1(Phatunique, xx(idxunique), U, 'spline');
        idx = find(X>xx(end));
        
        if ~isempty(idx)
                newX   = interp1(Phatunique, xx(idxunique), U(idx), 'makima');
                X(idx) = newX;
        end

    case 2 % linear
        
        [Phatunique, idxunique, ~] =unique(Phatsumres, 'stable'); 
        X  = interp1(Phatunique, xx(idxunique), U);
end

%% PDF
% density = @(x) a.*(exp(-a.*x)/pi).*h.*sum(real(exp(-1i.*([0:N-1]+0.5).*h.*x).*phi_s_t(([0:N-1]+0.5).*h-1i.*a)./(1i.*([0:N-1]+0.5).*h+a)))+...
%                        (exp(-a.*x)/pi).*h.*sum(real(1i*([0:N-1]+0.5).*h.*exp(-1i.*([0:N-1]+0.5).*h.*x).*phi_s_t(([0:N-1]+0.5).*h-1i.*a)./(1i.*([0:N-1]+0.5).*h+a)));
% 
% densityres = arrayfun(@(i) density(xx(i)), [1:length(xx)]');
% idx        = find(Phatsumres);
% densityres(1:idx(1)-1) = 0;
% 
% figure
% plot(xx, densityres, 'LineWidth', 2)
% title('PDF')

end



