function [X] = simulateATSFFT(M, Nsim, flagspline, t, s,sigma_t,sigma_s,eta_t,eta_s,k_t,k_s)
%
%
% Simulation of one forward increment of NIG ATS through CDF inversion via
% FFT. Provide as input parameters @corresponding time instants of the ATS
%
% INPUT
% M:          parameter that monitors the discretization of FFT grid
% Nsim:       number of simulations of U (for the inversion)
% flagspline: 1 for spline, 2 for linear interpolation (in the inversion)
% t:          yearfrac between settlement & second time instant
% s:          yearfrac between settlement & first time instant      
% sigma_t:    volatility of NIG ATS @ first time instant
% sigma_s:    volatility of NIG ATS @ second time instant
% eta_t:      skew of NIG ATS @ second time instant
% eta_s:      skew of NIG ATS @ second time instant
% k_t:        vol-of-vol of NIG ATS @ second time instant
% k_s:        vol-of-vol of NIG ATS @ first time instant
%
% OUTPUT
% X:          simulations of values st P(X)=U (U uniform rv simulated)
%
% NB: uncomment last rows to obtain a plot of the PDF, but this takes much
% time (No FFT to compute the PDF)
%
%

%% Assumption 2 parameters
alpha       = 0.5;
IBDaycount  = 3;
b           = (1-alpha)^(1-alpha)./(2^(alpha)*alpha).*(t./k_t.^(1-alpha).*sigma_t.^(2*alpha) - s./k_s.^(1-alpha).*sigma_s.^(2*alpha)) - 1e-4;
omega       = 2*alpha-1e-4;
g2          = -0.5 - eta_t - sqrt((0.5 + eta_t).^2 + 2*(1-alpha)./(sigma_t.^2.*k_t));
pplus       = -g2 - 1;  
a           = 0.5*(pplus +1);

%% characteristic function of the increment of the process (between s & t)
LaplaceExp_t   = @(x) (t./k_t) .* ((1-alpha)/alpha) .* (1 - (1 + (x.*k_t.*sigma_t.^2)./(1-alpha)).^alpha);
LaplaceExp_s   = @(x) (s./k_s) .* ((1-alpha)/alpha) .* (1 - (1 + (x.*k_s.*sigma_s.^2)./(1-alpha)).^alpha);
phi_s_t        = @(x) exp(-1i.*x.*(LaplaceExp_t(eta_t) - LaplaceExp_s(eta_s))).*exp(LaplaceExp_t((x.^2 + 1i.*(1+2.*eta_t).*x)./2)-LaplaceExp_s((x.^2 + 1i.*(1+2.*eta_s).*x)./2));


%% FFT parameters
N           = 2^M;
h           = (2.*pi.*a./(b*N.^(omega))).^(1/(omega+1));
delta       = t-s;
D           = 5;
gamma       = (2*pi)/(N*h);
X0          = -gamma*(N-1)*0.5;
Xk          = -X0;
z1          = -0.5*h*(N-1);
zn          = -z1;

% in a struct
Params.z1   = z1;
Params.zn   = zn;
Params.dz   = h;
Params.x1   = X0;
Params.xn   = Xk;
Params.dx   = gamma;
Params.M    = M;
xx          = X0:gamma:Xk;

%% CDF on whole grid
f       = @(u) phi_s_t(u-1i.*a)./(1i.*u+a);                               % f to be transformed via Fourier
Phat    = 1-exp(-a.*xx)/(2*pi).*computeIntegral(f, xx, [], Params, 1);    % evaluating the CDF on the whole grid
% Phatfun = @(x) 1-exp(-a.*x)/(2*pi).*computeIntegral(f, x, [], Params, 1); % as function handle 

% Uncomment to have a plot of the "not adjusted CDF"
% figure
% plot(xx,Phat,'LineWidth',2)
% grid on
% title('CDF')


%% CDF on truncated domain
point = D*sqrt(delta);
idx   = find(xx >= point);

% if last point of the grid is less than D*sqrt(delta) keep the last point
if isempty(idx)
    idx = length(xx);
end

% new grid
Xk   = xx(idx(1)).*(abs(xx(idx(1)) - point) < abs(xx(idx(1)-1) - point))+...
     xx(idx(1)-1).*(abs(xx(idx(1)) - point) > abs(xx(idx(1)-1) - point));
X0   = -Xk;
xx   = X0:gamma:Xk;
Phat = 1-exp(-a.*xx)/(2*pi).*computeIntegral(f, xx, [], Params, 1);

% put to 0 values before the CDF is increasing 
ind1       = find(Phat>1);
Phat(ind1) = 1;
diff       = Phat(2:end) - Phat(1:end-1);
for i = 1:length(Phat)-1
    Phat(i) = Phat(i).*(all(diff(i:end)>=0)).*(Phat(i)>=0);
end

% plot of the CDF
figure
plot(xx,Phat,'LineWidth',2)
grid on
title('CDF via FFT on yearfraction of ', num2str(t-s))
xlabel('x')
ylabel('Phat')
legend('Phat')

%% simulating & interpolating the grid values

% uncomment to simulate between max & min of the CDF (no akima needed then)
% X1 = min(Phat(find(Phat)));
% U  = X1 + (Phat(end) - X1).*rand(Nsim, 1);

U = rand(Nsim, 1); % uniform rv simulations

% interpolation
switch flagspline
    
    case 1 % spline
            [Phatunique, idxunique, ~] = unique(Phat, 'stable'); 
            X                          = interp1(Phatunique,   xx(idxunique), U, 'spline');
            idx                        = find(X>xx(end)); % indexes where spline does fail
            
            % better to use this method when spline fails
            if ~isempty(idx)
                newX   = interp1(Phatunique, xx(idxunique), U(idx), 'makima');
                X(idx) = newX;
            end
            
    case 2 % linear
        
        [Phatunique, idxunique, ~] =unique(Phat, 'stable'); 
        X  = interp1(Phatunique, xx(idxunique), U);
end

%% PDF (uncomment to have a plot NB: a bit slow)
% density = @(x) a.*(exp(-a.*x)/pi).*h.*sum(real(exp(-1i.*([0:N-1]+0.5).*h.*x).*phi_s_t(([0:N-1]+0.5).*h-1i.*a)./(1i.*([0:N-1]+0.5).*h+a)))+...
%                         (exp(-a.*x)/pi).*h.*sum(real(1i*([0:N-1]+0.5).*h.*exp(-1i.*([0:N-1]+0.5).*h.*x).*phi_s_t(([0:N-1]+0.5).*h-1i.*a)./(1i.*([0:N-1]+0.5).*h+a)));
%  
% % adjustment as before
% densityres             = arrayfun(@(i) density(xx(i)), [1:length(xx)]');
% idx                    = find(Phat);
% densityres(1:idx(1)-1) = 0;
% 
% figure
% plot(xx, densityres, 'LineWidth', 2)
% grid on
% axis([-1.5 1.5 0 10])
% xlabel('x')
% ylabel('p')
% title('PDF on yearfraction of ', num2str(t-s))
% legend('PDF')


end
