%% Computer Exercise 1 - Time Series analysis
%% Task 1
clear;
A1 = [ 1 -1.79 0.84 ] ;
C1 = [ 1 -0.18 -0.11 ] ;

A2 = [ 1 -1.79 ] ;
C2 = [ 1 -0.18 -0.11 ] ;
%%

ARMA_poly1 = idpoly(A1, [], C1);
ARMA_poly2 = idpoly(A2, [], C2);

% ARMA_poly1 = idpoly(A1, [], C1); % To see pzplot
% ARMA_poly1.a  % To access A-polynomial
% ARMA_poly1.C % To access C-polynomial 
%%
rng(0);
N = 300;
sigma2 = 1.5; % Setting samples and variance of noise e
 % Generate some random noise
%% Accomplished by the created function
% e = sqrt(sigma2)*randn(N, 1) ;
% y1 = filter(ARMA_poly1.c , ARMA_poly1.a, e);
% y2 = filter(ARMA_poly2.c , ARMA_poly2.a, e);
%%
y1 = simulateMyARMA(N, sigma2, ARMA_poly1);
y2 = simulateMyARMA(N, sigma2, ARMA_poly2);
y1 = y1(101:end); % omit initial samples due to unexpected behaviour
y2 = y2(101:end);
%% Plotting and studying behaviour 
subplot(211);
plot(y1) % remains stable
subplot(212)
plot(y2); % diverges, lets study pzplot
%%
% title("Pole-Zero plots for ARMA\_poly1 and ARMA\_poly2 respectively")
figure();
subplot(211);
pzplot(ARMA_poly1)
subplot(212);
pzplot(ARMA_poly2) % Pole (far) unit circle => the AR-part is non-stationary 
%% Compares theoretical and estimated covariance function for ARMA_poly1. 
% The estimated covariance seems to follow the theoretical result quite
% well but underestimates slightly.
% Why? estimated covariance is only asymptotically identical to
% theoretical. We use the biased estimator with 1/N in denominator which
% reduces the values more than the unbiased (but its better because the
% variance of the estimate is constant?)
figure();
m = 20; % Should never be higher than N/4 but thats obviously not a problem
r_theo = kovarians(ARMA_poly1.c, ARMA_poly1.a, m); % Computes theoretical covariance-function given C and A
stem(0:m, r_theo*sigma2);
hold on;
r_est = covf(y1, m+1); % Estimates covariance function given a realization;
stem(0:m, r_est, 'r');
%% Proceeding to re-estimate our model using the realization y1
% "Basic analysis using ACF, PACF and normplot"
[ACFest, PACFest] = ACFPACFNormplot(y1, m);
%%
data = iddata(y1);
%%
ar_model2 = arx(y1, 2);
ar_model3 = arx(y1, 3);
arma_model22= armax(y1, [2,2]);



