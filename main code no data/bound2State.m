% Two state bound comparison with kb = kd = k and y = alpha*x
clear all
clc
close all

% Set inputs
beta = 0:0.001:1000;

% Calculate slack bound
dslack = 1./(beta*log(2) + 2);

% Calculate tight bound
dtight = 1./(1 + sqrt(1 + 4*beta));

% Plot comparison
figure;
semilogx(beta, [dtight' dslack']);
legend('tight', 'slack', 'location', 'best');
xlabel('beta');
ylabel('bounds');
title('Comparison of tight and slack bound for 2 state CTMC');

