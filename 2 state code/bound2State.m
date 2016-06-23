%%
% Two state bound comparison with kb = kd = k and y = alpha*x
clear all
clc
close all

% Set inputs
beta = 0:0.01:100;

% Calculate slack bound
% dslack = 1./(beta*log(2) + 2);
dslack = 1./(beta*log(2) + 4); % <---------- corrected version

% Calculate tight bound
dtight = 1./(1 + sqrt(1 + 4*beta));

% Plot comparison
figure;
semilogx(beta, [dtight' dslack']);
legend('tight', 'slack', 'location', 'best');
xlabel('beta');
ylabel('bounds');
title('Comparison of tight and slack bound for 2 state CTMC');

%%
% A comparison of the behaviour of a maximum bound on the mse based on the
% q solution minimum

% Obtain the mse bound
a = beta + 0.5;
b = sqrt(beta.^2 + 0.25);
mseMax = (1 - (a - b)).^2;

% Obtain psi values
psi_tight = mseMax./dtight;
psi_slack = mseMax./dslack;

% Plot the comparisons
figure;
semilogx(beta, [dtight' dslack' mseMax']);
legend('tight', 'slack', 'max mse', 'location', 'best');
xlabel('beta');
ylabel('bounds and mse');
title('Comparison of LVP bounds for 2 state CTMC to a max mse bound');

figure;
loglog(beta, [psi_tight' psi_slack']);
legend('tight', 'slack', 'location', 'best');
xlabel('beta');
ylabel('psi values');
ylim([1 max([psi_tight psi_slack])])
title('Comparison of LVP psi for 2 state CTMC');