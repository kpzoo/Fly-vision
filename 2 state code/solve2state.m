% This first cell solves the single ODE that determines the q behaviour
% between jumps of the output

% Solution to 2 state ODE
clear all
clc
close all

% Assign parameters and solve coefficients
alpha = 20;
k = 10;
T = 2/k;
delt = 2/alpha;
duty = k/alpha;
a = duty + 0.5;
b = sqrt(duty^2 + 0.25);
B = atanh((a - 1)/b);

% Obtain ODE solution across time
t = 0:0.001:T;
q = a - b*tanh(b*alpha*t + B);

% Obtain mean of curve and mean error
qmean = a - (1/(alpha*T))*(log(cosh(b*alpha*T + B)) - log(cosh(B)));
emean = 0.5 - qmean;

% Calculate steady state q value
qSS = a - b;

% Estimate numerically the statistics of e
e = 0.5 - q;
eEst.mean = mean(e);
eEst.var = var(e);
eEst.mse = eEst.var + eEst.mean^2;

% Plot resulting curve and mean
figure;
plot(t, q);
hold on
plot(t, qmean*ones(size(t)), 'r');
xlabel('t');
ylabel('q');
title(['Solution to 2 state Markov chain at [alpha k] = ' [num2str(alpha) ' ' num2str(k)]]);

%% 
% This cell looks at some general behaviour of the q solutions

% Obtain how steady state q depends on alpha
alpha1 = 0.01:0.01:10000;
par = k./alpha1;
a1 = par + 0.5;
b1 = sqrt(par.^2 + 0.25);
qSS1 = a1 - b1;

% Plot the variation of qSS with k/alpha
figure;
semilogx(1./par, qSS1);
xlabel('alpha/k');
ylabel('steady state q');
title('Variation of steady state q from ODE with parameters');
