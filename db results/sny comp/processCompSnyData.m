% Simple script to process the Snyder compilation and create sufficient
% plots with necessary comparisons
clear all
clc
close all

% Obtain files --- assume code in same folder as data
files = dir('*.mat');
flen = length(files);

% Assume beta same across files and get theoretical bounds
load(files(1).name);
beta = beta_ratio;
dslack = 1./(beta*log(2) + 4);
dtight = 1./(1 + sqrt(1 + 4*beta));
dlen = length(beta);

% Declare output variables
MSE = zeros(flen, dlen);
MSEb = zeros(flen, dlen);
MSEc = zeros(flen, dlen);
kbset = zeros(1, flen);

% Loop through files and extract most important variables
for in = 1:flen
    load(files(in).name);
    MSE(in, :) = msex1a1;
    MSEb(in, :) = msex1a2;
    MSEc(in, :) = msex1a3;
    kbset(in) = kbirth;
    disp(['Processed ' files(in).name]);
end

% Calculate gamma assuming width of 100ms for QBs
gamma = 1./kbset/100;

% Plot comparison of all MSE values over beta for various gamma
figure;
plot(beta, dslack, 'k--', beta, dtight, 'r--');
hold on
plot(beta, MSE);
hold off
legend('dslack', 'dtight', 'MSE Snyder', 'location', 'best');
xlabel('beta');
ylabel('mse and dslack');
title(['Comparison of MSE Snyder across beta across gamma = [' [num2str(gamma(1)) ' ' num2str(gamma(end))] ']']);

% Similar plots with other mse estimates
figure;
plot(beta, dslack, 'k--', beta, dtight, 'r--');
hold on
plot(beta, MSEb);
hold off
legend('dslack', 'dtight', 'MSE Snyder', 'location', 'best');
xlabel('beta');
ylabel('mse and dslack');
title(['Comparison of MSE Snyder across beta across gamma = [' [num2str(gamma(1)) ' ' num2str(gamma(end))] ']']);

% Similar plots with other mse estimates
figure;
plot(beta, dslack, 'k--', beta, dtight, 'r--');
hold on
plot(beta, MSEc);
hold off
legend('dslack', 'dtight', 'MSE Snyder', 'location', 'best');
xlabel('beta');
ylabel('mse and dslack');
title(['Comparison of MSE Snyder across beta across gamma = [' [num2str(gamma(1)) ' ' num2str(gamma(end))] ']']);

% Plot of all MSE estimates
figure;
plot(beta, MSE, 'b', beta, MSEb, 'r', beta, MSEc, 'm');
xlabel('beta');
ylabel('mse estimates');
title('Comparison of several mse estimates');