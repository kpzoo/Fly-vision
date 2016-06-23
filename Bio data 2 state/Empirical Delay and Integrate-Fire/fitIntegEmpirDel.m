% Script to fit the empirically delayed MSE and that obtained from the
% integrate and fire method
clear all
clc
close all

% Load appriopriate data
load('empirAndInteg');

% Fit the empirical delay data
n1 = 6;
[p,S,mu] = polyfit(betaSet, mseD, n1);
[fD del] = polyval(p, betaSet, S, mu);
figure;
plot(betaSet, mseD, 'o', betaSet, fD);
xlabel('beta');
ylabel('mse estimate');
legend('real data', ['fitted curve (n = ' num2str(n1) ')'], 'location', 'best');
title('Fitting of empirical delay data');

% Fit the integrate-fire data
n2 = 6;
[p,S,mu] = polyfit(betaSet, mseI, n2);
[fI del] = polyval(p, betaSet, S, mu);
figure;
plot(betaSet, mseI, 'o', betaSet, fI);
xlabel('beta');
ylabel('mse estimate');
legend('real data', ['fitted curve (n = ' num2str(n2) ')'], 'location', 'best');
title('Fitting of integrate-fire data');

% Compare the fits
figure;
plot(betaSet, fD, betaSet, fI);
xlabel('beta');
ylabel('mse estimate');
legend('fit to empirical delay', 'fit to integrate-fire', 'location', 'best');
title(['Comparisons of fits at n1 = ' num2str(n1) ' and n2 = ' num2str(n2)]);

% Percentage difference between the curves
ed = abs(fD - fI);
p1 = 100*ed./fI;
p2 = 100*ed./fD;
mp1 = mean(p1);
mp2 = mean(p2);