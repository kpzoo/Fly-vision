% Simple script to check that a gamma distribution can be fit
clear all
clc
close all

% Load the empirical data, draw samples and obtain pdf and cdf
load('figdata.mat');
nSamps = 10000;
samples = drawEmpiricalDistr(xdata, ydata, nSamps);
emp_pdf = ydata./sum(ydata);
emp_cdf = cumsum(emp_pdf);

% Fit gamma distribution and obtain pdfs and cdfs
param = gamfit(samples);
fit_pdf = gampdf(xdata, param(1), param(2));
fit_cdf = gamcdf(xdata, param(1), param(2));
fitsamples = gamrnd(param(1), param(2), nSamps, 1);

% Renormalise the fitted pdf so that it has values only at the desired
% discrete points and compare cumulative distributions
fit_pdf = fit_pdf./(sum(fit_pdf));
fit_cdf2 = cumsum(fit_pdf);


% Plot comparison of the fitted and empirical distributions
figure;
stairs(xdata, fit_cdf);
hold on
stairs(xdata, fit_cdf2, 'r');
hold off
xlabel('delay (ms)');
ylabel('cumulative probability');
title('Comparison of gamma and empirical cumulative distributions with');
legend('directly fitted cdf', 'cdf from fitted, normalised pdf', 'location', 'best'); 

% Plot comparison of the fitted and empirical distributions
figure;
stairs(xdata, emp_pdf);
hold on
stairs(xdata, fit_pdf, 'r');
hold off
xlabel('delay (ms)');
ylabel('probability');
title('Comparison of gamma and empirical distributions with');
legend('empirical', ['gamma @ ' [num2str(param(1)) ' ' num2str(param(2))]], 'location', 'best'); 

% Plot comparison of the fitted and empirical cumulative distributions
figure;
stairs(xdata, emp_cdf);
hold on
stairs(xdata, fit_cdf, 'r');
hold off
xlabel('delay (ms)');
ylabel('cumulative probability');
title('Comparison of gamma and empirical cumulative distributions with');
legend('empirical', ['gamma @ ' [num2str(param(1)) ' ' num2str(param(2))]], 'location', 'best'); 