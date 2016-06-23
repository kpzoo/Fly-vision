% Simple code to compare GP, Snyder and Linear filter results
clear all
clc
close all

% Assume data file with GP and Snyder data in a compiled subfolder
cd('compiled');
load('GPandSnyData', 'MSEc', 'MSEdb', 'beta', 'gamma');
cd ..

% Sort data across beta if not already sorted
if ~all(beta == sort(beta))
    [betaSG isort] = sort(beta);
    MSEc = MSEc(:, isort);
    MSEdb = MSEdb(:, isort);
else
    betaSG = beta;
end
if ~all(gamma == sort(gamma))
    [gammaSG isort] = sort(gamma);
    MSEc = MSEc(isort, :);
    MSEdb = MSEdb(isort, :);
else
    gammaSG = gamma;
end

% Obtain linear filter files which are in current directory
linfiles = dir('*.mat');
linflen = length(linfiles);

% Load a single file to obtain filter orders which should be consistent
load(linfiles(1).name, 'filOrder', 'betaSet');
filOrder1 = filOrder;
nFil = length(filOrder1);
betaFil = betaSet;

% Declare loop variables
gammaFil = zeros(1, linflen);
MSELinSet = zeros(linflen, length(betaSet), nFil);
minTest = 0;
minTrain = 0;

% Loop across files and obtain important variables
for ilin = 1:linflen
    % Load file and find minimum training and test sets
    load(linfiles(ilin).name);
    minTest = min(minTest, min(min(testSet)));
    minTrain = min(minTrain, min(min(trainSet)));
    
    % Obtain simulation parameters and check filter orders
    if all(gammaSet == gammaSet(1))
        gammaFil(ilin) = gammaSet(1);
    else
        error('The gamma values are not consistent beween filter data');
    end
    if max(abs((betaFil - betaSet))) > 10^-8
        error('The beta values are not consistent beween filter data');
    end
    if max(abs((filOrder1 - filOrder))) > 10^-8
        error('The filter orders are not consistent beween filter data');
    end
    
    % Obtain the MSE from several filter orders
    for iord = 1:length(filOrder1)
        MSELinSet(ilin, :, iord) = MSESet(iord, :);
    end
end

% Sort data across beta and gamma if not already sorted
if ~all(betaFil == sort(betaFil))
    [betaFil isortFil1] = sort(betaFil);
    MSELinSet = MSELinSet(:, isortFil1, :);
end
if ~all(gammaFil == sort(gammaFil))
    [gammaFil isortFil2] = sort(gammaFil);
    MSELinSet = MSELinSet(isortFil2, :, :);
end

% Check consistency between linear filter and Snyder-GP parameters
if max(abs((betaFil - betaSG))) > 10^-8
    error('The beta values are not consistent beween filter and Snyder-GP data');
end
if max(abs((gammaFil - gammaSG))) > 10^-8
    error('The gamma values are not consistent beween filter and Snyder-GP data');
end

% Plot comparisons between MSE of Snyder and linear filter across all
% filter orders as well as ratios
for ilin = 1:nFil
    % Comparison of raw MSE    
    MSELin = MSELinSet(:, :, ilin);
    figure;
    plot(betaSG, MSEc);
    hold on
    plot(betaFil, MSELin, 's-');
    hold off
    xlabel('beta');
    ylabel('mse estimates');
    title(['Comparison of MSE from Snyder, GP and Linear filter order ' num2str(filOrder1(ilin))]);
    saveas(gcf, ['compA' num2str(ilin)], 'fig');
    
    % Comparison of MSE ratio
    psi = MSELin./MSEc;
    figure;
    plot(betaSG, psi, 'o-');
    xlabel('beta');
    ylabel('mse FIR / mse Sny');
    title(['Comparison of MSE ratio of Linear filter order ' num2str(filOrder1(ilin)) ' and Snyder']);    
    saveas(gcf, ['compB' num2str(ilin)], 'fig');
end

% Plot comparisons between MSE of GP and linear filter across all
% filter orders as well as ratios
for ilin = 1:nFil
    % Comparison of raw MSE
    MSELin = MSELinSet(:, :, ilin);
    figure;
    plot(betaSG, MSEdb);
    hold on
    plot(betaFil, MSELin, 's-');
    hold off
    xlabel('beta');
    ylabel('mse estimates');
    title(['Comparison of MSE from Snyder, GP and Linear filter order ' num2str(filOrder1(ilin))]);
    saveas(gcf, ['compC' num2str(ilin)], 'fig');
    
    % Comparison of MSE ratio
    psi = MSELin./MSEdb;
    figure;
    plot(betaSG, psi, 'o-');
    xlabel('beta');
    ylabel('mse FIR / mse GP');
    title(['Comparison of MSE ratio of Linear filter order ' num2str(filOrder1(ilin)) ' and GP']);
    saveas(gcf, ['compD' num2str(ilin)], 'fig');
end

% Simpler plot of the curve with highest filter order and gamma
gam = length(gammaSG) - 7;
mseFil = MSELinSet(gam, :, end);
mseGP = MSEdb(gam, :);
mseSny = MSEc(gam, :);
figure;
plot(betaSG, [mseGP' mseSny' mseFil']);
legend('GP estimate', 'Snyder filter', ['FIR filter ' num2str(filOrder1(end))], 'location', 'best');
xlabel('beta');
ylabel('mse estimates');
title(['Comparison of MSE from Snyder, GP and Linear filter for gamma = ' num2str(gammaSG(gam))]);
