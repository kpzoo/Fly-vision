% Simple plotting script for integrate-fire results with empirical delay -
% assumes all the required data is loaded into workspace
close all
clc

% Sort data from integrate-fire or thresholding code
[betaSet isort] = sort(betaSet);
mseAct = mseAct(isort);
mseEst = mseEst(isort);
distTrain = distTrain(isort);
distTest = distTest(isort);
lenTreal = lenTreal(isort);
lenTest = lenTest(isort);
meanAct = meanAct(isort);
meanEst = meanEst(isort);

% Plot the MSE comparisons
figure;
plot(betaSet, mseAct, betaSet, mseEst);
hold on
plot(beta, MSEest, 'k');
hold off
xlabel('beta');
ylabel('mse values');
legend('optimal', 'integrate-fire', 'empirical delay', 'location', 'best');
title(['Comparison of MSE estimates of 16 state bimodal MC at gamma = ' num2str(gammaSet(1))]);

% Plot the mean error comparisons
figure;
plot(betaSet, meanAct, betaSet, meanEst);
xlabel('beta');
ylabel('mean error values');
legend('optimal', 'integrate-fire', 'location', 'best');
title(['Comparison of mean error estimates of 16 state bimodal MC at gamma = ' num2str(gammaSet(1))]);

% Plot the estimated number of photons
figure;
plot(betaSet, lenTreal, betaSet, lenTest)
xlabel('beta');
ylabel('photon numbers');
title(['Estimated and actual photons for 16 state bimodal MC, gamma = ' num2str(gammaSet(1))]);
legend('real', 'estimated', 'location', 'best');

% Plot the spike metrics
figure;
plot(betaSet, distTrain, betaSet, distTest)
xlabel('beta');
ylabel('metric values');
title(['Comparison of spike metric for 16 state bimodal MC, gamma = ' num2str(gammaSet(1))]);
legend('training', 'testing', 'location', 'best');
