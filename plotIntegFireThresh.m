% Simple plotting script for integrate-fire results with empirical delay -
% assumes all the required data is loaded into workspace
close all
clc

% Set title script
titNo = 2;
switch(titNo)
    case 1
        str1 = ['Comparison of MSE estimates of 2 state symmetric MC at gamma = ' num2str(gammaSet(1))];
        str2 = ['Comparison of mean errors of 2 state symmetric MC at gamma = ' num2str(gammaSet(1))];
        str3 = ['Comparison of photon numbers for 2 state symmetric MC at gamma = ' num2str(gammaSet(1))];
        str4 = ['Comparison of metric values for 2 state symmetric MC at gamma = ' num2str(gammaSet(1))];
    case 2
        str1 = ['Comparison of MSE estimates of 16 state bimodal MC at gamma = ' num2str(gammaSet(1))];
        str2 = ['Comparison of mean errors of 16 state bimodal MC at gamma = ' num2str(gammaSet(1))];
        str3 = ['Comparison of photon numbers for 16 state bimodal MC at gamma = ' num2str(gammaSet(1))];
        str4 = ['Comparison of metric values for 16 state bimodal MC at gamma = ' num2str(gammaSet(1))];
end

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
title(str1);

% Plot the mean error comparisons
figure;
plot(betaSet, meanAct, betaSet, meanEst);
xlabel('beta');
ylabel('mean error values');
legend('optimal', 'integrate-fire', 'location', 'best');
title(str2);

% Plot the estimated number of photons
figure;
plot(betaSet, lenTreal, betaSet, lenTest)
xlabel('beta');
ylabel('photon numbers');
title(str3);
legend('real', 'estimated', 'location', 'best');

% Plot the spike metrics
figure;
plot(betaSet, distTrain, betaSet, distTest)
xlabel('beta');
ylabel('metric values');
title(str4);
legend('training', 'testing', 'location', 'best');