% Simple code to plot comparison between GP and the threshold methods -
% must be inserted into the complete folder <-----------------------------
clear all
clc
close all

% Set names of GP and threshold Snyder data
nameSny = 'threshBat10';
nameGP = 'dbatchCont10';

% Load the GP data
load(nameGP, 'mseErrx1');
load(nameGP, 'beta');
betaGP = beta;
mseGP = mseErrx1;
[betaGP isort] = sort(betaGP);
mseGP = mseGP(isort);

% Load all of the Snyder threshold data 
load(nameSny);

% Sort data
[betaSet isort] = sort(betaSet);
mseAct = mseAct(isort);
mseEst = mseEst(isort);
distTest = distTest(isort);
distTrain = distTrain(isort);

% Check for constant gamma
if ~all(gammaSet == gammaSet(1))
    error('Data not for a constant gamma');
else
    gam = gammaSet(1);
end

% Plot comparisons of MSE
figure;
plot(betaGP, mseGP, 'bo-');
hold all
plot(betaSet, mseAct, 'ko-', betaSet, mseEst, 'ro-');
hold off
xlabel('beta');
ylabel('mse estimates');
legend('GP', 'true Sny', 'est Sny', 'location', 'best');
title(['Comparison of GP and threshold Snyder at gamma = ' num2str(gam)]);

% Plot comparisons of metrics in training (optimal) and testing
figure;
plot(betaSet, distTrain, betaSet, distTest);
xlabel('beta');
ylabel('metric at optimal threshold');
legend('training', 'testing', 'location', 'best');
title(['Training and test metric values from optimal thresholds at gamma = ' num2str(gam)]);
