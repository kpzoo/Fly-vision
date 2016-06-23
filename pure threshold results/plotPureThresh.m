% Script to plot results - assumes data already sorted and compiled - first
% plots are of MSE and mean error
figure;
plot(beta, MSEact, beta, MSEest)
xlabel('beta');
ylabel('mse estimates');
legend('real photon', 'pure threshold estimates', 'location', 'best')
title(['Effect of pure thresholding on Snyder MSE at gamma = ' num2str(gammaSet(1))]);

figure;
plot(beta, Mact, beta, Mest)
xlabel('beta');
ylabel('mean error estimates');
legend('real photon', 'pure threshold estimates', 'location', 'best')
title(['Effect of pure thresholds on Snyder mean error at gamma = ' num2str(gammaSet(1))]);

% Plot the real and estimated number of photons
figure;
plot(beta, ltreal, beta, ltest)
xlabel('beta');
ylabel('no. photons');
legend('real', 'estimate', 'location', 'best')
title(['Photons counted via pure thresholding at gamma = ' num2str(gammaSet(1))]);

% Plot the training and testing metric values
figure;
plot(beta, dTrain, beta, dTest)
xlabel('beta');
ylabel('metric distance');
legend('training', 'testing', 'location', 'best')
title(['Training and testing metrics at [dtcost gamma] = ' [num2str(dtcost) ' ' num2str(gammaSet(1))]]);