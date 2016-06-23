% Simple plotting function for converting the photon time series into point
% process plots and stem plots
function plotPhotonStreamComparison(Treal, Test)

% Obtain births and spikes for each photon stream
xreal = 1:length(Treal);
xest = 1:length(Test);
deltaReal = ones(size(Treal));
deltaEst = ones(size(Test));

% Plot the point processes
figure;
stairs(Treal, xreal);
hold on
stairs(Test, xest, 'r');
hold off
xlabel('time');
ylabel('photon counts');
legend('original (real)', 'estimated (distorted)', 'location', 'best');
title('Comparison of photon streams in point process form');

% Plot the photon events as spikes
figure;
stem(Treal, deltaReal);
hold on
stem(Test, deltaEst, 'r');
hold off
xlabel('time');
ylabel('photon streams');
legend('original (real)', 'estimated (distorted)', 'location', 'best');
title('Comparison of photon streams in spike form');
ylim([0 1.1]);
xlim([min(min(Treal), min(Test)), max(max(Treal), max(Test))]);