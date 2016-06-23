% Modified to allow direct comparison of the deterministically delayed and
% normal Snyder results which ensures the random trajectory used in both
% cases is exactly the same

% Modified version of loopLiMod which allows for the application of a
% deterministic delay on the photons of all streams

% File to loop across the position-direction based light models with
% respect to the constant dot intensity parameter
clear all
clc
close all

% Set the intensity set of the dot implicitly and the explicit delay (ms)
intenPm = 1:2:50;
% intenPm = 1:5:20;
delayPm = 43.3*10;

% Output variable declaration
mseDirec = zeros(size(intenPm));
mseDirecSub = zeros(size(intenPm));
mseCalc = zeros(size(intenPm));
dStats = cell(size(intenPm));
dStatsSub = cell(size(intenPm));
x1Stats = cell(size(intenPm));

% Loop across the main function
for i = 1:length(intenPm)
    [mseDirec(i) mseDirecSub(i) mseCalc(i) dStats{i} dStatsSub{i} x1Stats{i} MCrates delayOn] = liMod3Fn(intenPm(i), delayPm);
    disp(['Completed iteration ' num2str(i) ' of ' num2str(length(intenPm))])
end

% Save data
save('delData10x.mat');

% Plot the main result
figure;
inten = MCrates(1)*intenPm*1000;
plot(inten, mseDirec, inten, mseDirecSub);
legend('optimal', 'suboptimal', 'location', 'best');
xlabel('intensity (s^{-1})');
ylabel('mse for directional estimates');
title('Comparison of optimal and suboptimal directional estimates with delay');
saveas(gcf, 'compDel10x.fig');