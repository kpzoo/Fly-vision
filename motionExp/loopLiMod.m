% Modified to allow for calculation of directional estimates by 2 different
% schemes with comparison of the optimal and suboptimal measures

% File to loop across the position-direction based light models with
% respect to the constant dot intensity parameter
clear all
clc
close all

% Set the intensity set of the dot implicitly
intenPm = 1:2:50;
mseDirec = zeros(size(intenPm));
mseDirecSub = zeros(size(intenPm));
mseCalc = zeros(size(intenPm));
dStats = cell(size(intenPm));
dStatsSub = cell(size(intenPm));
x1Stats = cell(size(intenPm));

% Loop across the main function
for i = 1:length(intenPm)
    [mseDirec(i) mseDirecSub(i) mseCalc(i) dStats{i} dStatsSub{i} x1Stats{i} MCrates] = liMod2Fn(intenPm(i));
    disp(['Completed iteration ' num2str(i) ' of ' num2str(length(intenPm))])
end

% Save data
save('noDelData.mat');

% Plot the main result of directional MSE
figure;
inten = MCrates(1)*intenPm*1000;
plot(inten, mseDirec, inten, mseDirecSub);
legend('optimal', 'suboptimal', 'location', 'best');
xlabel('intensity (s^{-1})');
ylabel('mse for directional estimates');
title('Comparison of optimal and suboptimal directional estimates without delay');
saveas(gcf, 'compNoDel.fig');
