% Code assumes that lam = lam(x1) is constant between x1 jumps while lamcap
% can change between x1 jumps and thus must be obtained via interpolation

% Function to calculate statistics of actual and estimated intensity
function [lamn lamcapn Tn stats] = getIntensityStats(lam, lamcap, T)

% Obtain exact samples (or zoh interpolation) of lam from its data
samp_inter = mean(diff(T))/20;
[lamn Tn] = getSamples(samp_inter, lam, T);

% Obtain linearly interpolated samples of lamcap from its timeseries
lamcapn = interp1(T, lamcap, Tn);
if ~all(size(lamcapn) == size(lamn))
    assignin('base', 'lamn', lamn);
    assignin('base', 'lamcapn', lamcapn);
    error('Mismatch in the sampled data sizes of lamcapn and lamn');
end

% Obtain relevant statistics and display
en = lamn - lamcapn;
stats.meanErr = mean(en);
stats.mseErr = mean(en.^2);
stats.varErr = var(en);
disp(['Stats are: [<e> var(e) <e^2>] = ' num2str(stats.meanErr) ' '...
    num2str(stats.varErr) ' ' num2str(stats.mseErr)]);