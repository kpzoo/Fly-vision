% Function to obtain statistics of Markov chains aimed at representing the
% fluctuations in light intensity (provided lam(x1) = kx1)
function [param m v stat distr] = getMarkovStats(T, x1, plot_on, distr)

% Sample data with equal time steps and obtain exact samples
samp_inter = mean(diff(T))/10;
ver = 0;
[x1n Tn] = getSamples2(ver, samp_inter, x1, T, []);

% Obtain normalised histogram of state distribution (assumed min:max)
len = length(x1n);
nbins = min(x1):max(x1); % <--------- assumption on state space of x1
nVals = hist(x1n, nbins);

% Fit data to a standard distribution type
[param m v stat] = fitdistr2(x1n, plot_on, distr, nbins);
