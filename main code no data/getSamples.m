% Simple function to re-sample the stochastic reaction timeseries
function [xn tn] = getSamples(samp_inter, x, t)

% Divide t into number of samples
tn = min(t):samp_inter:max(t);
len = length(tn);
disp(['Actual number of samples is ' num2str(len)]);

% Obtain relevant sample x values noting that x is stepwise
xn = zeros(1, len);
for i = 1:len
    id = find(t <= tn(i), 1, 'last');
    xn(i) = x(id);
end
    