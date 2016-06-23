% Modification to allow extrapolation with chosen values such as NaN or 0
% as indicators of the extrapolation

% Simple function to re-sample the stochastic reaction timeseries - updated
% to do uniform sampling with a specified interval or to use the actual
% samples of time provided
function [xn tn] = getSamplesExtrap(ver, samp_inter, x, t, tn, extrapID)

% Obtain appropriate output vector length and calculate tn if needed
if ~ver
    % Divide t into number of samples if none provided
    tn = min(t):samp_inter:max(t);
    len = length(tn);
    disp(['Actual number of samples is ' num2str(len)]);
else
    len = length(tn);
end

% Obtain relevant sample x values noting that x is stepwise continuous
xn = zeros(1, len);
maxt = max(t);
mint = min(t);
for i = 1:len
    % Extrapolate if tn is outside range of t
    if tn(i) < mint || tn(i) > maxt
        xn(i) = extrapID;
    else
        % Take value from last event as the function has not changed since
        id = find(t <= tn(i), 1, 'last');
        xn(i) = x(id);
%         xtemp = x(t <= tn(i));
%         xn(i) = xtemp(end);
    end
end
