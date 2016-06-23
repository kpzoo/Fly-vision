% Function to draw from an empirical delay distribution specified by
% histogram values given by xdata = bins, ydata = freq
function samples = drawEmpiricalDistr(xdata, ydata, nSamps)

% Remove zero entries from ydata
id = find(ydata ~= 0);
ydata = ydata(id);
xdata = xdata(id);

% Normalise ydata to probabilities
ydata = ydata/sum(ydata);

% Obtain specified number of weighted random samples
samples = randsample(xdata, nSamps, true, ydata);