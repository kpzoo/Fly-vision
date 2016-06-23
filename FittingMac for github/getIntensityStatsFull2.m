% This version evaluates the intensity at all the ODE points rather than at
% the jump points only - it is assumed that both lam and lamcap are
% constant between the ODE points

% Modification to remove a transient from the error curve which results
% from the initial distribution guess on the filter

% Function to calculate statistics of actual and estimated intensity
function [lamn lamcapn Tn stats qn] = getIntensityStatsFull2(lam, lamcapFull, Tset, T, Qset, rem_trans)

% Ensure first element of cells are empty (as no ODE solutions here)
if ~isempty(lamcapFull{1} + Tset{1})
    error('The first cell is not empty');
end

% Concatenate the data into a single stream
lenOuter = length(lamcapFull);
elemLen = zeros(1, lenOuter);
for i = 1:lenOuter
    elemLen(i) = length(lamcapFull{i});
end
lenFull = sum(elemLen);

% Initialise main data vectors and assign in loop - use -1 as indicator
Tn = -ones(lenFull, 1);
lamcapn = -ones(lenFull, 1);
stop = 0;
qSize = size(Qset{2});
qSize = qSize(2);
qn = -ones(lenFull, qSize);
for i = 2:lenOuter
    % Loop calculates the start and end points along the array successively
    % and then assigns the appropriate cell element
    start = stop + 1;
    stop = elemLen(i) + start - 1;
    Tn(start:stop) = Tset{i};
    lamcapn(start:stop) = lamcapFull{i};
    qn(start:stop, :) = Qset{i};
end

% Check correct length obtained and initialise lamn
if ~(length(Tn) == length(lamcapn) && length(Tn) == lenFull)
    error('Data vectors have inconsistent lengths');
end

% Obtain exact samples (or zoh interpolation) of lam from its data
ver = 1;
samp_inter = 'N/A';
[lamn Tn] = getSamples2(ver, samp_inter, lam, T, Tn);
lamn = lamn';

% assignin('base', 'lamn', lamn);
% assignin('base', 'Tn', Tn);
% assignin('base', 'lam', lam);
% assignin('base', 'T', T);

% Obtain relevant statistics and display - in this case the calculations
% are exact with assumption that lamcap does not change between ODE points
en = lamn - lamcapn;
if rem_trans
    stop = length(en);
    start = ceil(0.25*stop);
    enSS = en(start:stop);
    TnSS = en(start:stop);
else
    enSS = en;
    TnSS = Tn;
end
[stats.maxErr stats.meanErr stats.varErr stats.minErr] = calcStats2(enSS, TnSS);
stats.mseErr = stats.meanErr^2 + stats.varErr;

% This statistics calculation only holds if the points are equidistant in
% time instead of
% stats.meanErr = mean(en);
% stats.mseErr = mean(en.^2);
% stats.varErr = var(en);

disp(['Stats are: [<e> var(e) <e^2>] = ' num2str(stats.meanErr) ' '...
    num2str(stats.varErr) ' ' num2str(stats.mseErr)]);