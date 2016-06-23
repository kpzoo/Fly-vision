% Modified to allow lamn to be inserted in the case of non-homogeneous
% Poisson processes rather than interpolated
% Modified to remove the stats calculation and rem_trans input

% By removing the statistics calculation there is no need to assumes the
% lam function is constant between ODE points

% Modification to remove a transient from the error curve which results
% from the initial distribution guess on the filter

% Function to calculate statistics of actual and estimated intensity
function [lamn lamcapn Tn qn] = getIntensityStatsFull3(lam, lamcapFull, Tset, T, Qset, calclamn)

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
if calclamn
    ver = 1;
    samp_inter = 'N/A';
    [lamn Tn] = getSamples2(ver, samp_inter, lam, T, Tn);
    lamn = lamn';
else
    lamn = 0;
end