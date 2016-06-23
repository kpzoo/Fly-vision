% Modified from getAltStats2 to rempve the unnecessary T and x1 inputs and
% to use as a standalone function for which Tn must be inserted

% Modified to include faster interpolation via nakeinterp1 mex files and
% fixed issue in calculation of

% Function to calculate mse based on timeseries of x1n and Tn and a linear
% interpolation - the code includes iterative methods to handle large data
% sets upon fine interpolation
function statSet = getStatsAltMeths3(Tn, x1n, x1capn)

% Remove duplicate points - the +1 in id calculation chooses the latter of
% the duplicates <---------------------- remove if want to remove former
dTn = diff(Tn);
id = find(dTn == 0) + 1; % <--------- the +1 has a slight effect if removed
lenTn = length(Tn);

% Tn = Tn - Tn(1);

% Obtain altered data vectors with duplicate time indices removed
idset = 1:lenTn;
trunc = setdiff(idset', id);
Tm = Tn(trunc);
x1m = x1n(trunc);
x1capm = x1capn(trunc);

% Remove effect of transient time on Tm if it has not already been removed
Tm = Tm - Tm(1);

% Get sample interval and sample points and separate into sets
inter = mean(diff(Tm))/20;
nSamps = (max(Tm) - min(Tm))/inter;
% nSamps = floor(nSamps);
nSetSamps = floor(nSamps/20000);
tSetSamps = (max(Tm) - min(Tm))/nSetSamps;

% Case in which both x1m and x1capm are to be linearly interpolated -
% perform iterative calculation of the statistics of interest
sum_em = 0;
sum_emSq = 0;
sum_len = 0;
tlim = -ones(1, nSetSamps+1);
tlim(1) = 0;  %<---------------- assumes that transient on Tm removed

% Loop across sets iteratively calculating statistics
for i = 2:nSetSamps+1
    % Define time limits and obtain relevant section of data
    tlim(i) = tSetSamps*(i-1);
    idtemp = find(Tm >= tlim(i-1) & Tm < tlim(i));
    Ttemp = Tm(idtemp);
    x1temp = x1m(idtemp);
    x1captemp = x1capm(idtemp);
    
    % Obtain sample times and interpolate em data to the sampled times
    Tsamp = Ttemp(1):inter:Ttemp(end);
%     x1Samp = nakeinterp1(Ttemp, x1temp, Tsamp');
    x1Samp = getSamplesExtrap2(1, 0, x1temp, Ttemp, Tsamp', NaN);
    x1capSamp = nakeinterp1(Ttemp, x1captemp, Tsamp');
    eSamp = x1Samp - x1capSamp;
    
    lenSamp = length(eSamp);
    if lenSamp ~= length(Tsamp)
        assignin('base', 'Tsamp', Tsamp);
        assignin('base', 'eSamp', eSamp);
        error(['The sampling of the error curve failed at i = ' num2str(i)]);
    end
    
    % Obtain iterative sums for statistics
    eSampSq = eSamp.*eSamp;
    sum_em = sum_em + sum(eSamp);
    sum_emSq = sum_emSq + sum(eSampSq);
    sum_len = sum_len + lenSamp;
end

% Check that the correct number of points were obtained and calculate means
if sum_len > nSamps 
    assignin('base', 'tlim', tlim);
    error(['Inconsistent sample size: [nSamps sum_len] = ' [num2str(nSamps)...
        ' ' num2str(sum_len)]]);
else
    em_mean = sum_em/sum_len;
    em_mse = sum_emSq/sum_len;
    em_var = em_mse - em_mean^2;
end

% Assign stats set with indicator string to indicate data order
statSet.order = {'mean', 'var', 'mse'};
statSet.interpMeth = [em_mean em_var em_mse];