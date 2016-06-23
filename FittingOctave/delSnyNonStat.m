% Function to run the compensated Snyder filter which accounts for a
% delayed photon stream - based on cleanDistort2 from Snyder Test 2
function out = delSnyNonStat(locfolder1, locfolder2, adjust, noiseTraits, limPhoton, simname)

% Obtain suitable parameters for Snyder filter
inpSny = prepareSnyderInputs(simname, locfolder1, locfolder2, adjust, limPhoton);
x1 = inpSny.x1;
x2 = inpSny.x2;
T = inpSny.T;
params = inpSny.params;
assignin('base', 'paramsInp', params);
Q = inpSny.Q;
S = inpSny.S;

% Get real photon times (x2 birth times) and ensure no bound was calculated
[Treal perc] = getEventTimes(T, x2, 'birth');
lenTreal = length(Treal);
assignin('base', 'TeventReal', Treal);
    
% Run normal Snyder filter on real photon times
no_bnd = 1;
[lamstatsAct x1StatsAct x1nAct x1capnAct TnAct] = ...
    bioSnyderMod(T, [x1 x2], params.kgain, Q, S, params, no_bnd);
if ~isnan(x1StatsAct.meth3(3))
    meanAct = x1StatsAct.meth3(1);
    mseAct = x1StatsAct.meth3(3);
else
    meanAct = x1StatsAct.meth2(1);
    mseAct = x1StatsAct.meth2(3);
end
assignin('base', 'TnAct', TnAct);
assignin('base', 'x1nAct', x1nAct);
assignin('base', 'x1capnAct', x1capnAct);
assignin('base', 'x1StatsAct', x1StatsAct);
% clear x1nAct x1capnAct TnAct

% Declare loop variables
lenTh = length(noiseTraits);
mseEst = zeros(lenTh, 1);
meanEst = zeros(lenTh, 1);
Test = cell(lenTh, 1);
lenTest = zeros(lenTh, 1);
Lam = cell(lenTh, 1);
Qs = cell(lenTh, 1);

% Inputs for compensated Snyder formulation
Sx = diag(min(x1):max(x1));
rateParams.alpha = params.kgain;
rateParams.k = params.kdeath;

% Loop across several distortion settings
for iw = 1:lenTh
    % Obtain distorted photon times based on specified method
    Test{iw} = distortPhotons2(noiseTraits{iw}, Treal, T, params.kgain);
    lenTest(iw) = length(Test{iw});
    
    % Parameters controlling the delay distortion as applied to the new
    % Snyder formulation
    params.eta = noiseTraits{iw}.paramDistr;
    params.eta = params.eta(1);
    rateParams.eta = params.eta; % <---------------------- need to implement

    
    % Construct new data sets featuring the delayed photon time estimates
    [Tnew Xnew] = appendStream(T, [x1 x2], Test{iw});
    z = Xnew(:, end);
    assignin('base', 'z', z);
    assignin('base', 'Tnew', Tnew);
    
    % Run compensated Snyder filter with the delayed photons z with
    % appropriate calculation of remaining inputs
    % Tnew = Testnew; % <-------------- check this input as not sure if only events or full series
    % [x1StatsEst Lam{iw} Qs{iw}] = compSnyder(Tnew, Xnew, Sx, Sy, rateParams, x1, T);
    Xnew(:, end-1) = Xnew(:, end-1) + z;
    y = Xnew(:, end-1);
    Sy = diag(min(y):max(y));
%     assignin('base', 'y', y);
%     assignin('base', 'Xnew', Xnew);
    [x1StatsEst Lam{iw} Qs{iw}] = compSnyderNonStat2(Tnew, Xnew, Sx, Sy, rateParams, x1, T);
    assignin('base', 'x1StatsEst', x1StatsEst);
    % Assign compensated Snyder outputs
    meanEst(iw) = x1StatsEst.meth3(1);
    mseEst(iw) = x1StatsEst.meth3(3);
end

% Assign main output structure
out.meanAct = meanAct;
out.mseAct = mseAct;
out.meanEst = meanEst;
out.mseEst = mseEst;
out.Test = Test;
out.Treal = Treal;
out.lenTreal = lenTreal;
out.lenTest = lenTest;
out.beta = inpSny.beta;
out.gamma = inpSny.gamma;
out.x1StatsEst = x1StatsEst;
out.Qs = Qs;
out.Lam = Lam;
out.Sx = Sx;
out.Sy = Sy;