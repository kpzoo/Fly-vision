% Modified to work with distortPhotons2 which mixes empirical distortions
% with photon deletions

% Function to distort a Gillespie time series photon stream (x2 events)
% with clearly defined distortions based on insertions, deletions and
% delays drawn from specified distributions
function out = cleanDistort2(locfolder1, locfolder2, adjust, noiseTraits, limPhoton, simname)

% Obtain suitable parameters for Snyder filter
inpSny = prepareSnyderInputs(simname, locfolder1, locfolder2, adjust, limPhoton);
x1 = inpSny.x1;
x2 = inpSny.x2;
T = inpSny.T;
params = inpSny.params;
Q = inpSny.Q;
S = inpSny.S;

% Calculate dtcost in terms of the mean time between photons for an IPP
dtcost = 2/params.kgain;
metricPhotonLim = 5000;

% Get real photon times (x2 birth times) and ensure no bound was calculated
[Treal perc] = getEventTimes(T, x2, 'birth');
lenTreal = length(Treal);
if lenTreal > metricPhotonLim
    idend1 = metricPhotonLim;
else
    idend1 = lenTreal;
end
    
% Run Snyder filter on real photon times
no_bnd = 1;
[lamcapAct x1capAct lamstatsAct x1StatsAct lamnAct lamcapnAct TnAct] = ...
    bioSnyder(T, [x1 x2], params.kgain, Q, S, params, no_bnd);
meanAct = x1StatsAct.meth3(1);
mseAct = x1StatsAct.meth3(3);
clear lamnAct lamcapnAct TnAct

% Declare loop variables
lenTh = length(noiseTraits);
mseEst = zeros(lenTh, 1);
meanEst = zeros(lenTh, 1);
Test = cell(lenTh, 1);
distTest = zeros(lenTh, 1);
lenTest = zeros(lenTh, 1);

% Loop across several distortion settings
for iw = 1:lenTh
    % Obtain distorted photon times based on specified method
    Test{iw} = distortPhotons2(noiseTraits{iw}, Treal, T, params.kgain);
    lenTest(iw) = length(Test{iw});
    
    % Obtain spike metric between estimated and real photons with limit on
    % length to prevent out of memory errors
    if lenTest(iw) > metricPhotonLim
        idend2 = metricPhotonLim;
    else
        idend2 = lenTest(iw);
    end
    distTest(iw) = spkd(Treal(1:idend1), Test{iw}(1:idend2), dtcost);

    % Construct new data sets featuring the photon time estimates
    [Testnew Xestnew] = getXnewTnew(T, [x1 x2], Test{iw});

    % Run filter with estimated photon times as x2 events and assume that the
    % same Q matrix and params apply and then run with actual photon times
    [lamcapEst x1capEst lamstatsEst x1StatsEst lamnEst lamcapnEst TnEst] = ...
        bioSnyder(Testnew, Xestnew, params.kgain, Q, S, params, no_bnd);
    clear lamnEst lamcapnEst TnEst

    % Assign distorted Snyder outputs
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
out.distTest = distTest;
out.lenTreal = lenTreal;
out.lenTest = lenTest;
out.beta = inpSny.beta;
out.gamma = inpSny.gamma;