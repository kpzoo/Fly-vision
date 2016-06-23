% Problem: not sure if when using the estimated photon times if it makes
% sense to reuse Q or if it is possible to estimate Q

% Function that takes a Gillespie simulation and runs the Snyder filter
% which follows from the runDSPPMCQ2 code and is adapted for mapQBPh
function [lamcap x1cap lamstats x1Stats lamn lamcapn Tn] = bioSnyder(T, X, kgain, Q, S, params, no_bnd)

% Obtain actual intensity, proportional assumption <-------------------
lenT = length(T);
x1 = X(1:lenT, 1);
x2 = X(1:lenT, 2);
lam = kgain*x1;

% Calculate B matrix from functional form of lam = lam(x) with a suitable
% array of coefficients
birStr2 = 'cross'; % <------------------------ set dependence of lam on x1
coeff = [kgain 0]; % <----------------- set to correspond to birth types
B = getBMx(Q, S, coeff, birStr2);

% Start with an arbitrary q0 as uniform and calculate lamcap0
lenS = length(S);
q0 = (1/lenS)*ones(1, lenS);
[lamcap0 x1cap0] = calclamEst(birStr2, coeff, q0, S);

% Apply suitable inputs to the filter and extract outputs
disp('Calculating Snyder filter offline');
outFil = filterHybrid3(lamcap0, q0, B, T, S, x1cap0, coeff, lenS, birStr2, params, x2);
lamcap = outFil.lamcap;
q = outFil.q;
t = outFil.t;
x1cap = outFil.x1cap;
qODE = outFil.qODE;
clc
disp('Finished estimate of modulating intensity');

% Obtain qODE stats and display
sumqODE = sum(qODE, 2);
mqODE = mean(sumqODE);
vqODE = var(sumqODE);
disp(['Mean and var of qODE = ' [num2str(mqODE) ' ' num2str(vqODE)]]);

% Check that the probability densities and are sensible and that tcum
% matches the T vector
if max(abs(sum(q, 2) - 1)) > 10^-9
    error('The probability density sum is not close enough to 1');
end
if min(min(q)) < 0
    error('The probability density has negative values');
end

% Obtain parameters for plots and statistics
Tset = outFil.Tset;
Qset = outFil.Qset;
lamcapFull = outFil.lamcapFull;
rem_trans = 0;
x1Set = cell(1, 1);
for i = 1:length(lamcapFull)
    x1Set{i} = lamcapFull{i}/params.kgain;
end

% assignin('base', 'Tset', Tset);
% assignin('base', 'Qset', Qset);
% assignin('base', 'x1Set', x1Set);

% Obtain 3 estimates of x1 statistics and lamstats via 1 method
[lamn lamcapn Tn lamstats] = getIntensityStatsFull2(lam, lamcapFull, Tset, T, Qset, rem_trans);
if no_bnd
    [x1Stats x1n x1capn Tn] = getAllStats(x1, x1Set, Tset, T, Qset, rem_trans, -1);
else
    [x1Stats x1n x1capn Tn] = getAllStats(x1, x1Set, Tset, T, Qset, rem_trans, rawbnd);
    disp('********************************************************************');
    disp(['The ratio of mse to bound is ' num2str(x1Stats.psi(1))]);
    disp(['The ratio of N2 to N1 is ' num2str(Nratio)]);
    disp('********************************************************************');
end