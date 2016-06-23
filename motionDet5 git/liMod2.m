% Code to model light based on the moving train positional markov chain
% concept which aims to maintain a markov formulation overall
clear all
clc
close all

% Booleans to control code
visual1 = 0;
visual2 = 0;


%% Main code to run the stochastic simulations

% Set the space for possible position and intensity values (note that
% the position space is twice the actual no. positions as the duplication
% of states is used to represent different directions
posSpace = 0:5;
intenSpace = 1;
Nev = 20000;
nPos = length(posSpace)/2;
I = zeros(Nev, length(intenSpace), nPos);

% Set the transition rates of the infinitesimal generator
a = 0.1;
b = a/10;
c = a;
inten = 10*a;

% Obtain the Q matrix for a 6 state (3 position) stimulus and the state
% distribution at equilibirum
if length(posSpace) == 6
    Qnodiag = [0 a 0 0 0 0;
        0 0 a b 0 0;
        c 0 0 0 b 0;
        0 b 0 0 0 c;
        0 0 b a 0 0;
        0 0 0 0 a 0];
    Q = Qnodiag - diag(sum(Qnodiag, 2));
    Pi = null(Q');
    Pi = Pi/sum(Pi);
else
    error('No code accounting for more complex light models');
end

% Obtain the number of reactions from the max jump in state found in Q -
% assumes a state structure of [0, 1, ... , m]
dimQ = length(Q);
jumpMax = 0;
for i = 1:dimQ
    % Check each row for highest difference between position of non-zero
    % rates and the current state
    Qrow = Q(i, 1:dimQ);
    idjump = find(Qrow > 0);
    currjumpMax = max(abs(idjump - i));
    jumpMax = max(currjumpMax, jumpMax);
end

% Assume 3 (based on pixels) poisson process modulated by reactions in Q and
% obtain the transit matrix based on the jumps assuming 1:jumpMax exists
nReacs = 2*jumpMax + nPos;
bulk = [1:jumpMax 1:jumpMax];
bulk = sort(bulk);
transitx1 = [bulk zeros(1, nPos)];
transitx1(2:2:end) = -transitx1(2:2:end);

% Calculate the further transit rows - one for each pixel
transitxn = [zeros(nPos, nReacs - nPos) eye(nPos)];
transit = [transitx1; transitxn];

% Set the input parameters to the Gillespie code
inpGill.N = Nev;
inpGill.Nstart = 1;
inpGill.len = 4;
inpGill.Q = Q;
inpGill.transit = transit;
inpGill.nReacs = nReacs;
inpGill.alpha = inten;
inpGill.x0 = zeros(1, inpGill.len);
inpGill.nPos = nPos;

% Run the Gillespie algorithm for the Q matrix specified
[X, ~, T] = gillespieManyReacsQSimple(inpGill);
pos = X(:, 1);
phot = X(:, 2:end);
clear Qrow Qnodiag X

% Obtain statistics of markov position and compare to theoretical
dT = diff(T);
posStats.posMean = sum(dT.*pos(1:end-1))/sum(dT);
pos2 = pos.*pos;
posStats.posVar = sum(dT.*pos2(1:end-1))/sum(dT);
posStats.posMeanTheo = sum(posSpace*Pi);
posSpace2 = posSpace.*posSpace;
posStats.posVarTheo = sum(posSpace2*Pi);


%% Main code to perform Snyder filtering on a vector of observed processes

% Set various rate matrices based on assumed form of state MC <------------
lamDiagIndiv = zeros(nPos, length(posSpace));
posZeros = zeros(1, length(posSpace));
for i = 1:nPos
    % Assume the coding of position and direction via a state then a
    % counting process represents 2 states e.g. 0 and 3 for space 0:5
    lamDiagIndiv(i, [i nPos+i]) = inten;
end
lamDiag = sum(lamDiagIndiv, 1);

% Set the input structures to the filtering function
lamInp.lamDiagIndiv = lamDiagIndiv;
lamInp.lamDiag = lamDiag;


% Obtain the times of the various photons from different pixels
Tphot = cell(1, nPos);
idphot = cell(1, nPos);
nEvents = 0;
for i = 1:nPos
    [Tphot{i}, ~, idphot{i}] = getEventTimes2(T, phot(:, i), 'birth');
    nEvents = nEvents + length(Tphot{i});
end

% Obtain the full event set and ensure the times are increasing
Tevent = zeros(nEvents, 1);
offset = 0;
for i = 1:nPos
    if i > 1
        offset = length(idphot{i-1}) + offset;
    end
    range = 1:length(idphot{i});
    range = range + offset;
    Tevent(range) = T(idphot{i});
end
Tevent = sort(Tevent);

% Obtain position state space matrix and set initial conditions
Sx = diag(posSpace);
q0 = ones(1, length(posSpace))/length(posSpace);
x1cap0 = sum(q0*Sx, 2);

% Run Snyder filtering for multiply observed point processes
outFil = generalSnyFilterBasic(q0, x1cap0, Q, lamInp, nPos, Tphot, Sx, nEvents, Tevent);

% Obtain 3 estimates of x1 statistics (with the appended x1 stream applied)
rem_trans = 0;
[x1Stats x1n x1capn Tn] = getAllStats(pos, outFil.x1capset, outFil.Tset, T, outFil.Qset, rem_trans, -1);
if ~isnan(x1Stats.meth3(3))
    meanCalc = x1Stats.meth3(1);
    mseCalc = x1Stats.meth3(3);
else
    disp('Method 3 failed for MSE calculation');
    meanCalc = x1Stats.meth2(1);
    mseCalc = x1Stats.meth2(3);
end

% Obtain directional estimates from the positional schemes
[direc dis0] = getDirecfromPos(pos, posSpace);
[direcReal dis1] = getDirecfromPos(x1n, posSpace);
[direcEst dis2] = getDirecfromPos(x1capn, posSpace);

% Calculate a MSE index between the directional estimates
direcStats = getStatsAltMeths2(Tn, direcReal, direcEst, direc, T);



%% Code to visualise the system

if visual1
    % Visualisation of the photon stream (only 100 events)
    figure;
    for i = 1:nPos
        hold all
        stem(Tphot{i}(1:100), ones(size(Tphot{i}(1:100))));
        xlabel('photon times');
    end
end

if visual2
    % Visualise the estimated and actual state trajectories
    figure;
    plot(Tn, x1n, Tn, x1capn);
    xlabel('time');
    ylabel('state');
    
    % Visualise directional estimates
    figure;
    stairs(Tn, direcReal, 'b');
    hold on
    stairs(Tn, direcEst, 'r');
    hold off
    xlabel('time');
    ylabel('state');    
end
