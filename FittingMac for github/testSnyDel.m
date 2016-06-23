% Simple function to test the pure delayed Snyder on fresh data
clear all
clc
close all

% Specify location inputs
locfolder1 = 'doubleRange';
locfolder2 = 'test2';
savename = 'linSny2';

% Load simulation data
cd(locfolder1);
cd(locfolder2);
files = dir('*.mat');
if length(files) == 1
    load(files.name);
else
    error('Code can only handle 1 file');
end
cd ..
cd ..

% Assign data (with normalisation of non-stationary cases)
X = outGil.X;
T = outGil.T;
T = T - T(1);
X(:, 2) = X(:, 2) - X(1, 2);
X(:, 3) = X(:, 3) - X(1, 3);
x1 = X(:, 1);
x2 = X(:, 2);
x3 = X(:, 3);
S = diag(params.Slim);
k = params.kdeath;
Q1 = [-k k; k -k];

% Limit the photon stream based on limPhoton
limPhoton = 100;
if max(x2) < limPhoton
    warning('Mat:fewPh', ['There are less than ' num2str(limPhoton) ' photons present ']);
else
    disp(['Clipping Gillespie data to ' num2str(limPhoton) ' photons']);
    idlim = find(x2 == limPhoton, 1, 'last');
    x2 = x2(1:idlim);
    x1 = x1(1:idlim);
    x3 = x3(1:idlim);
    T = T(1:idlim);
    X = X(1:idlim, :);
end

% Run normal Snyder filter on real photon times
no_bnd = 1;
[lamstatsAct x1StatsAct x1nAct x1capnAct TnAct] = ...
    bioSnyderMod(T, [x1 x2], params.kgain, Q1, S, params, no_bnd);
if ~isnan(x1StatsAct.meth3(3))
    meanAct = x1StatsAct.meth3(1);
    mseAct = x1StatsAct.meth3(3);
else
    disp('Method 3 failed for MSE calculation');
    meanAct = x1StatsAct.meth2(1);
    mseAct = x1StatsAct.meth2(3);
end
assignin('base', 'TnAct', TnAct);
assignin('base', 'x1nAct', x1nAct);
assignin('base', 'x1capnAct', x1capnAct);
assignin('base', 'x1StatsAct', x1StatsAct);


% Extract rate parameters with a = alpha, b = eta for a single delay stage
% and obtain state space
a = params.kgain;
b = params.eta;
k = params.kdeath;
Sx = diag(min(x1):max(x1));
Sy = diag(min(x2):max(x2));

% Ensure only a 2 state x input is applied and set boolean to change
% between integration limits and solution type
if length(diag(Sx)) ~= 2
    error('Code cannot handle information processes (x) with more (or less) than 2 states');
end
int_sw = 0;
linearODE = 1;

% Initialise appropriate Q matrix with suitable dimensions
dimQ = length(diag(Sx))*length(diag(Sy));
Q = zeros(dimQ, dimQ);
rep = dimQ/2; % must be even by assumption on Sx <------------------------

% Construct the Q matrix diagonals
d1 = repmat([k 0], 1, rep);
d1 = d1(1:dimQ - 1);
d1neg = d1;
d2 = repmat([0 a], 1, rep);
d2 = d2(1:dimQ - 2);

% Obtain the Q matrix (lam matrix now in ode formulation)
Q = Q + diag(d1, 1) + diag(d2, 2) + diag(d1neg, -1);
d0 = -sum(Q, 2);
Q = Q + diag(d0);

% Start with an arbitrary q0 as uniform and calculate lamcap0, x1cap0
% q0 = (1/dimQ)*ones(1, dimQ);
q0 = [0 1 zeros(1, dimQ-2)];
lamDiag0 = zeros(1, dimQ);
iV = 3;
for iL = 2:dimQ/2
    lamDiag0(iV:iV+1) = b*((iL-1) - X(1, end));
    iV = iV + 2;
end
lam0 = diag(lamDiag0);
lamcap0 = q0*lam0*ones(dimQ, 1);
x1cap0 = sum(q0(2:2:end)); % <------------ assumes 2 state x MC


% Obtain inter-event times of point process observations
[Tevent perc] = getEventTimes(T, x3, 'birth');
[Tevent2 percy] = getEventTimes(T, x2, 'birth');
len = length(Tevent);
dt = [0; diff(Tevent)];

% Check closeness of event times of y and z
normEv = norm(Tevent2 - Tevent);
disp(['The event time norm is ' num2str(normEv)]);


% Declare variables and initialise
q = zeros(len, dimQ);
q(1, :) = q0;
t = zeros(len, 1);
lamcap = zeros(len, 1);
lamcap(1) = lamcap0;
x1cap = zeros(len, 1);
x1cap2 = zeros(len, 1);
x1cap(1) = x1cap0;
qODE = zeros(len, dimQ);
qODE(1, :) = q0;
qpert = zeros(len, dimQ);
qpert2 = zeros(len, dimQ);
betaSum0 = zeros(len, 1);
betaSum1 = zeros(len, 1);
betaSumRatio = zeros(len, 1);

% Cell to save output of ODE solver and set options
Qset = cell(1, 1);
DQset = cell(1, 1);
Tset = cell(1, 1);
x1Set = cell(1, 1);

% Loop across dt values integrating the Snyder equations with the ode solver
% and then applying jump corrections due to observed events
for i = 2:len
    
    % Obtain lam matrix
    lamDiag = zeros(1, dimQ);
    iV = 3;
    for iL = 2:dimQ/2
        lamDiag(iV:iV+1) = b*((iL-1) - X(i, end));
        iV = iV + 2;
    end
    lam = diag(lamDiag);
    
    if linearODE
        % Linear ODE solution between events with normalisation
        [tset qset] = ode15s(@(ts, y) odeCompSnyderNonStatLinear(ts, y, Q, lam),...
            [Tevent(i-1) Tevent(i)], q(i-1, :)');
        for j = 1:size(qset, 1)
            qset(j, :) = qset(j, :)/(sum(qset(j, :)));
        end
    else
        % Non-linear ODE solution between events
        if ~int_sw
            % Delayed numerical solution with integration limits in absolute time
            [tset qset] = ode23tb(@(ts, y) odeCompSnyderNonStat(ts, y, Q, X(i, end), dimQ, b, lam),...
                [Tevent(i-1) Tevent(i)], q(i-1, :)');
        else
            % Delayed numerical solution with integration limits in relative time
            [tset qset] = ode15s(@(ts, y) odeCompSnyderNonStat(ts, y, Q, X(i, end), dimQ, b, lam),...
                [0 dt(i)], q(i-1, :)');
            tset = tset + Tevent(i-1);
        end
    end
    
    % Obtain the probability rates - actual conditional probabilities for
    % non-linear case and unnormalised rates for the linear case
    dqset = zeros(size(qset));
    for j = 1:size(qset, 1)
        if linearODE
            dqset(j, :) = odeCompSnyderNonStatLinear(0, qset(j, :)', Q, lam);
        else
            dqset(j, :) = odeCompSnyderNonStat(0, qset(j, :)', Q, X(i, end), dimQ, b, lam);
        end
    end
    
    % Assign output value at event times
    q(i, :) = qset(end, :);
    t(i) = tset(end);
    qODE(i, :) = q(i, :);
    Qset{i} = qset;
    DQset{i} = dqset;
    Tset{i} = tset;
    x1Set{i} = sum(qset(:, 2:2:end), 2); % <----------- assumes that x is 2 state
    
    % Check raw ODE q output and derivatives
    if any(q(i, :) < -10^-8)
        assignin('base', 'qODErr', q(i, :));
        disp(['qODE distribution has negative entries at i =' num2str(i)]);
    end
    if max(abs(sum(q(i, :)) - 1)) > 10^-4
        assignin('base', 'qODErr', q(i, :));
        disp(['qODE distribution sums to ' num2str(sum(q(i, :))) ' at i = ' num2str(i-1)]);
        assignin('base', 'qset', qset);
    end
    
    if max(abs(sum(dqset(end, :)))) > 10^-4
        assignin('base', 'dqODErr', dqset(end, :));
        disp(['dqODE distribution sums to ' num2str(sum(q(i, :))) ' at i = ' num2str(i-1)]);
        assignin('base', 'dqset', dqset);
    end
    
    % Calculate perturbation on q due to jump and calculate lamcap
    q(i, :) = q(i, :)*lam/(q(i, :)*lam*ones(dimQ, 1));
    lamcap(i) = q(i, :)*lam*ones(dimQ, 1);
    qnew = q(i, :);
    x1cap(i) = sum(qnew(2:2:end));
    qpert(i, :) = qnew;
    
    
    % Calculate another qpert
    lamDiag = zeros(1, dimQ);
    iV = 3;
    for iL = 2:dimQ/2
        lamDiag(iV:iV+1) = b*((iL-1) - X(i+1, end));
        iV = iV + 2;
    end
    lam = diag(lamDiag);
    qnew2 = q(i, :)*lam/(q(i, :)*lam*ones(dimQ, 1));
    x1cap2(i) = sum(qnew2(2:2:end));
    qpert2(i, :) = qnew2;
    q(i, :) = qnew2; % <------- use an updated lam matrix to get the IC for the next step
    
    % Obtain the beta sum values
    Dlam = diag(lam);
    betaSum0(i) = qset(end, 1:2:end)*Dlam(1:2:end);
    betaSum1(i) = qset(end, 2:2:end)*Dlam(2:2:end);
    betaSumRatio(i) = betaSum1(i)/(betaSum1(i) + betaSum0(i));
    
    % Display progress and debug set of results from ODE solver
    disp(['Finished iteration: ' num2str(i-1) ' of ' num2str(len)]);
end

% Obtain qODE stats and display
disp('Finished estimate of modulating intensity');
sumqODE = sum(qODE, 2);
mqODE = mean(sumqODE);
vqODE = var(sumqODE);
disp(['Mean and var of qODE = ' [num2str(mqODE) ' ' num2str(vqODE)]]);

% Obtain 3 estimates of x1 statistics (with the appended x1 stream applied)
rem_trans = 0;
[x1Stats x1n x1capn Tn] = getAllStats(x1, x1Set, Tset, T, Qset, rem_trans, -1);


% Quick check of minimum solution if eta high (result should match normal
% Snyder for high enough eta or small enough delay)
if b >= 100
    beta = a/k;
    A = 1/beta + 0.5;
    B = sqrt(1/(beta^2) + 0.25);
    steady.minTheo = A - B;
    steady.minEmp = min(x1capn);
    assignin('base', 'steady', steady);
end

% Plot comparison of solution trajectories
figure;
plot(TnAct, x1capnAct, Tn, x1capn);
xlabel('time');
ylabel('x1 estimates');
title(['Comparison of pure Snyder and delayed Snyder when the eta = ' num2str(params.eta)]);
legend('pure Snyder', 'delayed Snyder');