% Reimagined version of runDSPP2 which is altered to accommodate general
% Markov chains with non-tridiagonal intensity matrices. Main changes are:
%   - Removed birth and death types in favour of molec and cross types
%   - Designed a more general Q matrix and included stationary solution
%   - Assumes proportional relation x2dot = alpha*x1
%   - All Gillespie simulations done via dsppMC code
%   - Altered the plotIntensity functions
%   - No bound calculated for general MCs

% Code to run doubly stochastic filters based on the concept of first
% obtaining the DSPP observations and then discretising
% clear all
clc
close all

% Set boolean to control the saving of data and figures and Markov chain
save_sim = 0;
markov = 1;

% Mode to skip SSA simulations if want to reuse data and to profile
skip_SSA = 0;
profile_on = 0;

% Set booleans for plots
plot_on1 = 1;
plot_on2 = 0;


%%
% Cell to perform Gillespie simulation to obtain observation of counts -
% this code is solely concerned with the observed process and provides the
% actual modulated intensity that is to be estimated

% Run SSA if boolean set else run filter with pre-existing simulation data
if ~skip_SSA
    
    % Assign state space limits - assumes min state is 0 <-----------------
    params.SlimSet.min = [0 0];
    if markov
        params.SlimSet.max = [10 inf];
    else
        params.SlimSet.max = [inf inf];
    end
    
    % Rate structure based inputs - bulk must be positive, all variables must
    % be same length and transit represents the incremental changes that bulk
    % reactions yield - odd reactions are births and even deaths - jumpMax
    % determines the setting for the maximum bulk size
    jumpMax = 3;
    switch(jumpMax)
        case 1
            params.reacType = [1 1 1 1 2];
            params.r_const = [10 10 5 5 1000];
            params.bulk = [1 1 2 2 1];
            params.transit = [1 -1 2 -2 0; 0 0 0 0 1];
            params.molecType = [1 1 1 1 2];
            params.crossType = [1 1 1 1 1];
        case 2
            params.reacType = [1 1 1 1 2];
            params.r_const = [10 10 5 5 1000];
            params.bulk = [1 1 2 2 1];
            params.transit = [1 -1 2 -2 0; 0 0 0 0 1];
            params.molecType = [1 1 1 1 2];
            params.crossType = [1 1 1 1 1];
            
        case 3
            params.reacType = [1 1 1 1 1 1 2];
            params.r_const = [20 20 10 10 5 5 1000];
            params.bulk = [1 1 2 2 3 3 1];
            params.transit = [1 -1 2 -2 3 -3 0; 0 0 0 0 0 0 1];
            params.molecType = [1 1 1 1 1 1 2];
            params.crossType = [1 1 1 1 1 1 1];
            
        otherwise
            error('Unsupported maximum jump size');
    end
    
    % Check that bulk matches the jumpMax setting and lengths are correct
    if max(params.bulk) ~= jumpMax
        error('Inconsistency between bulk and jumpMax');
    end
    if all([length(params.bulk) length(params.r_const) length(params.transit)...
            length(params.molecType) length(params.crossType)] ~= length(params.reacType))
        error('The lengths of the input simulation vectors do not match the number of reactions');
    else
        nReacs = length(params.reacType);
    end
    
    % Assign simulation control parameters and initial population
    params.plot_on = 1;
    params.Nstart = 30000;
    params.N = 50000;
    params.len = 2;
    params.x0 = [0 0];
    params.avgR = 100*ones(1, nReacs);
    params.maxR = 1000*ones(1, nReacs);
    params.minR = 0*ones(1, nReacs);
    simname = 'simMC';
    
    % Assuming proportional structure between x2dot and x1 <---------------
    params.kgain = params.r_const(end);
    params.coeff = [params.kgain 0];
    
    % Run the Gillespie simulations and extract the actual rate process and
    % counting observations and other SSA data
    tic;
    disp('Simulation started');
    outGil = dsppMC(params);
    disp('Simulation complete');
    
    % Log running time and display
    runtime = toc;
    normaliseTime(runtime);
    
    % Save data
    disp(['Saved file: ' simname]);
    save(simname);
end


%%
% Cell to process the observed process and extract the state space, prior
% probabilities and run the filter

% Clear workspace and load data file
% clear all
clc
load(simname);

% Set profiler
if profile_on
    profile on
end

% Obtain time and remove offset from the transient samples <---------------
T = outGil.T;
Tref = T(1);
T = T - Tref;

% Obtain molecular counts
lenT = length(T);
X = outGil.X;
x1 = X(1:lenT, 1);
x2 = X(1:lenT, 2);

% Obtain actual modulating intensity
Xdot = outGil.Xdot;
lam = Xdot(1:lenT, end); % <----------- assumes that the last rate is lam

% Obtain the state space of x1 and check dimensions
if ~ markov
    S = diag(0:2*max(x1)); %<---------- assumes a minimum space value of 0
else
    SlimSet = params.SlimSet;
    S = diag(SlimSet.min(1):SlimSet.max(1));
end
maxS = max(max(S));
lenS = length(S);

% Obtain transition matrix Q and stationary distribution with specification
% of type of CTMC and the x1 rate constants only <-------------------------
MCtype = 3;
xr_const = params.r_const(1:end-1);
[Q P Pi] = getPQGenMC(MCtype, xr_const, params.SlimSet, params.bulk);

% Calculate B matrix from functional form of lam = lam(x) with a suitable
% array of coefficients
birStr2 = 'cross'; % <------------------------ set dependence of lam on x1
coeff = params.coeff; % <----------------- set to correspone to birth types
B = getBMx(Q, S, coeff, birStr2);

% Start with an arbitrary q0 as uniform and calculate lamcap0
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

% Parameter to determine the relative speed of dynamics between x2 and x1
rateMean = outGil.rate.mean;
x1rateMean = max(rateMean(1:2));
x2rateMean = rateMean(3);
reldyn = x2rateMean/x1rateMean;
disp(['The relative speed of x2 to x1 is: ' num2str(reldyn)]);

% Calculate linear encoding LVP bound in special case of nearest neighbour
if length(xr_const) == 2
    umean = outGil.rate.mean(1);
    x1mean = outGil.molec.mean(1);
    alpha = params.kgain;
    kdeath_x1 = params.r_const(2);
    [relbnd rawbnd N1 N2] = getLVPLinearBound(alpha, x1mean, kdeath_x1, umean);
    Nratio = N2/N1;
    no_bnd = 0;
else
    no_bnd = 1;
    disp('No bound specified for general multi-bulk reaction based MCs');
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

% Obtain plots of posterior density, intensity estimates and KL divergence
plotPosteriorMC(Tn, q, lamn, lamcapn, x1n, x1capn, T, x2, S, params.kgain, reldyn, save_sim, plot_on1);

% Obtain Markov chain statistics
if markov
    distr = 'gauss';
    [paramfit mfit vfit stat distr] = getMarkovStats(T, x1, plot_on2, distr);
end

% Obtain event statistics - only for simple MCs
if MCtype ~= 3
    eventStats = getRelMeanEventTimes(T, x1, x2);
end

% Acquire and plot distribution of times in states
plotState = 1;
prob = getStateDistr(T, x1, plotState);

% A comparative plot of the holding time empirical distribution and the
% known MC stationary distribution
figure;
plot(prob.state, [Pi' prob.val'], 'o-');
xlabel('state values');
ylabel('state probabilities');
legend('theoretical', 'empirical', 'location', 'best');
title('Comparison of empirical holding with MC stationary distribution');

% Save processed data if boolean set
if save_sim
    if params.kgain < 1
        simname = 'filterSSAinv';
    else
        simname = ['filterSSA' num2str(params.kgain)];
    end
    save(simname);
end
disp('Simulation and post processing complete');

% End profiling
if profile_on
    profile off
    profile report
end