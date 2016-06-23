% This version calculates Q and uses Kolmogorov as dP = PQdt. Further this
% code uses a filter in which q postmultiplies the B matrix and transpsoses
% the dimensions of q0

% Major correction to use filterHybrid3 which ensures the ODE event updates
% only occur at the x2 event times rather than at all event times

% Added two further mse calculations for comparison
% Removed old mse calculation, only the ODE point method is applied
% Modification to include the full stats calculation including ODE points
% Modification to calculate the LVP linear encoding bound
% Modification to work with Markov chain case via boolean
% Modification to include other functional forms of lam = lam(x1), to store
% x1cap and to work with filterSnyderHybrid5

% Code to run doubly stochastic filters based on the concept of first
% obtaining the DSPP observations and then discretising
function outB = runDSPP2Fn(inpB)

% Set boolean to control the saving of data and figures and Markov chain
save_sim1 = inpB.save_sim(1);
save_sim2 = inpB.save_sim(2);
markov = inpB.markov;

% Mode to skip SSA simulations if want to reuse data and to profile
skip_SSA = inpB.skip_SSA;
profile_on = inpB.profile;

% Set booleans for plots
plot_on1 = inpB.plot_on(1);
plot_on2 = inpB.plot_on(2);
plot_on3 = inpB.plot_on(3);

%%
% Cell to perform Gillespie simulation to obtain observation of counts -
% this code is solely concerned with the observed process and provides the
% actual modulated intensity that is to be estimated

% Run SSA if boolean set else run filter with pre-existing simulation data
if ~skip_SSA
    % Assign simulation control parameters and initial population
    params.len = 2;
    params.N = inpB.Nset(2);
    params.plot_on = plot_on3;
    params.Nstart = inpB.Nset(1);
    params.x0 = [0 0];
    params.bulk_size = [1 1];
    simname = inpB.simname;

    % Assign state space limits (for Markov case)
    if markov
        params.Slim = inpB.Slim;
        params.x0 = [params.Slim(1) params.Slim(1)];
    end

    % Assign mean and var of rates and types with fano factor
    params.avgR = inpB.avgR;
    params.birType = inpB.birType;
    params.deaType = inpB.deaType;

    % Assign rate control parameters including x1 death rate and
    % coefficient set relating x2 intensity to x1 events
    coeff = inpB.coeff;
    params.r = 0;
    params.kdeath = inpB.kdeath;
    params.kgain = coeff(1);
    params.kbirth = inpB.kbirth; %<--- alter to avgR(1) for const birth rate

    % Run the Gillespie simulations and extract the actual rate process and
    % counting observations and other SSA data
    tic;
    disp('Gillespie simulation started');
    if ~markov
        outGil = dsppFilterSim4(params);
    else
        outGil = dsppFilterMarkov(params);
    end
    disp('Gillespie simulation complete');

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
tic;
T = outGil.T;
Tref = T(1);
T = T - Tref;

% Obtain molecular counts
disp('Filter calculation started');
lenT = length(T);
X = outGil.X;
x1 = X(1:lenT, 1);
x2 = X(1:lenT, 2);

% Obtain actual modulating intensity
Xdot = outGil.Xdot;
lam = Xdot(1:lenT, 3);

% Obtain the state space of x1 and check dimensions
if ~ markov
    S = diag(0:2*max(x1)); %<---------- assumes a minimum space value of 0
else
    Slim = outGil.Slim;
    S = diag(Slim(1):Slim(2));
end
lenS = length(S);

% Obtain transition matrix Q
birStr1 = outGil.birth{1};
[Q P] = getPQMx(birStr1, params.kbirth, params.kdeath, Slim);

% Calculate B matrix from functional form of lam = lam(x) with a suitable
% array of coefficients
birStr2 = outGil.birth{2};
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

% Calculate linear encoding LVP bound
umean = outGil.rate.mean(1);
x1mean = outGil.molec.mean(1);
alpha = params.kgain;
kdeath_x1 = outGil.r_const(2);
[relbnd rawbnd N1 N2] = getLVPLinearBound(alpha, x1mean, kdeath_x1, umean);
Nratio = N2/N1;

% Obtain parameters for plots and statistics
Tset = outFil.Tset;
Qset = outFil.Qset;
lamcapFull = outFil.lamcapFull;
rem_trans = 0;
kg = params.kgain;
kb = params.kbirth;
kd = params.kdeath;
x1Set = cell(1, 1);
for i = 1:length(lamcapFull)
    x1Set{i} = lamcapFull{i}/kg;
end

% Obtain 3 estimates of x1 statistics and lamstats via 1 method
[lamn lamcapn Tn lamstats] = getIntensityStatsFull2(lam, lamcapFull, Tset, T, Qset, rem_trans);
[x1Stats x1n x1capn Tn] = getAllStats(x1, x1Set, Tset, T, Qset, rem_trans, rawbnd);
disp('********************************************************************');
disp(['The ratio of mse to bound is ' num2str(x1Stats.psi(1))]);
disp(['The ratio of N2 to N1 is ' num2str(Nratio)]);
disp('********************************************************************');

% Obtain plots of posterior density, intensity estimates and KL divergence
plotPosteriorIntensity5(Tn, q, lamn, lamcapn, x1n, x1capn, T, x2, S, kg, kb, kd, reldyn, save_sim2, plot_on1);

% Obtain Markov chain statistics
if markov
    distr = 'gauss';
    [paramfit mfit vfit stat distr] = getMarkovStats(T, x1, plot_on2, distr);
    fitSet.param = paramfit;
    fitSet.meanVar = [mfit vfit];
    fitSet.stat = stat;
    fitSet.distr = distr;
end

% Obtain event statistics
eventStats = getRelMeanEventTimes(T, x1, x2);

% Log running time and display
disp('Filter calculation complete');
runtime = toc;
normaliseTime(runtime);

% Save processed data if boolean set
if save_sim1
    save(simname);
end
disp('Simulation and post processing complete');

% Assign output from filter
outB.lamstats = lamstats;
outB.x1Stats = x1Stats;
outB.lamn = lamn;
outB.lamcapn = lamcapn;
outB.Tn = Tn;
outB.Tevent = outFil.Tevent;
outB.x1cap = x1cap;
outB.lamcap = lamcap;
outB.reldyn = reldyn;
outB.birth = outGil.birth;
outB.death = outGil.death;
outB.relbnd = relbnd;
outB.rawbnd = rawbnd;
outB.N1N2 = [N1 N2];
outB.Nratio = Nratio;
outB.fitSet = fitSet;
outB.filterSettings = outFil.settings;
outB.eventStats = eventStats;

% End profiling
if profile_on
    profile off
    profile report
end