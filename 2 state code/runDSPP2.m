% This version calculates Q and uses Kolmogorov as dP = PQdt. Further this
% code uses a filter in which q postmultiplies the B matrix and transpsoses
% the dimensions of q0

% Major correction to use filterHybrid3 which ensures the ODE event updates
% only occur at the x2 event times rather than at all event times

% Allow for 3 separate x1 statistics calculations
% Removed old mse calculation, only the ODE point method is applied and an
% additional trapz method is also used
% Modification to include the full stats calculation including ODE points
% Modification to calculate the LVP linear encoding bound
% Modification to work with Markov chain case via boolean
% Modification to include other functional forms of lam = lam(x1), to store
% x1cap and to work with filterSnyderHybrid5

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
    % Assign simulation control parameters and initial population
    params.len = 2;
    params.N = 60000;
    params.plot_on = 0;
    params.Nstart = 30000;
    params.bulk_size = [1 1];
    simname = 'simF1';

    % Assign state space limits
    if markov
        params.Slim = [0 1];
        params.x0 = [params.Slim(1) params.Slim(1)];
    else
        params.x0 = [0 0];
    end

    % Assign mean and var of rates and types with fano factor
    params.avgR = [100 100];
    if ~markov
        params.birType = [1 3];
        params.deaType = [1 0];
    else
        params.birType = [7 3];
        params.deaType = [2 0];
    end

    % Assign rate control parameters including x1 death rate and proportional
    % gain for x2 rate dependence on x1 and x1 birth rate
    params.r = 0.1;  % <---- setting not of importance unless birfoll
    params.kdeath = 10;
    params.kgain = 15;
    params.kbirth = 10; %<--- alter to avgR(1) for const birth rate

    % Run the Gillespie simulations and extract the actual rate process and
    % counting observations and other SSA data
    tic;
    disp('Simulation started');
    if ~markov
        outGil = dsppFilterSim4(params);
    else
        outGil = dsppFilterMarkov(params);
    end
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
load('simF1');

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
lam = Xdot(1:lenT, 3);

% Obtain the state space of x1 and check dimensions
if ~ markov
    S = diag(0:2*max(x1)); %<---------- assumes a minimum space value of 0
else
    Slim = outGil.Slim;
    S = diag(Slim(1):Slim(2));
end
maxS = max(max(S));
lenS = length(S);

% Obtain transition matrix Q
birStr1 = outGil.birth{1};
[Q P] = getPQMx(birStr1, params.kbirth, params.kdeath, Slim);

% Calculate B matrix from functional form of lam = lam(x) with a suitable
% array of coefficients
birStr2 = outGil.birth{2};
coeff = [params.kgain 0]; % <---- MUST CORRESPOND TO SIMULATION
B = getBMx(Q, S, coeff, birStr2);

% Start with an arbitrary q0 as uniform and calculate lamcap0
q0 = (1/lenS)*ones(1, lenS);
[lamcap0 x1cap0] = calclamEst(birStr2, coeff, q0, S);

% Apply suitable inputs to the filter and extract outputs
disp('Calculating Snyder filter offline');
% outFil = filterHybrid2(lamcap0, q0, B, T, S, x1cap0, coeff, lenS, birStr2, params);
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
plotPosteriorIntensity5(Tn, q, lamn, lamcapn, x1n, x1capn, T, x2, S, kg, kb, kd, reldyn, save_sim, plot_on1);

% Obtain Markov chain statistics
if markov
    distr = 'gauss';
    [paramfit mfit vfit stat distr] = getMarkovStats(T, x1, plot_on2, distr);
end

% Obtain event statistics
eventStats = getRelMeanEventTimes(T, x1, x2);

% Save processed data if boolean set
if save_sim
    if any([kg kb kd] < 1)
        simname = 'filterSSAinv';
    else
        simname = ['filterSSA' num2str(kg) '_' num2str(kb) '_' num2str(kd)];
    end
    save(simname);
end
disp('Simulation and post processing complete');

% End profiling
if profile_on
    profile off
    profile report
end