% Simple script to run Snyder filter only on a set of Gillespie data
clear all
clc
close all

% Obtain data for Snyder processing - in main folder
files = dir('*.mat');
flen = length(files);

% Start timing
tic;

% Initialise necessary inputs
inpB.save_sim = [0 0];
inpB.markov = 1;
inpB.skip_SSA = 1;
inpB.skip_Snyder = 0;
inpB.profile = 0;
inpB.plot_on = [0 0 0];

% Declare output storage variables
N1 = zeros(1, flen);
N2 = zeros(1, flen);
x1STATS = cell(1, flen);
mseA = zeros(1, flen);
mseB = zeros(1, flen);
mseC = zeros(1, flen);
dtightemp = zeros(1, flen);
psi_tightemp = zeros(1, flen);
lamMean = zeros(1, flen);
kbirth = zeros(1, flen);
kdeath = zeros(1, flen);
kgain = zeros(1, flen);

% Run Synder filter on each file in turn
for i = 1:flen
    % Set file name
    inpB.simname = files(i).name;
    disp(['Processing' inpB.simname]);
    
    % Run filter and extract outputs
    outB = runDSPP2Fn(inpB);
    N1(i) = outB.N1N2(1);
    N2(i) = outB.N1N2(2);
    mseA(i) = outB.x1Stats.meth1(3);
    mseB(i) = outB.x1Stats.meth2(3);
    mseC(i) = outB.x1Stats.meth3(3);
    x1STATS{i} = outB.x1Stats;
    dtightemp(i) = outB.rawbnd;
    psi_tightemp(i) = outB.relbnd;
    lamMean(i) = outB.x2rateMean;
    kbirth(i) = outB.kbirth;
    kdeath(i) = outB.kdeath;
    kgain(i) = outB.kgain;
    
    clear outB
end

% Save data
save('snydBat.mat');

% Check for symmetric MC
if ~all(kbirth == kdeath)
    warning('Mat:notSymm', 'Non-symmetric MC so theoretical tight and all slack bounds invalid');
else
    % Calculate empirical version of slack bound and theoretical versions
    % of both bounds
    Nratio = N2./N1;
    dslackemp = 1./(Nratio*log(2) + 4);
    beta = kgain/kbirth(1);
    dslack = 1./(beta*log(2) + 4);
    dtight = 1./(1 + sqrt(1 + 4*beta));
    
    % Plot dslack against mseSnyder
    figure;
    plot(lamMean, mseA, lamMean, dslack);
    xlabel('mean intensity');
    ylabel('mse and bound');
    title('Comparison of mse and slack bound');
    legend('mse Snyder', 'theoretical slack bound', 'location', 'best');
end

% Log time
t_run = toc;
normaliseTime(t_run);