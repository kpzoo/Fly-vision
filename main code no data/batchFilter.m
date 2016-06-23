% Script to run a batch of filter simulations across various parameters and
% log the compound results
clear all
clc
close all


% Set batch numbers to determine simulations
batchNo = [0 2];

% Set relevant booleans and inputs to control code
inpB.save_sim = [0 0];
inpB.markov = 1;
inpB.skip_SSA = 0;
inpB.profile = 0;
inpB.plot_on = [1 1 0];
inpB.simname = 'simF1';

% Set some initial simulation parameters
inpB.kbirth = 10;
inpB.kdeath = 10;
inpB.Nset = [20000 21000];
inpB.Slim = [0 10];

% Set initial birth and death rates and types
inpB.coeff = [100 0];
inpB.avgR = [inpB.coeff(1) inpB.coeff(1)];
inpB.birType = [7 3];
inpB.deaType = [2 0];


%%
% Cell to batch simulate a set encoder for various gain values and
% then at various state space sizes

% Set encoding function and birth-death forms for specific case as well as
% loop controlling variable
inpB.birType = [7 3];
inpB.deaType = [2 0];
% kgainSet = logspace(-1, 4, 2);
kgainSet = [0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 750 1000];
lenB1 = length(kgainSet);
% N = ceil(linspace(1000, 10000, lenB1));
N = [1000 1000 1000 1000 1000 1500 1500 1500 1500 2000 2000 2000 2000 3000 3000 3000 6000 6000];

% Declare loop variables
msex1a = zeros(1, lenB1);

% Loop across the kgain = coeff(1) parameter defined as the set below and
% alter N value based on gain
if any(batchNo == 1)
    disp(['Starting batch simulation: ' 1]);
    for i = 1:lenB1
        % Alter coeff based on gain setting and run filter
        inpB.coeff(1) = kgainSet(i);
        inpB.Nset(2) = inpB.Nset(1) + N(i);
        outB = runDSPPfilterFn(inpB);
        
        % Assign outputs and disp progress
        msex1a(i) = outB.x1stats.mseErr;
        disp(['Finished iteration: ' num2str(i) ' of ' num2str(lenB1)]);
    end
    disp(['Completed batch simulation: ' 1]);
    
    % Plot results
    figure;
    plot(kgainSet, msex1a, 'bo-');
    xlabel('kgain');
    ylabel('mse of x1 or normalised lam mse');
    title(['Evolution of mse with kgain for ' outB.birth{2} ' at Smax = ' num2str(inpB.Slim(2))]);
    saveas(gcf, 'batch1.fig');
end

% Set changed booleans and control inputs for batch across state space
% SmaxSet = 1:20;
SmaxSet = 100;
lenB2 = length(SmaxSet);
inpB.Nset = [20000 80000];
inpB.coeff(1) = 1000;

% Declare loop variables
msex1b = zeros(1, lenB2);

% Loop across the kgain = coeff(1) parameter defined as the set below and
% alter N value based on gain
if any(batchNo == 2)
    disp(['Starting batch simulation: ' 2]);
    for i = 1:lenB2
        % Alter Smax based on space setting and run filter
        inpB.Slim(2) = SmaxSet(i);
        outB = runDSPPfilterFn(inpB);
        
        % Assign outputs and disp progress
        msex1b(i) = outB.x1stats.mseErr;
        disp(['Finished iteration: ' num2str(i) ' of ' num2str(lenB2)]);
    end
    disp(['Completed batch simulation: ' 2]);
    
    % Plot results
    figure;
    plot(SmaxSet, msex1b, 'bo-');
    xlabel('Smax');
    ylabel('mse of x1 or normalised lam mse');
    title(['Evolution of mse with Smax for ' outB.birth{2} ' at kgain = ' num2str(inpB.coeff(1))]);
    saveas(gcf, 'batch2.fig');
end