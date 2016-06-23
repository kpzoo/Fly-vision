% Script to batch across the multi-state MCs created in runDSPPMCQ3 with
% loop parameter of interest as alpha
clear all
clc
close all

% Define range of beta and a fixed gamma and obtain rate parameters
beta = [1 10:10:100];
gamma = 5;
k = 1/(100*gamma);
alpha = beta*k;
lenb = length(beta);
simroot = 'multiSim';

% Declare output variables
paramSet = cell(1, lenb);

% Loop across settings
for ib = 1:lenb
    % Set name for saved data
    simname = [simroot num2str(ib)];
    
    % Main code for running SSA
    paramSet{ib} = runDSPPMCQ3(alpha(ib), k, simname);
    
    % Display progress
    disp(['Completed run: ' num2str(ib) ' of ' num2str(lenb)]);
end
