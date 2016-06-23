% Script to double batch Gillespie simulations with this outer loop aimed
% on altering kbirth in a 2 state symmetric MC and the inner loop aimed at
% altering kgain but ensuring a fixed beta_ratio across different kbirth
clear all
clc
close all

% Set controlling inputs
kbirthinv = 100:200:2000;
kbirthSet = 1./kbirthinv;
lendb = length(kbirthSet);
simInit = 'simFil';

% Special input to control the beta range and start timer
beta_ratio = 1:2:50;
tic;

% Loop the Gillespie simulation set
for i = 2:lendb
    
    % Run Gillespie set with updated inputs
    kbirth = kbirthSet(i);
    simroot = [simInit num2str(i)];
    batchFilter2Fn(kbirth, simroot, beta_ratio, i);
    
    % Display progress
    disp('****************************************************************');
    disp(['Finished iteration: ' num2str(i) ' of ' num2str(lendb)]);
    disp('****************************************************************');
    
end

    

