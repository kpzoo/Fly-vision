% Input test code to general Gillespie simulation involving reactions with
% different bulk sizes - after debugging this code will be incorporated
% into a function similar to dsppFilterMarkov
clear all
clc
close all

% This reaction set includes 1 and 2 state changes on a finite x1 chain and
% has the inputs set as below
plot_on = 1;
Nstart = 30000;
N = 35000;
len = 2;
x0 = [0 0];
avgR = [100 100 100 100 100];
maxR = [1000 1000 1000 1000 1000];
minR = [0 0 0 0 0];
SlimSet.min = [0 0];
SlimSet.max = [10 inf];

% Rate structure based inputs - bulk must be positive, all variables must
% be same length and transit represents the incremental changes that bulk
% reactions yield - odd reactions are births and even deaths
reacType = [1 1 1 1 2];
r_const = [10 10 5 5 10];
bulk = [1 1 2 2 1];
transit = [1 -1 2 -2 0; 0 0 0 0 1];
molecType = [1 1 1 1 2];
crossType = [1 1 1 1 1];

% Inputs to Gillespie Markov algorithm
inpGill.len = len;
inpGill.N = N;
inpGill.Nstart = Nstart;
inpGill.reacType = reacType;
inpGill.molecType = molecType;
inpGill.crossType = crossType;
inpGill.avgR = avgR;
inpGill.maxR = maxR;
inpGill.minR = minR;
inpGill.x0 = x0;
inpGill.r_const = r_const;
inpGill.bulk = bulk;
inpGill.transit = transit;
inpGill.SlimSet = SlimSet;

% Ensure molecular types sensible
if ~all(ismember([molecType crossType], 1:len))
    assignin('base', 'molecType', molecType);
    assignin('base', 'len', len);
    error('Incorrect molecular identifiers specified');
end

% Run the actually Gillespie algorithm for specified Markov chain
disp('Gillespie simulation of Markov chain started');
[X Xdot T] = gillespieManyReacs(inpGill);
disp('Gillespie simulation of Markov chain complete');
lenEquil = 1:length(T);

% Obtain mse and stats across time - note assumes that only want mse
% between first 2 species <-----------------------------------------------
for i = 1:length(r_const)
    [rate.max(i) rate.mean(i) rate.var(i) rate.min(i)] = calcStats2(Xdot(lenEquil, i), T);
end
for i = 1:len
    [molec.max(i) molec.mean(i) molec.var(i) molec.min(i)] = calcStats2(X(lenEquil, i), T);
end

% Specify important rates, populations and their moments <-----------
% specific to type of Markov chain initiated
x1 = X(lenEquil, 1);
x2 = X(lenEquil, 2);
f = Xdot(lenEquil, end);
fmean = rate.mean(end);
fvar = rate.var(end);

% Plots of responses
if plot_on
    % Trajectory of counts
    figure;
    stairs(T, x1);
    title('Trajectory of x1');
    xlabel('time');
    ylabel('no. molecules');
    xlim([min(T) max(T)]);
    figure;
    stairs(T, x2);
    title('Trajectory of x2');
    xlabel('time');
    ylabel('no. molecules');
    xlim([min(T) max(T)]);
        
    % Plot the rate f across time
    figure;
    stairs(T, f);
    hold on
    stairs(T, fmean*ones(size(T)), 'r');
    hold off
    xlabel('time');
    ylabel('modulated rate');
    title(['Modulated intensity, [mean var] = ' [num2str(fmean) ' ' num2str(fvar)]]);
    xlim([min(T) max(T)]);
end
disp('Post processing of SSA output complete and data stored if specified');