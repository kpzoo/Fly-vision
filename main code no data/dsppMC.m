% Functional form of general markov chain Gillespie simulation
function output = dsppMC(params)

% Obtain input parameters - no. molecules, iterations, rate extrema and
% molecular ICs, start of equilibrium in terms of N
len = params.len;
N = params.N;
x0 = params.x0;
plot_on = params.plot_on;
avgR = params.avgR;
maxR = params.maxR;
minR = params.minR;
Nstart = params.Nstart;
SlimSet = params.SlimSet;

% Rate structure based inputs - bulk must be positive, all variables must
% be same length and transit represents the incremental changes that bulk
% reactions yield - odd reactions are births and even deaths
reacType = params.reacType;
molecType = params.molecType;
crossType = params.crossType;
bulk = params.bulk;
transit = params.transit;
r_const = params.r_const;

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

% Assign outputs to a single data structure
output.X = X;
output.Xdot = Xdot;
output.T = T;
output.molec = molec;
output.rate = rate;

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