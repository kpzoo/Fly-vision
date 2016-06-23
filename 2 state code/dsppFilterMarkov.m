% Mofified to work with gillespieMarkov in order to simulate finite state
% space Markov chains with state dependent rates

% Modified from just dsppFilter to account for non-linear lam(x1) in
% accordance with runDSPPfilter4

% Function to run a Gillespie simulation, plot results and save input - no
% x1 modulation allowed no calculation of mse and bound as not relevant
function output = dsppFilterMarkov(params)

% Obtain input parameters - no. molecules, iterations, utilisation, birth
% and death types, start of equilibrium in terms of N
len = params.len;
N = params.N;
r = params.r;
x0 = params.x0;
kdeath = params.kdeath;
kgain = params.kgain;
kbirth = params.kbirth;
plot_on = params.plot_on;
avgR = params.avgR;
birType = params.birType;
deaType = params.deaType;
Nstart = params.Nstart;
bulk_size = params.bulk_size;
Slim = params.Slim;

% Check on consistency of inputs
if length(birType) ~= length(deaType) || length(birType) ~= len
    error('Inconsistent input length');
end

% Determine if any deaths - currently not supported
if all(deaType == 0)
    death_on = 0;
else
    death_on = 1;
end

% Assign birth and death based on specified types - these are codes
birth = cell(1, len);
death = cell(1, len);
for i = 1:len
    switch(deaType(i))
        case 1
            % Normal mass kinetics
            death{i} = 'mass';
        case 2
            % Linear Markov death kinetics for x1
            death{i} = 'deaMarkLin';
        otherwise
            % Default of zero death rate
            death{i} = 'none';
    end
    switch(birType(i))
        case 1
            % Constant birth rate
            birth{i} = 'const';
            if kbirth ~= avgR(1)
                error('Const rate parameter set incorrectly');
            end
        case 2
            % Birth following (:|M|1)
            birth{i} = 'birfoll';
        case 3
            % Cross dependence - x2 rate proportional to x1 count
            birth{i} = 'cross';
        case 4
            % Linear rate - x1 replicates in proportion to its population
            birth{i} = 'linear';
        case 5
            % Square cross dependence - x2 rate is quadratic in x1 count
            birth{i} = 'square';
        case 6
            % Negative exponential dependence - rate(x2) = ae^(-bx1)
            birth{i} = 'negexp';
        case 7
            % Linear Markov birth kinetics for x1
            birth{i} = 'birMarkLin';
        otherwise
            % Births must be assigned
            error('Unsupported birth type input');
    end
end

% Set initial rate constants (for mass action and const rates only) and 
% assume 3 reaction model
maxR = avgR./r;
minR = 0*ones(1, len);
r_const = [kbirth kdeath kgain 0];

% Gillespie simulations - note rate params must be manually set in this
% function <--------------------------------------------------------------

% Set inputs to algorithm
inpGill.len = len;
inpGill.N = N;
inpGill.Nstart = Nstart;
inpGill.deaType = deaType;
inpGill.birType = birType;
inpGill.avgR = avgR;
inpGill.maxR = maxR;
inpGill.minR = minR;
inpGill.r = r;
inpGill.x0 = x0;
inpGill.r_const = r_const;
inpGill.bulk_size = bulk_size;
inpGill.Slim = Slim;

% Obtain rates and molecular populations as output of Gillespie simulations
[X Xdot T] = gillespieMarkov(inpGill);
lenEquil = 1:length(T);

% Obtain mse and stats across time - note assumes that only want mse
% between first 2 species <-----------------------------------------------
for i = 1:2*len
    [rate.max(i) rate.mean(i) rate.var(i) rate.min(i)] = calcStats2(Xdot(lenEquil, i), T);
end
for i = 1:len
    [molec.max(i) molec.mean(i) molec.var(i) molec.min(i)] = calcStats2(X(lenEquil, i), T);
end

% Specify important rates, populations and their moments
x1 = X(lenEquil, 1);
x2 = X(lenEquil, 2);
f = Xdot(lenEquil, 3);
fmean = rate.mean(3);
fvar = rate.var(3);


% Assign outputs to a single data structure
output.X = X;
output.Xdot = Xdot;
output.T = T;
output.molec = molec;
output.rate = rate;
output.r = r;
output.birth = birth;
output.death = death;
output.death_on = death_on;
output.r_const = r_const;
output.Slim = Slim;

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