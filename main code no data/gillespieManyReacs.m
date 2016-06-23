% Modified version of gillespieMarkov that is generalised to allow for
% extra non-nearest neighbour reactions - meant to simulate the effect of
% big and little jumps - now have to directly input the transit matrix
% which must correspond to r_const form [bir dea bir dea bir dea...] - also
% will now calculate reaction rates in this order directly i.e. deaType
% provides all the even indices and birType the odd ones

% Simplified Gillespie algorithm for use with filter algorithms - no
% modulation types on x1 is included as it is only required that x1
% modulate the intensity of x2
function [X Xdot T] = gillespieManyReacs(inpGill)

% Deconstruct inputs into individual variables
len = inpGill.len;
x0 = inpGill.x0;
N = inpGill.N;
Nstart = inpGill.Nstart;
avgR = inpGill.avgR;
maxR = inpGill.maxR;
minR = inpGill.minR;

% Altered inputs defining the multiple reaction types
reacType = inpGill.reacType;
molecType = inpGill.molecType;
crossType = inpGill.crossType;
r_const = inpGill.r_const;
SlimSet = inpGill.SlimSet;
bulk = inpGill.bulk;
transit = inpGill.transit;

% Set initial parameters for time and populations
t = zeros(N, 1);
z = zeros(N, len);
z(1, :) = x0;

% Extra state information required by birfoll
state.setmax = zeros(N, len);
state.numMiss = zeros(N, len);

% Set initial rates and reaction incrementation - based on r_const
nReacs = length(r_const);
alpha = zeros(N, nReacs);
rdot = zeros(1, nReacs);

% Loop across specified number of iterations and simulate
for i = 2:N

    % Update molecular numbers and other data sequentially
    x = z(i - 1, :);
    told = t(i-1);

    % Account for time history even when none is available
    if i > 2
        % History is up to i-2 as the current x value is i-1 (i not
        % calculated yet but will be at end of this iteration)
        xhist = z(1:i-2, :);
        rprev = rdot;
    else
        xhist = zeros(2, len);
        rprev = zeros(1, 2*len);
    end

    % Obtain the complete set of reaction rates
    rdot = getGenReacRates(xhist, rprev, x, r_const, nReacs, reacType,...
        molecType, crossType, SlimSet, bulk);

    % Obtain characteristics of the next reaction
    rdotsum = sum(rdot);
    tnex = told - log(rand)/rdotsum;
    rdot_ratio = rdot/rdotsum;
    reac = 1 + sum(rand > cumsum(rdot_ratio));
    try
        xnex = x + transit(:, reac)';
    catch
        assignin('base', 'reac', reac);
        return;
    end
    % Assign variables
    alpha(i-1, :) = rdot;
    t(i) = tnex;
    z(i, :) = xnex;

    % Catch a possible error
    if rdotsum == 0
        assignin('base', 'z', z(1:i, :));
        assignin('base', 'alpha', alpha(1:i-1, :));
        error(['All rates are zero at i = ' num2str(i)]);
    end

end

% Save data from simulation for post processing and account for the control
x = z;
xdot = alpha;

% Account for equilibrium and remove last value (in previous code would
% calculate the next step at i = N)
X = x(Nstart:N-1, :);
Xdot = xdot(Nstart:N-1, :);
T = t(Nstart:N-1, :);
