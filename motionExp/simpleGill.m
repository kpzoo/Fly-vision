% Caution: method here assumes a specific type of rate and only 1 species
% <------------------------------------------------------------------------
% Simple application of the Gillespie algorithm - based on the more complex
% gillespieManyReacsQ.m file found in folder general markov code
function [X Xdot T] = simpleGill(inpGill)

% Deconstruct inputs into individual variables
len = inpGill.len;
x0 = inpGill.x0;
N = inpGill.N;
Nstart = inpGill.Nstart;
r_const = inpGill.r_const;
transit = inpGill.transit;
stateMax = inpGill.stateMax; %<------------ assume only one maximum
stateMin = inpGill.stateMin; %<------------ assume only one minimum

% Ensure only 1 species considered at the moment
if length(r_const) ~= 2 || length(stateMax) ~= 1
    assignin('base', 'inpGill', inpGill);
    error('Code assumes only 1 species with constitutive rates');
end

% Set initial parameters for time and populations
t = zeros(N, 1);
z = zeros(N, len);
z(1, :) = x0;

% Set initial rates and reaction incrementation - based on r_const
nReacs = length(r_const);
alpha = zeros(N, nReacs);
% rdot = zeros(1, nReacs);


% Loop across specified number of iterations and simulate
for i = 2:N
    
    % Update molecular numbers and other data sequentially
    x = z(i - 1, :);
    told = t(i-1);    
    
    % Obtain reaction rates - assumes constitutive rates and forces a fixed
    % state space <--------------------------------------------------------
    rdot = r_const;
    if x == stateMax
        rdot = [0 r_const(2)];
    elseif x == stateMin
        rdot = [r_const(1) 0];
    end
    
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

% Ensure fixed state space not violated
if any(x > stateMax) || any(x < stateMin)
    assignin('base', 'xErr', x);
    error('The state space has been violated');
end

% Account for equilibrium and remove last value (in previous code would
% calculate the next step at i = N)
X = x(Nstart:N-1, :);
Xdot = xdot(Nstart:N-1, :);
T = t(Nstart:N-1, :);