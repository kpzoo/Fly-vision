% Function is based on gillespieManyReactionsQSimple but modified to
% include more light models of different forms

% Function simulates several light models provided the Q matrix is defined
% and the model type - current settings assume M state Markov chains that
% code M/2 positions and codes 2 directions
function [X Xdot T] = gillespieQGeneral(inpGill)

% Deconstruct inputs into individual variables
len = inpGill.len;
x0 = inpGill.x0;
N = inpGill.N;
Nstart = inpGill.Nstart;
nPos = inpGill.nPos;
modelType = inpGill.type;
    
% Inputs for the Q matrix rate method
Q = inpGill.Q;
transit = inpGill.transit;
nReacs = inpGill.nReacs;
intenParams.i0 = inpGill.alpha;

% Set initial parameters for time and populations
t = zeros(N, 1);
z = zeros(N, len);
z(1, :) = x0;
alp = zeros(N, nReacs);

% Loop across specified number of iterations and simulate
for i = 2:N

    % Update molecular numbers and other data sequentially
    x = z(i - 1, :);
    told = t(i-1);

    % Get reaction rates of the state process from Q matrix
    rdotx1 = getRatesQMx2(Q, x(1), nReacs-nPos);
    
    % Based on model type obtain the rate at various positions for the
    % observed counting processes
    R = getPositionalIntensity(x(1), modelType, intenParams, nPos);
    
    % Obtain characteristics of the next reaction based on rate set
    rdot = [rdotx1 R];
    rdotsum = sum(rdot);
    tnex = told - log(rand)/rdotsum;
    rdot_ratio = rdot/rdotsum;
    reac = 1 + sum(rand > cumsum(rdot_ratio));
    
    % Ensure the next reaction is sensible
    try
        xnex = x + transit(:, reac)';
    catch excep1
        assignin('base', 'reac', reac);
        assignin('base', 'transit', transit);
        assignin('base', 'excep1', excep1);
        error('Problems assinging xnex');
    end
    
    % Assign variables and ensure dimensions are consistent
    try
        alp(i-1, :) = rdot;
    catch excep2
        assignin('base', 'rdot', rdot);
        assignin('base', 'alp', alp(i-1, :));
        assignin('base', 'excep2', excep2);
        error('Dimension mismatch');
    end
    t(i) = tnex;
    z(i, :) = xnex;

    % Catch a possible error
    if rdotsum == 0
        assignin('base', 'z', z(1:i, :));
        assignin('base', 'alp', alp(1:i-1, :));
        error(['All rates are zero at i = ' num2str(i)]);
    end

end

% Save data from simulation for post processing and account for the control
x = z;
xdot = alp;

% Account for equilibrium and remove last value (in previous code would
% calculate the next step at i = N)
X = x(Nstart:N-1, :);
Xdot = xdot(Nstart:N-1, :);
T = t(Nstart:N-1, :);
