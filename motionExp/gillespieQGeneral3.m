% Modified to allow for full wave rectification inputs and to use
% getPositionalIntensity3 function

% Adjusted to ensure that the intensity switching code only runs when
% position changes rather than on all events

% Modified to handle the reverse phi type models

% Function is based on gillespieManyReactionsQSimple but modified to
% include more light models of different forms

% Function simulates several light models provided the Q matrix is defined
% and the model type - current settings assume M state Markov chains that
% code M/2 positions and codes 2 directions
function [X Xdot T] = gillespieQGeneral3(inpGill)

% Deconstruct inputs into individual variables
len = inpGill.len;
x0 = inpGill.x0;
N = inpGill.N;
Nstart = inpGill.Nstart;
nPos = inpGill.nPos;
nStates = inpGill.nStates;
modelType = inpGill.type;

% Set model type separately if full wave rectification is performed
fullWave = inpGill.fullWave;
if fullWave
    fullWaveType = 3;
else
    fullWaveType = 0;
end
    
% Inputs for the Q matrix rate method
Q = inpGill.Q;
transit = inpGill.transit;
nReacs = inpGill.nReacs;
intenParams.i0 = inpGill.i0;
intenParams.ic = inpGill.ic;

% Set initial parameters for time and populations
t = zeros(N, 1);
z = zeros(N, len);
z(1, :) = x0;
alp = zeros(N, nReacs);
Rset = zeros(N, nPos);

% Loop across specified number of iterations and simulate
for i = 2:N

    % Update molecular numbers and other data sequentially
    x = z(i - 1, :);
    told = t(i-1);

    % Get reaction rates of the state process from Q matrix
    rdotx1 = getRatesQMx2(Q, x(1), nReacs-nPos);
    
    % <------------------------FIX THIS----------------------------------
    % Get the historical data needed to implement an intensity toggle which
    % should allow reverse phi via flashing
    if i > 2
        xold = z(i-2, 1);
        Rold = R;
    else
        % Initial settings of being in position 1 and flashing bright
        xold = 0;
        Rold = intenParams.i0*ones(1, nPos);
        Rold(1) = Rold(1) + intenParams.ic;
    end
    Rset(i, :) = Rold;
    xcurr = x(1);
    % <------------------------FIX THIS----------------------------------
    
    % Based on model type obtain the rate at various positions for the
    % observed counting processes
    if xcurr ~= xold
        if ~fullWave
            R = getPositionalIntensity3(xcurr, modelType, intenParams, nPos, Rold, xold, nStates);
        else
            R = getPositionalIntensity3(xcurr, fullWaveType, intenParams, nPos, Rold, xold, nStates);
        end
    else
        R = Rold;
    end
    
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

% Obtain the distinct positions
pos = x(:, 1);
idchange = ~ismember(pos, 0:nPos-1);
pos(idchange) = pos(idchange) - nPos;

% Variables for diagnosis
assignin('base', 'pos', pos);
assignin('base', 'Rset', Rset);

% Account for equilibrium and remove last value (in previous code would
% calculate the next step at i = N)
X = x(Nstart:N-1, :);
Xdot = xdot(Nstart:N-1, :);
T = t(Nstart:N-1, :);
