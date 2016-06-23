% Modified to include extra molecules to code photon streams at various
% pixel locations - in a sense each is still 1D at the moment <------------
% assumes only 3 pixels and a constant intensity alpha = inten

% Simplified version of gillespieManyReacsQ that only works with the Q
% matrix and does not use the r_const setting or useQ method switch. Also
% removed rprev and xhist, included a clear alpha and added the nReacs as a
% function input to remove dependence on r_const

% This version allows for input of Q matrix and then calculates the x1 rdot
% entries via lookup on the infinitesimal generator. Further with this Q
% version the x1 is assumed to have only births with lam = alpha*x1

% Modified version of gillespieMarkov that is generalised to allow for
% extra non-nearest neighbour reactions - meant to simulate the effect of
% big and little jumps - now have to directly input the transit matrix
% which must correspond to r_const form [bir dea bir dea bir dea...] - also
% will now calculate reaction rates in this order directly i.e. deaType
% provides all the even indices and birType the odd ones

% Simplified Gillespie algorithm for use with filter algorithms - no
% modulation types on x1 is included as it is only required that x1
% modulate the intensity of x2
function [X Xdot T] = gillespieManyReacsQSimple(inpGill)

% Deconstruct inputs into individual variables
len = inpGill.len;
x0 = inpGill.x0;
N = inpGill.N;
Nstart = inpGill.Nstart;
nPos = inpGill.nPos;
if len ~= 4
    error('Code currently specific to 4 species with 3 pixels');
end
    
% Inputs for the Q matrix rate method
Q = inpGill.Q;
transit = inpGill.transit;
alpha = inpGill.alpha;
nReacs = inpGill.nReacs;

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
    % rdotx2 = alpha*x(1);
    
    % Obtain counting process rates where each molecule represents photons
    % at a certain location <------------------ specific to 3 positions and
    % assumes a constant intensity alpha = inten
    rdotx2 = 0;
    rdotx3 = 0;
    rdotx4 = 0;
    if x(1) == 0 || x(1) == 3
        rdotx2 = alpha;
    end
    if x(1) == 1 || x(1) == 4
        rdotx3 = alpha;
    end
    if x(1) == 2 || x(1) == 5
        rdotx4 = alpha;
    end
    
    % Obtain characteristics of the next reaction based on rate set
    rdot = [rdotx1 rdotx2 rdotx3 rdotx4];
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
