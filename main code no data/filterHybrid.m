% Rectified to evaluate C matrix within ODE as this is necessary since it
% evolves across time. Redundant code removed including symbolic ODE

% Modified to use calclamEst and to correspond with the runDSPPfilter5
% function which accommodates different functional lam(x1) - also includes
% calculation of the x1 estimate x1cap

% Modified to include improved error checking on the q matrices to account
% for sum(q) < 1 cases and negative entries

% Function implements filtering equations from Snyder in combination with
% use of ODE integration
function outFil = filterHybrid(lamcap0, q0, B, T, S, x1cap0, coeff, lenS, birStr2, params)

% Boolean for debugging, switch between ode modes and for constraint
debug_on = 0;
ode_sw = 0;

% Find interevent times (maybe normalise T to start at 0) <---------------
dt = [0; diff(T)];
len = length(dt);

% Check B and q matrix dimensions
if ~all(size(B) == lenS) || (length(q0) ~= lenS)
    error('Input matrices of incompatible sizes');
end

% Declare variables and initialise
q = zeros(len, lenS);
q(1, :) = q0;
t = zeros(len, 1);
lamcap = zeros(len, 1);
lamcap(1) = lamcap0;
x1cap = zeros(len, 1);
x1cap(1) = x1cap0;
qODE = zeros(len, lenS);
qODE(1, :) = q0;

% Cell to save output of ODE solver and set options
Yset = cell(1, 1);
Tset = cell(1, 1);
options = odeset('NonNegative', 1:lenS);
% options = odeset(options, 'Refine', 10);
options = odeset(options, 'RelTol', 0.001);

% Loop across dt values integrating the Snyder equations with the ode solver
% and then applying jump corrections due to x2 events
for i = 2:len
    
    % Solve ODEs continuously with setting of options
    % options = odeset(options, 'MaxStep', 0.01*dt(i));
    switch(ode_sw)
        case 0
            % Normal numerical solution which may have unconstrained
            [tset, qset] = ode113(@(ts, y) odeSnyder(ts, y, B, coeff, S, birStr2),...
                [T(i-1) T(i)], q(i-1, :)', options);
        case 1
            % Simple constrained numerical solution derived from symbolics
            % Assumes space <= 3 and that kbirth = kdeath
            [tset, qset] = ode113(@(ts, y) odeSimpleState(ts, y, params.kbirth,...
                params.kgain), [T(i-1) T(i)], q(i-1, :)', options);
    end
    
    % Assign output value at event times
    q(i, :) = qset(end, :);
    t(i) = tset(end);
    qODE(i, :) = q(i, :);
    Yset{i} = qset;
    Tset{i} = tset;
    
    % Check raw ODE q output
    if any(q(i, :) < -10^-5)
        assignin('base', 'qODE', q(i, :));
        error(['qODE distribution has negative entries at i =' num2str(i)]);
    end
    if max(abs(sum(q(i, :)) - 1)) > 10^-4
        assignin('base', 'qODE', q(i, :));
        error(['qODE distribution does not sum to 1 at i = ' num2str(i)]);
    end
    
    % Determine lamcap and x1cap estimates before jump and debug
    [lamcap(i) x1cap(i)] = calclamEst(birStr2, coeff, q(i, :), S);
    if debug_on
        disp(['Dim(qset) = ' num2str(size(qset)) ' at i = ' num2str(i)]);
        disp(['lamcap before jump = ' num2str(lamcap(i))]);
        disp(['x1cap before jump = ' num2str(x1cap(i))]);
    end
    
    % Calculate perturbation on q due to jump
    q(i, :) = calcGMx(birStr2, coeff, S, q(i, :), lamcap(i), i);

    % Determine lamcap and x1cap estimates after jump and debug
    [lamcap(i) x1cap(i)] = calclamEst(birStr2, coeff, q(i, :), S);
    if debug_on
        disp(['lamcap after jump = ' num2str(lamcap(i))]);
        disp(['x1cap after jump = ' num2str(x1cap(i))]);
    end
    
    % Display progress and debug set of results from ODE solver
    disp(['Finished iteration: ' num2str(i-1) ' of ' num2str(len)]);

end

% Assign output data
outFil.q = q;
outFil.lamcap = lamcap;
outFil.t = t;
outFil.x1cap = x1cap;
outFil.qODE = qODE;

% Assign ODE solver data directly to workspace
assignin('base', 'Yset', Yset);
assignin('base', 'Tset', Tset);