% Major correction to only evaluate loop at the x2 event times as opposed
% to all the event times which was the case before

% Rectified to evaluate C matrix within ODE as this is necessary since it
% evolves across time. Redundant code removed including symbolic ODE and
% debbugging variables and conditions lifted

% Modified to allow for proper switching of the integration limits from
% T(i-1) to T(i) to 0 to dt(i)

% Function implements filtering equations from Snyder in combination with
% use of ODE integration
function outFil = filterHybrid3(lamcap0, q0, B, T, S, x1cap0, coeff, lenS, birStr2, params, x2)

% Boolean to switch between ode modes and integration limits
ode_sw = 0;
int_sw = 0;

% Find interevent times of x2 births
[Tevent perc] = getEventTimes(T, x2, 'birth');
len = length(Tevent);
dt = [0; diff(Tevent)];

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
Qset = cell(1, 1);
Tset = cell(1, 1);
lamcapFull = cell(1, 1);
options = odeset('NonNegative', 1:lenS);
% options = odeset(options, 'Refine', 20); % <----- this option works well
% options = odeset(options, 'RelTol', 0.001);

% Loop across dt values integrating the Snyder equations with the ode solver
% and then applying jump corrections due to x2 events
for i = 2:len

    % Solve ODEs continuously with setting of options, integration limits
    % options = odeset(options, 'MaxStep', 0.01*dt(i));
    switch(ode_sw)
        case 0
            if ~int_sw
                % Normal numerical solution which may have unconstrained
                [tset qset] = ode113(@(ts, y) odeSnyder(ts, y, B, coeff, S, birStr2),...
                    [Tevent(i-1) Tevent(i)], q(i-1, :)', options);
            else
                [tset qset] = ode113(@(ts, y) odeSnyder(ts, y, B, coeff, S, birStr2),...
                    [0 dt(i)], q(i-1, :)', options);
                tset = tset + Tevent(i-1);
            end
        case 1
            if ~int_sw
                % Simple constrained numerical solution derived from symbolics
                % Assumes space <= 3 and that kbirth = kdeath
                [tset qset] = ode113(@(ts, y) odeSimpleState(ts, y, params.kbirth,...
                    params.kgain), [Tevent(i-1) Tevent(i)], q(i-1, :)', options);
            else
                [tset qset] = ode113(@(ts, y) odeSimpleState(ts, y, params.kbirth,...
                    params.kgain), [0 dt(i)], q(i-1, :)', options);
                tset = tset + Tevent(i-1);
            end
    end
                % Obtain lamcapset from qset output - assume linear cross case
            x1capset = qset*S;
            lamcapset = coeff(1)*sum(x1capset, 2);

    % Assign output value at event times
    q(i, :) = qset(end, :);
    t(i) = tset(end);
    qODE(i, :) = q(i, :);
    Qset{i} = qset;
    Tset{i} = tset;
    lamcapFull{i} = lamcapset;

    % Check raw ODE q output
    if any(q(i, :) < -10^-8)
        assignin('base', 'qODE', q(i, :));
        error(['qODE distribution has negative entries at i =' num2str(i)]);
    end
    if max(abs(sum(q(i, :)) - 1)) > 10^-4
        assignin('base', 'qODE', q(i, :));
        error(['qODE distribution does not sum to 1 at i = ' num2str(i)]);
    end

    % Determine lamcap and x1cap estimates before jump and check lamcapset
    [lamcap(i) x1cap(i)] = calclamEst(birStr2, coeff, q(i, :), S);
    if lamcap(i) ~= lamcapset(end)
        assignin('base', 'lamcapErr', lamcap(i));
        assignin('base', 'lamcapsetErr', lamcapset);
        error(['lamcap estimates do not match at i = ' num2str(i)]);
    end

    % Calculate perturbation on q due to jump
    q(i, :) = calcGMx(birStr2, coeff, S, q(i, :), lamcap(i), i);

    % Determine lamcap and x1cap estimates after jump and check results
    [lamcap(i) x1cap(i)] = calclamEst(birStr2, coeff, q(i, :), S);

    % Display progress and debug set of results from ODE solver
    disp(['Finished iteration: ' num2str(i-1) ' of ' num2str(len)]);

end

% Assign output data to a single structure
outFil.q = q;
outFil.lamcap = lamcap;
outFil.t = t;
outFil.x1cap = x1cap;
outFil.qODE = qODE;
outFil.Tevent = Tevent;
outFil.perc = perc;
outFil.settings = [ode_sw int_sw];
outFil.Qset = Qset;
outFil.Tset = Tset;
outFil.lamcapFull = lamcapFull;

% Assign ODE solver data directly to workspace
% assignin('base', 'Qset', Qset);
% assignin('base', 'Tset', Tset);
% assignin('base', 'lamcapFull', lamcapFull);