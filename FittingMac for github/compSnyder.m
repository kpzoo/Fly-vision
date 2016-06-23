% Function to calculate the compensated Snyder form - assumes 2 state x and
% a simple exponential delay
function [x1Stats lam Q] = compSnyder(T, X, Sx, Sy, rateParams, x1, T1)

% Ensure only a 2 state x input is applied and set boolean to change
% between integration limits
if length(diag(Sx)) ~= 2
    error('Code cannot handle information processes (x) with more than 2 states');
end
int_sw = 0;

% Extract rate parameters with a = alpha, b = eta for a single delay stage
a = rateParams.alpha;
b = rateParams.eta;
k = rateParams.k;

% Initialise appropriate Q matrix with suitable dimensions
dimQ = length(diag(Sx))*length(diag(Sy));
Q = zeros(dimQ, dimQ);
rep = dimQ/2; % must be even by assumption on Sx <------------------------

% Construct the Q matrix diagonals
d1 = repmat([k 0], 1, rep);
d1 = d1(1:dimQ - 1);
d1neg = d1;
d2 = a*ones(1, dimQ - 2);
d2(1) = 0; % no rate for x = 0
d2neg = b*ones(1, dimQ - 2);
d2negscale = sort(repmat(1:rep, 1, 2));
d2neg = d2neg.*d2negscale(1:dimQ - 2);

% Obtain the Q and lambda matrices
Q = Q + diag(d1, 1) + diag(d2, 2) + diag(d1neg, -1) + diag(d2neg, -2);
d0 = -sum(Q, 2);
Q = Q + diag(d0);
lamScale = sort(repmat(0:rep, 1, 2));
lam = diag(lamScale(1:dimQ));          
lam = b*lam;

% Start with an arbitrary q0 as uniform and calculate lamcap0, x1cap0
q0 = (1/dimQ)*ones(1, dimQ);
% Improved ICs for high eta
if b >= 100 && dimQ == 4
    q0 = [0.5 0.5 0 0];
end
lamcap0 = q0*lam*ones(dimQ, 1);
x1cap0 = sum(q0(2:2:end)); % <------------ assumes 2 state x MC

% Obtain inter-event times of point process observations
[Tevent perc] = getEventTimes(T, X(:, end), 'birth');
len = length(Tevent);
dt = [0; diff(Tevent)];

% Declare variables and initialise
q = zeros(len, dimQ);
q(1, :) = q0;
t = zeros(len, 1);
lamcap = zeros(len, 1);
lamcap(1) = lamcap0;
x1cap = zeros(len, 1);
x1cap(1) = x1cap0;
qODE = zeros(len, dimQ);
qODE(1, :) = q0;
qpert = zeros(len, dimQ);
qLimTest = zeros(len, 1);
ycap = zeros(len, 1);

% Cell to save output of ODE solver and set options
Qset = cell(1, 1);
Tset = cell(1, 1);
x1Set = cell(1, 1);
lamcapFull = cell(1, 1);
ySet = cell(1, 1);
% options = odeset('NonNegative', 1:dimQ); % <--- include options as input to odeCompSnyder
% options = odeset(options, 'Refine', 20); % <----- this option works well
% options = odeset('Refine', 20);
% options = odeset(options, 'RelTol', 0.001);
% options = odeset('RelTol', 0.00001);
% options = odeset(options, 'AbsTol', 0.0000001);

% Loop across dt values integrating the Snyder equations with the ode solver
% and then applying jump corrections due to observed events
for i = 2:len 
    if ~int_sw
        % Delayed numerical solution with integration limits in absolute time
        [tset qset] = ode23s(@(ts, y) odeCompSnyder(ts, y, Q, lam, dimQ),...
            [Tevent(i-1) Tevent(i)], q(i-1, :)');
    else
        % Delayed numerical solution with integration limits in relative time
        [tset qset] = ode23s(@(ts, y) odeCompSnyder(ts, y, Q, lam, dimQ),...
            [0 dt(i)], q(i-1, :)');
        tset = tset + Tevent(i-1);
    end
    
    % Assign output value at event times
    q(i, :) = qset(end, :);
    t(i) = tset(end);
    qODE(i, :) = q(i, :);
    Qset{i} = qset;
    Tset{i} = tset;
    x1Set{i} = sum(qset(:, 2:2:end), 2); % <----------- assumes that x is 2 state
    ySet{i} = sum(qset(:, 3:4), 2); % <----------- assumes that y is 2 state
    lamcapset = qset*lam*ones(dimQ, 1);
    lamcapFull{i} = lamcapset;
    
    % Check raw ODE q output
    if any(q(i, :) < -10^-8)
        assignin('base', 'qODE', q(i, :));
        error(['qODE distribution has negative entries at i =' num2str(i)]);
    end
    if max(abs(sum(q(i, :)) - 1)) > 10^-4
        assignin('base', 'qODErr', q(i, :));
        disp(['qODE distribution sums to ' num2str(sum(q(i, :))) ' at i = ' num2str(i)]);
    end
    
    % Test based on a comparison of these calculations with the normal
    % Snyder - only makes sense at high eta
    if dimQ == 4
        qLimTest(i) = b*q(i, 3)*(1 + q(i, 1) + q(i, 2)) - b*q(i, 4)*(q(i, 2) - q(i, 1));
    end
    
    % Calculate perturbation on q due to jump and calculate lamcap
    q(i, :) = q(i, :)*lam/(q(i, :)*lam*ones(dimQ, 1));
    lamcap(i) = q(i, :)*lam*ones(dimQ, 1);
    qnew = q(i, :);
    x1cap(i) = sum(qnew(2:2:end));
    qpert(i, :) = qnew;
    
    % Display progress and debug set of results from ODE solver
    disp(['Finished iteration: ' num2str(i-1) ' of ' num2str(len)]);
end

% Obtain qODE stats and display
disp('Finished estimate of modulating intensity');
sumqODE = sum(qODE, 2);
mqODE = mean(sumqODE);
vqODE = var(sumqODE);
disp(['Mean and var of qODE = ' [num2str(mqODE) ' ' num2str(vqODE)]]);

% Obtain 3 estimates of x1 statistics (with the appended x1 stream applied)
rem_trans = 0;
% [x1Stats x1n x1capn Tn] = getAllStats(x1, x1Set, Tset, T1, Qset, rem_trans, -1);
[x1Stats x1n x1capn Tn] = getAllStats(X(:, 1), x1Set, Tset, T, Qset, rem_trans, -1);
% [yStats yn ycapn Tn] = getAllStats(X(:, end-1), ySet, Tset, T, Qset, rem_trans, -1);
assignin('base', 'x1n', x1n);
assignin('base', 'x1capn', x1capn);
assignin('base', 'Tnn', Tn);
assignin('base', 'qODE', qODE);
assignin('base', 'qpert', qpert);
assignin('base', 'qLimTest', qLimTest);
% assignin('base', 'yn', yn);
% assignin('base', 'ycapn', ycapn);
% assignin('base', 'yStats', yStats);
assignin('base', 'X', X);
clear lamn

% Quick check of minimum solution if eta high (result should match normal
% Snyder for high enough eta or small enough delay)
if b >= 100
    beta = a/k;
    A = 1/beta + 0.5;
    B = sqrt(1/(beta^2) + 0.25);
    steady.minTheo = A - B;
    steady.minEmp = min(x1capn);
    assignin('base', 'steady', steady);
end