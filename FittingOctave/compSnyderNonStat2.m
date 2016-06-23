% Modified to include calculation of the linear ODE set between
% observations as an alternate Snyder calculation

% Note - this function is very different in the matrices constructed when
% compared to those in compSnyder and assumes X = [x y z] where y = ybir

% Modified version of compSnyder which handles marked point process case in
% which the information process is non-stationary

% Function to calculate the compensated Snyder form - assumes 2 state x and
% a simple exponential delay
function [x1Stats lam Q] = compSnyderNonStat2(T, X, Sx, Sy, rateParams, x1, T1)

% Ensure only a 2 state x input is applied and set boolean to change
% between integration limits and solution type
if length(diag(Sx)) ~= 2
    error('Code cannot handle information processes (x) with more (or less) than 2 states');
end
int_sw = 0;
linearODE = 1;

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
d2 = repmat([0 a], 1, rep);
d2 = d2(1:dimQ - 2);

% Obtain the Q matrix (lam matrix now in ode formulation)
Q = Q + diag(d1, 1) + diag(d2, 2) + diag(d1neg, -1);
d0 = -sum(Q, 2);
Q = Q + diag(d0);

% Start with an arbitrary q0 as uniform and calculate lamcap0, x1cap0
q0 = (1/dimQ)*ones(1, dimQ);
% q0 = [0 1 zeros(1, dimQ-2)];
lamDiag0 = zeros(1, dimQ);
iV = 3;
for iL = 2:dimQ/2
    lamDiag0(iV:iV+1) = b*((iL-1) - X(1, end));
    iV = iV + 2;
end
lam0 = diag(lamDiag0);
lamcap0 = q0*lam0*ones(dimQ, 1);
x1cap0 = sum(q0(2:2:end)); % <------------ assumes 2 state x MC

% Obtain inter-event times of point process observations
[Tevent perc] = getEventTimes(T, X(:, end), 'birth');
[Teventy percy] = getEventTimes(T, X(:, end-1), 'birth');
len = length(Tevent);
dt = [0; diff(Tevent)];

% Check closeness of event times of y and z
assignin('base', 'Teventz', Tevent);
assignin('base', 'Teventy', Teventy);
normEv = norm(Teventy - Tevent);
disp(['The event time norm is ' num2str(normEv)]);

% Declare variables and initialise
q = zeros(len, dimQ);
q(1, :) = q0;
t = zeros(len, 1);
lamcap = zeros(len, 1);
lamcap(1) = lamcap0;
x1cap = zeros(len, 1);
x1cap2 = zeros(len, 1);
x1cap(1) = x1cap0;
x1cap2(1) = x1cap0;
qODE = zeros(len, dimQ);
qODE(1, :) = q0;
qpert = zeros(len, dimQ);
qpert2 = zeros(len, dimQ);
betaSum0 = zeros(len, 1);
betaSum1 = zeros(len, 1);
betaSumRatio = zeros(len, 1);

% Cell to save output of ODE solver and set options
Qset = cell(1, 1);
DQset = cell(1, 1);
Tset = cell(1, 1);
x1Set = cell(1, 1);
% options = odeset('NonNegative', 1:dimQ); % <--- include options as input to odeCompSnyder
% options = odeset(options, 'Refine', 20); % <----- this option works well
% options = odeset('Refine', 20);
% options = odeset(options, 'RelTol', 0.001);
% options = odeset('RelTol', 0.00001);
% options = odeset(options, 'AbsTol', 0.0000001);

% Loop across dt values integrating the Snyder equations with the ode solver
% and then applying jump corrections due to observed events
for i = 2:len
    
    % Obtain lam matrix
    lamDiag = zeros(1, dimQ);
    iV = 3;
    for iL = 2:dimQ/2
        lamDiag(iV:iV+1) = b*((iL-1) - X(i, end));
        iV = iV + 2;
    end 
    lam = diag(lamDiag);
    
    if linearODE
        % Linear ODE solution between events with normalisation
        %[tset qset] = ode15s(@(ts, y) odeCompSnyderNonStatLinear(ts, y, Q, lam),...
        %        [Tevent(i-1) Tevent(i)], q(i-1, :)');
        %for j = 1:size(qset, 1)    
        %    qset(j, :) = qset(j, :)/(sum(qset(j, :)));
        %end
		
		% OCTAVE: modification for lsode function
		dtset = (Tevent(i) - Tevent(i-1))/100;
		tset = Tevent(i-1):dtset:Tevent(i);
		tset = tset';
		odeFN = @(qset, ts) odeCompSnyderNonStatLinear(qset, ts, Q, lam);
		lsode_options('integration method', 'stiff');
		qset = lsode(odeFN, q(i-1, :)', tset);
		
	else
        % Non-linear ODE solution between events
        if ~int_sw
            % Delayed numerical solution with integration limits in absolute time
        %    [tset qset] = ode23tb(@(ts, y) odeCompSnyderNonStat(ts, y, Q, X(i, end), dimQ, b, lam),...
        %        [Tevent(i-1) Tevent(i)], q(i-1, :)');
			
			% OCTAVE: modification for lsode function
			dtset = (Tevent(i) - Tevent(i-1))/100;
			tset = Tevent(i-1):dtset:Tevent(i);
			tset = tset';
			odeFN = @(qset, ts) odeCompSnyderNonStat(qset, ts, Q, X(i, end), dimQ, b, lam);
			lsode_options('integration method', 'stiff');
			qset = lsode(odeFN, q(i-1, :)', tset);
        else
            % Delayed numerical solution with integration limits in relative time
        %    [tset qset] = ode15s(@(ts, y) odeCompSnyderNonStat(ts, y, Q, X(i, end), dimQ, b, lam),...
        %        [0 dt(i)], q(i-1, :)');
        %    tset = tset + Tevent(i-1);
        end
    end
    
    % Obtain the probability rates - actual conditional probabilities for
    % non-linear case and unnormalised rates for the linear case
    dqset = zeros(size(qset));
    for j = 1:size(qset, 1)
        if linearODE
            % dqset(j, :) = odeCompSnyderNonStatLinear(0, qset(j, :)', Q, lam);
        else
            dqset(j, :) = odeCompSnyderNonStat(qset(j, :)', 0, Q, X(i, end), dimQ, b, lam);
        end
    end
    
    % Assign output value at event times
    q(i, :) = qset(end, :);
    t(i) = tset(end);
    qODE(i, :) = q(i, :);
    Qset{i} = qset;
    DQset{i} = dqset;
    Tset{i} = tset;
    x1Set{i} = sum(qset(:, 2:2:end), 2); % <----------- assumes that x is 2 state
    
    % Check raw ODE q output and derivatives
    if any(q(i, :) < -10^-8)
        assignin('base', 'qODErr', q(i, :));
        disp(['qODE distribution has negative entries at i =' num2str(i)]);
    end
    if max(abs(sum(q(i, :)) - 1)) > 10^-4
        assignin('base', 'qODErr', q(i, :));
        disp(['qODE distribution sums to ' num2str(sum(q(i, :))) ' at i = ' num2str(i-1)]);
        assignin('base', 'qset', qset);
    end

    if max(abs(sum(dqset(end, :)))) > 10^-4
        assignin('base', 'dqODErr', dqset(end, :));
        disp(['dqODE distribution sums to ' num2str(sum(q(i, :))) ' at i = ' num2str(i-1)]);
        assignin('base', 'dqset', dqset);
    end
    
    % Calculate perturbation on q due to jump and calculate lamcap
    q(i, :) = q(i, :)*lam/(q(i, :)*lam*ones(dimQ, 1));
    lamcap(i) = q(i, :)*lam*ones(dimQ, 1);
    qnew = q(i, :);
    x1cap(i) = sum(qnew(2:2:end));
    qpert(i, :) = qnew;
    
    
    % Calculate another qpert
    lamDiag = zeros(1, dimQ);
    iV = 3;
    for iL = 2:dimQ/2
        lamDiag(iV:iV+1) = b*((iL-1) - X(i+1, end));
        iV = iV + 2;
    end 
    lam = diag(lamDiag);
    qnew2 = q(i, :)*lam/(q(i, :)*lam*ones(dimQ, 1));
    x1cap2(i) = sum(qnew2(2:2:end));
    qpert2(i, :) = qnew2; 
    q(i, :) = qnew2; % <------- use an updated lam matrix to get the IC for the next step
    
    % Obtain the beta sum values
    Dlam = diag(lam);
    betaSum0(i) = qset(end, 1:2:end)*Dlam(1:2:end);
    betaSum1(i) = qset(end, 2:2:end)*Dlam(2:2:end);
    betaSumRatio(i) = betaSum1(i)/(betaSum1(i) + betaSum0(i));
    
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

% Assign variables to workspace for debugging
assignin('base', 'x1n', x1n);
assignin('base', 'x1capn', x1capn);
assignin('base', 'Tnn', Tn);
assignin('base', 'qODE', qODE);
assignin('base', 'qpert', qpert);
assignin('base', 'qpert2', qpert2);
assignin('base', 'X', X);
assignin('base', 'x1cap', x1cap);
assignin('base', 'x1cap2', x1cap2);
assignin('base', 'betaSum0', betaSum0);
assignin('base', 'betaSum1', betaSum1);
assignin('base', 'betaSumRatio', betaSumRatio);
assignin('base', 'DQset', DQset);
assignin('base', 'x1Set', x1Set);

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