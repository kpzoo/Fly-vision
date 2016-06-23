% Modified to include account for a deterministic delay in the estimation
% procedure which in this case simply involves delaying the photons only

% Modified to include a plot of the estimated parameter and to allow a
% scaling on the intensity to see if the MSE rises with higher gain

% Code assumes a certain form of SHM as (r+1)sin(wt) + max(rspace+1),
% changes to this form will invalidate the simulations<--------------------

% Implementation of a deterministic light model for SHM motion with random
% initial condition that is Snyder filtered
clear all
clc
close all

% SHM motion with random parameter r in [0 1] and Nev observations - note
% the motion used is of form rsin(wt) + max({r}) so rates are non-negative
rspace = 0:5;
lenq = length(rspace);
r = randsample(max(rspace)+1, 1);
r = r - 1;
w = 2;

% Set simulation and other control variables
Nev = 5000;
Nstart = 1;
alpha = 1;
delayON = 1;
delay = 5;

%% Run the stochastic simulation of an inhomogeneous Poisson process

% Set the parameters for the Gillespie simulation accounting for the form
% of the SHM chosen <-----------------------------------------------------
rmax = max(rspace + 1);
lam_max = alpha*(r + 1 + rmax);

% Simulate a homogeneous Poisson process with 3*Nev*rmax/r points at max rate
Nhpp = ceil(3*Nev*rmax/(r+1));
thpp = exprnd(1/lam_max, Nhpp, 1);
thpp = cumsum(thpp);

% Thin the homogeneous process to obtain the correct intensity
randNos = rand(Nhpp, 1);
rateNonHomo = alpha*((r + 1)*sin(w*thpp) + rmax);
rateRatio = rateNonHomo/lam_max;
tihpp = thpp(randNos < rateRatio);
xihpp = 1:length(tihpp);
xihpp = xihpp';

% Assign simulation vectors and account for delayed photons
x2 = xihpp(Nstart:Nev);
T = tihpp(Nstart:Nev);
if delayON
    T = T + delay;
end
avgRate = x2(end)/T(end);
nEvents = length(x2);
x1 = (r + 1)*sin(w*T) + rmax;

%% Snyder filtering for Poisson process modulated by a random variable

% Declare variables and initialise
q = zeros(nEvents, lenq);
q(1, :) = ones(1, lenq)/lenq;
t = zeros(nEvents, 1);

% In this case x1 refers to position which is continuous and initialised as
% 0 since at t = 0, sin(wt) = 0 <------------------------------------------
x1cap = zeros(nEvents, 1);
x1cap(1) = 0;
paramcap = zeros(nEvents, 1);
lamDiag = zeros(1, lenq);

% Cell to save output of ODE solver and set options
Qset = cell(1, 1);
Tset = cell(1, 1);
x1capset = cell(1, 1);
options = odeset('NonNegative', 1:lenq);
stiffODE = 0;

% Loop across dt values integrating the Snyder equations with the ode solver
% and then applying jump corrections due to x2 events
for i = 2:nEvents
    
    % Solve non-linear ODEs continuously with setting of options 
    if ~stiffODE
        [tset qset] = ode113(@(ts, y) odeSnyRVNonLin2(ts, y, rspace, w, alpha),...
            [T(i-1) T(i)], q(i-1, :)', options);
    else
        [tset qset] = ode23t(@(ts, y) odeSnyRVNonLin2(ts, y, rspace, w, alpha),...
            [T(i-1) T(i)], q(i-1, :)');
    end
    
    % Obtain the position vector from qset with account of SHM and the fact
    % that lam changes with time explicitly
    x1captemp = zeros(size(tset));
    lampert = zeros(lenq);
    for j = 1:length(tset)
        for k = 1:lenq
            % Construct the lam diagonal at every time
            lamDiag(k) = (rspace(k) + 1)*sin(w*tset(j)) + max(rspace + 1);
        end
        % Obtain the position estimate at each time
        lam = alpha*diag(lamDiag);
        x1captemp(j) = sum(qset(j, :)*lam, 2)/alpha;
        % Last lam calculation will be used in the perturbation calculation
        if j == length(tset)
            lampert = lam;
        end
    end
    x1capset{i} = x1captemp;
    
    % Assign output value at event times
    q(i, :) = qset(end, :);
    t(i) = tset(end);
    Qset{i} = qset;
    Tset{i} = tset;
    
    % Check raw ODE q output
    if any(q(i, :) < -10^-8)
        assignin('base', 'qODE', q(i, :));
        error(['qODE distribution has negative entries at i =' num2str(i)]);
    end
    if max(abs(sum(q(i, :)) - 1)) > 10^-4
        assignin('base', 'qODErr', q(i, :));
        disp(['qODE distribution does not sum to 1 at i = ' num2str(i)]);
    end
    
    % Calculate perturbation on q due to jump
    q(i, :) = q(i, :)*lampert/sum(q(i, :)*lampert);
    x1cap(i) = sum(q(i, :)*lampert, 2)/alpha;
    paramcap(i) = sum(q(i, :).*rspace, 2);
    
    % Display progress and debug set of results from ODE solver
    disp(['Finished iteration: ' num2str(i-1) ' of ' num2str(nEvents-1)]);
    
end

% Obtain the position statistics with calculations at Tn (ODE times)
[~, x1capn, Tn, ~] = getIntensityStatsFull3(x1, x1capset, Tset, T, Qset, 0);
x1n = (r + 1)*sin(w*Tn) + rmax;
x1Stats = getStatsAltMeths3(Tn, x1n, x1capn);
msePos = x1Stats.interpMeth(3);
disp(['MSE for position is ' num2str(msePos)]);

% Plot the actual and estimated trajectories
figure;
plot(Tn, x1n, Tn, x1capn);
xlabel('time');
ylabel('position');
legend('actual', 'estimated', 'location', 'best');

% Plot the actual and estimated parameter at event updates
figure;
plot(T, r*ones(size(T)), T, paramcap);
xlabel('time');
ylabel('parameter');
legend('actual', 'estimated', 'location', 'best');