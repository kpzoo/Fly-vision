% This version has both y and z as non-stationary processes while the
% previous version had y as stationary

% Gillespie code to simulate the effect of stages of exponential delay with
% x as the hidden MC, y the non-delayed photons and z the stages of delay
% Assumes a 2 state symmetrical MC at the moment
clear all
clc
close all


% Initialise main variables
N = 180000;
Nstart = 20000;
nMolec = 3;
T = zeros(N, 1);
Z = zeros(N, nMolec);
Zdot = zeros(N, 4); %<-------- modify as needed

% Initialise rate parameters
k = 5e-4;
alpha = 0.01;
eta = 100;

% Main Gillespie algorithm
for i = 2:N
    % Markov state and undelayed photons
    x = Z(i - 1, 1);
    y = Z(i - 1, 2);
    rdot3 = [k*(1 - x) k*x alpha*x];
    % Delayed molecules based on stages
    switch(nMolec - 2)
        case 0
            % The undelayed case
            rdot = [rdot3 0];
            R1 = [1 -1 0 0];
            R2 = [0 0 1 -1];
            R = [R1; R2];
        case 1
            z = Z(i - 1, 3);
            rdot = [rdot3 eta*(y - z)];
            R1 = [1 -1 0 0];
            R2 = [0 0 1 0];
            R3 = [0 0 0 1];
            R = [R1; R2; R3];                        
        otherwise
            error('Incorrect number of exponential stages specified');
    end
    % Calculate the next reaction time from rate sum and update molecular
    % numbers based on reaction chosen
    rdotsum = sum(rdot);
    T(i) = T(i - 1) - log(rand)/rdotsum;
    rdot = rdot/rdotsum;
    reac = 1 + sum(rand > cumsum(rdot));
    Z(i, :) = Z(i - 1, :) + R(:, reac)';
    Zdot(i, :) = rdot;
end

% Compare closeness of y and z births
[Teventz percz] = getEventTimes(T, Z(:, end), 'birth');
[Teventy percy] = getEventTimes(T, Z(:, end-1), 'birth');
normEv = norm(Teventy - Teventz);
disp(['The event time norm is ' num2str(normEv)]);

% Plot the main molecular trajectories
for i = 1:nMolec
    figure;
    stairs(T, Z(:, i));
    xlabel('time');
    switch(i)
        case 1
            ylabel('x');
        case 2
            ylabel('y');
        otherwise
            ylabel(['z_{' num2str(i - 2) '}']);
    end
end

% Extra plot comparing z with y births
if nMolec == 3
    figure;
    hold on
    stairs(T, Z(:, end-1));
    stairs(T, Z(:, end), 'r');
    hold off
    xlabel('time');
    ylabel('molecular numbers');
    legend('y births', 'z photons', 'location', 'best');
end

% Save simulation data with suitable parameter set
X = Z(Nstart:N, :);
T = T(Nstart:N, :);
Xdot = Zdot(Nstart:N, :);
Slim = [min(X(:, 1)) max(X(:, 1))];
outGil.X = X;
outGil.T = T;
outGil.Xdot = Xdot;
outGil.r_const = [k k alpha 0]; % based on x-y only
outGil.Slim = Slim;
clear Z T Zdot X
params.kdeath = k;
params.kbirth = k;
params.kgain = alpha;
params.eta = eta;
params.Slim = Slim;
params.Nstart = Nstart;
params.N = N;
params.len = 2;
params.x0 = [0 0];
markov = 1;
coeff = [alpha 0];
inpB.markov = markov;
inpB.coeff = coeff;
inpB.simname = 'sim1StageDel';
inpB.kbirth = k;
inpB.kdeath = k;
inpB.kgain = alpha;
inpB.Slim = Slim;
inpB.Nset = [Nstart N];
save(inpB.simname);

    