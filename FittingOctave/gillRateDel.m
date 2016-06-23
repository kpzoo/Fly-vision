% Gillespie code to simulate the effect of stages of exponential delay with
% x as the hidden MC, y the non-delayed photons and z an altered version of
% y with a reduced rate due to an introduced delay
% Assumes a 2 state symmetrical MC at the moment
clear all
clc
close all


% Initialise main variables
N = 300000;
nMolec = 3;
T = zeros(N, 1);
Z = zeros(N, nMolec);

% Initialise rate parameters
k = 10;
alpha = 100;

% Main Gillespie algorithm
for i = 2:N
    % Markov state and undelayed photons
    x = Z(i - 1, 1);
    y = Z(i - 1, 2);
    rdot3 = [k*(1 - x) k*x alpha*x];
    
    % Delayed molecules based on a reduced rate
    z = Z(i - 1, 3);
    b = exprnd(1/80);
    rdot = [rdot3 1/(1/(alpha*x) + b)];
    R1 = [1 -1 0 0];
    R2 = [0 0 1 0];
    R3 = [0 0 0 1];
    R = [R1; R2; R3];
    
    % Calculate the next reaction time from rate sum and update molecular
    % numbers based on reaction chosen
    rdotsum = sum(rdot);
    T(i) = T(i - 1) - log(rand)/rdotsum;
    rdot = rdot/rdotsum;
    reac = 1 + sum(rand > cumsum(rdot));
    Z(i, :) = Z(i - 1, :) + R(:, reac)';
end

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
    % Get birth process of y
    y = Z(:, 2);
    dy = diff(y);
    id_ydea = find(dy == -1);
    dy(id_ydea) = 0;
    ybir = cumsum([y(1); dy]);
    zbir = Z(:, 3);
    
    figure;
    hold on
    stairs(T, ybir);
    stairs(T, zbir, 'r');
    hold off
    xlabel('time');
    ylabel('molecular numbers');
    legend('y births', 'z photons', 'location', 'best');
end
    
    