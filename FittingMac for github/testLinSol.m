% Basic code to test the linear solution to the delayed Snyder scheme for
% the trajectory beteen events
clear all
clc
close all

% Set the initial parameters
k = 5e-5;
a = 0.01*(10^5);
beta = a/k;
b = 100000;
rho = a/b;
disp(['Utilisation is ' num2str(rho) ', beta is ' num2str(beta)]);

% Set the dimension of the state space and z input
m = 50;
z = 0;
dimQ = 2*m;
Q = zeros(dimQ, dimQ);
rep = dimQ/2;

% Obtain the Q matrix
d1 = repmat([k 0], 1, rep);
d1 = d1(1:dimQ - 1);
d1neg = d1;
d2 = repmat([0 a], 1, rep);
d2 = d2(1:dimQ - 2);
Q = Q + diag(d1, 1) + diag(d2, 2) + diag(d1neg, -1);
d0 = -sum(Q, 2);
Q = Q + diag(d0);

% Obtain the lam matrix and its updated form
lam = getLamZFn(z, dimQ, b);
lamUp = getLamZFn(z+1, dimQ, b);

% Set variables for linear solution
t = linspace(0, 50/a, 1000);
lent = length(t);
pstar = zeros(lent, dimQ);
q = zeros(lent, dimQ);
x = zeros(lent, 1);
qupdate = q; % assume updated to z+1
xupdate = x;

% Choose the initial conditions and initialise betaSums
% q0 = (1/dimQ)*ones(1, dimQ);
% q0 = [0 1 zeros(1, dimQ-2)];
q0 = zeros(1, dimQ);
q0(1:2:end) = 0;
q0(2:2:end) = 2/dimQ;
betaSum0 = zeros(lent, 1);
betaSum1 = zeros(lent, 1);
betaSumRatio = zeros(lent, 1);

% Implement the linear solution
Qlin = Q - lam;
eigen = eig(Qlin);
if all(eigen <= 0)
    disp('No positive eigenvalues');
end
Dlam = diag(lam);
for i = 1:lent
    % Obtain the normalised exponential x solution
    pstar(i, :) = q0*expm(Qlin*t(i));
    q(i, :) = pstar(i, :)/(sum(pstar(i, :)));
    x(i) = sum(q(i, 2:2:end));
    
    % Calculate updates for next event
    qupdate(i, :) = q(i, :)*lamUp/(q(i, :)*lamUp*ones(dimQ, 1));
    xupdate(i) = sum(qupdate(i, 2:2:end));
    
    % Obtain the beta sum values
    betaSum0(i) = q(i, 1:2:end)*Dlam(1:2:end);
    betaSum1(i) = q(i, 2:2:end)*Dlam(2:2:end);
    betaSumRatio(i) = betaSum1(i)/(betaSum1(i) + betaSum0(i));
end


% Obtain expected x1 estimate when delay is zero and eta high
A = k/a + 0.5;
B = sqrt((k^2)/(a^2) + 0.25);
C = atanh((A - 1)/B);
xexp = A - B*tanh(B*a*t + C);
xexp = xexp';
xr = x./xexp;
normDiff = norm(x - xexp);
disp(['The normed difference is ' num2str(normDiff)]);

% Plot the comparison 
figure;
plot(t, x, t, xexp);
xlabel('time');
ylabel('x trajectories');
legend('linear delayed solution', 'pure analytical solution', 'location', 'best');
title(['Comparison of x trajectory at eta = ' num2str(b) ' for ' num2str(m) ' y states']);
