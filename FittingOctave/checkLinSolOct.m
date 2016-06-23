% Test of linear solution for the delayed Snyder form between events
clear all
clc
close all

% Set variables for initialisation assuming z = 0
a = 0.01;
k = 5e-4;
beta = a/k;
t = linspace(0, 50/a, 1000);

% Set the state size, eta = b and the z input
z = 0;
m = 10;
b = 100;
dimQ = 2*m;
rep = dimQ/2;

% Construct the Q matrix
Q = zeros(dimQ, dimQ);
d1 = repmat([k 0], 1, rep);
d1 = d1(1:dimQ - 1);
d1neg = d1;
d2 = repmat([0 a], 1, rep);
d2 = d2(1:dimQ - 2);
Q = Q + diag(d1, 1) + diag(d2, 2) + diag(d1neg, -1);
d0 = -sum(Q, 2);
Q = Q + diag(d0);

% Obtain lam matrix
lamDiag = zeros(1, dimQ);
iV = 3;
for iL = 2:dimQ/2
   lamDiag(iV:iV+1) = b*((iL-1) - z);
   iV = iV + 2;
end 
lam = diag(lamDiag);

% Obtain the linearised Q matrix and set the initial q
Qlin = Q - lam;
lent = length(t);
q = zeros(dimQ, lent);
eigen = eig(Qlin);
% q0 = (1/dimQ)*ones(1, dimQ);
q0 = [0 1 zeros(1, dimQ-2)];

% Obtain trajectory across time
for i = 1:lent
	q(:, i) = q0*expm(Qlin*t(i));
	q(:, i) = q(:, i)/sum(q(:, i));
end
x = sum(q(2:2:end, :), 1);

% Analytical solution for high eta
A = k/a + 0.5;
B = sqrt((k^2)/(a^2) + 0.25);
C = atanh((A - 1)/B);
xexp = A - B*tanh(B*a*t + C);
xr = xexp./x;
normDiff = norm(x - xexp);
disp(['The normed difference between linear and analytic solutions is '...
 num2str(normDiff)]);


% Plot x trajectories
plot(t, x, t, xexp);
xlabel('time');
ylabel('x(t) between events');
legend('linear delayed', 'analytic pure', 'location', 'best');

	