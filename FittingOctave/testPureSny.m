% Check of linear solution for the pure Snyder form between events
clear all
clc
close all

% Set variables for initialisation 
q0 = [0 1];
t = 0:0.1:1000;
a = 0.1;
k = 5e-4;
beta = a/k;
Q = [-k k; k -k-a];
lent = length(t);
q = zeros(2, lent);
eigen = eig(Q);

% Obtain trajectory across time
for i = 1:lent
	q(:, i) = q0*expm(Q*t(i));
	q(:, i) = q(:, i)/sum(q(:, i));
end
x = q(2, :);

% Analytical solution
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
legend('linear', 'analytic', 'location', 'best');

	