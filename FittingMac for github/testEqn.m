% Some basic test code for the rate and rate estimate matrices of a Snyder
% type formulation for marked point processes
clear all
clc
close all

% Obtain rate and rate estimate matrices for truncated state space
syms q0 q1 q2 q3 q4 q5
syms e z
q = [q0 q1 q2 q3 q4 q5];
L = diag([0 0 e*(1 - z) e*(1 - z) e*(2 - z) e*(2 - z)]);
lhat = sum(q*L);
Lhat = lhat*eye(length(q));

% Calculate the differences of the rate matrices with q premutiplied
D = Lhat - L;
qD = q*D;
qx1D = sum(qD(2:2:end));
qx1D = simplify(qx1D);

% Obtain rate and rate estimate matrices for truncated state space
syms q00 q10 q01 q11
syms e z
q = [q00 q10 q01 q11];
L = diag([0 0 e*(1 - z) e*(1 - z)]);
lhat = sum(q*L);
Lhat = lhat*eye(length(q));

% Calculate the differences of the rate matrices with q premutiplied
D = Lhat - L;
qD = q*D;
qx1D = sum(qD(2:2:end));
qx1D = simplify(qx1D);