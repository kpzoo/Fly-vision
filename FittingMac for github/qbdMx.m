% Script to check QBD state and generator matrices for a simple 2 state
% exponential delay due to y with x as the usual symmetric MC
% Code assumes a total of 4 joint states due to y also being 2 state
clear all
clc
close all

% Set mode - numerical values or purely symbolic
sym_on = 0;

% Rates for MCs due to x and y and specify number of states as m, n
if sym_on
    syms k a b
else
    a = 0.01;
    b = 10^10;
    k = 5e-4;
end
m = 2;
n = 2;

% Block matrices of interest (A course on queueing models - Jain)
A = [-k k; k -k];
A0 = diag([0, a]);
A2 = diag([b b]);
A1 = A - (A0 + A2);

% Assemble generator matrix and check row sums and if Q is conservative -
% modification for 2 state made by introducing A0 in (4,4) block term
Q = [A1 + A2 A0; A2 A1 + A0];
rowSum = sum(Q');
cons1 = (A1 + A2 + A0)*ones(m, 1);
cons2 = (A2 + A1 + A0)*ones(m, 1);

% Calculate stationary distribution
Pi = null(Q');
Pi = Pi./sum(Pi);

% Set prior probabilities and then marginalise
p = Pi';
px0 = p(1) + p(3);
px1 = p(2) + p(4);
py0 = p(1) + p(2);
py1 = p(3) + p(4);
prior = [px0 px1 py0 py1];

% Set posterior probabilities
syms q00 q01 q10 q11 dq00 dq01 dq10 dq11
q = [q00 q10 q01 q11];
dq = [dq00 dq10 dq01 dq11];
qx0 = q(1) + q(3);
qx1 = q(2) + q(4);

% Produce lam and lamhat matrices
I = eye(m*n);
lam = diag([0 0 b b]);
lamhat = q*lam*ones(m*n, 1)*I;

% Obtain Snyder continuous time equation (between jumps) and marginalise
dq = q*(Q + lamhat - lam);
G = Q + lamhat - lam;
dqx0 = dq(1) + dq(3);
dqx1 = dq(2) + dq(4);
dqy0 = dq(1) + dq(2);
dqy1 = dq(3) + dq(4);

% Simplify x posterior equations and check probability conditions met
dqx0 = factor(dqx0);
dqx1= factor(dqx1);
sumqx = factor(dqx0 + dqx1);


% Attempt to solve the linear version of the posterior solution with
% appropriate transformations
Qlin = Q - lam;
if sym_on
    syms t
    pstar = [0 0 0 1]*expm(Qlin*t);
    qlin = pstar/(sum(pstar));
    qlin = simplify(qlin);
    qx1lin = qlin(2) + qlin(4);
    qx1lin = simplify(qx1lin);
    qratio = qlin(2)/qlin(4);
    qratio = simplify(qratio);
else
    % Initialise trajectory variables
    t = linspace(0, 50/a, 1000);
    lent = length(t);
    pstar = zeros(lent, m*n);
    qlin = zeros(lent, m*n);
    x1lin = zeros(lent, 1);
    qratio = zeros(lent, 1);
    
    % Obtain linear x1 estimated trajectories
    for i = 1:lent
        pstar(i, :) = [0 0 0 1]*expm(Qlin*t(i));
        qlin(i, :) = pstar(i, :)/(sum(pstar(i, :)));
        x1lin(i) = qlin(i, 2) + qlin(i, 4);
        qratio(i) = qlin(i, 2)/qlin(i, 4);
    end
    
    % Compare qratio value to desired ratio
    qratio = qratio(2:end);
    if all(qratio - qratio(end) < 10^-4)
        rr = qratio(end)/(b/a);
        disp(['Factor of difference between limits is ' num2str(rr)]);
    else
        disp('The factor varies with time');
    end
    
    % Obtain expected x1 estimate when delay is zero
    A = k/a + 0.5;
    B = sqrt((k^2)/(a^2) + 0.25);
    C = atanh((A - 1)/B);
    x1exp = A - B*tanh(B*a*t + C);
    x1exp = x1exp';
    xr = x1lin./x1exp;
    
    % Calculation of when delay is zero with linear equations
    Q_no_del = [-k k; k -k];
    lam_no_del = [0 0; 0 a];
    Qlin1 = Q_no_del - lam_no_del;
    pstar1 = zeros(lent, 2);
    qlin1 = zeros(lent, 2);
    x1exp1 = zeros(lent, 1);
    qratio1 = zeros(lent, 1);
    for i = 1:lent
        pstar1(i, :) = [0 1]*expm(Qlin1*t(i));
        qlin1(i, :) = pstar1(i, :)/(sum(pstar1(i, :)));
        x1exp1(i) = qlin1(i, 2);
        qratio1(i) = qlin1(i, 2)/qlin1(i, 1);
    end
    norm_diff = norm(x1exp - x1exp1);
    disp(['The norm difference between pure Snyder from theory and linear ODEs is ' num2str(norm_diff)]);
    
    % Compare curves for solution between events
    figure;
    plot(t, x1lin, t, x1exp);
    xlabel('t');
    ylabel('x1 estimates');
    legend('delayed Snyder', 'pure Snyder', 'location', 'best');
    title(['Comparison of x1 solutions between events, [alpha k eta] = ' [num2str(a) ' ' num2str(k) ' ' num2str(b)]]);
    
    figure;
    plot(t, xr);
    xlabel('t');
    ylabel('ratio of solutions');
    title(['Ratio of delayed Snyder to pure Snyder estimate, [alpha k eta] = ' [num2str(a) ' ' num2str(k) ' ' num2str(b)]]);
    
end


