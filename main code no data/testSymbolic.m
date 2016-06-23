% Test code to evaluate the explicit representation of the ODEs - assumes
% that the birth and death rate are according to the Markov chain case and
% that min state = 0
function [dq sum_dq Q C B dqnew sum_dqnew q a kb kd] = testSymbolic(space)

% Define symbolic variables and type and general q
syms q0 q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 k a kb kd
q = [q0 q1 q2 q3 q4 q5 q6 q7 q8 q9 q10];

% Check space not too large for defined symbolics and min not > 0
if length(q) < max(space) || min(space) > 0
    assignin('base', 'q', q);
    assignin('base', 'space', space);
    error('Either input space too large for defined symbolics or minimum > 0');
end

% Define initial matrices
S = diag(space);
q = q(1+space);
L = a*S;
len = length(diag(S));
Lcap = a*q*S*ones(len, 1)*eye(len);

% Compute important matrices
e = Lcap - L;
Smax = max(space);
Smin = min(space);
Q = getQSymbolic(kb, kd, Smax, Smin);
C = Q + e;
B = C - Lcap;

% Obtain differential equations
dq = q*C;
sum_dq = sum(dq);

% Factorise outputs
dq = factor(dq);
sum_dq = factor(sum_dq);
Q = factor(Q);
C = factor(C);
B = factor(B);

% Substitute expression for sum(q) = 1 and assign outputs
qSt = q(Smax+1);
qSub = 1 - sum(q(1:(Smax)));
dqnew = subs(dq, qSt, qSub);
dqnew = collect(dqnew);
sum_dqnew = simplify(sum(dqnew));
assignin('base', 'qSt', qSt);
assignin('base', 'qSub', qSub);

% Sub-method to calculate the Q matrix in symbolic form
function Q = getQSymbolic(kb, kd, Smax, Smin)

% Markov linear rate on x1 births
if Smax == 1
    birR = kb*(Smax - [Smin:(Smax-1)]);
    deaR = kd*([(Smin+1):Smax] - Smin);
    Q = [-birR birR; deaR -deaR];
else
    birR = kb*(Smax - [Smin:(Smax-1)]);
    deaR = kd*([(Smin+1):Smax] - Smin);
    diagR = -[birR 0] -[0 deaR];
    Q = diag(birR, 1) + diag(deaR, -1) + diag(diagR, 0);
end