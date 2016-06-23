% ODE function that obtains symbolic expressions for the differential
% equations and then substitutes the data into the expressions
function dy = odeSymbolic(ts, y, kb1, kd1, a1, constr, dq, dqnew, q, a, kb, kd)

% Obtain symbolic ODE functional forms and substitute for q vector
% [dq sum_dq Q C B dqnew sum_dqnew q a kb kd] = testSymbolic(space);
q1 = subs(q, q, y);
q1 = double(q1);
% assignin('base', 'q', q);

% Substitute the data into the equations with account for constraint
if constr
    dy = subs(dqnew, q, q1);
else
    dy = subs(dq, q, q1);
end
assignin('base', 'dy1', dy);

% Obtain numerical value of dy with column vector form
dy = subs(dy, a, a1);
dy = subs(dy, kb, kb1);
dy = subs(dy, kd, kd1);
dy = dy';
