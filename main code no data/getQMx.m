% Function to obtain the constant Q transition matrices for different birth
% forms and and exponential deaths <---------------------------------------
function [Q Qt] = getQMx(birStr, kbirth, kdeath, maxS)

% Birth rate options - must match other code <-----------------------------
birSet = {'const', 'linear'};
birType = strmatch(birStr, birSet, 'exact');
if isempty(birType)
    error('Input birth rate does not match any birth type');
end

% Obtain matrix based on birth type - assume deaths are proportional
switch(birType)
    case 1
        % Constant rate of x1 birth
        Q = diag([-kbirth*ones(1, maxS) 0], 0) + diag(kbirth*ones(1, maxS), 1);
        Q = Q + diag(-kdeath*(0:1:maxS), 0) + diag(kdeath*(1:maxS), -1);
        
    case 2  %<------ not debugged yet
        % Linear rate of x1 births - probably useless <--------------------
        Q = diag(-(kbirth + kdeath)*(0:-1:maxS), 0);
        Q = Q + diag(kbirth*ones(1, maxS), 1) + diag(kdeath*(1:maxS), -1);
end

% Check if rows all sum to 0
if any(sum(Q, 2))
    error('Q matrix incorrectly composed');
end

% Qt matrix is in agreement with Neil's form so obtained as transpose
Qt = Q';

        
        