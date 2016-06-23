% Modification to include the MArkov birth types on x1

% Function to obtain the constant Q transition matrices for different birth
% forms and and exponential deaths <---------------------------------------
function [Q Qt] = getQMxMarkov(birStr, kbirth, kdeath, Slim)

% Birth rate options - must match other code <-----------------------------
birSet = {'const', 'linear', 'birMarkLin'};
birType = strmatch(birStr, birSet, 'exact');
if isempty(birType)
    error('Input birth rate does not match any birth type');
end

% Decompose S into Smax and Smin
Smin = Slim(1);
Smax = Slim(2);

% Obtain matrix based on birth type - assume deaths are proportional
switch(birType)
    case 1
        % Constant rate of x1 birth - assume Slim(1) = 0
        if Smin ~= 0
            error('Code not modified to account for Smin > 0 in this case');
        end
        Q = diag([-kbirth*ones(1, maxS) 0], 0) + diag(kbirth*ones(1, maxS), 1);
        Q = Q + diag(-kdeath*(0:maxS), 0) + diag(kdeath*(1:maxS), -1);
        
    case 2  %<------ not debugged yet
        % Linear rate of x1 births - probably useless <--------------------
        Q = diag(-(kbirth + kdeath)*(0:-1:maxS), 0);
        Q = Q + diag(kbirth*ones(1, maxS), 1) + diag(kdeath*(1:maxS), -1);
        
    case 3
        % Markov linear rate on x1 births
        birR = kbirth*(Smax - [Smin:(Smax-1)]);
        deaR = kdeath*([(Smin+1):Smax] - Smin);
        diagR = -[birR 0] -[0 deaR];
        Q = diag(birR, 1) + diag(deaR, -1) + diag(diagR, 0);
end

% Check if rows all sum to 0 - a requirement for intensity matrices
if any(sum(Q, 2))
    assignin('base', 'Q', Q);
    error('Q matrix incorrectly composed');
end

% Check size of Q matrix and obtain the transition probability matrix at a
% high enough time (assumed) <-------------------------------------------
if all(size(Q) == size(diag(Smin:Smax)))
    Qt = Q';
else
    assignin('base', 'Q', Q);
    error('Q matrix of incorrect size');
end

        
        