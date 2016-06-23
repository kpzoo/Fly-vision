% Modified version of getPQMx aimed at handling more general CTMCs, though
% each has to be set manually in this function with appended MCType

% Removed the dependence on birth type of x1 to simply favour certain MC
% structures which are identified by MCtype

% Function to obtain the constant Q transition matrices for different birth
% forms and exponential deaths <---------------------------------------
function [Q P Pi] = getPQGenMC(MCtype, xr_const, SlimSet, bulk)

% Decompose S into Smax and Smin for x1
Smin = SlimSet.min(1);
Smax = SlimSet.max(1);

% Check only unity bulk sizes for nearest neighbour MCs
if ismember(MCtype, [1 2]) && ~all(bulk == 1)
    assignin('base', 'bulk', bulk);
    assignin('base', 'MCtype', MCtype);
    error('The bulk strcture is non-unity for a nearest neighbour MC');
end

% Obtain matrix based on birth type - assume deaths are proportional
switch(MCtype)
    case 1
        % Tridiagonal MC which has only nearest neighbour transitions with
        % constant birth rates and linear death rates
        if Smin ~= 0
            error('Code not modified to account for Smin > 0 in this case');
        end
        kbirth = xr_const(1);
        kdeath = xr_const(2);
        Q = diag([-kbirth*ones(1, maxS) 0], 0) + diag(kbirth*ones(1, maxS), 1);
        Q = Q + diag(-kdeath*(0:maxS), 0) + diag(kdeath*(1:maxS), -1);
                
    case 2
        % Tridiagonal MC which has only nearest neighbour transitions with
        % linear restricted space birth and death rates
        kbirth = xr_const(1);
        kdeath = xr_const(2);
        birR = kbirth*(Smax - (Smin:(Smax-1)));
        deaR = kdeath*(((Smin+1):Smax) - Smin);
        diagR = -[birR 0] -[0 deaR];
        Q = diag(birR, 1) + diag(deaR, -1) + diag(diagR, 0);
        
    case 3
        % General MC with linear rates of different bulk size with r_const
        % providing all the rates <-------- only provide the x values of it
        nxReacs = length(xr_const);
        if isfinite(Smax) && (max(bulk) > Smax - Smin)
            error('State space not large enough for all r_const reactions');
        else
            % Obtain successive off diagonal components corresponding to
            % different reactions and initialise Q and indices
            ieven = 1;
            iodd = 1;
            Qlen = Smax - Smin + 1;
            Q = zeros(Qlen, Qlen);
            for i = 1:nxReacs
                kreac = xr_const(i);
                if rem(i, 2) == 0
                    % Death reactions - lower diagonals
                    deaR = kreac*(((Smin+ieven):Smax) - Smin);
                    Q = Q + diag(deaR, -ieven);
                    ieven = ieven + 1;
                else
                    % Birth reactions - upper diagonals
                    birR = kreac*(Smax - (Smin:(Smax-iodd)));                    
                    Q = Q + diag(birR, iodd);
                    iodd = iodd + 1;
                end              
            end
            % Add main diagonal based on negative sum of row elements
            diagR = -sum(Q, 2);
            Q = Q + diag(diagR, 0);
        end
           
    otherwise
        assignin('base', 'MCtype', MCtype);
        error('Unsupported MC type specified');
end

% Check if rows all sum to 0 - a requirement for intensity matrices
if any(sum(Q, 2))
    assignin('base', 'Q', Q);
    error('Q matrix incorrectly composed');
end

% Check size of Q matrix and obtain the transition probability matrix at a
% high enough time
if all(size(Q) == size(diag(Smin:Smax)))
    % Kolmogorov solution
    t = 100;
    P = expm(t*Q);
    
    % Check that time is high enough
    Psq = P^2;
    eP = abs(Psq - P);
    if max(max(eP)) > 10^-8
        warning('Mat:Pcalc', 'Specified time not large enough in P calculation');
    end
else
    assignin('base', 'Q', Q);
    error('Q matrix of incorrect size');
end

% Solve the stationary distribution of the MC and force the sum = 1
Pi = null(Q');
Pi = Pi/sum(Pi);
Pi = Pi';        