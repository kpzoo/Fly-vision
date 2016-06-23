% Modified to further include a missing link of b between teh edge states
% to ensure symmetry in the Q matrix
% Modified to include extra b transitions in the opposite direction

% Function to construct Q matrix of several specific light models
function [Q Pi] = getLightModelQ(type, inpType)

% Determine the model form specified
switch(type)
    case 1
        % A constant intensity dot moves across several positions with the
        % Markov states coding both the direction and position. MArkov
        % rates are obtained below
        a = inpType.a;
        b = inpType.b;
        c = inpType.c;
        
        % Declare Q matrix size ensuring even no. of states
        posSpace = inpType.posSpace;
        M = length(posSpace);
        if rem(M, 2) ~= 0
            assignin('base', 'posSpaceErr', posSpace);
            error('Markov model for case 2 requires even no. states');
        end
        Q = zeros(M);
        
        % Set the non-diagonal Q elements in accordance with this model
        for i = 1:M
            % Basic motion in 1 direction
            if ismember(i, 1:(M/2 - 1))
                Q(i, i + 1) = a;
            end
            % Basic motion in the other direction
            if ismember(i, (M/2 + 2):M)
                Q(i, i - 1) = a;
            end
            % Transitions between directions
            if ismember(i, 2:M/2)
                Q(i, i + M/2 - 1) = b;
                Q(i + M/2 - 1, i) = b;
            end
            % This is probably the more sensible directional form <--------
            if ismember(i, 1:(M/2 - 1))
                Q(i, i + M/2 + 1) = b;
                Q(i + M/2 + 1, i) = b;
            end
        end
        % Continuity of motion in a given direction by cyclic transitions
        Q(M/2, 1) = c;
        Q(M/2 + 1, M) = c;
        
        % Further edge state transitions to ensure a symmetric Q
        Q(1, M) = b;
        Q(M, 1) = b;
        Q(M/2, M/2 + 1) = b;
        Q(M/2 + 1, M/2) = b;
        
        % Make Q conservative and obtain the stationary distribution
        Q = Q - diag(sum(Q, 2));
        Pi = null(Q');
        Pi = Pi/sum(Pi);
        
    otherwise
        assignin('base', 'typeErr', type);
        error('Model type specified is not valid');
end