% Setting nIntens is separate from that in inpType so that one can obtain a
% reduced model that has no intensity connections by inputting a 0 here -
% useful to obtain the model for Qphi

% Function to construct Q matrix of position-direction-intensity light models
function [Q Pi] = getQLightExp(inpType, testDir, nIntens)

% Obtain test forms
testDir1 = testDir(1);
testDir2 = testDir(2);

% A variable intensity dot moves across several positions with the
% Markov states coding direction, intensity and position. Markov
% rates are set below
a = inpType.a;
b = inpType.b;
c = inpType.c;
h = inpType.h;

% Declare Q matrix size ensuring even no. of positional states -
% this Q matrix will be constructed in block form from component Qs
posDirSpace = inpType.posDirSpace;
M = length(posDirSpace);
MI = nIntens;
if rem(M, 2) ~= 0
    assignin('base', 'posDirSpaceErr', posDirSpace);
    error('Markov model for case 2 requires even no. states');
end

% Construct the P1 matrix which controls the connections between the
% inner MCs which have generators of Q1
P1 = h*eye(M);
Q1 = zeros(M);

% Set the non-diagonal Q1 elements in accordance with this model -
% note that Q1 is not conservative so no diagonal elements set yet
for i = 1:M
    % Basic motion in 1 direction
    if ismember(i, 1:(M/2 - 1))
        Q1(i, i + 1) = a;
    end
    % Basic motion in the other direction
    if ismember(i, (M/2 + 2):M)
        Q1(i, i - 1) = a;
    end
    % Transitions between directions
    if ismember(i, 2:M/2)
        Q1(i, i + M/2 - 1) = b;
        Q1(i + M/2 - 1, i) = b;
    end
    % This is probably the more sensible directional form <--------
    if ismember(i, 1:(M/2 - 1)) && ~any(testDir)
        Q1(i, i + M/2 + 1) = b;
        Q1(i + M/2 + 1, i) = b;
    end
end
% Continuity of motion in a given direction by cyclic transitions
if testDir2
    Q1(M/2, 1) = c;
    Q1(M/2 + 1, M) = c;
end

% Further edge state transitions to ensure a symmetric Q
if testDir1
    Q1(1, M) = b;
    Q1(M, 1) = b;
    Q1(M/2, M/2 + 1) = b;
    Q1(M/2 + 1, M/2) = b;
end

% Construct the complete Q matrix which assumes only 2 dimensions
% to the MC and check dimensions
numQblks = MI;
Q = blktridiag(Q1, P1, P1, numQblks);
Q = full(Q);
if ~all(size(Q) == [M*MI M*MI])
    assignin('base', 'Q', Q);
    assignin('base', 'Q1', Q1);
    assignin('base', 'P1', P1);
    error('Block Q matrix has incorrect dimensions');
end

% Make Q conservative and obtain the stationary distribution
Q = Q - diag(sum(Q, 2));
Pi = null(Q');
Pi = Pi/sum(Pi);
