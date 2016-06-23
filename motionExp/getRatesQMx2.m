% Generalised to account for non-nearest neighbour reactions and set so
% that a birth death pair is obtained from the Q matrix for each possible
% reaction in order to match r_const

% Function that supplies rates based on the Q matrix of a Markov chain
% which represents x1's behaviour
function rdotx1 = getRatesQMx2(Qx1, x1, nx1Reacs)

% Essentially a lookup table on the Q matrix of x1 with the row
% corresponding to state x1 representing its reactions
nState = length(Qx1);
state = x1 + 1; % <----------------- x1 starts from 0 = 1st row
if state < 1
    assignin('base', 'x1', x1);
    error('Incorrect state value');
end
Qrow = Qx1(state, 1:nState);

% Construct the rdot variable based on the Q matrix and current state by
% finding values left and right of the diagonal and accounting for
% dimensions of matrix in determining how many values exist
diagLoc = state;
rdotx1 = zeros(1, nx1Reacs);
if diagLoc + nx1Reacs/2 <= nState
    Qbir = Qrow(diagLoc+1:diagLoc+nx1Reacs/2);
else
    Qbir = Qrow(diagLoc+1:end);
end
if diagLoc - nx1Reacs/2 >= 1
    Qdea = Qrow(diagLoc-1:-1:diagLoc-nx1Reacs/2);
else
    Qdea = Qrow(diagLoc-1:-1:1);
end
Qbirlen = length(Qbir);
Qdealen = length(Qdea);

% Append zeros to account for all possible reactions
desiredLen = nx1Reacs/2;
if Qbirlen < desiredLen
    nPad = desiredLen - Qbirlen;
    Qbir = [Qbir zeros(1, nPad)];
end
if Qdealen < desiredLen 
    nPad = desiredLen - Qdealen;
    Qdea = [Qdea zeros(1, nPad)];
end

% Check consistent lengths and assign some check variables
if length(Qdea) ~= length(Qbir)
    error('The birth and death lengths are not consistent');
end
% assignin('base', 'Qbir', Qbir);
% assignin('base', 'Qdea', Qdea);
% assignin('base', 'Qrow', Qrow);

% Assingn the rdot values with account for reactions which cannot occur at
% the given state and maintain [bir dea bir dea...] form
for i = 1:nx1Reacs
    if rem(i, 2) == 1
        % Odd reactions are births
        rdotx1(i) = Qbir((i+1)/2);
    else
        % Even reactions are deaths
        rdotx1(i) = Qdea(i/2);
    end
end