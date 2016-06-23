% Function to produce Q matrices representating the summation of
% independent and identical IPPs 
function [Q L] = getSumIPPQMx(k1, k2, lam, numIPP)

% Obtain individual generator and rate matrix
q = [-k1 k1; k2 -k2];
l = [0 0; 0 lam];

% Get the summed IPP matrices
Q = q;
L = l;
for i = 2:numIPP
    % Sequentially sum the Q and L matrices
    Q = getKronSum(Q, q);
    L = getKronSum(L, l);
end

% Check that results are of correct dimensions
if ~all(size(Q) == size(L))
    assignin('base', 'Q', Q);
    assignin('base', 'L', L);
    error('The Q and L matrices are dimensionally inconsistent');
end
[nr nc] = size(q);
NR = nr^numIPP;
NC = nc^numIPP;
if ~all(size(Q) == [NR NC])
    assignin('base', 'Q', Q);
    assignin('base', 'NR', NR);
    assignin('base', 'NC', NC);
    error('The Q (and L) matrix have incorrect dimension');
end


% Sub-function to perform Kronecker sum
function C = getKronSum(A, B)

% Get relevant identity matrices
IA = eye(size(A));
IB = eye(size(B));

% Perform Kronecker sum based on produce decomposition
C = kron(A, IB) + kron(IA, B);