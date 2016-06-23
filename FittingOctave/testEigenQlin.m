% Simple code to check the eigenvalues of the Qlin matrix at successive z
clear all
clc
close all

% Set the initial parameters
k = 5e-4;
a = 0.01;
beta = a/k;
b = 100;
rho = a/b;
disp(['Utilisation is ' num2str(rho) ', beta is ' num2str(beta)]);

% Set the dimension of the state space and z input
m = 10;
dimQ = 2*m + 2; % ensure max(z) < 0.5*dimQ <-----------------------
Q = zeros(dimQ, dimQ);
rep = dimQ/2;

% Obtain the Q matrix
d1 = repmat([k 0], 1, rep);
d1 = d1(1:dimQ - 1);
d1neg = d1;
d2 = repmat([0 a], 1, rep);
d2 = d2(1:dimQ - 2);
Q = Q + diag(d1, 1) + diag(d2, 2) + diag(d1neg, -1);
d0 = -sum(Q, 2);
Q = Q + diag(d0);

% Set the z values for consideration and initialise variables
z = 0:m;
lenz = length(z);
QLin = cell(1, lenz);
eigen = zeros(dimQ, lenz);
posEig = zeros(1, lenz);
eigenSub = cell(1, 1);
posEigSub = zeros(1, lenz);

% Loop through z values and obtain lam, Qlin and eigenvalues
for i = 1:lenz
    % Obtain lam, Qlin and the full eigenvalue set
    lam = getLamZFn2(z(i), dimQ, b); % <----- removed negative values
    Qlin = Q - lam;
    eigen(:, i) = eig(Qlin);
    posEig(i) = length(find(eigen(:, i) > 0));
    if all(eigen(:, i) <= 0)
        disp(['No positive eigenvalues at i = ' num2str(i)]);
    end
    % Attempt to remove lower states than z so use a Qlin submatrix
    % Qsub = Qlin(z(i)+1:end, z(i)+1:end);
	lamSub = diag(lam);
	lamSub(lamSub < 0) = 0;
	lamSub = diag(lamSub);
	Qsub = Q - lamSub;
    eigenSub{i} = eig(Qsub);
    posEigSub(i) = length(find(eigenSub{i} > 0));
end