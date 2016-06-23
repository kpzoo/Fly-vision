% Script to construct a Q matrix for a given stationary state distribution
% Pi, the code assumes that Pi has no zero entries <----------------------
clear all
clc
close all

% Set boolean for tridiagonal or general Q and type of Pi
struc = 4;
nState = 20;
statdistr = 2;
states = 0:nState-1;

% Construct different desired stationary distributions
switch(statdistr)
    case 1
        % Probabilities obtained from random numbers
        Pi = rand(1, nState);
    case 2
        % Bimodal distribution obtained from 2 Gaussian distributions
        if nState <= 5
            error('Not enough states for sensible bimodal form');
        end
        mu1 = nState/4;
        sigma1 = mu1/4;
        mu2 = 3*nState/4;
        sigma2 = mu2/12;
        Pi = normpdf(states, mu1, sigma1) + normpdf(states, mu2, sigma2);
    otherwise
        error(['No distribution available for statdistr = ' num2str(statdistr)]);
end

% Ensure Pi sums to 1 and choose the unconstrainted Q entries as a
Pi = Pi/sum(Pi);
a = 100*rand(1, nState-1);

% Ensure no zero Pi entries
if any(Pi == 0)
    error('No support for Pi with zero entries');
end

% Construct a Q matrix based on required structure
switch(struc)
    
    case 1
        % Solution for tridiagonal Q matrix
        r = Pi(1:end-1)./Pi(2:end);
        b = r.*a;
        Q = diag(a, 1) + diag(b, -1);
        Q = Q - diag(sum(Q, 2));
        
    case 2
        % Solution for tridiagonal with two corner entries - first
        % calculate homogeneous general solution
        r = Pi(1:end-1)./Pi(2:end);
        b = r.*a;
        
        % Calculate extra parameters from settings
        delta = 10;
        epsn = 100;
        eps1 = (delta + epsn*Pi(end))/Pi(1);
        
        % Obtain final b vector - particular solution
        b = b + delta./Pi(2:end);
        Q = diag(a, 1) + diag(b, -1);
        Q(1, nState) = eps1;
        Q(nState, 1) = epsn;
        Q = Q - diag(sum(Q, 2));
        
    case 3
        % Solution in which the b vector is given with kx death form and
        % the a vector are the unknowns and eps1 and epsn exist with epsn
        % also of a linear in x form
        kdeath = 1000;
        keps = 1;
        b = kdeath*states(2:end);
        epsn = keps*states(end);
        delta  = 0.01;
        eps1 = (delta + epsn*Pi(end))/Pi(1);
        
        % Check that the specified b vector will always yield positive a
        baug = b - delta./Pi(2:end);
        if any(baug < 0)
            error('Need to increase b values or decrease delta');
        end
        
        % Obtain the a vector by inverting the case 2 equation
        r = Pi(1:end-1)./Pi(2:end);
        rinv = 1./r;
        a = rinv.*baug;
        Q = diag(a, 1) + diag(b, -1);
        Q(1, nState) = eps1;
        Q(nState, 1) = epsn;
        Q = Q - diag(sum(Q, 2));
        
    case 4
        % General matrix solution - iterative with a as input and b unknown
        % and the extra parameters expressed in matrix that has zeros for
        % the unknown areas
        Qinit = 10*rand(nState);
        Qinit = Qinit - diag(diag(Qinit));
        
        % Having set diag to zero also set the b components to zero - these
        % also simplify the sumR calculation
        for i = 2:nState
            Qinit(i, i-1) = 0;
        end
        Qtemp = Qinit;
        
        % Having obtained input matrix solve for b vector iteratively
        b = zeros(1, nState-1);
        for i = 1:nState-1
            sumR = sum(Qtemp(i, :));
            prodC = Pi*Qtemp(:, i);
            b(i) = (Pi(i)*sumR - sum(prodC))/Pi(i+1);
            
            % Check for values
            if (Pi(i)*sumR - sum(prodC)) < 0
                error(['Bad parameters with negative b at i = ' num2str(i)]);
            end
            
            Qtemp(i+1, i) = b(i);
            Qtemp(i, i) = -sumR;
        end
        
        % Assign Q matrix and check b positive
        Qtemp(nState, nState) = -sum(Qtemp(nState, :));
        Q = Qtemp;
        if any(b <= 0)
            error('Output vector b does not have all positive values');
        end

                
    otherwise
        error(['Q structure of struc = ' num2str(struc) ' not supported']);
end

% Check that the Q matrix is properly conditioned
if any(diag(Q) > 0)
    error('Q matrix has a positive diagonal element');
end
sumRow = sum(sum(Q, 2));
if sumRow > 10^-10
    error('The row sums are not close enough to 0');
end
Qnodiag = Q - diag(diag(Q));
if min(min(Qnodiag)) < 0
    error('Off diagonal elements are not all non-negative');
end

% Check that the Q really has Pi as its null space
nu = Pi*Q;
sum_nu = sum(nu);
if sum_nu < 10^-10
    disp('Estimated Q matrix seems correct');
end

% Other checking method with plot
PiEst = null(Q');
PiEst = PiEst/sum(PiEst);
PiEst = PiEst';

% Plot estimate against original
figure;
plot(states, [Pi' PiEst']);
xlim([states(1) states(end)]);
xlabel('states');
ylabel('stationary distributions');
legend('Pi', 'Pi estimate', 'location', 'best');
title('Comparison of actual and estimated stationary solution');

% Plot the a and b vectors against states
figure;
subplot(2, 1, 1);
plot(states, [a 0], 'o-');
xlabel('states');
ylabel('birth rate with +1 transitions');
subplot(2, 1, 2);
plot(states, [0 b], 'o-');
xlabel('states');
ylabel('death rate with -1 transitions');