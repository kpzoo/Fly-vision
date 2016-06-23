% Script to construct a Q matrix for a given stationary state distribution
% Pi, the code assumes that Pi has no zero entries <----------------------
clear all
clc
close all

% Set boolean for tridiagonal or general Q and type of Pi
struc = 2;
nState = 10;
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

                
    otherwise
        error(['Q structure of struc = ' num2str(struc) ' not supported']);
end

% Check that the Q matrix is properly conditioned
if any(diag(Q) > 0)
    error('Q matrix has a positive diagonal element');
end
if sum(sum(Q, 2)) > 10^-10
    error('The row sums are not close enough to 0');
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