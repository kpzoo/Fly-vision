% Modified from just gillespieFilter to account for non-linear lam(x1) in
% accordance with runDSPPfilter4

% Simplified Gillespie algorithm for use with filter algorithms - no
% modulation types on x1 is included as it is only required that x1
% modulate the intensity of x2
function [X Xdot T] = gillespieFilter4(inpGill)

% Deconstruct inputs into individual variables
len = inpGill.len;
r = inpGill.r;
x0 = inpGill.x0; 
N = inpGill.N;
Nstart = inpGill.Nstart;
deaType = inpGill.deaType;
birType = inpGill.birType;
avgR = inpGill.avgR;
maxR = inpGill.maxR;
minR = inpGill.minR;
r_const = inpGill.r_const;
bulk_size = inpGill.bulk_size;

% Set initial parameters for time and populations
t = zeros(N, 1);
z = zeros(N, len);
z(1, :) = x0;

% Extra state information required by birfoll
state.setmax = zeros(N, len);
state.numMiss = zeros(N, len);

% Set initial rates and reaction incrementation
alpha = zeros(N, 2*len);
rdot = zeros(1, 2*len);

% Control reaction increments with bulk sizes
transit = zeros(len, 2*len);
for i = 1:len
    transit(i, 2*i-1:2*i) = [1*bulk_size(i) -1*bulk_size(i)];
end

% Loop across specified number of iterations and simulate
for i = 2:N
    
    % Update molecular numbers and other data sequentially
    x = z(i - 1, :);
    told = t(i-1);
    
    % Account for time history even when none is available
    if i > 2
        xhist = z(1:i-2, :);
        rprev = rdot;
    else
        xhist = zeros(2, len);
        rprev = zeros(1, 2*len);
    end
    
    % Obtain reaction rates for each species
    for j = 1:len
        % Assign death rates
        switch(deaType(j))
            case 1
                % Mass action
                rdot(2*j) = x(j)*r_const(2*j);
            otherwise
                % Zero rate
                rdot(2*j) = 0;
        end
        % Assign birth rates
        switch(birType(j))
            case 1
                % Constant rate --> set for input process
                rate_x1 = avgR(j);
                rdot(2*j - 1) = rate_x1;
                
            case 2
                % Birth following - assumes jth molecule follows j-1th and
                % uses last numMiss value ---> set for estimator only
                x1prev = xhist(end, j-1);
                x2prev = xhist(end, j);
                x1_del = x(j-1) - x1prev;
                x2_del = x(j) - x2prev;
                numMiss = state.numMiss(i-1, j);
                
                % Respond to births on x1
                if x1_del > 0
                    numMiss = numMiss + 1;
                end
                
                % Respond to births on x2
                if x2_del > 0
                    numMiss = numMiss - 1;
                end
                
                % Respond to numMiss - set rate of zero if numMiss <= 0
                if numMiss > 0
                    setmax = 1;
                else
                    setmax = 0;
                end
                
                % Set rates according to state and update numMiss
                if setmax
                    rate_x2 = maxR(j);
                else
                    rate_x2 = 0;
                end
                rdot(2*j - 1) = rate_x2;                
                state.numMiss(i, j) = numMiss;
                
            case 3
                % Proportional cross rate - set for estimator based on x1
                % (assumed to be (j-1)th molecule)
                rate_x2 = x(j - 1)*r_const(2*j - 1);
                rdot(2*j - 1) = rate_x2;
                
            case 4
                % Linear rate - set for input x1
                rate_x1 = x(j)*r_const(2*j - 1);
                rdot(2*j - 1) = rate_x1;
                
            case 5
                % Square cross rate - set for estimator based on x1
                % (assumed to be (j-1)th molecule)
                rate_x2 = (x(j - 1)^2)*r_const(2*j - 1);
                rdot(2*j - 1) = rate_x2;
                
            case 6
                % Negative exponential cross rate - set for estimator of x1
                rate_x2 = exp(-x(j - 1))*r_const(2*j - 1);
                rdot(2*j - 1) = rate_x2;
                
            otherwise
                % Zero rate
                rdot(2*j - 1) = 0;                
        end
        
        
    end
    
    % Obtain characteristics of the next reaction
    rdotsum = sum(rdot);
    tnex = told - log(rand)/rdotsum;
    rdot_ratio = rdot/rdotsum;
    reac = 1 + sum(rand > cumsum(rdot_ratio));
    xnex = x + transit(:, reac)';
        
    % Assign variables
    alpha(i-1, :) = rdot;
    t(i) = tnex;
    z(i, :) = xnex;
    
    % Catch a possible error
%     if rdotsum == 0
%         assignin('base', 'z', z(1:i, :));
%         assignin('base', 'alpha', alpha(1:i-1, :));
%         error(['All rates are zero at i = ' num2str(i)]);
%     end
           
end

% Save data from simulation for post processing and account for the control
x = z;
xdot = alpha;

% Account for equilibrium and remove last value (in previous code would
% calculate the next step at i = N)
X = x(Nstart:N-1, :);
Xdot = xdot(Nstart:N-1, :);
T = t(Nstart:N-1, :);
