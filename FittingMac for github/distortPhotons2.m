% Modified to include case that includes both the empirical delay and a
% chosen rate of event deletion
% Modified to include a case to allow for the empirical QB distribution

% Function to distort a photon stream by inserting, deleting or delaying
% photons with a specified distribution
function Test = distortPhotons2(noiseTraits, Treal, T, kg)

% Ensure that Treal is a column vector
if size(Treal, 2) ~= 1
    error('Input photon stream is not a column vector');
end

% Determine distortion form based on specified method
switch(noiseTraits.meth)
    case 1
        % Bernoulli deletion of photons using a set probability
        uniformDel = rand(size(Treal));
        probDel = noiseTraits.probDel;
        Test = Treal(uniformDel <= (1 - probDel));
        
    case 2
        % Exponential insertion of photons using a set probability based on
        % a Poisson number of extra events
        probIns = noiseTraits.probIns;
        rateIns = probIns*kg/2;
        nExtra = poissrnd(rateIns*max(T));
        expIns = exprnd(1/rateIns, 1, nExtra);
        TIns = cumsum(expIns);
        Test = [Treal; TIns'];
        Test = sort(Test);
        
    case 3
        % Additon of an extra delay based on a set distribution
        switch(noiseTraits.delayDistr)
            case 1
                % Exponential delays of given mean rate
                paramDistr = noiseTraits.paramDistr;
                if length(paramDistr) ~= 1
                    assignin('base', 'noiseTraits', noiseTraits);
                    error('Incorrect number of distribution parameters specified');
                end
                rateDelay = paramDistr(1);
                expDelay = exprnd(1/rateDelay, size(Treal));
                
                % Delays may rearrange photons so sort
                Test = Treal + expDelay;
                Test = sort(Test);
                
            otherwise
                disp('Specified delay distribution not supported');
        end
        
    case 4
        % Addition of an extra delay based on the empirical latency
        % distribution of the QBs
        empirData = load('figdata', 'xdata', 'ydata');
        
        % After loading data sample from the empirical distribution
        delay = drawEmpiricalDistr(empirData.xdata, empirData.ydata, length(Treal));
                
        % Delays may rearrange photons so sort
        Test = Treal + delay';
        Test = sort(Test);
        
    case 5
        % Obtain empirically delayed photons
        empirData = load('figdata', 'xdata', 'ydata');
        delay = drawEmpiricalDistr(empirData.xdata, empirData.ydata, length(Treal));
        Test = Treal + delay';
        Test = sort(Test);
        
        % Delete a fixed percent of the delayed photons
        uniformDel = rand(size(Test));
        probDel = noiseTraits.probDel;
        Test = Test(uniformDel <= (1 - probDel));
        
    otherwise
        error('Specified case not supported');
end