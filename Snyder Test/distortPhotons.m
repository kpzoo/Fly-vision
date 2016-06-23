% Function to distort a photon stream by inserting, deleting or delaying
% photons with a specified distribution
function Test = distortPhotons(noiseTraits, Treal, T)

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
        Test = Treal(uniformDel <= probDel);
        
    case 2
        % Exponential insertion of photons using a set probability based on
        % a Poisson number of extra events
        rateIns = noiseTraits.rateIns;
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
        
    otherwise
        error('Specified case not supported');
end