% Modified to work with multi-state MCs which are specified via the Q matrix

% Function to prepare the inputs for the Snyder filter using only the
% Gillespie data specified
function inpSny = prepareSnyderInputs(simname, locfolder1, locfolder2, adjust, limPhoton)

% Extract Gillespie simulation data based on location and adjust type
if ~adjust
    cd(locfolder1);
    cd(locfolder2);
    load(simname, 'outGil');
    load(simname, 'r_const');
    load(simname, 'markov');
    load(simname, 'kg');
    load(simname, 'kb');
    load(simname, 'kd');
    load(simname, 'params'); % <--- for bioSnyder
    cd ..
    cd ..
    T = outGil.T;
    x2 = outGil.X(1:length(outGil.T), 2);
    x1 = outGil.X(1:length(outGil.T), 1);
    kgain = r_const(length(r_const));
    lam = kgain*x1;
    lamMean = outGil.rate.mean(end);
    lamVar = outGil.rate.var(end);
else
    cd(locfolder1);
    cd(locfolder2);
    load(simname, 'outGil');
    load(simname, 'markov');
    load(simname, 'kg');
    load(simname, 'kb');
    load(simname, 'kd');
    load(simname, 'params'); % <--- for bioSnyder
    cd ..
    cd ..
    T = outGil.T;
    x2 = outGil.X(1:length(outGil.T), 2);
    x1 = outGil.X(1:length(outGil.T), 1);
    r_const = outGil.r_const;
    kgain = r_const(length(r_const) - 1);
    lam = kgain*x1;
    lamMean = outGil.rate.mean(end-1);
    lamVar = outGil.rate.var(end-1);
end

% Account for missing parameters
if ~exist('markov', 'var')
    markov = 1;
end
if ~exist('params', 'var')
    params.kbirth = kb;
    params.kgain = kg;
    params.kdeath = kd;
end

% Obtain the state space of x1 and check dimensions
if ~markov
    S = diag(0:2*max(x1)); %<---------- assumes minimum value of 0
else
	if adjust
		Slim = outGil.Slim;
	else
		SlimSet = params.SlimSet;
		Slim(1) = SlimSet.min(1);
		Slim(2) = SlimSet.max(1);
	end
    S = diag(Slim(1):Slim(2));
end

% Obtain transition matrix Q which is assumed to hold when apply Snyder to
% both actual and estimated photon times
if adjust
	birStr1 = outGil.birth{1};
	[Q P] = getPQMx(birStr1, params.kbirth, params.kdeath, Slim);
else
	Q = params.Q;
	P = params.P;
end
clear outGil

% Check if 2 state symmetric MC and obtain beta and gamma (gamma is more a 
% measure of the k value here as no QBs)
if max(x1) == 1 && params.kbirth == params.kdeath
    twoStateSymm = 1;
else
    twoStateSymm = 0;
end
avgQBwidth = 100;
if adjust
	beta = params.kgain/params.kbirth; % <------- note modification from GSout to params
	gamma = 1/params.kbirth/avgQBwidth;
else
	beta = params.beta;
	gamma = params.gamma;
end

% Remove offsets on x2 and T due to transient cutoff
x2 = x2 - x2(1);
T = T - T(1);

% Limit the photon stream based on limPhoton input
if max(x2) < limPhoton
    warning('Mat:fewPh', ['There are less than ' num2str(limPhoton) ' photons present ']);
else
    disp(['Clipping Gillespie data to ' num2str(limPhoton) ' photons']);
    idlim = find(x2 == limPhoton, 1, 'last');
    x2 = x2(1:idlim);
    x1 = x1(1:idlim);
    T = T(1:idlim);
end

% Assign output structure
inpSny.x1 = x1;
inpSny.x2 = x2;
inpSny.T = T;
inpSny.lam = lam;
inpSny.P = P;
inpSny.Q = Q;
inpSny.S = S;
inpSny.params = params;
inpSny.beta = beta;
inpSny.gamma = gamma;
inpSny.twoStateSymm = twoStateSymm;
