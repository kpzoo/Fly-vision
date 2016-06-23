% Modified to get the estimate of the positions in addition to the states
% and the direction and includes input calcBools over calcDirec only

% Modified to calculate the directional estimate differently depending on
% whether a reverse phi model or normal model is applied - uses a
% subfunction to perform these computations (extra input modType)

% Modified version of generalSnyFilterBasic to include the effect of a
% deterministic delay - assumes the Tphot and Tevent inputs are delayed
% <----------------------------------------------------------------------

% Removed the octave flexibility and the ode_sw and int_sw functionality
% and forced use of linear ODE scheme alone (no non-linear solution)

% Function used to implement the Snyder filtering when faced with M sensory
% cells that result in M observed processes - based on paper by Bobrowski
% in 2008 - Bayesian Filtering in Spiking Neural Networks: Noise,
% Adaptation and Multisensory Integration

% Modified from filterHybrid3Oct to allow a simpler Snyder filtering
% without requiring as many specific inputs or demanding a linear coupling
% between x1 and x2

% Function implements filtering equations from Snyder in combination with
% use of ODE integration
function outFil = generalSnyFilterBasicDetDel3(q0, x1cap0, Q, lamInp, nPos, Tphot, Sx, nEvents, Tevent, delay, calcBools, posDirSpace, intenSpace, modtype)

% Determine boolean settings
calcDirec = calcBools(1);
calcPos = calcBools(2);

% Obtain lam matrix and lenQ
lam = diag(lamInp.lamDiag);
lenQ = length(Q);

% Obtain state space of the positions based on nPos as well as initialise
% other needed variables <-------------------------------------------------
Spos = diag(0:(nPos-1));
poscapset = cell(1, 1);
poscap = zeros(nEvents, 1);
diagPos = diag(Spos);

% Declare variables and initialise
q = zeros(nEvents, lenQ);
q(1, :) = q0;
t = zeros(nEvents, 1);
x1cap = zeros(nEvents, 1);
x1cap(1) = x1cap0;
qODE = zeros(nEvents, lenQ);
qODE(1, :) = q0;

% Initialise directional estimate variables
dcap = zeros(nEvents, 1);
dcapset = cell(1, 1);

% Cell to save output of ODE solver and set options
Qset = cell(1, 1);
Tset = cell(1, 1);
x1capset = cell(1, 1);
options = odeset('NonNegative', 1:lenQ);

% Check that the position space (which codes for left and right) is even
if rem(length(posDirSpace), 2) == 1
    assignin('base', 'posDirSpace', posDirSpace);
    error('Code assumes an even no. states');
end
    
% Obtain the initial estimate of the direction
if calcDirec
    Pleft_ind = getPleftIndices(modtype, intenSpace, posDirSpace);
    dcap(1) = 1 - 2*sum(q0(Pleft_ind));
end

% Obtain the initial estimate of the position by first obtaining the
% initial position probability vector Ppos
if calcPos
    pos_ind = getPosIndices(lenQ, nPos);
    Ppos = zeros(1, nPos);
    for i = 1:nPos
        Ppos(i) = sum(q0(pos_ind{i}));
    end
    poscap(1) = sum(Ppos*Spos, 2);
end

% Loop across dt values integrating the Snyder equations with the ode solver
% and then applying jump corrections due to x2 events
for i = 2:nEvents
    
    % Solve linear ODEs continuously with setting of options <-----------
    [tset qset] = ode113(@(ts, y) odeSnyLinearBasic(ts, y, Q, lam),...
        [Tevent(i-1) Tevent(i)], q(i-1, :)', options);
    
    % Normalise the posterior probabilities
    for j = 1:size(qset, 1)
        qset(j, :) = qset(j, :)/(sum(qset(j, :)));
    end
            
    % Obtain state, position and directional estimates from qset output
    x1capset{i} = sum(qset*Sx, 2);
    if calcDirec
        % This is the expected direction based on a [-1 1] configuration
        Pleft = sum(qset(:, Pleft_ind), 2);
        dcapset{i} = ones(size(Pleft)) - 2*Pleft;
    end
    if calcPos
        % This is the expected position based on a [0,...,nPos-1] form
        Ppos = cell(1, nPos);
        poscapset{i} = zeros(size(qset, 1), 1);
        for jj = 1:nPos
            Ppos{jj} = sum(qset(:, pos_ind{jj}), 2);
            poscapset{i} = poscapset{i} + diagPos(jj)*Ppos{jj};
        end
    end
    
    % Assign output value at event times
    q(i, :) = qset(end, :);
    t(i) = tset(end);
    qODE(i, :) = q(i, :);
    Qset{i} = qset;
    Tset{i} = tset;
    
    % Check raw ODE q output
    if any(q(i, :) < -10^-8)
        assignin('base', 'qODE', q(i, :));
        error(['qODE distribution has negative entries at i =' num2str(i)]);
    end
    if max(abs(sum(q(i, :)) - 1)) > 10^-4
        assignin('base', 'qODErr', q(i, :));
        disp(['qODE distribution does not sum to 1 at i = ' num2str(i)]);
    end
    
    % Calculate perturbation on q due to jump
    q(i, :) = calcPertEvent(q(i, :), lamInp.lamDiagIndiv, Tevent(i), Tphot, nPos);
    x1cap(i) = sum(q(i, :)*Sx, 2);
    if calcDirec
        Pleft = sum(q(i, Pleft_ind), 2);
        dcap(i) = ones(size(Pleft)) - 2*Pleft;
    end
    if calcPos
        % This is the expected position based on a [0,...,nPos-1] form
        Ppos = zeros(1, nPos);
        for jj = 1:nPos
            Ppos(jj) = sum(q(i, pos_ind{jj}), 2);
        end
        poscap(i) = sum(Ppos*Spos, 2);
    end
    
    % Display progress and debug set of results from ODE solver
    disp(['Finished iteration: ' num2str(i-1) ' of ' num2str(nEvents-1)]);
    
end

% Apply Kolmogorov exponential solution for posterior evolution -
% assuming a delay has already been added to the event times
expSol = expm(delay*Q);
QsetOld = Qset;
x1capsetOld = x1capset;
dcapsetOld = dcapset;
poscapsetOld = poscapset;

% Update the estimates according to the open loop segment from Kolmogorov
for i = 2:length(x1capset)
    qtemp = Qset{i};
    qtemp = qtemp*expSol;
    Qset{i} = qtemp;
    x1capset{i} = sum(qtemp*Sx, 2);
    if calcDirec
        Pleft = sum(qtemp(:, Pleft_ind), 2);
        dcapset{i} = ones(size(Pleft)) - 2*Pleft;
    end
    if calcPos
        Ppos = cell(1, nPos);
        poscapset{i} = zeros(size(qtemp, 1), 1);
        for jj = 1:nPos
            Ppos{jj} = sum(qtemp(:, pos_ind{jj}), 2);
            poscapset{i} = poscapset{i} + diagPos(jj)*Ppos{jj};
        end
    end
end

% Assign output data to a single structure
outFil.q = q;
outFil.t = t;
outFil.x1cap = x1cap;
outFil.dcap = dcap;
outFil.poscap = poscap;
outFil.qODE = qODE;
outFil.Tevent = Tevent;
outFil.Qset = Qset;
outFil.QsetOld = QsetOld;
outFil.Tset = Tset; 
outFil.x1capset = x1capset;
outFil.x1capsetOld = x1capsetOld;
outFil.dcapset = dcapset;
outFil.dcapsetOld = dcapsetOld;
outFil.poscapset = poscapset;
outFil.poscapsetOld = poscapsetOld;


% Sub-function to calculate the normalised linear solution at the event
% times (note this is done after the general normalisation but in fact does
% not matter as the normalisations cancel out)
function qplus = calcPertEvent(qminus, lamDiagIndiv, TevCurr, Tphot, nPos)

% Determine which stream the current event came from (this will act as the
% initial condition for the next iteration)
count = 0;
for i = 1:nPos
    if ismember(TevCurr, Tphot{i})
        currStream = i;
        count = count + 1;
    end
end
lamIndiv = diag(lamDiagIndiv(currStream, :));

% Check that the current event only belongs to one stream
if count > 1
    assignin('base', 'TevCurr', TevCurr);
    assignin('base', 'Tphot', Tphot);
    error('The current event belongs to multiple streams');
end

% Event update based on which stream the observed event originates as this
% determines the rate matrix that is applied
qplus = qminus*lamIndiv./sum(qminus*lamIndiv);


% Sub-function calculates the indices of the probability vector that need
% to be summed to caculate the posterior left direction probability
function Pleft_ind = getPleftIndices(modtype, intenSpace, posDirSpace)

% Get state sizes
nPos = length(posDirSpace)/2;
nIntens = length(intenSpace);

switch(modtype)
    case 1
        % Normal model of constant intensity moving dot with directional
        % trains - first half of all states code the same direction
        discrimin = posDirSpace(nPos);
        Pleft_ind = 1:discrimin+1;
    case 2
        % 2D model that is the filter optimised case to be tested with the
        % reverse phi photon data - state redundancy depends on the number
        % of intensities with every other train coding the same direction
        discrimin = posDirSpace(nPos);
        init_ind = 1:discrimin+1;
        len_init = length(init_ind);
        
        % Get the recursive set of states that code the same position
        % noting that they repeat every 2*nPos values
        Pleft_ind = repmat(init_ind, 1, nIntens);
        idalter = 1:len_init;
        for i = 2:nIntens
            idalter = idalter + nPos*ones(size(idalter));
            Pleft_ind(idalter) = Pleft_ind(idalter) + (i-1)*2*nPos*ones(size(idalter));
        end
        
    otherwise
        assignin('base', 'modtype', modtype);
        error('Incorrect model type specified');
end

% Sub function to calculate the positional estimate indices
function pos_ind = getPosIndices(lenQ, nPos)

% Find number of states coding same position and initialise variables
numRepeat = lenQ/nPos;
pos_ind = cell(1, nPos);

% Determine the repeat and positional ids
repeat_id = 0:nPos:(numRepeat-1)*nPos;
temp = ones(size(repeat_id));
for i = 1:nPos
    pos_ind{i} = i*temp + repeat_id;
end