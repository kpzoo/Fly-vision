% Modified to calculate the directional estimate directly from the q
% vectors using the ecxtra boolean input - calcDirec. This assumes the
% first half states code 1 direction <*************************************

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
function outFil = generalSnyFilterBasic(q0, x1cap0, Q, lamInp, nPos, Tphot, Sx, nEvents, Tevent, calcDirec, posSpace)


% outFil = generalSnyFilterBasic(lamcap0, q0, B, T, S, x1cap0, coeff, lenS,
% birStr2, params, x2, octave)

% Obtain lam matrix and lenQ
lam = diag(lamInp.lamDiag);
lenQ = length(Q);

% Declare variables and initialise
q = zeros(nEvents, lenQ);
q(1, :) = q0;
t = zeros(nEvents, 1);
x1cap = zeros(nEvents, 1);
x1cap(1) = x1cap0;
dcap = zeros(nEvents, 1);
qODE = zeros(nEvents, lenQ);
qODE(1, :) = q0;

% Cell to save output of ODE solver and set options
Qset = cell(1, 1);
Tset = cell(1, 1);
x1capset = cell(1, 1);
dcapset = cell(1, 1);
options = odeset('NonNegative', 1:lenQ);

% Obtain state which gives the threshold between the directions and ensure
% that the number of states is even
if calcDirec
    if rem(length(posSpace), 2) == 1
        assignin('base', 'posSpace', posSpace);
        error('Code assumes an even no. states');
    end
    discrimin = posSpace(length(posSpace)/2);
    dcap(1) = 1 - 2*sum(q0(1:discrimin+1));
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
            
    % Obtain state and directional estimates from qset output
    x1capset{i} = sum(qset*Sx, 2);
    if calcDirec
        % This is the expected direction based on a [-1 1] configuration -
        % assumes the first half states represent left <*******************
         Pleft = sum(qset(:, 1:discrimin+1), 2);
         dcapset{i} = ones(size(Pleft)) - 2*Pleft;
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
        Pleft = sum(q(i, 1:discrimin+1), 2);
        dcap(i) = ones(size(Pleft)) - 2*Pleft;
    end
    
    % Display progress and debug set of results from ODE solver
    disp(['Finished iteration: ' num2str(i-1) ' of ' num2str(nEvents-1)]);
    
end

% Assign output data to a single structure
outFil.q = q;
outFil.t = t;
outFil.x1cap = x1cap;
outFil.dcap = dcap;
outFil.qODE = qODE;
outFil.Tevent = Tevent;
outFil.Qset = Qset;
outFil.Tset = Tset; 
outFil.x1capset = x1capset;
outFil.dcapset = dcapset;


% Sub function to calculate the normalised linear solution at the event
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

