% This code requires zeros for missing reactions in rdot e.g. if have
% increments of 1 and 5 need zeros in positions for 2, 3 and 4 to account
% for missing indices <---------------------------------------------------

% Simple function to check that the output rates are all members of the
% values specified in the Q matrix
% Assumes only the last reaction is for x2 <------------------------------
function check = checkRateQ(Q, Xdot)

% Extract the x1 reactions and initialise check variable
rdot = Xdot(:, 1:end-1);
nReacs = size(rdot, 2);
check = zeros(1, nReacs);

% Noting the structure [bir dea bir dea]
for i = 1:nReacs
    
    % Obtain correct diagonal of Q matrix
    rdiv = rem(i, 2);
    if rdiv == 1
        % Odd reactions are births
        reac = (i+1)/2;
        Qdiag = diag(Q, reac);      
    else
        % Even reactions are deaths
        reac = i/2;
        Qdiag = diag(Q, -reac);
    end
    
    % Check the Q entries with the rate values (include 0 rate which is not
    % in the Q matrix)
    rateVal = unique(rdot(:, i));
    if ismember(rateVal, [Qdiag; 0])
        check(i) = 1;
    end
end

% Determine overall correctness
sumCheck = sum(check);
if sumCheck == nReacs
    disp('All reactions have the correct values from the Q matrix');
else
    disp(['There are ' num2str(nReacs - sumCheck) ' erroneous reactions']);
end
    
    