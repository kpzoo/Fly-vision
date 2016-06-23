% Method assumes that the deleted events involve increments of only +1 so
% there is no account for +n increments yet

% Function to delete data points from a photon stream and create a data set
% for use with the Snyder filter
function [Tdel Xdel eventNo] = deleteDataPts(T, X, probDel, del_meth)

% Simple dimensional check on inputs
if ~all(size(T) == size(X(:, 1)))
    assignin('base', 'Terr', T);
    assignin('base', 'x2err', x2);
    error('The dimensions of T and x2 are inconsistent');
end

% Find the x2 event indices and times
x1 = X(:, 1);
x2 = X(:, 2);
idjump = find(diff([0; x2]) == 1);
Tjump = T(idjump);

% For every x2 event draw from a distribution whether that event should
% remain or be deleted
for i = 1:length(Tjump)
    % Check that the x2 event exists
    idevent = idjump(i);
    if idevent ~= 0
        if x2(idevent) - x2(idevent - 1) ~= 1
            error(['Incorrect event index at i = ' num2str(i)]);
        end
    end
    % Perform deletion based on specified method and probability
    isDel = isDeleted(del_meth, probDel);
    if isDel
        % For deletion to be sensible need to account for surrounding
        % values i.e. x2 does not increment
        x2(idevent+1:end) = x2(idevent+1:end) - (x2(idevent) - x2(idevent-1)); 
        
        % Deletion means no change from previous x2 value which is
        % equivalent to removing that row in X and T
        x2(idevent) = NaN;
        x1(idevent) = NaN;
        T(idevent) = NaN;
    end
end

% Actually delete the rows corresponding to the NaNs and assign outputs
assignin('base', 'x2NaN', x2);
x2del = x2(~isnan(x2));
x1del = x1(~isnan(x1));
Xdel = [x1del x2del];
Tdel = T(~isnan(T));


% Comparison of number of events after deletion
eventNo.orig = length(idjump);
idnewjump = find(diff([0; x2del]) == 1);
assignin('base', 'idnewjump', idnewjump);
eventNo.del = length(idnewjump);
disp(['Number of x2 events before deletion = ' num2str(eventNo.orig)]);
disp(['Number of x2 events after deletion = ' num2str(eventNo.del)]);
if eventNo.orig == eventNo.del && probDel ~= 0
    warning('Mat:noDel', 'No events were deleted');
end


%%
% Method to perform deletion based on a probability and distribution type
function isDel = isDeleted(del_meth, probDel)

% Specified method controls the distribution
switch(del_meth)
    case 0
        % Bernoulli deletion of x2 events (Poisson thinning)
        test = rand;
        if test <= probDel
            isDel = 1;
        else
            isDel = 0;
        end
        
    otherwise
        error('Specified deletion distribution not implemented');
end