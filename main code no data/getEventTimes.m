% Modification to the checking condition on the yid to dyid as this
% accounts for x1 which has births and deaths that can lead to diff(y(id))
% counts which fail the condition ~all(diff(y(id)) == event)

% Note no event = 0 inv = 0 settings as would need to account for the first
% 0 term in dy <------------------------------- needs to be done

% Simple function to separate a specific type of event - either births or
% deaths and to return the event times assuming that the input is the event
% times of all events as is the case for Gillespie simulations
function [Ty perc] = getEventTimes(T, y, type)

% Obtain transitions specific to y and pad with a zero to keep indices
dy = [0; diff(y)];

% Extract event indices based on type
list = {'birth', 'death', 'both'};
idtype = strmatch(type , list);

% Assign event indicators
switch(idtype)
    case 1
        % Search for +1 increments
        event = 1;
        inv = 0;
    case 2
        % Search for -1 increments
        event = -1;
        inv = 0;
    case 3
        % Search for +1 and -1 increments (not 0 increments)
        event = 0;
        inv = 1;
    otherwise
        assignin('base', 'type', eventType);
        error('Unsupported event type specified');
end

% Find index points corresponding to event indicator and invert if set
id = find(dy == event);
if inv
    % Remove the event = 0 indices from total set
    idAll = 1:length(T);
    id = setdiff(idAll, id);
end

% Check that indices are correct and assign outputs
dyid = dy(id);
if isempty(id)
    % Default settings for empty case
    Ty = [0 -1];
    perc = -1;
    warning('Mat:eventType',['No events of type = ' type ' found']);
else
    % If non-empty check validity of events based on inv
    if ~inv
        if ~all(dyid == event)
            assignin('base', 'y', y);
            assignin('base', 'id', id);
            error('Incorrect extraction of events');
        end
    else
        cond = (dyid == 1) + (dyid == -1);
        if ~all(cond == 1)
            assignin('base', 'y', y);
            assignin('base', 'id', id);
            error('Incorrect extraction of events');
        end
    end
    
    % Obtain corresponding times and %frequency of specified event
    Ty = T(id);
    perc = 100*length(Ty)/length(T);
end