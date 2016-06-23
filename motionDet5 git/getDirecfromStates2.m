% IMPORTANT - renamed direction from states as this is more sensible

% IMPORTANT - case 2 has not been tested yet <-----------------------------

% Modified to allow calculation of direction even if a 2D MC were used -
% note that for reverse phi a 1D MC is used as usual (i.e. set string input
% as 1D)

% Warning - only use this function to obtain the direction from the true
% state data - do not apply to xcap data as this would not be optimal -
% instead use the posterior probabilities

% Simple function to produce directional estimates which accounts for both
% normal constant intensity dots and reverse phi type stimuli - in this
% case the 
function [direc discrimin] = getDirecfromStates2(x1, posDirSpace, strType)

% Check that the position is within the posDirSpace
maxPos = max(x1);
minPos = min(x1);
if maxPos > posDirSpace(end) || minPos < posDirSpace(1)
    assignin('base', 'x1', x1);
    assignin('base', 'posDirSpace', posDirSpace);
    error('The position has values outside the limits of its space');
end

% Ensure that the number of states is even
if rem(length(posDirSpace), 2) == 1
    assignin('base', 'posDirSpace', posDirSpace);
    error('Code assumes an even no. states');
end

% Determine the type of data model used
strModels = {'1D', '2D'};
type = find(strcmp(strType, strModels));

% Obtain the mapping to the directional [-1 1] space
switch(type)
    case 1
        % 1D models - obtain state which gives the threshold between directions
        discrimin = posDirSpace(length(posDirSpace)/2);
        
        % Convert position to direction under assumption of states representing
        % both position and direction
        direc = x1;
        direc(x1 > discrimin) = 1;
        direc(x1 <= discrimin) = -1;
    case 2
        % 2D models - no simple threshold state but rather a recursive set
        % of states all represent 1 direction
        discrimin = posDirSpace(length(posDirSpace)/2);
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
        
        % Obtain the direction from the appropriate states - note posDirSpace
        % here includes duplicate states due to the intensity dimension of
        % the 2D MC <-----------------------------------------------------
        direc = x1;
        direc(ismember(x1, posDirSpace(Pleft_id))) = 1;
        direc(~ismember(x1, posDirSpace(Pleft_id))) = -1;
    otherwise
        assiginin('base', 'strType', strType);
        error('Model form applied is not supported');
end