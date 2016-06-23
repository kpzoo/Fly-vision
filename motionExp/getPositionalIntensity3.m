% Modified to include a 3rd case which full wave rectifies (in some sense)
% the intensity so both black and white are white

% Modified to include the flashing light intensity to allow for reverse phi
% modelling with contrast flipping on a grey background

% Function to produce the intensity at each of several position based on
% different light models
function R = getPositionalIntensity3(xcurr, modelType, intenParams, nPos, Rold, xold, nStates)

% Define the observed process rate stream - only a birth rate for each
R = ones(1, nPos);
redPos = nStates/nPos;
if round(redPos) ~= redPos
    assignin('base', 'redPos', redPos);
    error('State cardinality is not an integer multiple of position cardinality');
end

switch(modelType)
    case 1
        % Check the expected ratio for the model
        if redPos ~= 2
            assignin('base', 'redPos', redPos);
            error('Incorrect state to position ratio');
        end
        
        % Light model assigns a constant intensity to the dot and has no
        % background light so light only exists at the dot position. There
        % are M/2 positions for a M state Markov chain
        i0 = intenParams.i0;
        if ismember(xcurr, 0:(nPos - 1))
            R(xcurr+1) = i0;
        else
            % The second state that code the same position for a different
            % direction is assigned here
            R(xcurr+1 - nPos) = i0;
        end
        
    case 2
        % Model is similar to case 1 but now there is a background
        % intensity and flashing occurs based on the last intensity
        i0 = intenParams.i0;
        ic = intenParams.ic;
        
        % Apply the background intensity and obtain the positional ids
        R = i0*R;
        posidOld = getPosfromState(xold, nPos, redPos);
        posid = getPosfromState(xcurr, nPos, redPos);
        
        % Apply the constrast to the already grey R (i.e. it has the
        % background intensities) in a toggle switch manner
        if Rold(posidOld) == i0 + ic
            R(posid) = R(posid) - ic;
            % Modification to switch the last positional intensity
            % <------------------------------------------------------------
            R(posidOld) = R(posidOld) + ic;
%             R(posidOld) = R(posidOld) - ic;
        elseif Rold(posidOld) == i0 - ic
            R(posid) = R(posid) + ic;
            % Modification to switch the last positional intensity
            % <------------------------------------------------------------
            R(posidOld) = R(posidOld) - ic;
%             R(posidOld) = R(posidOld) + ic;
        else
            assiginin('base', 'Rold', Rold);
            error('Method can only handle 2 intensities when the dot is in a position');
        end
        
    case 3
        % Full wave rectified version of the flashing reverse phi model
        i0 = intenParams.i0;
        ic = intenParams.ic;
        
        % Apply the background intensity and obtain the positional ids
        R = i0*R;
        posidOld = getPosfromState(xold, nPos, redPos);
        posid = getPosfromState(xcurr, nPos, redPos);
        
        % Full wave rectify the intensities
        if Rold(posidOld) == i0 + ic || Rold(posidOld) == i0 - ic
            R(posid) = R(posid) + ic;
%             R(posidOld) = i0 - ic;
        else
            assiginin('base', 'Rold', Rold);
            error('Method can only handle 2 intensities when the dot is in a position');
        end
        
    otherwise
        assignin('base', 'modelType', modelType);
        error('Light model id specified is invalid');
end

% Sub function that converts the x values to position - assumes a
% representation as in case 2 where redPos states code for the same
% position in space
function posid = getPosfromState(x, nPos, redPos)

% Account for the different states that code the same position
for i = 1:redPos
    if ismember(x, ((i-1)*nPos):(i*nPos - 1))
        % The positional id is the state within a set of nPos values +1
        posid = x - (i-1)*nPos + 1;
    end
end
        


