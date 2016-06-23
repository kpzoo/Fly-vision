% Function to produce the intensity at each of several position based on
% different light models
function R = getPositionalIntensity(xcurr, modelType, intenParams, nPos)

% Define the observed process rate stream - only a birth rate for each
R = zeros(1, nPos);

switch(modelType)
    case 1
        % Light model assigns a constant intensity to the dot and has no
        % background light so light only exists at the dot position. There
        % are M/2 positions for a M state Markov chain
        i0 = intenParams.i0;
        if ismember(xcurr, 0:(nPos - 1))
            R(xcurr+1) = i0;
        else
            % The second state that codes the same position for a different
            % direction is assigned here
            R(xcurr+1 - nPos) = i0;
        end
        
    case 2
        % Model is exactly the same as in case 1 but now there is a
        % constant intensity in the background on all channels
        i0 = intenParams.i0;
        ic = intenParams.ic;
        R = ic*ones(size(R));
        if ismember(xcurr, 0:(nPos - 1))
            R(xcurr+1) = i0 + ic;
        else
            % The second state that codes the same position for a different
            % direction is assigned here
            R(xcurr+1 - nPos) = i0 + ic;
        end     
        
    otherwise
        assignin('base', 'modelType', modelType);
        error('Light model id specified is invalid');
end
