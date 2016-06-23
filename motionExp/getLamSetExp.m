% Note - if an even no. intensities are supplied then a gray dot or bar
% moving on a gray background is allowed else it is as in motionDet 5

% Function to derive the rate matrices for the streams being produced from
% the moving dot or bar positions - allows Tuthill stimuli
function lamInp = getLamSetExp(stimType, intMod, intenInp)

% Extract the inputs
cset = intenInp.cset;
ibgnd = intenInp.ibgnd;
nPos = intMod.nPos;
nStates = intMod.nStates;
nIntens = intMod.nIntens;

% Initialise the main rate matrices diagonals
lamDiagIndiv = zeros(nPos, nStates);

% Account for various stimulation models
switch(stimType)
    case 0
        % A dot moving across a background of middle intensity (if even no.
        % of intensities then the gray one is not among possible dot
        % intensities - if supply intenSpace = [0 1] - motionDet 5 case
        for i = 1:nPos
            for j = 1:nIntens
                % Apply appropriate contrasts to MC sections with earlier
                % states coding higher intensities <---------------------
                pos = 2*(j-1)*nPos + i;
                lamDiagIndiv(i, [pos pos + nPos]) = cset(nIntens - (j-1));
            end
            
            % Add background intensity
            lamDiagIndiv(i, :) = lamDiagIndiv(i, :) + ibgnd*ones(size(lamDiagIndiv(i, :)));
        end
        
    case 1
        % A series of bars that behave like the dot in case 0 but now there
        % are alternating intensities - depends on odd and even positions
        for i = 1:nPos
            % Separate if even or odd positions
            if rem(i, 2) == 1
                % Positions 1, 3, 5,...
                pos = 1;
            else
                % Positions 2, 4, 6,...
                pos = 2;
            end
            
            % Assign the alternating intensities
            for j = 1:nIntens
                rpiece = pos:2:nPos;
                range = [rpiece rpiece + nPos];
                range = range + 2*(j-1)*nPos;
                lamDiagIndiv(i, range) = cset(nIntens - (j-1));
            end
            
            % Add background intensity
            lamDiagIndiv(i, :) = lamDiagIndiv(i, :) + ibgnd*ones(size(lamDiagIndiv(i, :)));
        end      
        
    otherwise
        error('Improper type of rate matrix model specified');
end


% Set the function outputs (inputs to the filtering function) after
% calculating the overall intensity matrix
lamDiag = sum(lamDiagIndiv, 1);
lamInp.lamDiagIndiv = lamDiagIndiv;
lamInp.lamDiag = lamDiag;