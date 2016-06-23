% Function to derive the rate matrices for the streams being produced from
% the moving dot positions
function lamInp = getLamManyStreams(type, nPos, nStates, intenInp)

switch(type)
    case 1
        % Assumption of 2 train MC means nStates should be 2*nPos
        if nStates ~= 2*nPos
            assiginin('base', 'nPos', nPos);
            assiginin('base', 'nStates', nStates);
            error('Incorrect light model form specified');
        end
        
        % Model involves a constant intensity dot moving across the
        % Markov chain of positions where each state codes both position
        % and direction of motion
        lamDiagIndiv = zeros(nPos, nStates);
        for i = 1:nPos
            % Assume the coding of position and direction via a state then a
            % counting process represents 2 states e.g. 0 and 3 for space 0:5
            lamDiagIndiv(i, [i nPos+i]) = intenInp;
        end
        lamDiag = sum(lamDiagIndiv, 1);
        
        % Set the input structures to the filtering function
        lamInp.lamDiagIndiv = lamDiagIndiv;
        lamInp.lamDiag = lamDiag;
        
    case 2
        % Assumption of a 2D MC means nStates should be 4*nPos
        if nStates ~= 4*nPos
            assiginin('base', 'nPos', nPos);
            assiginin('base', 'nStates', nStates);
            error('Incorrect light model form specified');
        end
        
        % Model involves a 2D MC in which the states code position,
        % direction and intensity
        lamDiagIndiv = zeros(nPos, nStates);
        i0 = intenInp.i0;
        ic = intenInp.ic;
        
        % At each position account for the background (i0) and constrast
        % intensity (ic) expected
        for i = 1:nPos
            % Positive contrast section of MC
            lamDiagIndiv(i, [i nPos+i]) = ic;
            % Negative contrast section of MC
            lamDiagIndiv(i, [2*nPos + i 3*nPos+i]) = -ic;
            % Background intensity
            lamDiagIndiv(i, :) = lamDiagIndiv(i, :) + i0*ones(size(lamDiagIndiv(i, :))); 
        end
        lamDiag = sum(lamDiagIndiv, 1);
        
        % Set the input structures to the filtering function
        lamInp.lamDiagIndiv = lamDiagIndiv;
        lamInp.lamDiag = lamDiag;        
        
    otherwise
        error('Improper type of rate matrix model specified');
end