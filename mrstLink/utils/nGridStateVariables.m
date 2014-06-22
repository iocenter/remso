function [ nSG] = nGridStateVariables( comp )
% Returns the number of state variables per grid-block.  This is valid
% implemented for OW and OWG


if ~comp.gas && ~comp.polymer && ~(comp.T || comp.MI)
    nSG = 2; % OW
elseif comp.gas && comp.oil && comp.water
    nSG = 3; % OWG
else
    error('not supported')
end


end

