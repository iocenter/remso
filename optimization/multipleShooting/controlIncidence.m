function [i] = controlIncidence(stepControl,k)
%
%  PIECEWISE CONTINUOUS CONTROL!
%
%  this function map the simulation step k to the active control index k
%  given the piecewise control enconded in stepControl
%


if isempty(stepControl)
    i = k;
else
    i = stepControl(k);
end

end

