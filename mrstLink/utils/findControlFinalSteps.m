function [ steps ] = findControlFinalSteps( stepControl )
%
%  find the steps where the controls are switching!
% stepControl =  schedule.step.control

steps = find([arrayfun(@(x,xn)x~=xn,stepControl(1:end-1),stepControl(2:end));true]);


end

