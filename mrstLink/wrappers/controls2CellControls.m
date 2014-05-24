function [vals] = controls2CellControls(u,schedule)
%
%  divide the controls in cells according to the number of wells
%

% assuming same number of wells during the schedule
nW  = getnW( schedule );
nC = numel(schedule.control);

vals = mat2cell(u,nW*ones(nC,1),1);



end