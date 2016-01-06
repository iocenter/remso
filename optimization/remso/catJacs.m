function [Jac,m] = catJacs(J,xDims,vDims,uDims)

Jac.Jx = cell2mat(vertcat(J(:).Jx));
m = size(Jac.Jx,1);

Jac.Jx = mat2cell(Jac.Jx,m,xDims);

if isfield(J,'Ju')
    Jac.Ju = mat2cell(cell2mat(vertcat(J(:).Ju)),m,uDims);
end

if isfield(J,'Jv')
    Jac.Jv = mat2cell(cell2mat(vertcat(J(:).Jv)),m,vDims);
end


end