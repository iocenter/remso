function [Jac] = catJacs(J,xDims,vDims,uDims,withAlgs)

Jac.Jx = cell2mat(vertcat(J(:).Jx));
m = size(Jac.Jx,1);

Jac.Jx = mat2cell(Jac.Jx,m,xDims);
Jac.Ju = mat2cell(cell2mat(vertcat(J(:).Ju)),m,uDims);

if withAlgs
    Jac.Jv = mat2cell(cell2mat(vertcat(J(:).Jv)),m,vDims);
end


end