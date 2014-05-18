function [ vals, type, control ] = wells2Values(W)
%
% From the well structure, extract the control value, the type of control,
% and the sign
%
vals = vertcat(W.val);

if nargout > 1
    control = cell(numel(W),1);
    
    for k = 1:numel(W)
        control{k} = upper(W(k).type);
    end
    
    type = {W.sign}';
end



end