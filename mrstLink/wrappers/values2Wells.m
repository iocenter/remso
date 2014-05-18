function [ W ] = values2Wells( vals,W )
% set the values to the well control values

for k = 1:numel(vals)
    W(k).val = vals(k);
end


end

